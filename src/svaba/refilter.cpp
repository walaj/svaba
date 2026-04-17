// refilter.cpp — re-apply LOD / PASS filtering to an existing bps.txt.gz
//
// The goal here is to avoid re-running assembly when the user only wants to
// tweak filtering thresholds. We read bps.txt.gz line-by-line, reconstruct
// BreakPoint objects via BreakPoint(line, sc*), rescore, redump, and
// regenerate the VCFs. The `dbsnp` column is round-tripped verbatim from
// the dump; DBSnp re-querying has been sunset.
//
// Requirements on the input file:
//   - Header line present (tab-delimited). Per-sample columns are expected to
//     begin with 't' (tumor) or 'n' (normal). Core column count is either 41
//     (legacy pre-SvABA2.0) or 51 (SvABA2.0 refilter-aware dump). Both are
//     accepted.
//   - Column formatting matches BreakPoint::toFileString(). The BreakPoint
//     parse constructor tolerates sentinel values ("x", "NA", "nan", empty)
//     for ints/doubles.
//
// Caveats:
//   - Without contigs/reads, assembly-derived gates that re-examine sequence
//     cannot run. We rely on the dumped split/read counts and the scoring
//     primitives in svabaModels.
//   - Older dumps may leave per-end LocalAlignment, cpos, match, split_cov
//     bounds, and contig_len at defaults; those filters will degrade
//     gracefully rather than erroring out.

#include "refilter.h"

#include <getopt.h>
#include <sstream>
#include <iostream>
#include <memory>

#include "gzstream.h"
#include "SeqLib/BamReader.h"

#include "vcf.h"
#include "BreakPoint.h"
#include "SvabaUtils.h"
#include "SvabaLogger.h"
#include "SvabaOptions.h"
#include "SvabaOutputWriter.h"
#include "SvabaSharedConfig.h"


namespace opt {

  static std::string input_file;
  static std::string analysis_id = "refilter";

  static std::string bam;  // any BAM; used only for the sequence header

  static int verbose = 1;

  // Indel probability cutoffs (will be pushed into SvabaOptions before
  // scoring). Defaults mirror the in-tree SvabaOptions defaults so refilter
  // is non-destructive when no --lod* flags are passed.
  static double lod            = 1.0;
  static double lod_db         = 1.0;
  static double lod_somatic    = 0.0;
  static double lod_somatic_db = 2.0;

  // Assumed mean read length when the dump lacks it. modelSelection uses this
  // for per-read error scaling; in practice most short-read BAMs are 101-151.
  static int readlen = 150;
}

enum {
  OPT_LOD,
  OPT_LOD_DB,
  OPT_LOD_SOMATIC,
  OPT_LOD_SOMATIC_DB,
  OPT_READLEN
};


static const char* shortopts = "hi:a:v:G:b:";
static const struct option longopts[] = {
  { "help",                    no_argument,       NULL, 'h' },
  { "input-bps",               required_argument, NULL, 'i' },
  { "bam",                     required_argument, NULL, 'b' },
  { "analysis-id",             required_argument, NULL, 'a' },
  { "verbose",                 required_argument, NULL, 'v' },
  { "lod",                     required_argument, NULL, OPT_LOD },
  { "lod-dbsnp",               required_argument, NULL, OPT_LOD_DB },
  { "lod-somatic",             required_argument, NULL, OPT_LOD_SOMATIC },
  { "lod-somatic-dbsnp",       required_argument, NULL, OPT_LOD_SOMATIC_DB },
  { "readlen",                 required_argument, NULL, OPT_READLEN },
  { NULL, 0, NULL, 0 }
};

static const char *BP_USAGE_MESSAGE =
"Usage: svaba refilter [OPTION] -i bps.txt.gz -b <bam>\n\n"
"  Description: re-run LOD/PASS filtering against an existing bps.txt.gz\n"
"  without re-doing assembly. Samples, coverage, split counts, contig\n"
"  coordinates, and local-alignment flags are read back from the dump.\n"
"\n"
"  General options\n"
"  -v, --verbose                Verbosity level (0-4). Default: 1\n"
"  -h, --help                   Display this help and exit\n"
"  -a, --analysis-id            Analysis ID for output filenames. [refilter]\n"
"  Required input\n"
"  -i, --input-bps              bps.txt.gz produced by `svaba run`.\n"
"  -b, --bam                    Any svaba-input BAM (used for sequence header only).\n"
"  Optional\n"
"      --lod                    LOD for non-REF PASS. [1.0]\n"
"      --lod-dbsnp              LOD for non-REF PASS at DBSnp sites. [1.0]\n"
"      --lod-somatic            LOD for SOMATIC. [0.0]\n"
"      --lod-somatic-dbsnp      LOD for SOMATIC at DBSnp sites. [2.0]\n"
"      --readlen                Assumed read length when dump lacks it. [150]\n"
"\n";

// parse the command line options
static void parseBreakOptions(int argc, char** argv) {
  bool die = false;

  if (argc <= 2)
    die = true;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
      case 'h': die = true; break;
      case 'i': arg >> opt::input_file; break;
      case 'v': arg >> opt::verbose; break;
      case 'a': arg >> opt::analysis_id; break;
      case 'b': arg >> opt::bam; break;
      case OPT_LOD:             arg >> opt::lod;            break;
      case OPT_LOD_DB:          arg >> opt::lod_db;         break;
      case OPT_LOD_SOMATIC:     arg >> opt::lod_somatic;    break;
      case OPT_LOD_SOMATIC_DB:  arg >> opt::lod_somatic_db; break;
      case OPT_READLEN:         arg >> opt::readlen;        break;
    }
  }

  if (opt::input_file.empty()) {
    std::cerr << "ERROR: --input-bps is required\n";
    die = true;
  }
  if (opt::bam.empty()) {
    std::cerr << "ERROR: --bam is required (for the sequence header)\n";
    die = true;
  }

  if (die) {
    std::cerr << "\n" << BP_USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }
}

// Parse a bps.txt sample-column header token like "t001_/path/to.bam" into
// (prefix, path). If there is no underscore, the whole token is treated as
// the prefix.
static std::pair<std::string,std::string>
splitSampleHeader(const std::string& tok) {
  const auto us = tok.find('_');
  if (us == std::string::npos) return { tok, "" };
  return { tok.substr(0, us), tok.substr(us + 1) };
}

void runRefilterBreakpoints(int argc, char** argv) {

  parseBreakOptions(argc, argv);

  if (opt::verbose > 0) {
    std::cerr << "Input bps file:  " << opt::input_file << std::endl
              << "Output analysis id: " << opt::analysis_id << std::endl
              << "    LOD cutoff (non-REF):           " << opt::lod << std::endl
              << "    LOD cutoff (non-REF, at DBSNP): " << opt::lod_db << std::endl
              << "    LOD somatic cutoff:             " << opt::lod_somatic << std::endl
              << "    LOD somatic cutoff (at DBSNP):  " << opt::lod_somatic_db << std::endl;
  }

  if (!SeqLib::read_access_test(opt::input_file)) {
    std::cerr << "ERROR: Cannot read " << opt::input_file << std::endl;
    exit(EXIT_FAILURE);
  }

  // Open the BAM only to grab its header; refilter does not read any reads.
  SeqLib::BamReader bwalker;
  if (!bwalker.Open(opt::bam)) {
    std::cerr << "ERROR: could not open BAM " << opt::bam << std::endl;
    exit(EXIT_FAILURE);
  }
  const SeqLib::BamHeader hdr = bwalker.Header();

  // ------------------------------------------------------------------
  // Read the bps.txt.gz header and recover sample prefixes.
  //
  // Layouts supported:
  //   legacy (41 core cols): ..., contig_conf1, contig_conf2, <samples...>
  //   svaba2 (51 core cols): ..., contig_conf1, contig_conf2,
  //                              cpos1, cpos2, lmatch, rmatch, scov1, scov2,
  //                              local1, local2, ctglen, flipped,
  //                              <samples...>
  //
  // Sample header tokens are of the form "{prefix}_{bampath}" where prefix
  // begins with 't' or 'n'. We locate the sample block by scanning for the
  // first column whose name starts with 't' or 'n' past the fixed-width
  // prefix (col 34+ is safe: no core column from there starts t/n).
  // ------------------------------------------------------------------
  igzstream infile(opt::input_file.c_str(), std::ios::in);

  std::string headerLine;
  if (!std::getline(infile, headerLine)) {
    std::cerr << "ERROR: bps.txt file is empty or unreadable" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::vector<std::string> headerv =
      svabaUtils::tokenize_delimited(headerLine, '\t');

  if (headerv.size() < 41) {
    std::cerr << "ERROR: bps.txt header has only " << headerv.size()
              << " columns; expected at least 41 (legacy) or 51 (svaba2)."
              << " File may be truncated or from an unsupported svaba version."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  size_t sample_start = 0;
  for (size_t i = 34; i < headerv.size(); ++i) {
    if (!headerv[i].empty() &&
        (headerv[i].at(0) == 't' || headerv[i].at(0) == 'n')) {
      sample_start = i;
      break;
    }
  }
  if (sample_start == 0) {
    std::cerr << "ERROR: could not locate sample-name columns in bps.txt header."
              << " No column starting with 't' or 'n' found after the core block."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  // Build a minimal SvabaSharedConfig sufficient to drive BreakPoint parsing
  // and rescoring. SvabaOutputWriter is constructed but not `init`-ed; no
  // files are opened by it. We open our own bps output stream below.
  SvabaLogger       logger;
  SvabaOptions      opts;
  SvabaOutputWriter writer(logger, opts);

  // Push CLI-provided cutoffs into opts; score_indel() reads these.
  opts.lod          = opt::lod;
  opts.lodDb        = opt::lod_db;
  opts.lodSomatic   = opt::lod_somatic;
  opts.lodSomaticDb = opt::lod_somatic_db;

  // Populate opts.bams from the header. The BreakPoint parse constructor
  // iterates sc->opts.bams (a std::map) in sorted key order, so as long as
  // the dump also iterated opts.bams in sorted order (which it does), the
  // per-sample token order is round-trip-safe. Key is the sample prefix
  // (e.g. "t001"); the value (bam path from the header token) is only used
  // here for the reconstructed header line.
  std::vector<std::string> allele_names;  // raw header tokens, for info/debug
  std::vector<std::string> sample_prefixes;
  for (size_t i = sample_start; i < headerv.size(); ++i) {
    allele_names.push_back(headerv[i]);
    auto [pref, path] = splitSampleHeader(headerv[i]);
    if (pref.empty() ||
        (pref.at(0) != 't' && pref.at(0) != 'n')) {
      std::cerr << "ERROR: unexpected sample header token '"
                << headerv[i] << "' (col " << i << "); "
                << "expected to start with 't' or 'n'." << std::endl;
      exit(EXIT_FAILURE);
    }
    // Preserve the original bam path from the header (may be empty if the
    // header token was just a bare prefix). Used only for the re-emitted
    // header line and (incidentally) as the value in the map iterated by
    // the BreakPoint parse constructor.
    opts.bams[pref] = path;
    sample_prefixes.push_back(pref);
  }

  if (opt::verbose > 0) {
    std::cerr << "...bps.txt has " << headerv.size() << " columns; "
              << sample_start << " core, "
              << (headerv.size() - sample_start) << " sample(s): ";
    for (const auto& p : sample_prefixes) std::cerr << p << " ";
    std::cerr << std::endl;
  }

  // Build the SvabaSharedConfig
  SvabaSharedConfig sc(logger, opts, writer);
  sc.header  = hdr;
  sc.readlen = opt::readlen;

  // (DBSnp filtering has been sunset — refilter no longer supports it. The
  // `dbsnp` column in bps.txt.gz is round-tripped verbatim from the dump.)

  // VCF header (minimal)
  VCFHeader vheader;
  vheader.filedate  = svabaUtils::fileDateString();
  vheader.source    = "svaba refilter";
  vheader.reference = "";

  // Open the rewritten bps file and write its header. We re-emit the current
  // BreakPoint::header() layout so downstream tooling always sees the
  // newest column set (the reader above was backward-compatible with the
  // older 41-col dumps, but we prefer to always emit the modern format).
  const std::string new_bps_file = opt::analysis_id + ".bps.txt.gz";
  ogzstream os_bps;
  svabaUtils::fopen(new_bps_file, os_bps);
  os_bps << BreakPoint::header();
  for (auto& p : opts.bams)
    os_bps << "\t" << p.first << "_" << p.second;
  os_bps << "\n";

  // --------------------------- main loop ---------------------------
  size_t line_count = 0;
  size_t parse_errors = 0;
  std::string line;
  while (std::getline(infile, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;  // defensive

    if (opt::verbose > 0 && line_count && (line_count % 100000) == 0)
      std::cerr << "...read " << opt::input_file
                << " at line " << SeqLib::AddCommas(line_count) << std::endl;

    std::unique_ptr<BreakPoint> bp;
    try {
      bp.reset(new BreakPoint(line, &sc));
    } catch (const std::exception& e) {
      ++parse_errors;
      if (parse_errors <= 5) {
        std::cerr << "WARN: parse error at line " << (line_count + 1)
                  << ": " << e.what() << std::endl;
      }
      ++line_count;
      continue;
    }

    // Roll discordant counts from per-sample disc into dc.{tcount,ncount}.
    // The original scorer sets these during assembly; they may not be
    // present in the dump's dc fields, so derive from per-sample `disc`.
    bp->dc.tcount = 0;
    bp->dc.ncount = 0;
    for (auto& [pref, si] : bp->allele) {
      if (!pref.empty() && pref.at(0) == 't') bp->dc.tcount += si.disc;
      else                                    bp->dc.ncount += si.disc;
    }

    // Clear derived fields so scoreBreakpoint() can repopulate from scratch.
    // The parse sets `confidence` from the dump; scoreBreakpoint() asserts
    // it's empty (except for BLACKLIST), so blank it here unless already
    // BLACKLIST (which is a terminal state we want to preserve).
    if (bp->confidence != "BLACKLIST")
      bp->confidence.clear();

    // Re-score
    try {
      bp->scoreBreakpoint();
    } catch (const std::exception& e) {
      std::cerr << "WARN: scoreBreakpoint failed at line " << (line_count + 1)
                << ": " << e.what() << std::endl;
      ++line_count;
      continue;
    }

    // toFileString() needs the BAM header to render chromosome names;
    // read-tracking is no longer a per-call toggle (supporting-read qnames
    // are stored or stripped elsewhere in the write path).
    os_bps << bp->toFileString(hdr) << "\n";
    ++line_count;
  }
  os_bps.close();

  if (opt::verbose > 0) {
    std::cerr << "...refilter read " << line_count << " lines"
              << " (" << parse_errors << " parse error"
              << (parse_errors == 1 ? "" : "s") << ")" << std::endl;
  }

  // ---- regenerate VCFs from the refiltered bps ----
  if (SeqLib::read_access_test(new_bps_file)) {

    if (opt::verbose > 0)
      std::cerr << "...making the primary VCFs (unfiltered and filtered) from file "
                << new_bps_file << std::endl;

    VCFFile snowvcf(new_bps_file, opt::analysis_id, sc, vheader, true,
                    opt::verbose > 0);

    const bool onefile = (allele_names.size() == 1);

    std::string basename = opt::analysis_id + ".svaba.unfiltered.";
    snowvcf.include_nonpass = true;
    snowvcf.writeIndels(basename, false, onefile, hdr);
    snowvcf.writeSVs   (basename, false, onefile, hdr);

    basename = opt::analysis_id + ".svaba.";
    snowvcf.include_nonpass = false;
    snowvcf.writeIndels(basename, false, onefile, hdr);
    snowvcf.writeSVs   (basename, false, onefile, hdr);

  } else {
    std::cerr << "Failed to make VCF. Could not find rewritten bps file "
              << new_bps_file << std::endl;
  }
}
