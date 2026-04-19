// tovcf.cpp
//
// `svaba tovcf` — convert a deduplicated bps.txt.gz into a pair of VCF
// files (one for SVs, one for indels). Target spec: VCFv4.5.
//
// This is the standalone format-conversion entry point. It assumes the
// input bps.txt.gz has already been sorted and deduplicated by
// scripts/svaba_postprocess.sh — so the VCFFile internal dedup pass is
// skipped. All calls (somatic + germline) land in one SV file and one
// indel file, with the SOMATIC INFO flag distinguishing somatic rows
// for downstream `bcftools filter` / grep workflows.
//
// Knobs (with defaults chosen for the new pipeline):
//   --sv-out PATH        output SV VCF (default ${id}.sv.vcf.gz)
//   --indel-out PATH     output indel VCF (default ${id}.indel.vcf.gz)
//   --always-bnd         force paired BND records for every SV
//                        (otherwise intrachrom events collapse into
//                         <DEL>/<DUP>/<INV> symbolic records).
//   --qual MODE          missing | maxlod | sum (default: missing).
//                        Controls the VCF QUAL column. "missing" ('.')
//                        pushes users toward INFO/MAXLOD / INFO/SOMLOD
//                        for filtering; "maxlod" puts 10*MAXLOD there;
//                        "sum" is the legacy behavior.
//   --include-nonpass    also emit records with FILTER != PASS.
//   --dedup              opt back in to the legacy interval-tree dedup
//                        (off by default; the postprocess pipeline
//                         already dedupes upstream).
//   --plain              write plain .vcf instead of bgzip'd .vcf.gz.
//                        (Appends .vcf; by default outputs .vcf.gz.)
//   -v / --verbose       chatty diagnostics.
//   -h / --help          show help.
//
// Required:
//   -i, --input-bps      deduped bps.txt.gz (from svaba_postprocess)
//   -b, --bam            BAM used to recover the chromosome name/length
//                        table; no reads are actually read.
//   -a, --id-string      analysis id for output names.

#include <getopt.h>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "gzstream.h"
#include "SeqLib/BamReader.h"

#include "BreakPoint.h"
#include "SvabaLogger.h"
#include "SvabaOptions.h"
#include "SvabaOutputWriter.h"
#include "SvabaSharedConfig.h"
#include "SvabaUtils.h"
#include "vcf.h"

namespace {

struct ToVcfOpts {
  std::string input_file;
  std::string bam;
  std::string analysis_id;
  std::string sv_out;       // explicit path; filled from analysis_id if empty
  std::string indel_out;
  bool        always_bnd   = false;
  bool        include_nonpass = false;
  bool        dedup         = false;
  bool        gzip          = true;
  int         verbose       = 0;
  QualMode    qual_mode     = QualMode::MISSING;
};

constexpr const char* TOVCF_USAGE =
"Usage: svaba tovcf -i BPS.txt.gz -b BAM -a ID [options]\n"
"\n"
"  Convert a deduplicated bps.txt.gz into VCFv4.5 output. Emits a SV\n"
"  VCF and an indel VCF; all calls (somatic+germline) go into one file\n"
"  each, with the SOMATIC INFO flag distinguishing somatic rows.\n"
"\n"
"  Required:\n"
"    -i, --input-bps FILE        deduplicated bps.txt.gz\n"
"    -b, --bam FILE              BAM used to source the chrom name/length table\n"
"    -a, --id-string STR         analysis id for default output names\n"
"\n"
"  Output:\n"
"        --sv-out FILE           override SV output path\n"
"        --indel-out FILE        override indel output path\n"
"        --plain                 write plain .vcf instead of .vcf.gz\n"
"\n"
"  Format knobs:\n"
"        --always-bnd            force paired BND for every SV\n"
"        --qual MODE             missing | maxlod | sum (default: missing)\n"
"        --include-nonpass       include records with FILTER != PASS\n"
"        --dedup                 re-run legacy interval-tree dedup on input\n"
"\n"
"  Misc:\n"
"    -v, --verbose\n"
"    -h, --help\n";

constexpr const char* SHORTOPTS = "hi:b:va:";

// Long-only flags are dispatched via the `val` column of longopts; keep
// these unique across the option surface.
constexpr int OPT_SV_OUT          = 2001;
constexpr int OPT_INDEL_OUT       = 2002;
constexpr int OPT_ALWAYS_BND      = 2003;
constexpr int OPT_QUAL            = 2004;
constexpr int OPT_INCLUDE_NONPASS = 2005;
constexpr int OPT_DEDUP           = 2006;
constexpr int OPT_PLAIN           = 2007;

const struct option LONGOPTS[] = {
  { "help",            no_argument,       nullptr, 'h' },
  { "input-bps",       required_argument, nullptr, 'i' },
  { "bam",             required_argument, nullptr, 'b' },
  { "id-string",       required_argument, nullptr, 'a' },
  { "verbose",         no_argument,       nullptr, 'v' },
  { "sv-out",          required_argument, nullptr, OPT_SV_OUT },
  { "indel-out",       required_argument, nullptr, OPT_INDEL_OUT },
  { "always-bnd",      no_argument,       nullptr, OPT_ALWAYS_BND },
  { "qual",            required_argument, nullptr, OPT_QUAL },
  { "include-nonpass", no_argument,       nullptr, OPT_INCLUDE_NONPASS },
  { "dedup",           no_argument,       nullptr, OPT_DEDUP },
  { "plain",           no_argument,       nullptr, OPT_PLAIN },
  { nullptr, 0, nullptr, 0 }
};

// Parse --qual {missing|maxlod|sum}. Case-insensitive prefixes accepted.
QualMode parse_qual_mode(const std::string& s) {
  if (s == "missing" || s == "."   || s == "none" ) return QualMode::MISSING;
  if (s == "maxlod"  || s == "max" || s == "mlod") return QualMode::MAXLOD_PHRED;
  if (s == "sum"     || s == "legacy")              return QualMode::SUM_LO_PHRED;
  std::cerr << "ERROR: --qual must be one of: missing | maxlod | sum (got '"
            << s << "')\n";
  std::exit(EXIT_FAILURE);
}

ToVcfOpts parse_cli(int argc, char** argv) {
  ToVcfOpts o;
  bool die = (argc <= 2);

  optind = 1;  // reset getopt state (svaba.cpp consumes argv[1] before dispatch)
  for (int c; (c = getopt_long(argc, argv, SHORTOPTS, LONGOPTS, nullptr)) != -1;) {
    switch (c) {
      case 'h': die = true; break;
      case 'i': o.input_file  = optarg ? optarg : ""; break;
      case 'b': o.bam          = optarg ? optarg : ""; break;
      case 'a': o.analysis_id  = optarg ? optarg : ""; break;
      case 'v': o.verbose      = 1; break;
      case OPT_SV_OUT:          o.sv_out    = optarg ? optarg : ""; break;
      case OPT_INDEL_OUT:       o.indel_out = optarg ? optarg : ""; break;
      case OPT_ALWAYS_BND:      o.always_bnd = true; break;
      case OPT_QUAL:            o.qual_mode = parse_qual_mode(optarg ? optarg : ""); break;
      case OPT_INCLUDE_NONPASS: o.include_nonpass = true; break;
      case OPT_DEDUP:           o.dedup = true; break;
      case OPT_PLAIN:           o.gzip = false; break;
      default: die = true; break;
    }
  }

  if (!die) {
    if (o.input_file.empty()) {
      std::cerr << "ERROR: -i / --input-bps is required\n"; die = true;
    }
    if (o.bam.empty()) {
      std::cerr << "ERROR: -b / --bam is required (for chrom name/length table)\n";
      die = true;
    }
    if (o.analysis_id.empty()) {
      std::cerr << "ERROR: -a / --id-string is required\n"; die = true;
    }
  }

  if (die) {
    std::cerr << "\n" << TOVCF_USAGE;
    std::exit(1);
  }

  const std::string ext = o.gzip ? ".vcf.gz" : ".vcf";
  if (o.sv_out.empty())    o.sv_out    = o.analysis_id + ".sv"    + ext;
  if (o.indel_out.empty()) o.indel_out = o.analysis_id + ".indel" + ext;

  return o;
}

// Copy of the refilter.cpp helper — splits "prefix_/path/to.bam" header
// tokens into (prefix, path).
std::pair<std::string,std::string> splitSampleHeader(const std::string& tok) {
  const auto us = tok.find('_');
  if (us == std::string::npos) return { tok, "" };
  return { tok.substr(0, us), tok.substr(us + 1) };
}

} // namespace

// Forward declared in svaba.cpp dispatch.
void runToVCF(int argc, char** argv) {

  const ToVcfOpts o = parse_cli(argc, argv);

  if (o.verbose > 0) {
    std::cerr << "svaba tovcf:\n"
              << "  input:       " << o.input_file  << "\n"
              << "  bam:         " << o.bam         << "\n"
              << "  sv-out:      " << o.sv_out      << "\n"
              << "  indel-out:   " << o.indel_out   << "\n"
              << "  qual:        "
              << (o.qual_mode == QualMode::MISSING      ? "missing"
                : o.qual_mode == QualMode::MAXLOD_PHRED ? "maxlod"
                :                                          "sum")
              << "\n"
              << "  sv-format:   "
              << (o.always_bnd ? "BND_ALWAYS" : "SYMBOLIC_WHEN_OBVIOUS") << "\n"
              << "  dedup:       " << (o.dedup ? "yes" : "no (input assumed pre-deduped)") << "\n"
              << "  nonpass:     " << (o.include_nonpass ? "yes" : "no") << "\n";
  }

  if (!SeqLib::read_access_test(o.input_file)) {
    std::cerr << "ERROR: cannot read " << o.input_file << "\n";
    std::exit(EXIT_FAILURE);
  }

  // ---- 1. Open BAM only to grab the reference chromosome header ----------
  SeqLib::BamReader bwalker;
  if (!bwalker.Open(o.bam)) {
    std::cerr << "ERROR: cannot open BAM " << o.bam << "\n";
    std::exit(EXIT_FAILURE);
  }
  const SeqLib::BamHeader hdr = bwalker.Header();

  // ---- 2. Read bps.txt.gz header and discover sample prefixes ----------
  //
  // Matches the logic in refilter.cpp:runRefilterBreakpoints. The sample
  // block starts at the first header token whose first char is 't' or 'n'
  // past the fixed-width prefix (scan from col 34 to be safe on all v1-v3
  // schemas). Once we know the prefixes, we stuff them into opts.bams so
  // BreakPoint's line-parser iterates sample tokens in the right order.
  igzstream infile(o.input_file.c_str(), std::ios::in);
  std::string headerLine;
  if (!std::getline(infile, headerLine)) {
    std::cerr << "ERROR: " << o.input_file << " is empty or unreadable\n";
    std::exit(EXIT_FAILURE);
  }
  const auto headerv = svabaUtils::tokenize_delimited(headerLine, '\t');
  if (headerv.size() < 41) {
    std::cerr << "ERROR: bps.txt header has " << headerv.size()
              << " columns; expected >=41 (legacy) or >=52 (v3).\n";
    std::exit(EXIT_FAILURE);
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
    std::cerr << "ERROR: could not locate sample columns in bps.txt header.\n";
    std::exit(EXIT_FAILURE);
  }

  // ---- 3. Build a minimal SvabaSharedConfig ------------------------------
  //
  // Same pattern refilter.cpp uses: logger/opts/writer on the stack (they
  // wrap a no-op output writer; we never init it or open any BAMs through
  // it). sc.header carries the BamHeader so BreakPoint parsing can turn
  // chrom names into IDs.
  SvabaLogger       logger;
  SvabaOptions      opts;
  SvabaOutputWriter writer(logger, opts);
  SvabaSharedConfig sc(logger, opts, writer);
  sc.header = hdr;

  // Populate opts.bams (std::map<prefix,path>) and the VCFHeader's
  // sample + colname state in lockstep, so the per-sample VCF column
  // order matches the sample block order in bps.txt.
  VCFHeader base_header;
  base_header.filedate  = svabaUtils::fileDateString();
  base_header.source    = "svaba tovcf " + o.input_file;
  base_header.reference = "";

  // contig lines: chrom + length from the BAM
  for (int i = 0; i < hdr.NumSequences(); ++i) {
    base_header.addContigField(hdr.IDtoName(i), hdr.GetSequenceLength(i));
  }

  for (size_t i = sample_start; i < headerv.size(); ++i) {
    auto [pref, path] = splitSampleHeader(headerv[i]);
    if (pref.empty() || (pref.at(0) != 't' && pref.at(0) != 'n')) {
      std::cerr << "ERROR: unexpected sample header token '" << headerv[i]
                << "' (col " << i << "); expected t***/n***.\n";
      std::exit(EXIT_FAILURE);
    }
    opts.bams[pref] = path;
    base_header.addSampleField(pref);
    base_header.colnames += "\t" + pref;
  }

  if (o.verbose > 0) {
    std::cerr << "...bps.txt has " << headerv.size() << " columns; "
              << sample_start << " core, "
              << (headerv.size() - sample_start) << " sample(s)\n";
  }

  // ---- 4. Construct VCFFile ---------------------------------------------
  //
  // `nopass=true` in the ctor tells VCFFile to RETAIN non-PASS rows while
  // parsing; whether they're emitted is a separate decision controlled by
  // `include_nonpass`. We always ingest everything from disk so the knob
  // surface below can choose later.
  //
  // `skip_dedup_in_ctor` is set from the --dedup flag (default on): the
  // ctor runs deduplicate() internally, and we need skip_dedup=true to
  // be in effect BEFORE that call fires. Passing it through the ctor
  // param gets that ordering right in one place.
  const bool skip_internal_dedup = !o.dedup;
  VCFFile vcf(o.input_file, o.analysis_id, sc, base_header,
              /*nopass=*/true,
              /*verbose=*/o.verbose > 0,
              /*skip_dedup_in_ctor=*/skip_internal_dedup);

  vcf.qual_mode       = o.qual_mode;
  vcf.sv_format       = o.always_bnd ? SvFormat::BND_ALWAYS
                                     : SvFormat::SYMBOLIC_WHEN_OBVIOUS;
  vcf.include_nonpass = o.include_nonpass;
  // vcf.skip_dedup is already set from the ctor param above.

  // ---- 5. Emit --------------------------------------------------------
  if (o.verbose > 0) {
    std::cerr << "...writing SV VCF -> "    << o.sv_out    << "\n"
              << "...writing indel VCF -> " << o.indel_out << "\n";
  }
  vcf.writeSvsSingleFile   (o.sv_out,    o.gzip, hdr);
  vcf.writeIndelsSingleFile(o.indel_out, o.gzip, hdr);

  if (o.verbose > 0) {
    std::cerr << "svaba tovcf: done.\n";
  }
}
