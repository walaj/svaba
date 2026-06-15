// SvabaExtractPairs.cpp — see SvabaExtractPairs.h for design rationale.
//
// Replaces scripts/extract_pairs_by_seq.sh. Two-pass BAM-native extract:
//
//   Pass 1: stream BAM, scan the 4-bit-packed SEQ field of each record
//           through an Aho-Corasick automaton built over the union of
//           query sequences and (unless --no-rc) their reverse
//           complements. Collect QNAMEs that hit at least one pattern.
//   Pass 2: re-stream BAM, write every record whose QNAME is in the set.
//
// The output preserves input record order, so coord-sorted in →
// coord-sorted out. Sort is skipped iff the input header declares
// @HD SO:coordinate; otherwise we shell out to `samtools sort`, matching
// the fallback used in SvabaPostprocess.cpp. Output is indexed via
// sam_index_build directly (no samtools index fork).

#include "SvabaExtractPairs.h"

#include <fcntl.h>
#include <getopt.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "htslib/sam.h"
#include "zlib.h"

#include "SeqLib/BamHeader.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/SeqLibUtils.h"  // AddCommas

#include "SvabaOptions.h"  // SVABA_VERSION

namespace {

// ---------- options ----------

struct Opts {
  std::string in_bam;
  std::string out_bam;
  std::vector<std::string> seqs;
  std::string seq_file;
  std::string counts_file;   // optional; only meaningful when seq_file is a bps
  int  threads  = 4;
  int  verbose  = 1;
  bool include_rc    = true;
  bool include_pairs = true;  // false → single-pass; emit only matched records
  bool whole_bam     = false; // true → scan entire BAM (skip breakpoint windows)
  int  pad           = 1000;  // breakend window half-width (bp) for region mode
};

enum { OPT_NO_RC = 1000, OPT_NO_PAIRS, OPT_COUNTS, OPT_WHOLE_BAM, OPT_PAD };

const char* kShortOpts = "hi:o:s:f:t:v:";
const struct option kLongOpts[] = {
  { "help",     no_argument,       nullptr, 'h' },
  { "input",    required_argument, nullptr, 'i' },
  { "output",   required_argument, nullptr, 'o' },
  { "seq",      required_argument, nullptr, 's' },
  { "seq-file", required_argument, nullptr, 'f' },
  { "threads",  required_argument, nullptr, 't' },
  { "verbose",  required_argument, nullptr, 'v' },
  { "no-rc",     no_argument,       nullptr, OPT_NO_RC },
  { "no-pairs",  no_argument,       nullptr, OPT_NO_PAIRS },
  { "counts",    required_argument, nullptr, OPT_COUNTS  },
  { "whole-bam", no_argument,       nullptr, OPT_WHOLE_BAM },
  { "pad",       required_argument, nullptr, OPT_PAD },
  { nullptr, 0, nullptr, 0 }
};

void printUsage() {
  std::cerr <<
    "Usage: svaba extract-pairs -i IN.bam -o OUT.bam (-s SEQ ... | -f FILE)\n"
    "                           [-t THREADS] [--no-rc] [--no-pairs]\n"
    "                           [--whole-bam] [--pad N] [-v V]\n"
    "\n"
    "  Extract all read pairs from IN.bam where either mate's SEQ contains\n"
    "  any of the given query sequences or their reverse complements.\n"
    "  Output is sorted (or input-order-preserved if input is sorted) and\n"
    "  indexed (.bai). Replaces scripts/extract_pairs_by_seq.sh.\n"
    "\n"
    "  Region targeting (default when -f is a bps.txt[.gz]): the kmer scan\n"
    "  (pass 1) is restricted to reads within +/-PAD bp of any breakend in\n"
    "  the bps file, using the BAM index. This avoids matching a kmer at a\n"
    "  distant, unrelated reference site that just happens to share the\n"
    "  20-mer. ALL kmers are searched within that merged window set, so\n"
    "  cross-breakpoint connections are still found. Pass --whole-bam to\n"
    "  scan every read (the old behavior); required if you used -s or a\n"
    "  plain seq-list (no breakend coordinates are available then).\n"
    "\n"
    "  -i, --input  <bam>   Input BAM (required).\n"
    "  -o, --output <bam>   Output BAM (required).\n"
    "  -s, --seq    <SEQ>   Query sequence; repeatable.\n"
    "  -f, --seq-file <p>   File of query sequences. Two formats accepted,\n"
    "                       auto-detected from the first non-empty line:\n"
    "                         * one ACGTN sequence per line ('#' and blank\n"
    "                           lines ignored), or\n"
    "                         * a svaba bps.txt[.gz] dump — every row's\n"
    "                           jxn_kmer column (col 53) is used as a query;\n"
    "                           rows with kmer == \".\" are skipped.\n"
    "                       Plain or gzip is fine for either format.\n"
    "  -t, --threads <n>    BGZF reader+writer threads. [4]\n"
    "      --whole-bam      Scan the entire BAM in pass 1 instead of just\n"
    "                       the breakpoint windows. Old default behavior;\n"
    "                       also the only option when queries come from -s\n"
    "                       or a plain seq-list (no coordinates).\n"
    "      --pad <N>        Half-width (bp) of the window placed around each\n"
    "                       breakend in region-targeted mode. [1000]\n"
    "      --no-rc          Skip reverse-complement augmentation.\n"
    "      --no-pairs       Emit only records whose own SEQ matched. Skips\n"
    "                       the pair-mate / supplementary pickup, runs in a\n"
    "                       single BAM pass (~2x faster). Useful when you\n"
    "                       just want to inspect which reads contain a motif.\n"
    "      --counts <FILE>  Emit a per-bp_id table of reads carrying each\n"
    "                       bp's junction kmer. Four columns with a header:\n"
    "                       bp_id<TAB>jxn_kmer<TAB>n_total_hits<TAB>n_unique_reads,\n"
    "                       sorted by n_total_hits DESC (tiebreak bp_id),\n"
    "                       so the worst over-matchers are at the top.\n"
    "                       n_total_hits = every matching alignment (incl.\n"
    "                       dup/secondary/supplementary); n_unique_reads =\n"
    "                       primary non-dup reads dedup'd by (bp_id, qname,\n"
    "                       mate). The jxn_kmer column makes low-complexity\n"
    "                       kmers (poly-A/T, simple repeats) obvious. Only\n"
    "                       valid when -f is a bps.txt[.gz] file (the source\n"
    "                       of the bp_id<->kmer mapping); rejected otherwise.\n"
    "  -v, --verbose <0-3>  Verbosity. [1]\n"
    "  -h, --help           This message.\n"
    "\n"
    "  At least one query sequence must be supplied (via -s or -f).\n"
    "  Match alphabet is ACGTN (case-insensitive); other IUPAC bases are\n"
    "  treated as mismatches.\n";
}

Opts parseOpts(int argc, char** argv) {
  Opts o;
  if (argc <= 1) { printUsage(); std::exit(EXIT_FAILURE); }

  for (int c; (c = getopt_long(argc, argv, kShortOpts, kLongOpts, nullptr)) != -1;) {
    std::istringstream arg(optarg ? optarg : "");
    switch (c) {
      case 'h': printUsage(); std::exit(EXIT_SUCCESS);
      case 'i': arg >> o.in_bam;   break;
      case 'o': arg >> o.out_bam;  break;
      case 's': o.seqs.emplace_back(optarg ? optarg : ""); break;
      case 'f': arg >> o.seq_file; break;
      case 't': arg >> o.threads;  break;
      case 'v': arg >> o.verbose;  break;
      case OPT_NO_RC:     o.include_rc    = false; break;
      case OPT_NO_PAIRS:  o.include_pairs = false; break;
      case OPT_COUNTS:    arg >> o.counts_file;    break;
      case OPT_WHOLE_BAM: o.whole_bam     = true;  break;
      case OPT_PAD:       arg >> o.pad;            break;
      default: printUsage(); std::exit(EXIT_FAILURE);
    }
  }

  if (o.in_bam.empty() || o.out_bam.empty()) {
    std::cerr << "ERROR: --input and --output are required\n";
    printUsage();
    std::exit(EXIT_FAILURE);
  }
  if (o.seqs.empty() && o.seq_file.empty()) {
    std::cerr << "ERROR: at least one of --seq or --seq-file must be given\n";
    printUsage();
    std::exit(EXIT_FAILURE);
  }
  if (o.threads < 1) o.threads = 1;
  if (o.pad < 0)     o.pad = 0;
  return o;
}

// ---------- small utilities ----------

bool fileExists(const std::string& p) {
  struct stat st{};
  return ::stat(p.c_str(), &st) == 0 && S_ISREG(st.st_mode);
}

// Wrap a std::string in single quotes for /bin/sh, escaping embedded quotes.
std::string shQuote(const std::string& s) {
  std::string out;
  out.reserve(s.size() + 2);
  out.push_back('\'');
  for (char c : s) {
    if (c == '\'') out += "'\\''";
    else out.push_back(c);
  }
  out.push_back('\'');
  return out;
}

// Uppercase + drop whitespace; sequences from CLI/files come in any case
// and may have stray whitespace from copy-paste.
std::string normalizeSeq(std::string_view s) {
  std::string out;
  out.reserve(s.size());
  for (char c : s) {
    if (std::isspace(static_cast<unsigned char>(c))) continue;
    out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
  }
  return out;
}

// Reverse-complement an ACGTN string (any other char passes through as-is —
// validation happens at AC insertion time).
std::string reverseComplement(std::string_view s) {
  std::string out;
  out.resize(s.size());
  for (std::size_t i = 0; i < s.size(); ++i) {
    char c = s[s.size() - 1 - i];
    switch (c) {
      case 'A': c = 'T'; break;
      case 'C': c = 'G'; break;
      case 'G': c = 'C'; break;
      case 'T': c = 'A'; break;
      case 'N': c = 'N'; break;
      default: break;  // leave invalid char so AC::add() rejects it cleanly
    }
    out[i] = c;
  }
  return out;
}

// Slurp a (gzip-aware) text file into a vector of lines. zlib's gzopen
// transparently handles plain and gzip inputs by sniffing the magic bytes,
// so we can hand it either .txt or .txt.gz without branching at the call
// site. Returns one vector entry per line, newline stripped.
std::vector<std::string> slurpLines(const std::string& path) {
  gzFile gz = ::gzopen(path.c_str(), "rb");
  if (!gz) throw std::runtime_error("cannot open file: " + path);
  std::vector<std::string> out;
  // Per-call line buffer big enough for any realistic bps.txt row
  // (which run ~kilobyte-scale once per-sample blocks are appended).
  // gzgets reads one line up to size-1 chars + null terminator; if a
  // line overflows, gzgets returns a partial fragment and we'd miss
  // the rest. 64K is comfortable.
  constexpr size_t LINE_BUF = 65536;
  std::vector<char> buf(LINE_BUF);
  while (::gzgets(gz, buf.data(), static_cast<int>(buf.size())) != nullptr) {
    std::string line(buf.data());
    while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
      line.pop_back();
    out.push_back(std::move(line));
  }
  ::gzclose(gz);
  return out;
}

// Decide whether a file is a bps.txt[.gz] (svaba breakpoint dump) or a
// plain one-seq-per-line list. Detection is content-based: bps files
// always start with the BreakPoint::header() comment line, which begins
// with `#chr1\t`. Anything else is treated as one-seq-per-line.
//
// (We could detect on extension too, but content-sniffing is robust to
// .txt.gz files that have been renamed and to user habits like `head -1`
// piping an uncompressed bps file through some intermediate tool.)
bool looksLikeBpsFile(const std::vector<std::string>& lines) {
  for (const auto& l : lines) {
    if (l.empty()) continue;
    return l.rfind("#chr1\t", 0) == 0;
  }
  return false;
}

// What `readSeqsFromFile` returns. For bps inputs we also carry the
// kmer→bp_id mapping so the counts pass can attribute hits back to the
// breakpoint they came from. For plain inputs `bp_ids_by_kmer` is empty
// and `from_bps` is false.
//
// Note `bp_ids_by_kmer` is keyed by the FORWARD-strand kmer string (as
// it appears in col 53 of bps.txt). The caller is responsible for
// looking up reverse-complement hits under the forward key — equivalent,
// since fwd and rc both descend from the same bp_id.
struct LoadedSeqs {
  std::vector<std::string> kmers;   // raw, may contain duplicates; caller de-dups
  std::unordered_map<std::string, std::vector<std::string>> bp_ids_by_kmer;
  // Both breakends (chr name, 1-based pos) of every row that contributed a
  // kmer. Drives region-targeted scanning. Empty for plain seq-list inputs.
  std::vector<std::pair<std::string, int>> breakends;
  bool from_bps = false;
};

// Locate a named column in a bps header line. Header lines begin with '#'.
int findBpsCol(const std::string& header_line, const std::string& name) {
  std::string h = header_line;
  if (!h.empty() && h.front() == '#') h.erase(0, 1);
  std::istringstream iss(h);
  std::string col;
  int idx = 0;
  while (std::getline(iss, col, '\t')) {
    if (col == name) return idx;
    ++idx;
  }
  return -1;
}

// Read junction-kmer queries (and their owning bp_ids) out of a
// bps.txt[.gz] dump. Skips rows where the kmer is "." (no precise
// junction — DSCRD-only, etc.) or empty. Caller is responsible for
// de-duplicating the returned kmer list; the bp_ids_by_kmer map already
// folds duplicates (same kmer string from multiple rows → all their
// bp_ids accumulated under one key).
LoadedSeqs readJxnKmersFromBps(const std::vector<std::string>& lines) {
  LoadedSeqs out;
  out.from_bps = true;
  if (lines.empty()) return out;

  // First non-empty line is the header — find the columns we need. Besides
  // jxn_kmer + bp_id we also pull both breakend coordinates (chr1/pos1 and
  // chr2/pos2) so the caller can target the BAM scan to breakpoint windows.
  int kmer_col = -1;
  int bp_col   = -1;
  int chr1_col = -1, pos1_col = -1, chr2_col = -1, pos2_col = -1;
  std::string header;
  for (const auto& l : lines) {
    if (l.empty()) continue;
    header   = l;
    kmer_col = findBpsCol(l, "jxn_kmer");
    bp_col   = findBpsCol(l, "bp_id");
    chr1_col = findBpsCol(l, "chr1");
    pos1_col = findBpsCol(l, "pos1");
    chr2_col = findBpsCol(l, "chr2");
    pos2_col = findBpsCol(l, "pos2");
    break;
  }
  if (kmer_col < 0) {
    throw std::runtime_error(
      "bps file has no `jxn_kmer` column — too old to use as a query "
      "source for extract-pairs (need v4 or later; re-run `svaba run` "
      "to regenerate)");
  }
  // bp_col may be -1 on a very early v4 dump that lacked bp_id; that's
  // fine — kmers still work, only the counts attribution will be empty.
  // chr/pos cols are standard (cols 0/1/3/4); if any are missing we simply
  // emit no breakends and the caller falls back to a whole-BAM scan.
  const bool have_coords =
      chr1_col >= 0 && pos1_col >= 0 && chr2_col >= 0 && pos2_col >= 0;

  // Walk data rows, splitting on tabs to grab all needed columns in one pass.
  int last_col = std::max(kmer_col, bp_col);
  if (have_coords)
    last_col = std::max(last_col,
                        std::max(std::max(chr1_col, pos1_col),
                                 std::max(chr2_col, pos2_col)));
  for (const auto& l : lines) {
    if (l.empty() || l[0] == '#') continue;
    std::string kmer_norm;
    std::string bp_id_val;
    std::string chr1_v, pos1_v, chr2_v, pos2_v;
    int col = 0;
    std::size_t s = 0;
    while (s <= l.size()) {
      const std::size_t e = l.find('\t', s);
      const std::size_t end = (e == std::string::npos) ? l.size() : e;
      std::string_view tok(&l[s], end - s);
      if (col == kmer_col) {
        if (!tok.empty() && tok != ".") {
          kmer_norm = normalizeSeq(tok);
        }
      } else if (col == bp_col) {
        if (!tok.empty() && tok != ".") bp_id_val.assign(tok.data(), tok.size());
      } else if (have_coords && col == chr1_col) {
        chr1_v.assign(tok.data(), tok.size());
      } else if (have_coords && col == pos1_col) {
        pos1_v.assign(tok.data(), tok.size());
      } else if (have_coords && col == chr2_col) {
        chr2_v.assign(tok.data(), tok.size());
      } else if (have_coords && col == pos2_col) {
        pos2_v.assign(tok.data(), tok.size());
      }
      if (e == std::string::npos || col >= last_col) break;
      s = e + 1;
      ++col;
    }
    if (kmer_norm.empty()) continue;
    out.kmers.push_back(kmer_norm);
    if (!bp_id_val.empty())
      out.bp_ids_by_kmer[kmer_norm].push_back(std::move(bp_id_val));

    // Record both breakends of this (kmer-bearing) row for region targeting.
    if (have_coords) {
      auto add_be = [&](const std::string& chr, const std::string& pos) {
        if (chr.empty() || chr == "." || pos.empty() || pos == ".") return;
        errno = 0;
        char* endp = nullptr;
        const long p = std::strtol(pos.c_str(), &endp, 10);
        if (endp == pos.c_str() || errno != 0 || p <= 0) return;
        out.breakends.emplace_back(chr, static_cast<int>(p));
      };
      add_be(chr1_v, pos1_v);
      add_be(chr2_v, pos2_v);
    }
  }
  // Dedup the bp_id list per kmer — same bp_id might legitimately appear
  // twice across rows (shouldn't, but defensive). Cheap and bounded.
  for (auto& kv : out.bp_ids_by_kmer) {
    auto& v = kv.second;
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
  }
  return out;
}

// Public seq-file reader. Auto-detects bps.txt[.gz] vs plain one-per-line
// list by sniffing the first non-empty line. gzip is handled transparently
// for both formats.
LoadedSeqs readSeqsFromFile(const std::string& path) {
  const auto lines = slurpLines(path);
  if (looksLikeBpsFile(lines))
    return readJxnKmersFromBps(lines);

  // Plain seq-per-line file (the original behavior).
  LoadedSeqs out;
  out.from_bps = false;
  for (const auto& l : lines) {
    if (l.empty() || l[0] == '#') continue;
    std::string s = normalizeSeq(l);
    if (!s.empty()) out.kmers.push_back(std::move(s));
  }
  return out;
}

// Header inspection: does the input BAM declare itself coord-sorted?
// Mirrors SvabaPostprocess.cpp::isCoordinateSorted (kept local rather than
// pulling in the postprocess header — it's a 10-line helper).
bool isCoordinateSorted(const SeqLib::BamHeader& h) {
  const std::string text = h.AsString();
  std::istringstream iss(text);
  std::string line;
  while (std::getline(iss, line)) {
    if (line.rfind("@HD\t", 0) == 0)
      return line.find("SO:coordinate") != std::string::npos;
    if (!line.empty() && line[0] == '@') break;  // @HD must be first
  }
  return false;
}

// Build the merged set of breakpoint windows to restrict pass-1 scanning to.
// Each breakend (chr name, 1-based pos) becomes a [pos-pad, pos+pad] window
// (0-based, clipped to the contig). Names absent from the BAM header are
// skipped (counted in n_skipped). Overlapping/touching windows are merged so
// the sequential region iterator in BamReader::SetRegions doesn't return the
// same read from two windows (a residual cross-window duplicate is still
// possible for a read spanning a sub-read-length gap; the scan path guards
// against that separately).
SeqLib::GRC buildBreakpointRegions(
    const std::vector<std::pair<std::string, int>>& breakends,
    const SeqLib::BamHeader& hdr, int pad, std::size_t& n_skipped) {
  SeqLib::GRC grc;
  n_skipped = 0;
  for (const auto& be : breakends) {
    const int tid = hdr.Name2ID(be.first);
    if (tid < 0) { ++n_skipped; continue; }
    const int seqlen = hdr.GetSequenceLength(tid);   // -1 if unknown
    const int pos0   = be.second - 1;                // bps pos is 1-based
    int start = pos0 - pad; if (start < 0) start = 0;
    int end   = pos0 + pad;
    if (seqlen > 0 && end > seqlen - 1) end = seqlen - 1;
    if (end < start) end = start;
    grc.add(SeqLib::GenomicRegion(tid, start, end));
  }
  // Guard: MergeOverlappingIntervals walks a std::list and increments past
  // begin() unconditionally, which is UB on an empty list. Only merge when
  // we actually added windows (all-skipped → empty GRC → whole-BAM fallback).
  if (grc.size())
    grc.MergeOverlappingIntervals();  // sorts + collapses overlapping/touching
  return grc;
}

// ---------- @PG stamp ----------
//
// Append a single @PG line tagged "svaba_extract_pairs" onto the input
// header so derived BAMs trace back to the exact extract command. We do
// NOT bother with the full postprocess-style PgChainInfo / uniquify
// logic: this tool produces a fresh output BAM, so collisions with a
// pre-existing svaba_extract_pairs ID would only arise if the user piped
// the same extract output through extract again — at which point a
// duplicate ID is the correct, informative trace.

std::string sanitizeHeaderValue(const std::string& s) {
  std::string out;
  out.reserve(s.size());
  for (char c : s) {
    if      (c == '\t')                 out.push_back(' ');
    else if (c == '\n' || c == '\r')    continue;
    else                                 out.push_back(c);
  }
  return out;
}

std::string appendExtractPairsPg(const std::string& hdr_text,
                                 const std::string& cl) {
  std::string line = std::string("@PG\tID:svaba_extract_pairs\tPN:svaba\tVN:") +
                     SVABA_VERSION + "\tCL:" + sanitizeHeaderValue(cl) + "\n";
  std::string out = hdr_text;
  if (!out.empty() && out.back() != '\n') out.push_back('\n');
  out += line;
  return out;
}

// Stitch the original CLI back together for the @PG CL: field.
std::string buildCommandLine(int argc, char** argv) {
  std::string s;
  for (int i = 0; i < argc; ++i) {
    if (i) s.push_back(' ');
    s += argv[i];
  }
  return s;
}

// ---------- Aho-Corasick over {A,C,G,T,N} ----------
//
// 5-letter alphabet keeps the goto table tight (5 * sizeof(int) per node).
// Patterns containing any non-ACGTN base are silently rejected at add() —
// in practice this only kicks in if the user passes a seq with IUPAC
// ambiguity codes like R/Y/M, which would never match an htslib-stored
// nibble anyway (regex-engine-wise R != [AG], it's a literal R).
//
// We precompute a goto-function table during BFS so the search inner loop
// is one indexed read per character with no failure-link walking. The
// tradeoff (vs. the textbook trie+fail-link variant) is a 5x larger goto
// table; for 100 patterns of length 30 that's 1500 nodes * 5 ints ≈ 30 KB,
// trivially L1-resident.
//
// Each node carries:
//   - pattern_id (-1 unless this node is the exact endpoint of a
//     user-added pattern). Multiple user patterns mapping to the same
//     trie endpoint (i.e. duplicate strings) keep the FIRST id assigned.
//   - has_match: propagated boolean along the failure chain so the fast
//     `searchNibbles()` bool API can short-circuit on first match without
//     walking output links.
//   - output_link: nearest ancestor along the failure chain whose
//     pattern_id is set. Used by `searchNibblesCollect()` to enumerate
//     every pattern_id that ended at the current input position,
//     including shorter suffix patterns that the longest match subsumes.

class AhoCorasick {
 public:
  AhoCorasick() : nodes_(1) {}  // node 0 = root

  // Add one pattern with an explicit pattern_id (>= 0). Pattern is
  // normalized (upper-case, whitespace-free) by the caller. Returns true
  // if added; false if rejected for containing non-ACGTN characters or
  // being empty. Calling `add` multiple times with the same string is
  // fine — it just leaves the first pattern_id in place (e.g. when fwd
  // and rc of the same kmer happen to be identical, like a palindrome).
  bool add(std::string_view pat, int pattern_id) {
    if (pat.empty()) return false;
    if (pattern_id < 0) return false;
    int cur = 0;
    for (char c : pat) {
      int idx = baseIdx(c);
      if (idx < 0) return false;  // reject pattern entirely
      // IMPORTANT: read by value, not by reference. A reference into
      // nodes_[cur].child[idx] would dangle across emplace_back if the
      // vector reallocates — UB that has occasionally produced spurious
      // matches in the wild. Resolve the index, then look up again post-
      // realloc to assign.
      int next = nodes_[cur].child[idx];
      if (next < 0) {
        next = static_cast<int>(nodes_.size());
        nodes_.emplace_back();              // may reallocate the storage
        nodes_[cur].child[idx] = next;      // safe: index lookup post-realloc
      }
      cur = next;
    }
    nodes_[cur].has_match = true;
    if (nodes_[cur].pattern_id < 0) nodes_[cur].pattern_id = pattern_id;
    ++n_added_;
    return true;
  }

  // Build failure links and the precomputed goto table. Must be called
  // exactly once after all patterns are added and before any search.
  void build() {
    // Flatten goto[u][i]: existing trie children stay; missing children at
    // root point back to root; missing children at deeper nodes are
    // resolved during BFS using the failure link of the parent.
    std::queue<int> q;
    for (int i = 0; i < kAlphabet; ++i) {
      int c = nodes_[0].child[i];
      if (c >= 0) {
        nodes_[c].fail = 0;
        q.push(c);
      } else {
        nodes_[0].child[i] = 0;  // root self-loop on misses
      }
    }
    while (!q.empty()) {
      const int u = q.front(); q.pop();
      const int uf = nodes_[u].fail;
      // output_link: first true endpoint reachable along the failure
      // chain. Used only by searchNibblesCollect to enumerate all
      // suffix-pattern matches at the current input position.
      nodes_[u].output_link = (nodes_[uf].pattern_id >= 0)
                            ? uf
                            : nodes_[uf].output_link;
      // has_match propagates the bool "some pattern ends here or via
      // failure chain" — keeps searchNibbles a single branch per char.
      if (nodes_[uf].has_match) nodes_[u].has_match = true;
      for (int i = 0; i < kAlphabet; ++i) {
        const int v = nodes_[u].child[i];
        if (v >= 0) {
          nodes_[v].fail = nodes_[uf].child[i];
          q.push(v);
        } else {
          // Missing child at u: redirect to where uf's child for i lands.
          nodes_[u].child[i] = nodes_[uf].child[i];
        }
      }
    }
  }

  // Scan a 4-bit-packed BAM SEQ (htslib bam_get_seq layout) of length qlen
  // for any pattern hit. Returns true on first match.
  //
  // We translate htslib's nibble encoding (A=1,C=2,G=4,T=8,N=15, others
  // = ambiguity codes) directly into our 0..4 alphabet via a 16-entry
  // lookup table. Ambiguity nibbles map to -1, which resets the AC state
  // to root (no match can span them).
  bool searchNibbles(const uint8_t* seq4, int qlen) const {
    int cur = 0;
    for (int i = 0; i < qlen; ++i) {
      const uint8_t nib = bam_seqi(seq4, i);
      const int idx = nib16_to_idx_[nib];
      if (idx < 0) { cur = 0; continue; }
      cur = nodes_[cur].child[idx];
      if (nodes_[cur].has_match) return true;
    }
    return false;
  }

  // Scan as above, but collect EVERY pattern_id that matched anywhere in
  // the input. Result is sorted+deduplicated. Slower than `searchNibbles`
  // because we don't early-exit (and we walk the output_link chain on
  // matches), but the per-char cost on non-matching inputs is identical.
  // Used by the per-bp_id counting path to attribute each hit back to
  // the right user-supplied kmer.
  void searchNibblesCollect(const uint8_t* seq4, int qlen,
                            std::vector<int>& out) const {
    out.clear();
    int cur = 0;
    for (int i = 0; i < qlen; ++i) {
      const uint8_t nib = bam_seqi(seq4, i);
      const int idx = nib16_to_idx_[nib];
      if (idx < 0) { cur = 0; continue; }
      cur = nodes_[cur].child[idx];
      if (nodes_[cur].pattern_id >= 0) out.push_back(nodes_[cur].pattern_id);
      for (int v = nodes_[cur].output_link; v >= 0;
           v = nodes_[v].output_link) {
        out.push_back(nodes_[v].pattern_id);
      }
    }
    if (out.size() > 1) {
      std::sort(out.begin(), out.end());
      out.erase(std::unique(out.begin(), out.end()), out.end());
    }
  }

  std::size_t pattern_count() const { return n_added_; }
  std::size_t node_count()    const { return nodes_.size(); }

 private:
  static constexpr int kAlphabet = 5;  // A,C,G,T,N

  struct Node {
    std::array<int, kAlphabet> child;
    int  fail        = 0;
    int  output_link = -1;   // next true-endpoint ancestor along failure chain
    int  pattern_id  = -1;   // -1 = not a user-pattern endpoint
    bool has_match   = false;
    Node() { child.fill(-1); }
  };

  static int baseIdx(char c) {
    switch (c) {
      case 'A': return 0;
      case 'C': return 1;
      case 'G': return 2;
      case 'T': return 3;
      case 'N': return 4;
      default:  return -1;
    }
  }

  // bam_seqi nibble -> AC alphabet index. -1 means "any match-breaking
  // base" (=, M, R, S, V, W, Y, H, K, D, B). Indices: A=1, C=2, G=4,
  // T=8, N=15 (per htslib's seq_nt16_str ordering).
  static constexpr std::array<int, 16> nib16_to_idx_ = {
    -1,  // 0  '='
     0,  // 1  'A'
     1,  // 2  'C'
    -1,  // 3  'M'
     2,  // 4  'G'
    -1,  // 5  'R'
    -1,  // 6  'S'
    -1,  // 7  'V'
     3,  // 8  'T'
    -1,  // 9  'W'
    -1,  // 10 'Y'
    -1,  // 11 'H'
    -1,  // 12 'K'
    -1,  // 13 'D'
    -1,  // 14 'B'
     4,  // 15 'N'
  };

  std::vector<Node> nodes_;
  std::size_t       n_added_ = 0;
};

// ---------- pass 1: collect QNAMEs ----------
//
// Stream the input BAM and run AC on each record's SEQ. We use SeqLib's
// pull-style Next() and access bam_get_seq through BamRecord::raw() so the
// inner loop never allocates a std::string for the decoded sequence. The
// QNAME hash set holds one std::string per matching read pair (plus
// supplementary/secondary aliases under the same QNAME).
//
// In counts mode (caller passes non-null pattern_to_bp_ids), we additionally
// switch the AC scan to `searchNibblesCollect` to enumerate every pattern_id
// that matched in this read's SEQ, and attribute each hit to the bp_id(s)
// that own the kmer. Counts are restricted to primary, non-duplicate,
// non-secondary, non-supplementary alignments, dedup'd by (bp_id, qname,
// mate). QNAME collection itself does NOT apply the flag filter — we still
// want every read whose SEQ matched (including duplicates / supplementaries
// of a real read pair) to flow into the pair-pickup pass.

struct Pass1Stats {
  std::size_t reads        = 0;
  std::size_t matched_reads = 0;
  std::size_t unique_qnames = 0;
  std::size_t counted_reads = 0;   // primary non-dup matches actually tallied
  double      seconds      = 0;
};

Pass1Stats collectQnames(const std::string& in_bam,
                         const AhoCorasick& ac,
                         int threads,
                         int verbose,
                         std::unordered_set<std::string>& out_qnames,
                         const std::vector<std::vector<std::string>>*
                             pattern_to_bp_ids,
                         std::unordered_map<std::string, std::size_t>*
                             bp_id_counts,
                         std::unordered_map<std::string, std::size_t>*
                             bp_id_total,
                         const SeqLib::GRC* regions) {
  using clock = std::chrono::steady_clock;
  Pass1Stats st;
  const auto t0 = clock::now();

  SeqLib::BamReader r;
  if (!r.Open(in_bam))
    throw std::runtime_error("pass1: cannot open " + in_bam);
  r.SetThreads(std::max(1, threads));

  // Region targeting: restrict the iterator to the breakpoint windows. If the
  // BAM has no index, SetRegions fails — fall back to a whole-BAM scan rather
  // than aborting, so the command still produces a (less targeted) result.
  bool region_dedup_active = false;
  if (regions && regions->size()) {
    if (r.SetRegions(*regions)) {
      region_dedup_active = true;
    } else {
      std::cerr << "[pass1] WARNING: could not restrict to breakpoint regions "
                   "(missing .bai?); scanning the whole BAM instead\n";
    }
  }
  // Collapses a read returned by two adjacent windows (see is_region_dup).
  std::unordered_set<std::string> region_seen;

  // Cached chr name for the progress line. ChrName() hits the header every
  // call, so we resolve it only on (rare) chromosome change, not per record.
  // Sentinel -2 won't match any valid tid (-1 = unmapped, 0..N-1 = mapped).
  int32_t     cur_chr      = -2;
  int32_t     cur_pos      = -1;
  std::string cur_chr_name = "*";

  const bool counts_mode = (pattern_to_bp_ids != nullptr) &&
                           (bp_id_counts != nullptr) &&
                           (bp_id_total != nullptr);
  // Dedup set for per-bp_id counts. Key is bp_id + '\0' + qname + '\0' + mate.
  // Lives for the whole pass — bounded by the number of (bp, read) pairs that
  // actually match a kmer, which is small even on WGS (1e4–1e6 entries).
  std::unordered_set<std::string> bp_seen;
  std::vector<int> hit_ids;        // reused per matching record

  while (auto opt = r.Next()) {
    SeqLib::BamRecord& rec = *opt;
    ++st.reads;

    const bam1_t* b = rec.raw();
    if (!b) continue;

    // Update cached locus before any continue — we want progress to keep
    // ticking past records without SEQ (they're rare but they shouldn't
    // freeze the displayed position).
    const int32_t chr = b->core.tid;
    const int32_t pos = b->core.pos;
    if (chr != cur_chr) {
      cur_chr_name = (chr >= 0) ? rec.ChrName(r.Header()) : std::string("*");
      cur_chr      = chr;
    }
    cur_pos = pos;

    const int qlen = b->core.l_qseq;
    if (qlen <= 0) continue;
    const uint8_t* seq4 = bam_get_seq(b);

    // In region mode, the sequential per-window iterator can hand back a read
    // that straddles the gap between two adjacent (post-merge) windows twice.
    // Collapse those by (tid, pos, qname, flag) so matches aren't double-
    // counted and (single-pass) output isn't duplicated. Only matched reads
    // ever enter the set, so it stays small.
    auto is_region_dup = [&]() -> bool {
      if (!region_dedup_active) return false;
      std::string key;
      key.reserve(64);
      key.append(std::to_string(b->core.tid)); key.push_back(':');
      key.append(std::to_string(b->core.pos)); key.push_back(':');
      key.append(rec.Qname());                 key.push_back(':');
      key.append(std::to_string(b->core.flag));
      return !region_seen.insert(std::move(key)).second;
    };

    if (counts_mode) {
      ac.searchNibblesCollect(seq4, qlen, hit_ids);
      if (!hit_ids.empty()) {
        if (is_region_dup()) continue;
        ++st.matched_reads;
        out_qnames.insert(rec.Qname());

        // RAW total: credit every matching alignment (incl. duplicate,
        // secondary, supplementary) to each owning bp_id. hit_ids is
        // already deduplicated per read, so a bp_id is counted at most
        // once per alignment. This is the magnitude that explains the
        // "matches" number in the progress log — the right signal for
        // spotting low-complexity kmers that match far more than expected.
        for (int pid : hit_ids) {
          if (pid < 0 || pid >= static_cast<int>(pattern_to_bp_ids->size()))
            continue;
          for (const auto& bp_id : (*pattern_to_bp_ids)[pid])
            ++(*bp_id_total)[bp_id];
        }

        // UNIQUE: only primary non-dup non-sec non-sup reads, deduplicated
        // by (bp_id, qname, mate) — the conservative "real read pairs"
        // count, unaffected by PCR duplication.
        const uint16_t flag = b->core.flag;
        const bool primary = (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY
                                       | BAM_FDUP)) == 0;
        if (primary) {
          const char* mate = (flag & BAM_FREAD1) ? "1"
                            : (flag & BAM_FREAD2) ? "2" : "0";
          // Build dedup key on the stack to avoid per-bp string churn.
          std::string qn = rec.Qname();
          std::string key;
          key.reserve(48);
          for (int pid : hit_ids) {
            if (pid < 0 ||
                pid >= static_cast<int>(pattern_to_bp_ids->size()))
              continue;
            const auto& bps = (*pattern_to_bp_ids)[pid];
            for (const auto& bp_id : bps) {
              key.clear();
              key.append(bp_id);
              key.push_back('\0');
              key.append(qn);
              key.push_back('\0');
              key.append(mate);
              if (bp_seen.insert(key).second) {
                ++(*bp_id_counts)[bp_id];
              }
            }
          }
          ++st.counted_reads;
        }
      }
    } else if (ac.searchNibbles(seq4, qlen)) {
      if (is_region_dup()) continue;
      ++st.matched_reads;
      out_qnames.insert(rec.Qname());
    }

    constexpr std::size_t PROGRESS_EVERY = 25'000'000ULL;
    if (verbose >= 1 && (st.reads % PROGRESS_EVERY) == 0) {
      const double elapsed =
          std::chrono::duration<double>(clock::now() - t0).count();
      // BAM positions are 0-based; humans expect 1-based. Unmapped (-1)
      // displayed as "*:0" via the chr_name sentinel and the >=0 guard.
      std::cerr << "[pass1] at " << cur_chr_name << ":"
                << SeqLib::AddCommas<int32_t>(cur_pos >= 0 ? cur_pos + 1 : 0)
                << " | "
                << SeqLib::AddCommas(st.reads) << " reads scanned, "
                << SeqLib::AddCommas(st.matched_reads) << " matches, "
                << SeqLib::AddCommas(out_qnames.size()) << " unique qnames, "
                << std::fixed << std::setprecision(1) << elapsed
                << "s elapsed\n";
    }
  }

  st.unique_qnames = out_qnames.size();
  st.seconds = std::chrono::duration<double>(clock::now() - t0).count();
  return st;
}

// ---------- pass 2: extract all alignments for matched QNAMEs ----------

struct Pass2Stats {
  std::size_t reads_in   = 0;
  std::size_t reads_out  = 0;
  double      seconds    = 0;
};

Pass2Stats extractByQname(const std::string& in_bam,
                          const std::string& out_bam,
                          const std::unordered_set<std::string>& qnames,
                          int threads,
                          int verbose,
                          const std::string& cl) {
  using clock = std::chrono::steady_clock;
  Pass2Stats st;
  const auto t0 = clock::now();

  SeqLib::BamReader r;
  if (!r.Open(in_bam))
    throw std::runtime_error("pass2: cannot open " + in_bam);
  r.SetThreads(std::max(1, threads));

  // Build the stamped header.
  const SeqLib::BamHeader src_hdr = r.Header();
  const SeqLib::BamHeader stamped(appendExtractPairsPg(src_hdr.AsString(), cl));

  SeqLib::BamWriter w;
  if (!w.Open(out_bam))
    throw std::runtime_error("pass2: cannot open output " + out_bam);
  w.SetThreads(std::max(1, threads));
  w.SetHeader(stamped);
  if (!w.WriteHeader())
    throw std::runtime_error("pass2: WriteHeader failed on " + out_bam);

  // Cached chr name for the progress line — same pattern as pass 1.
  int32_t     cur_chr      = -2;
  int32_t     cur_pos      = -1;
  std::string cur_chr_name = "*";

  while (auto opt = r.Next()) {
    SeqLib::BamRecord& rec = *opt;
    ++st.reads_in;

    if (const bam1_t* b = rec.raw()) {
      const int32_t chr = b->core.tid;
      if (chr != cur_chr) {
        cur_chr_name = (chr >= 0) ? rec.ChrName(r.Header()) : std::string("*");
        cur_chr      = chr;
      }
      cur_pos = b->core.pos;
    }

    // Look up by std::string view of bam_get_qname to avoid allocating
    // a new std::string when the qname isn't in the set. unordered_set's
    // transparent lookup is C++20-only, so we compare via temporary.
    // (Allocation cost on a miss is bounded by qname length; on dense
    //  BAMs with rare matches, the alternative — a custom heterogeneous
    //  hasher — is a marginal optimization we can revisit if profiling
    //  ever flags pass 2 as the bottleneck. It's currently not.)
    std::string qn = rec.Qname();
    if (qnames.find(qn) == qnames.end()) continue;

    if (!w.WriteRecord(rec))
      throw std::runtime_error("pass2: WriteRecord failed on " + out_bam);
    ++st.reads_out;

    constexpr std::size_t PROGRESS_EVERY = 25'000'000ULL;
    if (verbose >= 2 && (st.reads_in % PROGRESS_EVERY) == 0) {
      const double elapsed =
          std::chrono::duration<double>(clock::now() - t0).count();
      std::cerr << "[pass2] at " << cur_chr_name << ":"
                << SeqLib::AddCommas<int32_t>(cur_pos >= 0 ? cur_pos + 1 : 0)
                << " | "
                << SeqLib::AddCommas(st.reads_in) << " in, "
                << SeqLib::AddCommas(st.reads_out) << " written, "
                << std::fixed << std::setprecision(1) << elapsed
                << "s elapsed\n";
    }
  }

  if (!w.Close())
    throw std::runtime_error("pass2: writer Close failed on " + out_bam);

  st.seconds = std::chrono::duration<double>(clock::now() - t0).count();
  return st;
}

// ---------- single-pass mode: emit only matched records ----------
//
// Used when --no-pairs is passed. Skips the QNAME hash-set bookkeeping
// and the second BAM pass entirely: every record's SEQ is scanned, and
// records that hit are written directly to the output. The mate of a
// matched read is NOT pulled in, and supplementary/secondary alignments
// are emitted only if their own SEQ also matched (which usually it
// won't, since supplementaries often carry SEQ="*").
//
// Roughly 2x faster than the two-pass mode: one BGZF decompress pass
// instead of two, no hash set, no qname allocation per record.

struct SinglePassStats {
  std::size_t reads_in  = 0;
  std::size_t reads_out = 0;
  double      seconds   = 0;
};

SinglePassStats extractMatchedOnly(const std::string& in_bam,
                                   const std::string& out_bam,
                                   const AhoCorasick& ac,
                                   int threads,
                                   int verbose,
                                   const std::string& cl,
                                   const SeqLib::GRC* regions) {
  using clock = std::chrono::steady_clock;
  SinglePassStats st;
  const auto t0 = clock::now();

  SeqLib::BamReader r;
  if (!r.Open(in_bam))
    throw std::runtime_error("single-pass: cannot open " + in_bam);
  r.SetThreads(std::max(1, threads));

  // Region targeting (see collectQnames). Fall back to whole-BAM if the
  // index is missing so the command still works.
  bool region_dedup_active = false;
  if (regions && regions->size()) {
    if (r.SetRegions(*regions)) {
      region_dedup_active = true;
    } else {
      std::cerr << "[single-pass] WARNING: could not restrict to breakpoint "
                   "regions (missing .bai?); scanning the whole BAM instead\n";
    }
  }
  std::unordered_set<std::string> region_seen;

  // Stamped output header (same shape as the two-pass mode's pass-2 header).
  const SeqLib::BamHeader src_hdr = r.Header();
  const SeqLib::BamHeader stamped(appendExtractPairsPg(src_hdr.AsString(), cl));

  SeqLib::BamWriter w;
  if (!w.Open(out_bam))
    throw std::runtime_error("single-pass: cannot open output " + out_bam);
  w.SetThreads(std::max(1, threads));
  w.SetHeader(stamped);
  if (!w.WriteHeader())
    throw std::runtime_error("single-pass: WriteHeader failed on " + out_bam);

  // Cached chr name for the progress line — same pattern as pass 1.
  int32_t     cur_chr      = -2;
  int32_t     cur_pos      = -1;
  std::string cur_chr_name = "*";

  while (auto opt = r.Next()) {
    SeqLib::BamRecord& rec = *opt;
    ++st.reads_in;

    const bam1_t* b = rec.raw();
    if (!b) continue;

    const int32_t chr = b->core.tid;
    if (chr != cur_chr) {
      cur_chr_name = (chr >= 0) ? rec.ChrName(r.Header()) : std::string("*");
      cur_chr      = chr;
    }
    cur_pos = b->core.pos;

    const int qlen = b->core.l_qseq;
    if (qlen <= 0) continue;
    const uint8_t* seq4 = bam_get_seq(b);

    if (ac.searchNibbles(seq4, qlen)) {
      // Drop a read handed back by two adjacent windows (region mode only).
      if (region_dedup_active) {
        std::string key;
        key.reserve(64);
        key.append(std::to_string(b->core.tid)); key.push_back(':');
        key.append(std::to_string(b->core.pos)); key.push_back(':');
        key.append(rec.Qname());                 key.push_back(':');
        key.append(std::to_string(b->core.flag));
        if (!region_seen.insert(std::move(key)).second) continue;
      }
      if (!w.WriteRecord(rec))
        throw std::runtime_error("single-pass: WriteRecord failed on " + out_bam);
      ++st.reads_out;
    }

    constexpr std::size_t PROGRESS_EVERY = 25'000'000ULL;
    if (verbose >= 1 && (st.reads_in % PROGRESS_EVERY) == 0) {
      const double elapsed =
          std::chrono::duration<double>(clock::now() - t0).count();
      std::cerr << "[single-pass] at " << cur_chr_name << ":"
                << SeqLib::AddCommas<int32_t>(cur_pos >= 0 ? cur_pos + 1 : 0)
                << " | "
                << SeqLib::AddCommas(st.reads_in) << " reads scanned, "
                << SeqLib::AddCommas(st.reads_out) << " written, "
                << std::fixed << std::setprecision(1) << elapsed
                << "s elapsed\n";
    }
  }

  if (!w.Close())
    throw std::runtime_error("single-pass: writer Close failed on " + out_bam);

  st.seconds = std::chrono::duration<double>(clock::now() - t0).count();
  return st;
}

// ---------- post-pass2 sort fallback ----------
//
// Only invoked when the input header does NOT declare SO:coordinate. The
// common case is sorted input, where we skip this entirely and the output
// inherits the input's coord-sort.

int sortInPlace(const std::string& bam, int threads, int verbose) {
  const std::string tmp = bam + ".extract.sort.tmp.bam";
  std::string cmd = "samtools sort -@ " + std::to_string(std::max(1, threads)) +
                    " -o " + shQuote(tmp) + " " + shQuote(bam);
  if (verbose >= 2) std::cerr << "[sort] exec: " << cmd << "\n";
  const auto t0 = std::chrono::steady_clock::now();
  const int rc = std::system(cmd.c_str());
  const double elapsed =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  if (rc != 0) {
    ::unlink(tmp.c_str());
    std::cerr << "[sort] samtools sort FAILED (exit " << rc << ") after "
              << std::fixed << std::setprecision(1) << elapsed << "s\n";
    return rc;
  }
  if (std::rename(tmp.c_str(), bam.c_str()) != 0) {
    ::unlink(tmp.c_str());
    std::cerr << "[sort] rename failed: " << std::strerror(errno) << "\n";
    return -1;
  }
  if (verbose >= 1)
    std::cerr << "[sort] done in " << std::fixed << std::setprecision(1)
              << elapsed << "s\n";
  return 0;
}

}  // namespace

void runExtractPairs(int argc, char** argv) {
  const Opts o = parseOpts(argc, argv);

  if (!fileExists(o.in_bam)) {
    std::cerr << "ERROR: input BAM not found: " << o.in_bam << "\n";
    std::exit(EXIT_FAILURE);
  }

  // ---- gather + normalize patterns ----
  // bp_ids_by_kmer is populated only when a bps.txt was the input. It maps
  // forward-strand kmer string → list of bp_ids that emitted it. We carry
  // it through to the AC-pattern_id → bp_ids vector below.
  std::vector<std::string> seqs = o.seqs;
  for (auto& s : seqs) s = normalizeSeq(s);
  std::unordered_map<std::string, std::vector<std::string>> bp_ids_by_kmer;
  std::vector<std::pair<std::string, int>> breakends;  // (chr, 1-based pos)
  bool seq_file_is_bps = false;
  if (!o.seq_file.empty()) {
    LoadedSeqs loaded = readSeqsFromFile(o.seq_file);
    seq_file_is_bps   = loaded.from_bps;
    bp_ids_by_kmer    = std::move(loaded.bp_ids_by_kmer);
    breakends         = std::move(loaded.breakends);
    seqs.insert(seqs.end(),
                std::make_move_iterator(loaded.kmers.begin()),
                std::make_move_iterator(loaded.kmers.end()));
  }
  // De-dup pattern list.
  {
    std::sort(seqs.begin(), seqs.end());
    seqs.erase(std::unique(seqs.begin(), seqs.end()), seqs.end());
  }

  if (seqs.empty()) {
    std::cerr << "ERROR: no usable query sequences\n";
    std::exit(EXIT_FAILURE);
  }

  // ---- --counts validity check ----
  // Counts only make sense when the source of the kmer list is a bps file
  // (the only place we have the kmer↔bp_id mapping). Reject otherwise so
  // the user doesn't silently get an empty TSV. Also incompatible with
  // --no-pairs: the count tally lives on pass 1 of the two-pass path,
  // and the single-pass writer doesn't carry that state.
  const bool want_counts = !o.counts_file.empty();
  if (want_counts && !seq_file_is_bps) {
    std::cerr << "ERROR: --counts requires -f to be a svaba bps.txt[.gz] file; "
                 "got "
              << (o.seq_file.empty() ? "no seq-file (only -s)"
                                     : "a plain seq-per-line file")
              << "\n";
    std::exit(EXIT_FAILURE);
  }
  if (want_counts && !o.include_pairs) {
    std::cerr << "ERROR: --counts is incompatible with --no-pairs (counting "
                 "happens in the two-pass mode's pass 1)\n";
    std::exit(EXIT_FAILURE);
  }

  // ---- build Aho-Corasick over patterns + reverse complements ----
  // pattern_to_bp_ids[i] is the list of bp_ids that emitted seqs[i]
  // (empty when -s was used or when the bps file had no bp_id column).
  // Both the forward kmer and its reverse complement share the same
  // pattern_id (== index into `seqs`) because they refer to the same
  // breakpoint kmer from the user's POV.
  AhoCorasick ac;
  std::vector<std::vector<std::string>> pattern_to_bp_ids;
  pattern_to_bp_ids.reserve(seqs.size());
  std::size_t added = 0, rejected = 0;
  for (std::size_t i = 0; i < seqs.size(); ++i) {
    const auto& s   = seqs[i];
    const int   pid = static_cast<int>(i);
    if (ac.add(s, pid)) { ++added; } else { ++rejected; }
    if (o.include_rc) {
      const std::string rc = reverseComplement(s);
      // Add only if RC is a different string from the forward kmer
      // (palindromes naturally collapse). The pattern_id stays == i so
      // that hits via either direction attribute to the same bp_ids.
      if (rc != s) {
        if (ac.add(rc, pid)) { ++added; } else { ++rejected; }
      }
    }
    // Look the forward kmer up in bp_ids_by_kmer to build pattern_to_bp_ids.
    auto it = bp_ids_by_kmer.find(s);
    if (it != bp_ids_by_kmer.end())
      pattern_to_bp_ids.emplace_back(it->second);
    else
      pattern_to_bp_ids.emplace_back();
  }
  ac.build();

  if (o.verbose >= 1) {
    std::cerr << "[extract-pairs] " << seqs.size() << " query sequences ("
              << (o.include_rc ? "with" : "without")
              << " reverse complements); " << added << " added to AC, "
              << rejected << " rejected (non-ACGTN)\n";
    if (added == 0) {
      std::cerr << "ERROR: no patterns survived ACGTN normalization\n";
      std::exit(EXIT_FAILURE);
    }
    if (o.verbose >= 2) {
      std::cerr << "[extract-pairs] AC node count: " << ac.node_count() << "\n";
    }
  }

  const std::string cl = buildCommandLine(argc, argv);

  // ---- detect input sort order + build breakpoint regions ----
  // One header open serves both: the sort check (cheap) and resolving
  // breakend chrom names to tids for the region windows.
  bool input_sorted = false;
  SeqLib::GRC regions;       // breakpoint windows; empty => whole-BAM scan
  {
    SeqLib::BamReader hdr_only;
    if (hdr_only.Open(o.in_bam)) {
      input_sorted = isCoordinateSorted(hdr_only.Header());

      // Region targeting is the default whenever breakend coordinates are
      // available (i.e. -f was a bps.txt) and --whole-bam wasn't requested.
      if (!o.whole_bam && seq_file_is_bps && !breakends.empty()) {
        std::size_t skipped = 0;
        regions = buildBreakpointRegions(breakends, hdr_only.Header(),
                                         o.pad, skipped);
        if (o.verbose >= 1 && skipped)
          std::cerr << "[extract-pairs] " << SeqLib::AddCommas(skipped)
                    << " breakend(s) skipped (chrom not in BAM header)\n";
      }
    }
  }
  const bool use_regions = regions.size() > 0;
  const SeqLib::GRC* regions_ptr = use_regions ? &regions : nullptr;

  if (o.verbose >= 1) {
    std::cerr << "[extract-pairs] input is "
              << (input_sorted ? "coord-sorted (sort step will be skipped)"
                               : "not coord-sorted (samtools sort will run)")
              << "\n";
    if (use_regions)
      std::cerr << "[extract-pairs] region-targeted scan: "
                << SeqLib::AddCommas(regions.size())
                << " merged windows (+/-" << o.pad << " bp around "
                << SeqLib::AddCommas(breakends.size())
                << " breakends); pass 1 reads only these regions "
                   "(--whole-bam to scan everything)\n";
    else if (o.whole_bam)
      std::cerr << "[extract-pairs] whole-BAM scan (--whole-bam)\n";
    else if (!seq_file_is_bps)
      std::cerr << "[extract-pairs] whole-BAM scan (no breakend coordinates; "
                   "pass -f a bps.txt[.gz] to enable region targeting)\n";
  }

  std::size_t out_records = 0;
  std::size_t qname_count = 0;  // unused in single-pass mode

  if (!o.include_pairs) {
    // ---- single-pass mode ----
    if (o.verbose >= 1)
      std::cerr << "[extract-pairs] single-pass mode (--no-pairs): emitting "
                   "only records whose own SEQ matched\n";
    const SinglePassStats sp = extractMatchedOnly(
        o.in_bam, o.out_bam, ac, o.threads, o.verbose, cl, regions_ptr);
    if (o.verbose >= 1)
      std::cerr << "[extract-pairs] single-pass: "
                << SeqLib::AddCommas(sp.reads_in) << " reads scanned, "
                << SeqLib::AddCommas(sp.reads_out) << " written in "
                << std::fixed << std::setprecision(1) << sp.seconds << "s\n";
    if (sp.reads_out == 0) {
      // The output BAM exists but contains zero records. Fine for indexing,
      // but flag it.
      std::cerr << "[extract-pairs] WARNING: no records matched; output BAM is empty\n";
    }
    out_records = sp.reads_out;
  } else {
    // ---- two-pass mode (default) ----
    std::unordered_set<std::string> qnames;
    qnames.reserve(1 << 14);
    if (o.verbose >= 1)
      std::cerr << "[extract-pairs] pass 1/2: scanning " << o.in_bam
                << " for matching SEQ"
                << (want_counts ? " (counts mode: tracking per-bp_id support)"
                                : "")
                << "...\n";
    std::unordered_map<std::string, std::size_t> bp_id_counts;  // unique
    std::unordered_map<std::string, std::size_t> bp_id_total;   // raw hits
    const Pass1Stats p1 = collectQnames(
        o.in_bam, ac, o.threads, o.verbose, qnames,
        want_counts ? &pattern_to_bp_ids : nullptr,
        want_counts ? &bp_id_counts      : nullptr,
        want_counts ? &bp_id_total       : nullptr,
        regions_ptr);
    if (o.verbose >= 1) {
      std::cerr << "[extract-pairs] pass 1: "
                << SeqLib::AddCommas(p1.reads) << " reads, "
                << SeqLib::AddCommas(p1.matched_reads) << " matched, "
                << SeqLib::AddCommas(p1.unique_qnames) << " unique qnames";
      if (want_counts)
        std::cerr << ", " << SeqLib::AddCommas(p1.counted_reads)
                  << " primary non-dup reads counted toward bp_ids";
      std::cerr << " in " << std::fixed << std::setprecision(1)
                << p1.seconds << "s\n";
    }

    // Emit the counts TSV before any early exit on empty matches — the user
    // asked for a per-bp_id table; an empty result is still a result.
    //
    // Schema (4 cols): bp_id, jxn_kmer, n_total_hits, n_unique_reads.
    //   - jxn_kmer:       the forward-strand 20-mer (bps col 53) that drove
    //                     the matches. Low-complexity kmers (poly-A/T, simple
    //                     repeats) jump out here as the over-match culprits.
    //   - n_total_hits:   every matching alignment, incl. dup/secondary/supp.
    //   - n_unique_reads: primary non-dup reads, dedup'd by bp_id+qname+mate.
    // Rows are sorted by n_total_hits DESC (tiebreak bp_id ASC) so the worst
    // offenders are at the top — that's what you want when debugging.
    if (want_counts) {
      // Invert the pattern→bp_id map to bp_id→kmer (1:1 per bps row; if two
      // bp_ids share a kmer, each still resolves to that same kmer string).
      std::unordered_map<std::string, std::string> kmer_by_bp_id;
      for (std::size_t i = 0; i < seqs.size(); ++i)
        for (const auto& bp : pattern_to_bp_ids[i])
          kmer_by_bp_id.emplace(bp, seqs[i]);

      std::ofstream f(o.counts_file);
      if (!f) {
        std::cerr << "ERROR: cannot open --counts file for write: "
                  << o.counts_file << "\n";
        std::exit(EXIT_FAILURE);
      }
      f << "bp_id\tjxn_kmer\tn_total_hits\tn_unique_reads\n";

      // Build sortable rows. n_total_hits is a superset of n_unique_reads
      // (every counted read is also a raw hit), so iterating bp_id_total
      // covers every bp_id with >= 1 match.
      struct Row { const std::string* id; std::size_t total; std::size_t uniq; };
      std::vector<Row> rows;
      rows.reserve(bp_id_total.size());
      for (const auto& kv : bp_id_total) {
        const auto uit = bp_id_counts.find(kv.first);
        rows.push_back({ &kv.first, kv.second,
                         uit == bp_id_counts.end() ? 0 : uit->second });
      }
      std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) {
        if (a.total != b.total) return a.total > b.total;  // desc
        return *a.id < *b.id;                              // tiebreak asc
      });

      for (const auto& r : rows) {
        const auto kit = kmer_by_bp_id.find(*r.id);
        f << *r.id << '\t'
          << (kit == kmer_by_bp_id.end() ? std::string(".") : kit->second)
          << '\t' << r.total << '\t' << r.uniq << '\n';
      }
      if (o.verbose >= 1)
        std::cerr << "[extract-pairs] counts: "
                  << SeqLib::AddCommas(rows.size())
                  << " bp_ids with >= 1 matching alignment -> "
                  << o.counts_file << "\n";
    }

    if (qnames.empty()) {
      std::cerr << "[extract-pairs] no matches; output not created\n";
      std::exit(EXIT_SUCCESS);
    }

    if (o.verbose >= 1)
      std::cerr << "[extract-pairs] pass 2/2: extracting alignments for "
                << SeqLib::AddCommas(qnames.size()) << " qnames -> "
                << o.out_bam << "\n";
    const Pass2Stats p2 = extractByQname(o.in_bam, o.out_bam, qnames,
                                         o.threads, o.verbose, cl);
    if (o.verbose >= 1)
      std::cerr << "[extract-pairs] pass 2: "
                << SeqLib::AddCommas(p2.reads_in) << " in, "
                << SeqLib::AddCommas(p2.reads_out) << " written in "
                << std::fixed << std::setprecision(1) << p2.seconds << "s\n";
    out_records = p2.reads_out;
    qname_count = qnames.size();
  }

  // ---- conditional sort + index (shared between both modes) ----
  if (!input_sorted) {
    if (sortInPlace(o.out_bam, o.threads, o.verbose) != 0)
      std::exit(EXIT_FAILURE);
  }
  const int idx_rc = sam_index_build(o.out_bam.c_str(), 0);
  if (idx_rc < 0) {
    std::cerr << "ERROR: sam_index_build failed (rc=" << idx_rc << ") on "
              << o.out_bam << "\n";
    std::exit(EXIT_FAILURE);
  }

  if (o.verbose >= 1) {
    std::cerr << "[extract-pairs] done.\n";
    if (o.include_pairs)
      std::cerr << "  Matching qnames:  " << SeqLib::AddCommas(qname_count) << "\n";
    std::cerr << "  Output records:   " << SeqLib::AddCommas(out_records) << "\n"
              << "  Output:           " << o.out_bam
              << " (+ " << o.out_bam << ".bai)\n";
  }
}
