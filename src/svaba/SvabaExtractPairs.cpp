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
#include <unordered_set>
#include <vector>

#include "htslib/sam.h"
#include "zlib.h"

#include "SeqLib/BamHeader.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/SeqLibUtils.h"  // AddCommas

#include "SvabaOptions.h"  // SVABA_VERSION

namespace {

// ---------- options ----------

struct Opts {
  std::string in_bam;
  std::string out_bam;
  std::vector<std::string> seqs;
  std::string seq_file;
  int  threads  = 4;
  int  verbose  = 1;
  bool include_rc    = true;
  bool include_pairs = true;  // false → single-pass; emit only matched records
};

enum { OPT_NO_RC = 1000, OPT_NO_PAIRS };

const char* kShortOpts = "hi:o:s:f:t:v:";
const struct option kLongOpts[] = {
  { "help",     no_argument,       nullptr, 'h' },
  { "input",    required_argument, nullptr, 'i' },
  { "output",   required_argument, nullptr, 'o' },
  { "seq",      required_argument, nullptr, 's' },
  { "seq-file", required_argument, nullptr, 'f' },
  { "threads",  required_argument, nullptr, 't' },
  { "verbose",  required_argument, nullptr, 'v' },
  { "no-rc",    no_argument,       nullptr, OPT_NO_RC },
  { "no-pairs", no_argument,       nullptr, OPT_NO_PAIRS },
  { nullptr, 0, nullptr, 0 }
};

void printUsage() {
  std::cerr <<
    "Usage: svaba extract-pairs -i IN.bam -o OUT.bam (-s SEQ ... | -f FILE)\n"
    "                           [-t THREADS] [--no-rc] [--no-pairs] [-v V]\n"
    "\n"
    "  Extract all read pairs from IN.bam where either mate's SEQ contains\n"
    "  any of the given query sequences or their reverse complements.\n"
    "  Output is sorted (or input-order-preserved if input is sorted) and\n"
    "  indexed (.bai). Replaces scripts/extract_pairs_by_seq.sh.\n"
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
    "      --no-rc          Skip reverse-complement augmentation.\n"
    "      --no-pairs       Emit only records whose own SEQ matched. Skips\n"
    "                       the pair-mate / supplementary pickup, runs in a\n"
    "                       single BAM pass (~2x faster). Useful when you\n"
    "                       just want to inspect which reads contain a motif.\n"
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
      case OPT_NO_RC:    o.include_rc    = false; break;
      case OPT_NO_PAIRS: o.include_pairs = false; break;
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

// Locate the index of the `jxn_kmer` column in a bps header. Header lines
// always begin with '#'; columns are tab-separated. Returns -1 if the
// column isn't present (older v3 / v2 / legacy dumps without a junction
// kmer column).
int findJxnKmerCol(const std::string& header_line) {
  std::string h = header_line;
  if (!h.empty() && h.front() == '#') h.erase(0, 1);
  std::istringstream iss(h);
  std::string col;
  int idx = 0;
  while (std::getline(iss, col, '\t')) {
    if (col == "jxn_kmer") return idx;
    ++idx;
  }
  return -1;
}

// Read junction-kmer queries out of a bps.txt[.gz] dump. Skips rows
// where the kmer is "." (no precise junction — DSCRD-only, etc.) or
// empty. Caller is responsible for de-duplicating the returned list.
std::vector<std::string> readJxnKmersFromBps(const std::vector<std::string>& lines) {
  std::vector<std::string> out;
  if (lines.empty()) return out;

  // First non-empty line is the header — find jxn_kmer's position.
  int kmer_col = -1;
  std::string header;
  for (const auto& l : lines) {
    if (l.empty()) continue;
    header = l;
    kmer_col = findJxnKmerCol(l);
    break;
  }
  if (kmer_col < 0) {
    throw std::runtime_error(
      "bps file has no `jxn_kmer` column — too old to use as a query "
      "source for extract-pairs (need v4 or later; re-run `svaba run` "
      "to regenerate)");
  }

  // Walk data rows, splitting on tabs to grab the kmer column.
  // Tab-split is hot — std::find is fine here (rows are kilobyte-scale,
  // we touch one column per row, and this isn't on the BAM-stream path).
  for (const auto& l : lines) {
    if (l.empty() || l[0] == '#') continue;
    int col = 0;
    std::size_t s = 0;
    while (s < l.size()) {
      const std::size_t e = l.find('\t', s);
      const std::size_t end = (e == std::string::npos) ? l.size() : e;
      if (col == kmer_col) {
        std::string_view tok(&l[s], end - s);
        if (!tok.empty() && tok != ".") {
          // Normalize: upper-case + strip whitespace, just like the -s/-f
          // path. The contig is already upper-case in practice but be
          // defensive.
          std::string norm = normalizeSeq(tok);
          if (!norm.empty()) out.push_back(std::move(norm));
        }
        break;
      }
      if (e == std::string::npos) break;
      s = e + 1;
      ++col;
    }
  }
  return out;
}

// Public seq-file reader. Auto-detects bps.txt[.gz] vs plain one-per-line
// list by sniffing the first non-empty line. gzip is handled transparently
// for both formats.
std::vector<std::string> readSeqsFromFile(const std::string& path) {
  const auto lines = slurpLines(path);
  if (looksLikeBpsFile(lines))
    return readJxnKmersFromBps(lines);

  // Plain seq-per-line file (the original behavior).
  std::vector<std::string> out;
  for (const auto& l : lines) {
    if (l.empty() || l[0] == '#') continue;
    std::string s = normalizeSeq(l);
    if (!s.empty()) out.push_back(std::move(s));
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
// ambiguity codes like R/Y/M, which extract_pairs_by_seq.sh wouldn't have
// matched anyway (regex-engine-wise R != [AG], it's a literal R).
//
// We precompute a goto-function table during BFS so the search inner loop
// is one indexed read per character with no failure-link walking. The
// tradeoff (vs. the textbook trie+fail-link variant) is a 5x larger goto
// table; for 100 patterns of length 30 that's 1500 nodes * 5 ints ≈ 30 KB,
// trivially L1-resident.

class AhoCorasick {
 public:
  AhoCorasick() : nodes_(1) {}  // node 0 = root

  // Add one pattern. Pattern is normalized (upper-case, whitespace-free).
  // Returns true if added; false if rejected for containing non-ACGTN
  // characters or being empty.
  bool add(std::string_view pat) {
    if (pat.empty()) return false;
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
    nodes_[cur].terminal = true;
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
      // Propagate terminal: a state is terminal if any state on its
      // failure chain is terminal (a longer pattern ending here implies
      // any of its suffix patterns also matched).
      if (nodes_[uf].terminal) nodes_[u].terminal = true;
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
      if (nodes_[cur].terminal) return true;
    }
    return false;
  }

  std::size_t pattern_count() const { return n_added_; }
  std::size_t node_count()    const { return nodes_.size(); }
  void mark_added() { ++n_added_; }

 private:
  static constexpr int kAlphabet = 5;  // A,C,G,T,N

  struct Node {
    std::array<int, kAlphabet> child;
    int  fail     = 0;
    bool terminal = false;
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

struct Pass1Stats {
  std::size_t reads        = 0;
  std::size_t matched_reads = 0;
  std::size_t unique_qnames = 0;
  double      seconds      = 0;
};

Pass1Stats collectQnames(const std::string& in_bam,
                         const AhoCorasick& ac,
                         int threads,
                         int verbose,
                         std::unordered_set<std::string>& out_qnames) {
  using clock = std::chrono::steady_clock;
  Pass1Stats st;
  const auto t0 = clock::now();

  SeqLib::BamReader r;
  if (!r.Open(in_bam))
    throw std::runtime_error("pass1: cannot open " + in_bam);
  r.SetThreads(std::max(1, threads));

  // Cached chr name for the progress line. ChrName() hits the header every
  // call, so we resolve it only on (rare) chromosome change, not per record.
  // Sentinel -2 won't match any valid tid (-1 = unmapped, 0..N-1 = mapped).
  int32_t     cur_chr      = -2;
  int32_t     cur_pos      = -1;
  std::string cur_chr_name = "*";

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

    if (ac.searchNibbles(seq4, qlen)) {
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
                                   const std::string& cl) {
  using clock = std::chrono::steady_clock;
  SinglePassStats st;
  const auto t0 = clock::now();

  SeqLib::BamReader r;
  if (!r.Open(in_bam))
    throw std::runtime_error("single-pass: cannot open " + in_bam);
  r.SetThreads(std::max(1, threads));

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
  std::vector<std::string> seqs = o.seqs;
  for (auto& s : seqs) s = normalizeSeq(s);
  if (!o.seq_file.empty()) {
    auto from_file = readSeqsFromFile(o.seq_file);
    seqs.insert(seqs.end(),
                std::make_move_iterator(from_file.begin()),
                std::make_move_iterator(from_file.end()));
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

  // ---- build Aho-Corasick over patterns + reverse complements ----
  AhoCorasick ac;
  std::size_t added = 0, rejected = 0;
  for (const auto& s : seqs) {
    if (ac.add(s)) { ac.mark_added(); ++added; }
    else           { ++rejected; }
    if (o.include_rc) {
      const std::string rc = reverseComplement(s);
      if (ac.add(rc)) { ac.mark_added(); ++added; }
      else            { ++rejected; }
    }
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

  // ---- detect input sort order (cheap; just header inspection) ----
  // Done up-front so both modes (two-pass and single-pass) report it
  // consistently before they start the heavy work.
  bool input_sorted = false;
  {
    SeqLib::BamReader hdr_only;
    if (hdr_only.Open(o.in_bam))
      input_sorted = isCoordinateSorted(hdr_only.Header());
  }
  if (o.verbose >= 1)
    std::cerr << "[extract-pairs] input is "
              << (input_sorted ? "coord-sorted (sort step will be skipped)"
                               : "not coord-sorted (samtools sort will run)")
              << "\n";

  std::size_t out_records = 0;
  std::size_t qname_count = 0;  // unused in single-pass mode

  if (!o.include_pairs) {
    // ---- single-pass mode ----
    if (o.verbose >= 1)
      std::cerr << "[extract-pairs] single-pass mode (--no-pairs): emitting "
                   "only records whose own SEQ matched\n";
    const SinglePassStats sp = extractMatchedOnly(
        o.in_bam, o.out_bam, ac, o.threads, o.verbose, cl);
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
                << " for matching SEQ...\n";
    const Pass1Stats p1 = collectQnames(o.in_bam, ac, o.threads, o.verbose, qnames);
    if (o.verbose >= 1)
      std::cerr << "[extract-pairs] pass 1: "
                << SeqLib::AddCommas(p1.reads) << " reads, "
                << SeqLib::AddCommas(p1.matched_reads) << " matched, "
                << SeqLib::AddCommas(p1.unique_qnames) << " unique qnames in "
                << std::fixed << std::setprecision(1) << p1.seconds << "s\n";

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
