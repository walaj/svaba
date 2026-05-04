// SvabaPostprocess.cpp — see SvabaPostprocess.h for the design rationale.
//
// The flow for each suffix runs in its own std::thread:
//
//   [${ID}.${suffix}.bam]
//        |
//        |  (unless --dedup-only)
//        v
//   samtools sort -@ jthreads -m MEM -o *.sort.tmp.bam ${bam}    # shell out
//   rename *.sort.tmp.bam -> ${bam}
//        |
//        |  (if suffix in DEDUP set, unless --sort-only)
//        v
//   native streaming dedup+merge  -> *.dedup.tmp.bam             # htslib
//     (writer uses a header with our @PG line appended)
//   rename *.dedup.tmp.bam -> ${bam}                             # PG stamped
//        |
//        |  (only if dedup didn't run for this suffix)
//        v
//   samtools reheader <hdr+PG> ${bam} > *.reheader.tmp.bam       # shell out
//   rename *.reheader.tmp.bam -> ${bam}                          # PG stamped
//        |
//        v
//   sam_index_build(${bam}, 0)                                   # htslib
//     → ${bam}.bai
//
// Every final ${ID}.${suffix}.bam therefore ends up sorted, deduped (where
// applicable), carrying an @PG "svaba_postprocess" line chained onto the
// existing PG chain, and accompanied by a .bai index. The user never has to
// know the intermediate .postprocess.*.tmp.bam files existed.
//
// The dedup step is NOT a straight drop of duplicate (qname, flag) records
// at a locus: overlapping assembly windows can emit the same read twice
// with different bi:Z / bz:Z tags naming different supporting contigs.
// Dropping one would lose that support evidence. Instead we buffer per
// locus and union the comma-token bi/bz lists into the first record, using
// the same boundary-aware merge that SvabaOutputWriter.cpp::stamp_tag()
// performs during `svaba run`.
//
// Progress is reported per 1M processed reads, per suffix, serialized through
// a shared stderr mutex so lines don't interleave. We avoid any per-record
// system call — stat() is only used for file-existence checks and size-based
// ETA estimation.

#include "SvabaPostprocess.h"

#include <fcntl.h>
#include <getopt.h>
#include <glob.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "htslib/sam.h"
#include "zlib.h"

#include "SeqLib/BamHeader.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BamWriter.h"

#include "SvabaOptions.h"  // SVABA_VERSION

namespace {

// RAII helper to suppress stderr (htslib's "[E::idx_find_and_load]" messages
// when opening BAMs without a .bai — expected during merge/dedup phases).
// Thread-safe: uses a mutex so concurrent Open() calls from parallel workers
// don't race on the fd-level redirect.
std::mutex g_suppress_mu;

struct SuppressStderr {
  int saved_fd_, devnull_;
  SuppressStderr() : saved_fd_(-1), devnull_(-1) {
    g_suppress_mu.lock();
    saved_fd_ = dup(STDERR_FILENO);
    devnull_  = open("/dev/null", O_WRONLY);
    if (devnull_ >= 0) dup2(devnull_, STDERR_FILENO);
  }
  ~SuppressStderr() {
    if (devnull_ >= 0) { dup2(saved_fd_, STDERR_FILENO); close(devnull_); }
    if (saved_fd_ >= 0) close(saved_fd_);
    g_suppress_mu.unlock();
  }
  SuppressStderr(const SuppressStderr&) = delete;
  SuppressStderr& operator=(const SuppressStderr&) = delete;
};

// Serialize progress/log output across the per-suffix threads. We intentionally
// keep the critical section tiny — just the single ostringstream flush — so
// workers don't serialize on their own progress prints.
std::mutex g_stderr_mu;

// Suffixes eligible for dedup. `contigs` is sorted but never deduped because
// each contig is emitted exactly once by a single window and there is no
// overlap-induced duplication to clean up.
const std::vector<std::string> kDedupSuffixes = {
  "weird", "corrected", "discordant"
};

// Full suffix list to consider, in the same order as sort_output.sh. The
// contigs suffix goes last so that any shared stderr verbosity from it doesn't
// bury the (more interesting) dedup stats.
const std::vector<std::string> kAllSuffixes = {
  "weird", "corrected", "discordant", "contigs"
};

// ---------- options ----------
struct Opts {
  std::string id;
  int         threads   = 4;
  std::string mem       = "2G";  // per samtools sort thread; matches -m flag
  int         verbose   = 1;
  bool        sort_only  = false;  // skip dedup
  bool        dedup_only = false;  // skip sort (assume sorted)
  bool        split      = false;  // split BAMs by QNAME source prefix
};

enum { OPT_SORT_ONLY = 1000, OPT_DEDUP_ONLY, OPT_MEM, OPT_SPLIT };

const char* kShortOpts = "hi:t:m:v:";
const struct option kLongOpts[] = {
  { "help",       no_argument,       nullptr, 'h' },
  { "id",         required_argument, nullptr, 'i' },
  { "threads",    required_argument, nullptr, 't' },
  { "mem",        required_argument, nullptr, 'm' },
  { "verbose",    required_argument, nullptr, 'v' },
  { "sort-only",  no_argument,       nullptr, OPT_SORT_ONLY  },
  { "dedup-only", no_argument,       nullptr, OPT_DEDUP_ONLY },
  { "split",      no_argument,       nullptr, OPT_SPLIT      },
  { nullptr, 0, nullptr, 0 }
};

void printUsage() {
  std::cerr <<
    "Usage: svaba postprocess -i <ID> [options]\n"
    "\n"
    "  Coord-sort (via samtools sort) and natively streaming-dedup\n"
    "  svaba's per-suffix output BAMs. Much faster than the legacy\n"
    "  samtools-view | awk | samtools-view pipeline.\n"
    "\n"
    "  Expects files at ${ID}.${suffix}.bam for suffix in:\n"
    "    weird corrected discordant contigs\n"
    "\n"
    "  -i, --id <str>            Analysis ID (required).\n"
    "  -t, --threads <n>         Total threads budget, split across concurrent\n"
    "                            per-suffix jobs. [4]\n"
    "  -m, --mem <str>           Memory per samtools-sort thread (e.g. 2G). [2G]\n"
    "      --sort-only           Only sort; skip dedup.\n"
    "      --dedup-only          Only dedup; skip sort. Assumes input is already\n"
    "                            coord-sorted; output is undefined otherwise.\n"
    "      --split               Split dedup-eligible BAMs by QNAME source prefix\n"
    "                            (first 4 chars). Produces ${ID}.${suffix}.${PREFIX}.bam\n"
    "                            for each observed prefix. Runs after dedup.\n"
    "  -v, --verbose <0-3>       Verbosity. [1]\n"
    "  -h, --help                This message.\n";
}

Opts parseOpts(int argc, char** argv) {
  Opts o;
  if (argc <= 2) { printUsage(); std::exit(EXIT_FAILURE); }

  for (int c; (c = getopt_long(argc, argv, kShortOpts, kLongOpts, nullptr)) != -1;) {
    std::istringstream arg(optarg ? optarg : "");
    switch (c) {
      case 'h': printUsage(); std::exit(EXIT_SUCCESS);
      case 'i': arg >> o.id; break;
      case 't': arg >> o.threads; break;
      case 'm': arg >> o.mem; break;
      case 'v': arg >> o.verbose; break;
      case OPT_SORT_ONLY:  o.sort_only  = true; break;
      case OPT_DEDUP_ONLY: o.dedup_only = true; break;
      case OPT_SPLIT:      o.split      = true; break;
      default: printUsage(); std::exit(EXIT_FAILURE);
    }
  }
  if (o.id.empty()) {
    std::cerr << "ERROR: --id is required\n";
    printUsage();
    std::exit(EXIT_FAILURE);
  }
  if (o.sort_only && o.dedup_only) {
    std::cerr << "ERROR: --sort-only and --dedup-only are mutually exclusive\n";
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

off_t fileSize(const std::string& p) {
  struct stat st{};
  if (::stat(p.c_str(), &st) != 0) return 0;
  return st.st_size;
}

// ---------- glob helper ----------

// Return all paths matching a glob pattern, sorted. Uses POSIX glob(3).
// Returns empty vector on no match or error. Caller owns the result.
std::vector<std::string> globPaths(const std::string& pattern) {
  std::vector<std::string> result;
  glob_t g{};
  const int rc = ::glob(pattern.c_str(), GLOB_NOSORT, nullptr, &g);
  if (rc == 0) {
    result.reserve(g.gl_pathc);
    for (size_t i = 0; i < g.gl_pathc; ++i)
      result.emplace_back(g.gl_pathv[i]);
    std::sort(result.begin(), result.end());
  }
  ::globfree(&g);
  return result;
}

// Extract the numeric thread index from a filename like
// "${ID}.thread${N}.${suffix}.bam" or "${ID}.thread${N}.r2c.txt.gz".
// Returns -1 if parsing fails.
int extractThreadIndex(const std::string& path) {
  // Find ".thread" then parse the number after it.
  auto pos = path.find(".thread");
  if (pos == std::string::npos) return -1;
  pos += 7;  // skip ".thread"
  int idx = 0;
  bool found_digit = false;
  while (pos < path.size() && std::isdigit(static_cast<unsigned char>(path[pos]))) {
    idx = idx * 10 + (path[pos] - '0');
    found_digit = true;
    ++pos;
  }
  return found_digit ? idx : -1;
}

// Sort paths by their embedded thread index (ascending). Thread 1 first.
void sortByThreadIndex(std::vector<std::string>& paths) {
  std::sort(paths.begin(), paths.end(),
    [](const std::string& a, const std::string& b) {
      return extractThreadIndex(a) < extractThreadIndex(b);
    });
}

// Thread-safe stderr log. Accepts anything streamable into an ostringstream.
template <typename... Args>
void logLine(Args&&... args) {
  std::ostringstream oss;
  (oss << ... << std::forward<Args>(args));
  std::lock_guard<std::mutex> lg(g_stderr_mu);
  std::cerr << oss.str() << std::endl;
}

// Quote a path for use inside a /bin/sh command string. Wraps in single
// quotes and escapes embedded single quotes as '\''. This is safe against
// arbitrary IDs containing spaces or shell metacharacters.
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

// Forward decl: defined further down in the "per-suffix pipeline" section.
// Declared here so the reheader helper below can use std::rename via the same
// diagnostics-rich wrapper the dedup path uses.
void renameOrThrow(const std::string& from, const std::string& to);

// ---------- Phase 0: merge per-thread outputs ----------
//
// svaba run with -p N produces per-thread BAMs (${ID}.thread${N}.${suffix}.bam)
// and per-thread r2c TSVs (${ID}.thread${N}.r2c.txt.gz). This phase coalesces
// each family into a single ${ID}.${suffix}.bam / ${ID}.r2c.txt.gz.
//
// BAM merge: open all per-thread BAMs, concatenate records into one writer.
// The output is NOT coordinate-sorted (Phase 1 handles that), so we do a
// simple sequential dump — no priority queue needed. We use the header from
// the first input (all per-thread BAMs share the same ref/RG/PG header
// because they're all opened from the same SvabaSharedConfig).
//
// r2c merge: binary byte-copy of gzip streams. Gzip is concatenation-safe
// per RFC 1952, so cat(a.gz, b.gz) → valid .gz whose decompressed content
// is decompressed(a) + decompressed(b). Thread 1's file is placed first
// (it contains the TSV column-header line).

// Merge per-thread BAMs for a single suffix. Returns number of records merged,
// or 0 if nothing to merge (no per-thread files found). Removes per-thread
// inputs on success.
size_t mergeThreadBams(const std::string& id,
                       const std::string& suffix,
                       int threads,
                       int verbose) {
  const std::string pattern = id + ".thread*." + suffix + ".bam";
  std::vector<std::string> inputs = globPaths(pattern);
  if (inputs.empty()) return 0;

  const std::string target = id + "." + suffix + ".bam";

  // Single file — just rename.
  if (inputs.size() == 1) {
    logLine("[merge] mv ", inputs[0], " -> ", target);
    renameOrThrow(inputs[0], target);
    return 0;  // no records actually streamed
  }

  sortByThreadIndex(inputs);

  logLine("[merge] merging ", inputs.size(), " thread BAMs for '", suffix,
          "' -> ", target);

  const auto t0 = std::chrono::steady_clock::now();

  // Open all readers; use header from the first.
  // Suppress htslib "[E::idx_find_and_load]" — per-thread BAMs have no .bai.
  std::vector<SeqLib::BamReader> readers(inputs.size());
  {
    SuppressStderr quiet;
    for (size_t i = 0; i < inputs.size(); ++i) {
      if (!readers[i].Open(inputs[i]))
        throw std::runtime_error("merge: cannot open " + inputs[i]);
    }
  }

  SeqLib::BamWriter writer;
  const std::string tmp = target + ".postprocess.merge.tmp.bam";
  if (!writer.Open(tmp))
    throw std::runtime_error("merge: cannot open output " + tmp);
  writer.SetHeader(readers[0].Header());
  if (!writer.WriteHeader())
    throw std::runtime_error("merge: WriteHeader failed on " + tmp);
  // Enable BGZF thread pool on the writer for throughput.
  writer.SetThreads(std::max(1, threads));

  size_t total = 0;
  for (size_t i = 0; i < readers.size(); ++i) {
    while (auto rec = readers[i].Next()) {
      if (!writer.WriteRecord(*rec))
        throw std::runtime_error("merge: WriteRecord failed on " + tmp);
      ++total;
    }
  }

  if (!writer.Close())
    throw std::runtime_error("merge: Close failed on " + tmp);

  renameOrThrow(tmp, target);

  // Remove per-thread inputs.
  for (const auto& f : inputs)
    ::unlink(f.c_str());

  const double elapsed = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - t0).count();
  logLine("[merge] ", suffix, ": ", total, " records from ", inputs.size(),
          " files in ", std::fixed, std::setprecision(1), elapsed, "s");
  return total;
}

// Merge per-thread r2c.txt.gz files via binary concatenation.
// Returns number of files merged, or 0 if nothing to merge.
// Removes per-thread inputs on success.
size_t mergeThreadR2c(const std::string& id, int verbose) {
  const std::string pattern = id + ".thread*.r2c.txt.gz";
  std::vector<std::string> inputs = globPaths(pattern);
  if (inputs.empty()) return 0;

  const std::string target = id + ".r2c.txt.gz";

  // Single file — just rename.
  if (inputs.size() == 1) {
    logLine("[merge] mv ", inputs[0], " -> ", target);
    renameOrThrow(inputs[0], target);
    return 1;
  }

  // Sort by thread index so thread 1 (which has the header line) is first.
  sortByThreadIndex(inputs);

  logLine("[merge] merging ", inputs.size(), " thread r2c files -> ", target);

  const auto t0 = std::chrono::steady_clock::now();
  const std::string tmp = target + ".postprocess.merge.tmp.gz";

  // Binary byte-copy — gzip concatenation.
  {
    std::ofstream out(tmp, std::ios::binary);
    if (!out)
      throw std::runtime_error("merge: cannot open " + tmp + " for writing");

    constexpr size_t BUF_SIZE = 256 * 1024;  // 256 KB copy buffer
    std::vector<char> buf(BUF_SIZE);

    for (const auto& f : inputs) {
      std::ifstream in(f, std::ios::binary);
      if (!in)
        throw std::runtime_error("merge: cannot open " + f + " for reading");
      while (in.read(buf.data(), BUF_SIZE) || in.gcount() > 0) {
        out.write(buf.data(), in.gcount());
        if (!out)
          throw std::runtime_error("merge: write error on " + tmp);
      }
    }
  }

  renameOrThrow(tmp, target);

  // Remove per-thread inputs.
  for (const auto& f : inputs)
    ::unlink(f.c_str());

  const double elapsed = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - t0).count();
  logLine("[merge] r2c: ", inputs.size(), " files -> ", target,
          " in ", std::fixed, std::setprecision(1), elapsed, "s");
  return inputs.size();
}

// Run Phase 0: merge all per-thread outputs for all suffixes + r2c.
void runMergePhase(const std::string& id, int threads, int verbose) {
  // BAM suffixes that svaba emits per-thread (contigs is NOT per-thread).
  const std::vector<std::string> merge_suffixes = {
    "discordant", "weird", "corrected"
  };

  bool did_anything = false;
  for (const auto& suffix : merge_suffixes) {
    // Check if per-thread files exist BEFORE merging (since merge removes them).
    const std::string pat = id + ".thread*." + suffix + ".bam";
    if (!globPaths(pat).empty()) {
      mergeThreadBams(id, suffix, threads, verbose);
      did_anything = true;
    }
  }

  {
    const std::string pat = id + ".thread*.r2c.txt.gz";
    if (!globPaths(pat).empty()) {
      mergeThreadR2c(id, verbose);
      did_anything = true;
    }
  }

  if (!did_anything && verbose >= 1)
    logLine("[merge] no per-thread files found; nothing to merge");
}

// ---------- bps.txt.gz sort / dedup / PASS filter ----------
//
// Replaces the GNU sort + awk pipeline from svaba_postprocess.sh step 3.
//
// Approach:
//   1. Decompress bps.txt.gz into a contiguous slab (one big vector<char>).
//   2. Build a lightweight index: one BpsSortEntry per line with parsed sort
//      keys (chr1_id, pos1, strand1, chr2_id, pos2, strand2, maxlod) plus
//      offset/length into the slab. ~28 bytes per entry.
//   3. std::sort the index by the same key order as the shell script:
//      chr1(V), pos1(n), strand1, chr2(V), pos2(n), strand2, maxlod(desc).
//   4. Single pass over the sorted index to write three gzipped outputs:
//      - .bps.sorted.txt.gz         (all rows, sorted)
//      - .bps.sorted.dedup.txt.gz   (first row per unique breakpoint pair)
//      - .bps.sorted.dedup.pass.txt.gz (PASS only from the dedup set)
//      - .bps.sorted.dedup.pass.somatic.txt.gz (PASS + somatic)
//
// Memory: slab ≈ uncompressed size of bps.txt.gz (typically 25-500 MB for WGS),
//         index ≈ 28 bytes × N_lines (typically 1-30 MB). Well within modern
//         machine RAM for even the noisiest runs.
//
// Version-sort (V) for chromosome names: we parse chrN → numeric ID using the
// standard human convention (chr1=1..chr22=22, chrX=23, chrY=24, chrM=25,
// other=1000+lexicographic). This matches GNU sort -V on chr-prefixed names
// and sorts non-standard contigs at the end.

// Column indices (0-based) in bps.txt.gz. From BreakPoint::header():
//   0=chr1, 1=pos1, 2=strand1, 3=chr2, 4=pos2, 5=strand2
//   29=contig_and_region (cname), 31=conf, 35=somatic, 36=somlod, 37=maxlod
// Note: the shell script uses 1-based cols (30,32,36,37,38).
constexpr int kBpsCol_Chr1   = 0;
constexpr int kBpsCol_Pos1   = 1;
constexpr int kBpsCol_Strand1= 2;
constexpr int kBpsCol_Chr2   = 3;
constexpr int kBpsCol_Pos2   = 4;
constexpr int kBpsCol_Strand2= 5;
constexpr int kBpsCol_Cname  = 29;  // contig_and_region (join key for r2c)
constexpr int kBpsCol_Conf   = 31;  // "PASS" / "LOWMAPQ" / etc
constexpr int kBpsCol_Somatic= 35;  // "0" or "1"
constexpr int kBpsCol_Somlod = 36;  // somatic log-odds (lower = more conservative)
constexpr int kBpsCol_Maxlod = 37;

// Parse a chromosome name to a numeric sort key matching GNU sort -V behavior:
// chr1..chr22 → 1..22, chrX → 23, chrY → 24, chrM → 25.
// Anything else → 1000 + first 8 chars as a stable hash-like key (for
// lexicographic ordering among unknowns).
int chrSortKey(std::string_view name) {
  // Strip "chr" prefix if present.
  std::string_view base = name;
  if (base.size() > 3 && base[0] == 'c' && base[1] == 'h' && base[2] == 'r')
    base = base.substr(3);

  if (base == "X") return 23;
  if (base == "Y") return 24;
  if (base == "M" || base == "MT") return 25;

  // Try parsing as integer.
  int val = 0;
  bool is_num = true;
  for (char c : base) {
    if (c >= '0' && c <= '9') val = val * 10 + (c - '0');
    else { is_num = false; break; }
  }
  if (is_num && val >= 1 && val <= 22) return val;

  // Unknown chromosome — use a stable ordering based on the string.
  // We want deterministic sort so use a simple djb2-ish hash that
  // preserves lexicographic ordering for short strings.
  // Actually just use 1000 + lexicographic comparison in the comparator.
  return 1000;  // all unknowns tie here; comparator breaks tie lexicographically
}

struct BpsSortEntry {
  int32_t  chr1_key;    // chrSortKey(chr1)
  int32_t  pos1;
  int32_t  chr2_key;
  int32_t  pos2;
  char     strand1;
  char     strand2;
  float    maxlod;      // for descending sort; float precision is fine
  float    somlod;      // somatic LOD; used for dedup tie-breaking (keep lowest)
  uint32_t offset;      // byte offset into slab
  uint32_t length;      // line length (excl newline)
  // For tie-breaking unknown chromosomes lexicographically, we store
  // short views. These point into the slab, so they're valid as long as
  // the slab is alive.
  std::string_view chr1_name;
  std::string_view chr2_name;
};

// Comparator: chr1(V asc), pos1(n asc), strand1(asc), chr2(V asc), pos2(n asc),
//             strand2(asc), maxlod(desc).
struct BpsSortCmp {
  bool operator()(const BpsSortEntry& a, const BpsSortEntry& b) const {
    if (a.chr1_key != b.chr1_key) return a.chr1_key < b.chr1_key;
    // Tie-break unknowns lexicographically.
    if (a.chr1_key >= 1000 && a.chr1_name != b.chr1_name)
      return a.chr1_name < b.chr1_name;
    if (a.pos1 != b.pos1) return a.pos1 < b.pos1;
    if (a.strand1 != b.strand1) return a.strand1 < b.strand1;
    if (a.chr2_key != b.chr2_key) return a.chr2_key < b.chr2_key;
    if (a.chr2_key >= 1000 && a.chr2_name != b.chr2_name)
      return a.chr2_name < b.chr2_name;
    if (a.pos2 != b.pos2) return a.pos2 < b.pos2;
    if (a.strand2 != b.strand2) return a.strand2 < b.strand2;
    // maxlod descending (so highest comes first within same breakpoint pair).
    return a.maxlod > b.maxlod;
  }
};

// Get the Nth tab-separated field from a line (0-based). Returns empty view
// if col is out of range.
std::string_view getField(const char* line, uint32_t len, int col) {
  const char* p = line;
  const char* end = line + len;
  int cur = 0;
  while (cur < col && p < end) {
    if (*p == '\t') ++cur;
    ++p;
  }
  if (cur < col) return {};
  const char* start = p;
  while (p < end && *p != '\t') ++p;
  return std::string_view(start, p - start);
}

// Parse an integer from a string_view. Returns 0 on failure.
int32_t parseInt(std::string_view sv) {
  int32_t val = 0;
  bool neg = false;
  size_t i = 0;
  if (!sv.empty() && sv[0] == '-') { neg = true; i = 1; }
  for (; i < sv.size(); ++i) {
    char c = sv[i];
    if (c < '0' || c > '9') break;
    val = val * 10 + (c - '0');
  }
  return neg ? -val : val;
}

// Parse a float from a string_view. Returns 0 on failure or "NA".
float parseFloat(std::string_view sv) {
  if (sv.empty() || sv == "NA" || sv == "na") return -1e30f;  // sorts last
  // Use strtof on a null-terminated copy (string_view isn't guaranteed null-term).
  char buf[32];
  size_t n = std::min(sv.size(), size_t(31));
  std::memcpy(buf, sv.data(), n);
  buf[n] = '\0';
  return std::strtof(buf, nullptr);
}

// Main bps sort/dedup/filter. Returns number of data lines processed.
// Writes up to 4 output files:
//   ${id}.bps.sorted.txt.gz
//   ${id}.bps.sorted.dedup.txt.gz
//   ${id}.bps.sorted.dedup.pass.txt.gz
//   ${id}.bps.sorted.dedup.pass.somatic.txt.gz
// Populates pass_cnames / som_cnames with cnames of PASS / PASS+somatic
// winners (for downstream r2c filtering).
size_t sortDedupFilterBps(const std::string& id, int verbose,
                          std::unordered_set<std::string>& pass_cnames,
                          std::unordered_set<std::string>& som_cnames) {
  // Find the input: try .bps.txt.gz first, then .bps.txt
  const std::string gz_path  = id + ".bps.txt.gz";
  const std::string txt_path = id + ".bps.txt";

  std::string in_path;
  bool is_gzipped = false;
  if (fileExists(gz_path))       { in_path = gz_path;  is_gzipped = true; }
  else if (fileExists(txt_path)) { in_path = txt_path; is_gzipped = false; }
  else {
    if (verbose >= 1)
      logLine("[bps] skipping: neither ", gz_path, " nor ", txt_path, " found");
    return 0;
  }

  logLine("[bps] sort + dedup + filter: ", in_path);
  const auto t0 = std::chrono::steady_clock::now();

  // --- Step 1: decompress into slab + build index --------------------------
  // Read the file line by line. We use gzFile (from zlib/htslib) for .gz,
  // or plain ifstream for .txt.
  std::vector<char> slab;
  std::vector<BpsSortEntry> index;
  std::string header_line;

  // Reserve a rough estimate. Typical bps: 200k lines × 500 bytes = 100 MB.
  slab.reserve(64 * 1024 * 1024);  // 64 MB initial
  index.reserve(200000);

  // Use zlib's gzFile for reading (available via htslib linkage).
  gzFile gz = nullptr;
  if (is_gzipped) {
    gz = gzopen(in_path.c_str(), "rb");
    if (!gz) throw std::runtime_error("bps: cannot gzopen " + in_path);
    // Set a large internal buffer for better decompression throughput.
    gzbuffer(gz, 256 * 1024);
  }
  std::ifstream plain_in;
  if (!is_gzipped) {
    plain_in.open(in_path, std::ios::binary);
    if (!plain_in) throw std::runtime_error("bps: cannot open " + in_path);
  }

  // Line reading buffer.
  std::vector<char> linebuf(8192);
  size_t n_lines = 0;

  auto readLine = [&](std::string& out) -> bool {
    out.clear();
    if (is_gzipped) {
      while (true) {
        if (!gzgets(gz, linebuf.data(), static_cast<int>(linebuf.size())))
          return !out.empty();
        out.append(linebuf.data());
        if (!out.empty() && out.back() == '\n') { out.pop_back(); return true; }
        // Line longer than buffer — keep reading.
      }
    } else {
      if (!std::getline(plain_in, out)) return false;
      // Strip trailing \r if present (Windows line endings).
      if (!out.empty() && out.back() == '\r') out.pop_back();
      return true;
    }
  };

  std::string line;
  while (readLine(line)) {
    // Header line (starts with #) — save and skip.
    if (n_lines == 0 && !line.empty() && line[0] == '#') {
      header_line = line;
      continue;
    }
    ++n_lines;

    // Append line to slab.
    uint32_t offset = static_cast<uint32_t>(slab.size());
    uint32_t length = static_cast<uint32_t>(line.size());
    slab.insert(slab.end(), line.begin(), line.end());

    // Parse sort key fields from the slab (points into slab are stable because
    // we reserved upfront — but actually vector can realloc. We'll fix views
    // after loading is done. For now store offset/length and parse at the end.)
    // Actually: parse from `line` directly for key extraction, then store
    // offset/length. The string_view fields (chr names) will be fixed up
    // to point into the slab after all lines are loaded.
    BpsSortEntry entry{};
    entry.offset = offset;
    entry.length = length;

    // Parse fields from `line` (still alive).
    auto f_chr1   = getField(line.data(), length, kBpsCol_Chr1);
    auto f_pos1   = getField(line.data(), length, kBpsCol_Pos1);
    auto f_strand1= getField(line.data(), length, kBpsCol_Strand1);
    auto f_chr2   = getField(line.data(), length, kBpsCol_Chr2);
    auto f_pos2   = getField(line.data(), length, kBpsCol_Pos2);
    auto f_strand2= getField(line.data(), length, kBpsCol_Strand2);
    auto f_maxlod = getField(line.data(), length, kBpsCol_Maxlod);
    auto f_somlod = getField(line.data(), length, kBpsCol_Somlod);

    entry.chr1_key = chrSortKey(f_chr1);
    entry.pos1     = parseInt(f_pos1);
    entry.strand1  = f_strand1.empty() ? '+' : f_strand1[0];
    entry.chr2_key = chrSortKey(f_chr2);
    entry.pos2     = parseInt(f_pos2);
    entry.strand2  = f_strand2.empty() ? '+' : f_strand2[0];
    entry.maxlod   = parseFloat(f_maxlod);
    entry.somlod   = parseFloat(f_somlod);
    // chr_name views will be set after slab is finalized (no more reallocs).
    entry.chr1_name = {};
    entry.chr2_name = {};

    index.push_back(entry);
  }

  if (gz) gzclose(gz);

  // Fix up string_view fields to point into the finalized slab.
  for (auto& e : index) {
    const char* base = slab.data() + e.offset;
    e.chr1_name = getField(base, e.length, kBpsCol_Chr1);
    e.chr2_name = getField(base, e.length, kBpsCol_Chr2);
  }

  if (verbose >= 1)
    logLine("[bps] loaded ", n_lines, " data lines (",
            slab.size() / (1024 * 1024), " MB slab, ",
            index.size() * sizeof(BpsSortEntry) / (1024 * 1024), " MB index)");

  // --- Step 2: sort --------------------------------------------------------
  {
    const auto ts = std::chrono::steady_clock::now();
    std::sort(index.begin(), index.end(), BpsSortCmp{});
    const double sort_sec = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - ts).count();
    if (verbose >= 1)
      logLine("[bps] sorted in ", std::fixed, std::setprecision(2), sort_sec, "s");
  }

  // --- Step 3: write outputs (sorted, dedup, pass, somatic) ----------------
  const std::string out_sorted  = id + ".bps.sorted.txt.gz";
  const std::string out_dedup   = id + ".bps.sorted.dedup.txt.gz";
  const std::string out_pass    = id + ".bps.sorted.dedup.pass.txt.gz";
  const std::string out_som     = id + ".bps.sorted.dedup.pass.somatic.txt.gz";

  gzFile gz_sorted = gzopen(out_sorted.c_str(), "wb");
  gzFile gz_dedup  = gzopen(out_dedup.c_str(),  "wb");
  gzFile gz_pass   = gzopen(out_pass.c_str(),   "wb");
  gzFile gz_som    = gzopen(out_som.c_str(),    "wb");

  if (!gz_sorted || !gz_dedup || !gz_pass || !gz_som)
    throw std::runtime_error("bps: cannot open one or more output gz files");

  // Write header to all outputs.
  auto writeHeader = [&](gzFile f) {
    if (!header_line.empty()) {
      gzwrite(f, header_line.data(), static_cast<unsigned>(header_line.size()));
      gzwrite(f, "\n", 1);
    }
  };
  writeHeader(gz_sorted);
  writeHeader(gz_dedup);
  writeHeader(gz_pass);
  writeHeader(gz_som);

  // Dedup: keep one representative per breakpoint PAIR, canonicalized so
  // that AB == BA. An SV at chr1:100+ ↔ chr5:900- can be emitted from two
  // overlapping assembly windows with the breakends in either order:
  //   row A: chr1 100 + chr5 900 -
  //   row B: chr5 900 - chr1 100 +
  // The old shell-script dedup used adjacent-comparison on the raw key and
  // missed this. We canonicalize: always put the "lesser" breakend first
  // (compare chr_key, then pos, then strand), then pick the representative
  // with the LOWEST somlod — conservative against calling germline events
  // as somatic. If somlod ties, prefer higher maxlod.
  //
  // Two-pass approach:
  //   Pass 1: for each canonical pair, find the index of the winner (lowest
  //           somlod) and store it in a set.
  //   Pass 2: walk sorted index, write all rows to sorted output; write
  //           only winners to dedup/pass/somatic outputs.

  struct CanonicalKey {
    int32_t chr_a_key, pos_a, chr_b_key, pos_b;
    char strand_a, strand_b;

    bool operator==(const CanonicalKey& o) const {
      return chr_a_key == o.chr_a_key && pos_a == o.pos_a &&
             strand_a == o.strand_a &&
             chr_b_key == o.chr_b_key && pos_b == o.pos_b &&
             strand_b == o.strand_b;
    }
  };

  struct CanonicalKeyHash {
    size_t operator()(const CanonicalKey& k) const {
      size_t h = std::hash<int32_t>{}(k.chr_a_key);
      h ^= std::hash<int32_t>{}(k.pos_a) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= std::hash<int32_t>{}(k.chr_b_key) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= std::hash<int32_t>{}(k.pos_b) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= std::hash<char>{}(k.strand_a) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= std::hash<char>{}(k.strand_b) + 0x9e3779b9 + (h << 6) + (h >> 2);
      return h;
    }
  };

  // Canonicalize: put the "lesser" breakend as A.
  auto makeCanonical = [](const BpsSortEntry& e) -> CanonicalKey {
    bool swap = false;
    if (e.chr1_key != e.chr2_key)
      swap = e.chr1_key > e.chr2_key;
    else if (e.pos1 != e.pos2)
      swap = e.pos1 > e.pos2;
    else
      swap = e.strand1 > e.strand2;

    if (!swap)
      return {e.chr1_key, e.pos1, e.chr2_key, e.pos2, e.strand1, e.strand2};
    else
      return {e.chr2_key, e.pos2, e.chr1_key, e.pos1, e.strand2, e.strand1};
  };

  // --- Pass 1: select winners (lowest somlod per canonical pair) -----------
  // Map canonical key → index into `index` vector of the best row.
  std::unordered_map<CanonicalKey, size_t, CanonicalKeyHash> winners;
  winners.reserve(index.size() / 2);

  for (size_t i = 0; i < index.size(); ++i) {
    CanonicalKey ckey = makeCanonical(index[i]);
    auto [it, inserted] = winners.try_emplace(ckey, i);
    if (!inserted) {
      // Already have a candidate — replace if this row has lower somlod,
      // or same somlod but higher maxlod (better evidence it's real).
      const auto& cur_best = index[it->second];
      const auto& challenger = index[i];
      if (challenger.somlod < cur_best.somlod ||
          (challenger.somlod == cur_best.somlod &&
           challenger.maxlod > cur_best.maxlod)) {
        it->second = i;
      }
    }
  }

  // Build a fast lookup: is index[i] a winner?
  std::vector<bool> is_winner(index.size(), false);
  for (const auto& [key, idx] : winners)
    is_winner[idx] = true;

  if (verbose >= 1)
    logLine("[bps] dedup: ", index.size(), " rows -> ", winners.size(),
            " unique canonical pairs (lowest-somlod selection)");

  // --- Pass 2: write outputs -----------------------------------------------
  size_t n_sorted = 0, n_dedup = 0, n_pass = 0, n_som = 0;

  for (size_t i = 0; i < index.size(); ++i) {
    const auto& e = index[i];
    const char* line_ptr = slab.data() + e.offset;
    uint32_t line_len = e.length;

    // Write to sorted output (all rows, no dedup).
    gzwrite(gz_sorted, line_ptr, static_cast<unsigned>(line_len));
    gzwrite(gz_sorted, "\n", 1);
    ++n_sorted;

    // Only winners go to dedup/pass/somatic.
    if (!is_winner[i])
      continue;

    // Write to dedup output.
    gzwrite(gz_dedup, line_ptr, static_cast<unsigned>(line_len));
    gzwrite(gz_dedup, "\n", 1);
    ++n_dedup;

    // Check PASS.
    auto conf = getField(line_ptr, line_len, kBpsCol_Conf);
    if (conf != "PASS")
      continue;

    gzwrite(gz_pass, line_ptr, static_cast<unsigned>(line_len));
    gzwrite(gz_pass, "\n", 1);
    ++n_pass;

    // Collect cname for r2c filtering.
    auto cname = getField(line_ptr, line_len, kBpsCol_Cname);
    if (!cname.empty())
      pass_cnames.emplace(cname);

    // Check somatic.
    auto somatic = getField(line_ptr, line_len, kBpsCol_Somatic);
    if (somatic == "1") {
      gzwrite(gz_som, line_ptr, static_cast<unsigned>(line_len));
      gzwrite(gz_som, "\n", 1);
      ++n_som;
      if (!cname.empty())
        som_cnames.emplace(cname);
    }
  }

  gzclose(gz_sorted);
  gzclose(gz_dedup);
  gzclose(gz_pass);
  gzclose(gz_som);

  const double elapsed = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - t0).count();
  logLine("[bps] done in ", std::fixed, std::setprecision(1), elapsed, "s: ",
          n_sorted, " sorted, ", n_dedup, " dedup, ",
          n_pass, " pass, ", n_som, " somatic");

  return n_lines;
}

// ---------- r2c.txt.gz PASS / somatic filter ----------
//
// Streams ${id}.r2c.txt.gz and writes two filtered outputs:
//   ${id}.r2c.pass.txt.gz          — rows for PASS contigs
//   ${id}.r2c.pass.somatic.txt.gz  — rows for PASS+somatic contigs
//
// The r2c TSV has record_type in col 0 and contig_name in col 1.
// Both "contig" and "read" rows share the same contig_name (col 1),
// so filtering on col 1 membership in the cname set keeps contig+read
// rows together. The header row (first line) goes to both outputs.

constexpr int kR2cCol_RecordType = 0;
constexpr int kR2cCol_ContigName = 1;

size_t filterR2c(const std::string& id,
                 const std::unordered_set<std::string>& pass_cnames,
                 const std::unordered_set<std::string>& som_cnames,
                 int verbose) {
  const std::string in_path = id + ".r2c.txt.gz";
  if (!fileExists(in_path)) {
    if (verbose >= 1)
      logLine("[r2c] skipping: ", in_path, " not found");
    return 0;
  }

  if (pass_cnames.empty()) {
    if (verbose >= 1)
      logLine("[r2c] skipping: no PASS cnames to filter on");
    return 0;
  }

  logLine("[r2c] filtering ", in_path, " (",
          pass_cnames.size(), " PASS cnames, ",
          som_cnames.size(), " somatic cnames)");
  const auto t0 = std::chrono::steady_clock::now();

  const std::string out_pass = id + ".r2c.pass.txt.gz";
  const std::string out_som  = id + ".r2c.pass.somatic.txt.gz";

  gzFile gz_in   = gzopen(in_path.c_str(), "rb");
  gzFile gz_pass = gzopen(out_pass.c_str(), "wb1");  // level 1: ~3x faster compress
  gzFile gz_som  = gzopen(out_som.c_str(),  "wb1");

  if (!gz_in || !gz_pass || !gz_som)
    throw std::runtime_error("r2c filter: cannot open input/output gz files");

  gzbuffer(gz_in,   256 * 1024);
  gzbuffer(gz_pass, 256 * 1024);
  gzbuffer(gz_som,  256 * 1024);

  // Line reading with arbitrary-length support.
  std::vector<char> linebuf(8192);
  std::string line;

  auto readLine = [&]() -> bool {
    line.clear();
    while (true) {
      if (!gzgets(gz_in, linebuf.data(), static_cast<int>(linebuf.size())))
        return !line.empty();
      line.append(linebuf.data());
      if (!line.empty() && line.back() == '\n') { line.pop_back(); return true; }
    }
  };

  size_t n_in = 0, n_pass = 0, n_som = 0;
  bool first_line = true;

  // Cache: r2c rows are grouped by contig (one contig row then its read rows).
  // Avoid repeated hash lookups for the same cname by caching the last result.
  std::string last_cname;
  bool last_in_pass = false;
  bool last_in_som  = false;

  while (readLine()) {
    ++n_in;

    // Header line — write to both outputs.
    if (first_line) {
      first_line = false;
      // Detect header: starts with "record_type" or "#"
      if (!line.empty() && (line[0] == '#' || line.rfind("record_type", 0) == 0)) {
        gzwrite(gz_pass, line.data(), static_cast<unsigned>(line.size()));
        gzwrite(gz_pass, "\n", 1);
        gzwrite(gz_som, line.data(), static_cast<unsigned>(line.size()));
        gzwrite(gz_som, "\n", 1);
        continue;
      }
      // Not a header — fall through to normal processing.
    }

    // Get contig_name (col 1).
    auto cname_sv = getField(line.data(), static_cast<uint32_t>(line.size()),
                             kR2cCol_ContigName);
    if (cname_sv.empty())
      continue;

    // Check cache — avoid hash lookup when consecutive rows share a cname
    // (which they almost always do: one contig row + N read rows).
    if (cname_sv != last_cname) {
      last_cname.assign(cname_sv);
      last_in_pass = pass_cnames.count(last_cname) > 0;
      last_in_som  = last_in_pass && som_cnames.count(last_cname) > 0;
    }

    if (last_in_pass) {
      gzwrite(gz_pass, line.data(), static_cast<unsigned>(line.size()));
      gzwrite(gz_pass, "\n", 1);
      ++n_pass;

      if (last_in_som) {
        gzwrite(gz_som, line.data(), static_cast<unsigned>(line.size()));
        gzwrite(gz_som, "\n", 1);
        ++n_som;
      }
    }
  }

  gzclose(gz_in);
  gzclose(gz_pass);
  gzclose(gz_som);

  const double elapsed = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - t0).count();
  logLine("[r2c] done in ", std::fixed, std::setprecision(1), elapsed, "s: ",
          n_in, " rows in, ", n_pass, " pass, ", n_som, " somatic");

  return n_in;
}

// ---------- @PG stamping ----------
//
// Every final BAM produced by `svaba postprocess` gets a new @PG line stamped
// onto its header so downstream tools (and humans) can tell at a glance that
// postprocess ran and with what arguments. The line chains onto the existing
// @PG chain via PP: so samtools' `--pg` resolution keeps working.
//
// Conventions:
//   ID: "svaba_postprocess" (or ".1", ".2", ... if already present — IDs must
//       be unique within a header)
//   PN: "svaba"
//   VN: SVABA_VERSION
//   CL: the original argv reassembled, with "svaba " prepended so the CL
//       matches what the user actually typed
//   PP: the tail of the existing PG chain (the ID that is nobody's PP pointer)

struct PgChainInfo {
  std::vector<std::string> ids;   // all @PG IDs, in file order
  std::string              tail;  // ID that is nobody else's PP (or "" if none)
};

PgChainInfo scanPgChain(const std::string& hdr_text) {
  PgChainInfo info;
  std::set<std::string> is_pp_target;
  std::istringstream iss(hdr_text);
  std::string line;
  while (std::getline(iss, line)) {
    if (line.rfind("@PG\t", 0) != 0) continue;
    std::string id, pp;
    std::istringstream ls(line);
    std::string field;
    while (std::getline(ls, field, '\t')) {
      if      (field.rfind("ID:", 0) == 0) id = field.substr(3);
      else if (field.rfind("PP:", 0) == 0) pp = field.substr(3);
    }
    if (!id.empty()) info.ids.push_back(id);
    if (!pp.empty()) is_pp_target.insert(pp);
  }
  // Tail = an ID that no one else's PP points at. If multiple, take the last
  // one in file order — matches what samtools reheader/markdup do in practice.
  for (const auto& id : info.ids) {
    if (is_pp_target.find(id) == is_pp_target.end()) info.tail = id;
  }
  return info;
}

// Idempotency check: has *any* svaba_postprocess @PG line (including
// uniquified variants svaba_postprocess.1, .2, ...) already been stamped
// onto this header? Used to auto-skip the dedup/reheader phase on
// re-runs so `scripts/svaba_postprocess.sh` is safely rerunnable.
//
// Matches both the bare "svaba_postprocess" ID emitted on the first run
// and the "svaba_postprocess.<n>" variants uniquifyId() would produce if
// the user deliberately runs postprocess twice in a row. Either form
// means "dedup+PG-stamp has already happened at least once."
bool hasSvabaPostprocessPg(const std::string& hdr_text) {
  const PgChainInfo info = scanPgChain(hdr_text);
  for (const auto& id : info.ids) {
    if (id == "svaba_postprocess") return true;
    if (id.rfind("svaba_postprocess.", 0) == 0) return true;
  }
  return false;
}

std::string uniquifyId(const std::string& base,
                       const std::vector<std::string>& existing) {
  auto has = [&](const std::string& s) {
    return std::find(existing.begin(), existing.end(), s) != existing.end();
  };
  if (!has(base)) return base;
  for (int i = 1; i < 10000; ++i) {
    std::string cand = base + "." + std::to_string(i);
    if (!has(cand)) return cand;
  }
  return base + ".new";  // extremely unreachable fallback
}

// Strip tab/newline/CR from a value so it can't break @PG line structure.
std::string sanitizeHeaderValue(const std::string& s) {
  std::string out;
  out.reserve(s.size());
  for (char c : s) {
    if (c == '\t') out.push_back(' ');
    else if (c == '\n' || c == '\r') continue;
    else out.push_back(c);
  }
  return out;
}

// Append a single @PG line to hdr_text. Does not otherwise touch the header.
std::string appendSvabaPostprocessPg(const std::string& hdr_text,
                                     const std::string& cl) {
  const PgChainInfo info = scanPgChain(hdr_text);
  const std::string id   = uniquifyId("svaba_postprocess", info.ids);
  const std::string vn   = SVABA_VERSION;
  const std::string cl_s = sanitizeHeaderValue(cl);

  std::string line = "@PG\tID:" + id + "\tPN:svaba\tVN:" + vn + "\tCL:" + cl_s;
  if (!info.tail.empty()) line += "\tPP:" + info.tail;
  line.push_back('\n');

  std::string out = hdr_text;
  if (!out.empty() && out.back() != '\n') out.push_back('\n');
  out += line;
  return out;
}

// Build a new SeqLib::BamHeader with our @PG line appended onto `src`.
SeqLib::BamHeader stampedHeader(const SeqLib::BamHeader& src,
                                const std::string& cl) {
  return SeqLib::BamHeader(appendSvabaPostprocessPg(src.AsString(), cl));
}

// Read the header of `bam` (without walking records) and return it as a
// SeqLib::BamHeader. The BamReader is closed on scope exit.
SeqLib::BamHeader readHeaderOnly(const std::string& bam) {
  SeqLib::BamReader r;
  { SuppressStderr quiet;
    if (!r.Open(bam))
      throw std::runtime_error("reheader: cannot open " + bam);
  }
  return r.Header();
}

// Return true iff the BAM's header declares it coordinate-sorted via an
// @HD line with SO:coordinate. Samtools, SeqLib, and bwa all emit that
// form when producing a sorted BAM, so a simple substring check on the
// first @HD line is sufficient.
//
// Used by processSuffix() to short-circuit the sort step when it would be
// a no-op — lets the user rerun postprocess on already-sorted inputs
// without paying the sort cost again (e.g. after a --skip-dedup rerun).
//
// Errors (can't open, no @HD at all, @HD without SO:) conservatively
// return false so we fall through to running sort, which will produce a
// correct result regardless.
bool isCoordinateSorted(const std::string& bam) {
  try {
    const SeqLib::BamHeader h = readHeaderOnly(bam);
    const std::string text = h.AsString();
    std::istringstream iss(text);
    std::string line;
    while (std::getline(iss, line)) {
      if (line.rfind("@HD\t", 0) == 0) {
        // @HD line found; check for the SO:coordinate tag. Tags are
        // tab-separated, any order — a plain substring search on the
        // line is accurate enough because "SO:coordinate" can't appear
        // as a substring of any other valid @HD tag value.
        return line.find("SO:coordinate") != std::string::npos;
      }
      // @HD, if present, is required to be the first header line. If
      // the first non-empty line isn't @HD we can stop looking.
      if (!line.empty() && line[0] == '@') break;
    }
  } catch (...) {
    // header unreadable — let the sort step run and produce its own
    // error with better context.
  }
  return false;
}

// Use `samtools reheader` to swap in a new header on an existing BAM. Fast
// because samtools reheader streams the BGZF body as opaque blocks. We use
// this for the suffixes whose pipeline doesn't already rewrite the file
// (i.e. `contigs`, or any suffix when --sort-only is passed).
int reheaderBamWithPg(const std::string& bam,
                      const std::string& cl,
                      const std::string& tag,
                      int verbose) {
  const std::string tmp_hdr = bam + ".postprocess.hdr.tmp.sam";
  const std::string tmp_bam = bam + ".postprocess.reheader.tmp.bam";

  try {
    const SeqLib::BamHeader src = readHeaderOnly(bam);
    const std::string new_text = appendSvabaPostprocessPg(src.AsString(), cl);

    {
      std::ofstream ofs(tmp_hdr, std::ios::binary);
      if (!ofs)
        throw std::runtime_error("cannot write temp header " + tmp_hdr);
      ofs << new_text;
      if (!ofs)
        throw std::runtime_error("short write on " + tmp_hdr);
    }

    const std::string cmd = "samtools reheader " + shQuote(tmp_hdr) + " " +
                            shQuote(bam) + " > " + shQuote(tmp_bam);
    if (verbose >= 2) logLine("[", tag, "] exec: ", cmd);

    const auto t0 = std::chrono::steady_clock::now();
    const int rc = std::system(cmd.c_str());
    const double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    ::unlink(tmp_hdr.c_str());

    if (rc != 0) {
      ::unlink(tmp_bam.c_str());
      logLine("[", tag, "] samtools reheader FAILED (exit ", rc, ") after ",
              std::fixed, std::setprecision(1), elapsed, "s");
      return rc;
    }

    renameOrThrow(tmp_bam, bam);
    if (verbose >= 1)
      logLine("[", tag, "] PG stamped via reheader in ",
              std::fixed, std::setprecision(1), elapsed, "s");
    return 0;
  } catch (const std::exception& e) {
    ::unlink(tmp_hdr.c_str());
    ::unlink(tmp_bam.c_str());
    logLine("[", tag, "] reheader error: ", e.what());
    return -1;
  }
}

// ---------- BAI indexing ----------

// Build a .bai alongside `bam`. Uses htslib directly rather than shelling out
// to `samtools index` — one less fork, and we're already linked against
// htslib through SeqLib. `sam_index_build(fn, 0)` emits BAI; min_shift > 0
// would emit CSI (which svaba's consumers don't expect).
int indexBam(const std::string& bam, const std::string& tag, int verbose) {
  const auto t0 = std::chrono::steady_clock::now();
  const int rc = sam_index_build(bam.c_str(), 0);
  const double elapsed =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  if (rc < 0) {
    logLine("[", tag, "] indexing FAILED (sam_index_build rc=", rc, ") after ",
            std::fixed, std::setprecision(2), elapsed, "s");
    return rc;
  }
  if (verbose >= 1)
    logLine("[", tag, "] indexed (.bai) in ",
            std::fixed, std::setprecision(2), elapsed, "s");
  return 0;
}

// ---------- sort step ----------

// Run samtools sort on in_bam producing out_bam. Returns process exit code.
// We don't parse samtools' progress output; samtools -v isn't universally
// respected and its stderr format changes. Caller emits start/end log lines
// with timing.
int runSort(const std::string& in_bam,
            const std::string& out_bam,
            int sort_threads,
            const std::string& mem,
            const std::string& tag,
            int verbose) {
  std::string cmd =
      "samtools sort -@ " + std::to_string(std::max(1, sort_threads)) +
      " -m " + shQuote(mem) +
      " -o " + shQuote(out_bam) +
      " "    + shQuote(in_bam);

  if (verbose >= 2)
    logLine("[", tag, "] exec: ", cmd);

  const auto t0 = std::chrono::steady_clock::now();
  const int rc = std::system(cmd.c_str());
  const double elapsed =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  if (rc != 0) {
    logLine("[", tag, "] samtools sort FAILED (exit ", rc, ") after ",
            std::fixed, std::setprecision(1), elapsed, "s");
    return rc;
  }
  if (verbose >= 1)
    logLine("[", tag, "] sort done in ",
            std::fixed, std::setprecision(1), elapsed, "s");
  return 0;
}

// ---------- dedup step ----------

struct DedupStats {
  std::size_t in_records = 0;  // records read
  std::size_t kept       = 0;  // unique (qname, flag) per locus, written out
  std::size_t merged     = 0;  // duplicates folded into an earlier record's tags
  std::size_t tag_edits  = 0;  // number of bi/bz merges that actually changed tag
  double      seconds    = 0;
};

// Boundary-aware comma-token union: adds any tokens in `incoming` that aren't
// already in `cur` as comma-separated atoms. Mirrors stamp_tag() in
// SvabaOutputWriter.cpp so the postprocess dedup produces exactly the same
// bi:Z / bz:Z strings svaba run would have written if the two contributing
// windows had been merged in-memory.
//
// Returns true iff cur was modified.
bool mergeCommaTokens(std::string& cur, std::string_view incoming) {
  if (incoming.empty()) return false;
  if (cur.empty()) { cur.assign(incoming); return true; }

  bool changed = false;
  std::size_t s = 0;
  while (s <= incoming.size()) {
    std::size_t e = incoming.find(',', s);
    if (e == std::string_view::npos) e = incoming.size();
    std::string_view tok = incoming.substr(s, e - s);
    if (!tok.empty()) {
      bool present = false;
      std::size_t pos = 0;
      while ((pos = cur.find(tok, pos)) != std::string::npos) {
        const bool left_ok  = (pos == 0) || cur[pos - 1] == ',';
        const bool right_ok = (pos + tok.size() == cur.size()) ||
                              cur[pos + tok.size()] == ',';
        if (left_ok && right_ok) { present = true; break; }
        pos += tok.size();
      }
      if (!present) {
        cur.push_back(',');
        cur.append(tok);
        changed = true;
      }
    }
    if (e == incoming.size()) break;
    s = e + 1;
  }
  return changed;
}

// Pull `tag` off `incoming` and union it into `existing`'s `tag`. Returns
// true iff the existing record's tag string actually changed.
bool mergeZTagInto(SeqLib::BamRecord& existing,
                   const SeqLib::BamRecord& incoming,
                   const char* tag) {
  std::string incoming_val;
  if (!incoming.GetZTag(tag, incoming_val) || incoming_val.empty())
    return false;

  std::string cur;
  existing.GetZTag(tag, cur);
  if (!mergeCommaTokens(cur, incoming_val))
    return false;

  existing.RemoveTag(tag);
  existing.AddZTag(tag, cur);
  return true;
}

// Streaming dedup + tag-merge. Reads in_bam, writes out_bam, collapses exact
// (qname, flag) duplicates at the same (chr, pos). Duplicates are not simply
// discarded: their bi:Z and bz:Z tags (the comma-joined alt-supporting and
// all-supporting contig lists) are unioned into the first record's tags
// using the same boundary-aware merge that stamp_tag() in
// SvabaOutputWriter.cpp uses during `svaba run`. Since a duplicate record
// typically comes from an overlapping assembly window, its bi/bz may name
// a different contig than the first record — dropping it would lose that
// support evidence. Requires input to be coord-sorted.
//
// Memory: O(reads at a single locus). The per-locus index and record buffer
// are cleared whenever (chr, pos) advances, so even pathological pileups
// only blow up transiently.
DedupStats streamDedup(const std::string& in_bam,
                       const std::string& out_bam,
                       const SeqLib::BamHeader& out_hdr,
                       const std::string& tag,
                       int threads,
                       int verbose) {
  using clock = std::chrono::steady_clock;
  DedupStats st;
  const auto t0 = clock::now();
  auto last = t0;
  const off_t in_size = fileSize(in_bam);

  SeqLib::BamReader r;
  { SuppressStderr quiet;
    if (!r.Open(in_bam))
      throw std::runtime_error("dedup: cannot open input " + in_bam);
  }
  // Enable BGZF decompression thread pool on the reader. Must happen
  // after Open() (where fp_ is populated) and before the first Next().
  // Typical 3–5x end-to-end speedup at threads=4..8 on dense BAMs,
  // because without it the main thread does every BGZF block inflate
  // sequentially.
  const int io_threads = std::max(1, threads);
  r.SetThreads(io_threads);

  SeqLib::BamWriter w;
  if (!w.Open(out_bam))
    throw std::runtime_error("dedup: cannot open output " + out_bam);
  // Matching BGZF compression thread pool on the write side. The
  // output is the expensive half of this function — every kept record
  // goes through deflate — so pooling here matters most.
  w.SetThreads(io_threads);
  // Use the caller-supplied header (which already carries our @PG stamp)
  // rather than blindly mirroring the reader's — that's how postprocess
  // marks its output without a separate reheader pass for dedup suffixes.
  w.SetHeader(out_hdr);
  if (!w.WriteHeader())
    throw std::runtime_error("dedup: cannot write header " + out_bam);

  // Per-locus buffer. We must hold records in insertion order so flush
  // preserves the original coord-sort: same (chr, pos) records keep whatever
  // secondary ordering samtools sort produced (qname-stable tiebreak).
  //
  //   keep_buf[i] holds the record we're emitting for the i'th unique key.
  //   idx_by_key[key] is that buffer index, used to fold duplicate tags into
  //     the already-buffered record.
  std::vector<SeqLib::BamRecord> keep_buf;
  std::unordered_map<std::string, std::size_t> idx_by_key;
  keep_buf.reserve(64);
  idx_by_key.reserve(64);

  int32_t cur_chr = -2;  // -2 is a sentinel that won't match any valid ChrID
  int32_t cur_pos = -2;
  // Cached display name for cur_chr. ChrName() hits the header every call, so
  // we resolve it only on (rare) chromosome change, not per record.
  std::string cur_chr_name = "*";

  // When bucket_count() drifts beyond this threshold (after a pileup locus
  // inflates the map), flushLocus() swaps in a fresh small map instead of
  // calling clear(). See comment below for why this matters.
  constexpr std::size_t kBucketResetThreshold = 256;

  // Flush buffered records for the current locus in insertion order.
  //
  // Perf note: std::unordered_map::clear() is O(bucket_count), NOT O(size)
  // — it walks the entire bucket array to zero out each slot, even when
  // the map currently has 0 elements. The bucket array grows with the
  // max number of entries ever held and NEVER shrinks. So after one
  // pileup locus (centromere, simple repeat, HLA, etc.) bloats buckets
  // to ~130k, every subsequent locus transition — hundreds of millions
  // of them — pays a full memset of that inflated array. A perf profile
  // showed 95% of main-thread CPU going to this single memset before
  // this fix landed, explaining why the BGZF thread pools were starving.
  //
  // Fix: when the bucket count has grown past kBucketResetThreshold,
  // swap the inflated map with a fresh small one. The inflated map's
  // destructor pays the O(bucket_count) cost ONCE (on pileup exit)
  // instead of once per locus transition. Subsequent small loci see a
  // small bucket array and their clear() is fast.
  auto flushLocus = [&]() {
    for (auto& kept : keep_buf) {
      if (!w.WriteRecord(kept))
        throw std::runtime_error("dedup: WriteRecord failed on " + out_bam);
    }
    keep_buf.clear();
    if (idx_by_key.bucket_count() > kBucketResetThreshold) {
      std::unordered_map<std::string, std::size_t> fresh;
      fresh.reserve(64);
      idx_by_key.swap(fresh);
      // `fresh` (now holding the inflated map) destructs at end of this
      // lambda — that's where we pay the one-time big memset.
    } else {
      idx_by_key.clear();
    }
  };

  while (auto opt = r.Next()) {
    SeqLib::BamRecord& rec = *opt;

    const int32_t chr = rec.ChrID();
    const int32_t pos = rec.Position();
    if (chr != cur_chr || pos != cur_pos) {
      flushLocus();
      if (chr != cur_chr) {
        // Resolve name once per chromosome. Unmapped (chr == -1) stays "*".
        cur_chr_name = (chr >= 0) ? rec.ChrName(r.Header()) : std::string("*");
      }
      cur_chr = chr;
      cur_pos = pos;
    }

    // Build the dedup key: qname + '\t' + decimal flag. Decimal keeps keys
    // ASCII / comparable across platforms; '\t' as the separator can never
    // appear inside a BAM qname.
    std::string key = rec.Qname();
    key.push_back('\t');
    {
      char buf[16];
      const int n = std::snprintf(buf, sizeof(buf), "%u", rec.AlignmentFlag());
      key.append(buf, n);
    }

    auto [it, inserted] = idx_by_key.try_emplace(std::move(key), keep_buf.size());
    if (inserted) {
      // BamRecord is move-only (shared_ptr<bam1_t> under the hood); move
      // it into the buffer so we transfer the htslib pointer without
      // ref-bumping.
      keep_buf.emplace_back(std::move(rec));
      ++st.kept;
    } else {
      // Duplicate at same (chr, pos). Fold its bi/bz tags into the kept copy
      // so we don't drop alt/contig support from the overlapping window. We
      // read from rec here (still alive — only moved in the other branch).
      SeqLib::BamRecord& existing = keep_buf[it->second];
      if (mergeZTagInto(existing, rec, "bi")) ++st.tag_edits;
      if (mergeZTagInto(existing, rec, "bz")) ++st.tag_edits;
      ++st.merged;
    }
    ++st.in_records;

    // Progress print every 25M reads. Was 1M originally, then 5M; post
    // BGZF-threading + bucket-clear fix the dedup runs fast enough that
    // 5M was firing too often on deep BAMs. 25M gives ~a screenful of
    // lines for a typical 100–500M-read WGS BAM.
    constexpr std::size_t PROGRESS_EVERY = 25'000'000ULL;
    if (verbose >= 1 && (st.in_records % PROGRESS_EVERY) == 0) {
      const auto now = clock::now();
      const double elapsed = std::chrono::duration<double>(now - t0).count();
      const double dt      = std::chrono::duration<double>(now - last).count();
      const double inst_rps = dt > 0 ? static_cast<double>(PROGRESS_EVERY) / dt : 0.0;
      const double avg_rps  = elapsed > 0 ? st.in_records / elapsed : 0.0;
      last = now;

      // Column widths chosen so the line never shifts as counters grow:
      //   chr:    5 chars right-padded   ("chr22" max)
      //   :
      //   pos:    left-justified, with thousands separators
      //             (flush against colon; human chr1 max "248,956,422" = 11 chars)
      //   locus total padded to 17 so subsequent columns stay anchored.
      //   counts: 11 chars right-aligned with thousands separators
      //             (fits "999,999,999" = 9 figures + 2 commas)
      //   rps:    6 chars "%6.2f"        (up to 999.99 M/s)
      //   elapsed 8 chars "%8.1f"        (up to 99999.9 s ≈ 27 hr)
      //   MiB:    7 chars                (up to 9,999,999 MiB ≈ 10 TB)
      auto fmtN = [](std::size_t n) {
        std::ostringstream oss;
        oss << std::setw(11) << std::right << SeqLib::AddCommas<std::size_t>(n);
        return oss.str();
      };

      // Locus = right-padded chr + ':' + left-justified comma'd pos.
      // Flush to 17 chars total so " | " lands at a fixed column regardless
      // of chr length or position magnitude.
      std::ostringstream locus_oss;
      locus_oss << std::setw(5) << std::right << cur_chr_name
                << ":"
                << SeqLib::AddCommas<int32_t>(cur_pos >= 0 ? cur_pos + 1 : 0);
      std::string locus_str = locus_oss.str();
      if (locus_str.size() < 17) locus_str.append(17 - locus_str.size(), ' ');

      // Input-size hint for a coarse "how far into the job" feel. A true
      // ETA would need bytes-read from htslib, which SeqLib doesn't expose.
      std::string size_hint;
      if (in_size > 0) {
        std::ostringstream oss;
        oss << " (input " << std::setw(7) << std::right
            << (in_size / (1024 * 1024)) << " MiB)";
        size_hint = oss.str();
      }

      logLine("[", tag, "] dedup: at ", locus_str, " | ",
              fmtN(st.in_records), " reads | ",
              fmtN(st.kept),       " kept, ",
              fmtN(st.merged),     " merged (",
              fmtN(st.tag_edits),  " tag-edits) | ",
              std::fixed, std::setprecision(2),
              std::setw(6), inst_rps / 1e6, "M/s last-25M, ",
              std::setw(6), avg_rps  / 1e6, "M/s avg, ",
              std::setprecision(1), std::setw(8), elapsed, "s elapsed",
              size_hint);
    }
  }

  // Final locus flush.
  flushLocus();

  if (!w.Close())
    throw std::runtime_error("dedup: writer Close failed on " + out_bam);

  st.seconds = std::chrono::duration<double>(clock::now() - t0).count();
  return st;
}

// ---------- per-suffix pipeline ----------

// Move one file to another atomically within the same filesystem, via
// std::rename. Throws on failure so the caller can decide what to do.
void renameOrThrow(const std::string& from, const std::string& to) {
  if (std::rename(from.c_str(), to.c_str()) != 0) {
    throw std::runtime_error("rename " + from + " -> " + to +
                             " failed: " + std::strerror(errno));
  }
}

bool isDedupSuffix(const std::string& s) {
  for (const auto& d : kDedupSuffixes)
    if (d == s) return true;
  return false;
}

// Run the full sort + dedup + reheader + index pipeline for a single suffix.
// Safe to invoke from its own thread. Never throws out — caller joins; any
// failure is logged and the pipeline aborts for this suffix only.
//
// End state on success: the file at ${id}.${suffix}.bam is the final output,
// coord-sorted, dedup-collapsed (for dedup-eligible suffixes), @PG-stamped,
// and accompanied by a ${id}.${suffix}.bam.bai index. All intermediate
// filenames use a .postprocess.*.tmp.bam suffix so anything left behind on
// failure is obviously transient and trivially cleanable.
// do_finalize controls whether the dedup + reheader + index phase runs.
// Phase-splitting the driver lets sort run in parallel across suffixes
// (narrow thread budget each) and dedup run serially with the full thread
// budget, which is what the user wants since htslib BGZF pools are
// per-file — sharing them across concurrent workers oversubscribes.
void processSuffix(const std::string& id,
                   const std::string& suffix,
                   int  per_job_threads,
                   const std::string& mem,
                   bool do_sort,
                   bool do_dedup,
                   bool do_finalize,
                   const std::string& cl,
                   int  verbose) {
  const std::string bam       = id + "." + suffix + ".bam";
  // Intermediate temp files. Named so there's no ambiguity about the step
  // that produced them — anything matching .postprocess.*.tmp.bam is safe
  // to remove after the run is done.
  const std::string sort_tmp  = id + "." + suffix + ".postprocess.sort.tmp.bam";
  const std::string dedup_tmp = id + "." + suffix + ".postprocess.dedup.tmp.bam";

  auto cleanup_tmps = [&]() {
    ::unlink(sort_tmp.c_str());
    ::unlink(dedup_tmp.c_str());
  };

  if (!fileExists(bam)) {
    logLine("[", suffix, "] skipping: ", bam, " not found");
    return;
  }

  bool pg_stamped = false;

  try {
    // --- Sort -----------------------------------------------------------
    // Skip the sort step if the BAM already declares itself coordinate-
    // sorted. Cheap O(header) check; turns a repeat postprocess run on
    // an already-postprocessed BAM into an effectively instant no-op
    // for this suffix. Lets users pass --skip-dedup at the shell layer
    // (or --sort-only from the C++ CLI) and rerun without paying 30+s
    // per suffix for a redundant sort.
    if (do_sort && isCoordinateSorted(bam)) {
      logLine("[", suffix, "] already coordinate-sorted (per @HD SO:coordinate); "
              "skipping sort");
    } else if (do_sort) {
      const int rc = runSort(bam, sort_tmp, per_job_threads, mem, suffix, verbose);
      if (rc != 0)
        throw std::runtime_error("samtools sort exited " + std::to_string(rc));
      renameOrThrow(sort_tmp, bam);
    }

    if (!do_finalize) {
      // Caller wants sort-only (the parallel Phase 1 of the driver); leave
      // dedup + reheader + index to the caller's Phase 2 invocation. This
      // split exists so Phase 2 can run serially with the full thread
      // budget dedicated to BGZF read/write pooling per BAM, instead of
      // fragmenting threads across concurrent workers.
      return;
    }

    // Read the current header ONCE and reuse it for both the PG-idempotency
    // check and (if we do run dedup) the stamped output header. Avoids a
    // second BAM open just to look up PG state.
    const SeqLib::BamHeader existing_hdr = readHeaderOnly(bam);
    const std::string       existing_txt = existing_hdr.AsString();
    const bool already_postprocessed     = hasSvabaPostprocessPg(existing_txt);

    // --- Dedup (also stamps PG on the way) ------------------------------
    if (do_dedup && isDedupSuffix(suffix) && !already_postprocessed) {
      if (verbose >= 1)
        logLine("[", suffix, "] starting dedup on ", bam);

      // Append our @PG line to the existing header text, then hand the
      // stamped header to the writer. This is free — dedup is already
      // rewriting the file, so we avoid a separate reheader pass.
      const SeqLib::BamHeader stamped =
          SeqLib::BamHeader(appendSvabaPostprocessPg(existing_txt, cl));

      const DedupStats ds =
          streamDedup(bam, dedup_tmp, stamped, suffix,
                      per_job_threads,   // BGZF read/write pool size
                      verbose);
      renameOrThrow(dedup_tmp, bam);
      pg_stamped = true;

      logLine("[", suffix, "] dedup complete: ",
              ds.in_records, " in, ",
              ds.kept, " kept, ",
              ds.merged, " merged (",
              (ds.in_records ?
                100.0 * ds.merged / static_cast<double>(ds.in_records) : 0.0),
              "% dup, ",
              ds.tag_edits, " bi/bz tag-edits) in ",
              std::fixed, std::setprecision(1), ds.seconds, "s");
    } else if (do_dedup && isDedupSuffix(suffix) && already_postprocessed) {
      // BAM has already been through svaba_postprocess. Skip the expensive
      // dedup rewrite AND the reheader below — PG is already stamped.
      logLine("[", suffix, "] already has svaba_postprocess @PG in header; "
              "skipping dedup (rerun-safe no-op)");
      pg_stamped = true;
    } else if (do_dedup && !isDedupSuffix(suffix)) {
      if (verbose >= 2)
        logLine("[", suffix, "] dedup skipped (not a dedup suffix)");
    }

    // --- Reheader to stamp @PG for any suffix that didn't go through -----
    // --- dedup (e.g. `contigs`, or anything under --sort-only). Also      -----
    // --- skipped when `already_postprocessed` set `pg_stamped` above      -----
    // --- (nothing to add to an already-PG-stamped header).                -----
    if (!pg_stamped && !already_postprocessed) {
      const int rc = reheaderBamWithPg(bam, cl, suffix, verbose);
      if (rc != 0)
        throw std::runtime_error("reheader failed with rc " + std::to_string(rc));
    }

    // --- Index ----------------------------------------------------------
    // Index last so the .bai matches the BGZF offsets of the final file.
    // Any earlier index would be invalidated by the dedup/reheader rewrite.
    // Always runs — even on a rerun, the .bai may be missing or stale
    // (e.g. if the user renamed the BAM after a previous run). Cheap op.
    (void) indexBam(bam, suffix, verbose);

  } catch (const std::exception& e) {
    logLine("[", suffix, "] ERROR: ", e.what());
    cleanup_tmps();
  }
}

// ==========================================================================
// Phase 5: split-by-source — demux each dedup-eligible BAM by the first 4
// characters of QNAME into per-prefix BAMs, then sort + index each.
// ==========================================================================

void splitBySource(const std::string& id, const std::string& suffix,
                   int threads, const std::string& mem, int verbose) {
  const std::string input_bam = id + "." + suffix + ".bam";
  if (!fileExists(input_bam)) {
    if (verbose >= 1)
      std::cerr << "  [split/" << suffix << "] no input BAM, skipping"
                << std::endl;
    return;
  }

  if (verbose >= 1)
    std::cerr << "  [split/" << suffix << "] reading " << input_bam
              << std::endl;

  // Open input (suppress missing-index warning — we only need sequential access)
  SeqLib::BamReader reader;
  { SuppressStderr quiet;
    if (!reader.Open(input_bam)) {
      std::cerr << "  [split/" << suffix << "] ERROR: cannot open " << input_bam
                << std::endl;
      return;
    }
  }
  const SeqLib::BamHeader& hdr = reader.Header();

  // Single-pass approach: open writers on-demand as new prefixes are
  // encountered. Each prefix gets its own unsorted tmp BAM.
  std::unordered_map<std::string, std::unique_ptr<SeqLib::BamWriter>> writers;
  std::unordered_map<std::string, std::string> tmp_paths;  // prefix -> tmp path

  size_t n_records = 0;
  while (auto rec = reader.Next()) {
    std::string qname = rec->Qname();
    std::string prefix = qname.substr(0, std::min<size_t>(4, qname.size()));

    auto it = writers.find(prefix);
    if (it == writers.end()) {
      // Open a new writer for this prefix
      std::string out_path = id + "." + suffix + "." + prefix
                             + ".postprocess.split.tmp.bam";
      auto w = std::make_unique<SeqLib::BamWriter>();
      w->Open(out_path);
      w->SetHeader(hdr);
      w->WriteHeader();
      writers[prefix] = std::move(w);
      tmp_paths[prefix] = out_path;
    }
    writers[prefix]->WriteRecord(*rec);
    ++n_records;

    if (verbose >= 2 && (n_records % 5000000) == 0)
      std::cerr << "  [split/" << suffix << "] " << n_records
                << " records routed to " << writers.size()
                << " prefixes" << std::endl;
  }
  reader.Close();

  if (verbose >= 1)
    std::cerr << "  [split/" << suffix << "] " << n_records
              << " records -> " << writers.size() << " prefixes" << std::endl;

  // Close all writers
  for (auto& [pfx, w] : writers)
    w->Close();
  writers.clear();

  // Sort + index each prefix BAM
  int n_sorted = 0;
  for (auto& [prefix, tmp_path] : tmp_paths) {
    const std::string final_path = id + "." + suffix + "." + prefix + ".bam";

    // Sort via samtools
    std::ostringstream cmd;
    cmd << "samtools sort"
        << " -@ " << threads
        << " -m " << mem
        << " -o " << final_path
        << " " << tmp_path;

    if (verbose >= 2)
      std::cerr << "  [split/" << suffix << "/" << prefix << "] "
                << cmd.str() << std::endl;

    int rc = std::system(cmd.str().c_str());
    if (rc != 0) {
      std::cerr << "  [split/" << suffix << "/" << prefix
                << "] ERROR: samtools sort failed (rc=" << rc << ")"
                << std::endl;
      // Clean up tmp and skip indexing
      std::remove(tmp_path.c_str());
      continue;
    }

    // Remove tmp
    std::remove(tmp_path.c_str());

    // Index
    int idx_rc = sam_index_build(final_path.c_str(), 0);
    if (idx_rc != 0) {
      std::cerr << "  [split/" << suffix << "/" << prefix
                << "] WARNING: index build failed for " << final_path
                << std::endl;
    }

    ++n_sorted;
    if (verbose >= 2)
      std::cerr << "  [split/" << suffix << "/" << prefix << "] done -> "
                << final_path << std::endl;
  }

  if (verbose >= 1)
    std::cerr << "  [split/" << suffix << "] sorted+indexed "
              << n_sorted << "/" << tmp_paths.size() << " prefix BAMs"
              << std::endl;
}

}  // namespace

void runPostprocess(int argc, char** argv) {
  // Reconstruct the invocation string for the @PG CL: tag BEFORE parseOpts
  // consumes argv via getopt_long. We prepend "svaba" so the CL in the
  // final BAM headers reads as the user actually typed it (the dispatch in
  // svaba.cpp strips the "svaba" token before calling us).
  std::string cl = "svaba";
  for (int i = 0; i < argc; ++i) { cl.push_back(' '); cl += argv[i]; }

  const Opts o = parseOpts(argc, argv);

  // --- Phase 0: merge per-thread outputs ---------------------------------
  // Must run BEFORE active-suffix detection, since merge creates the single
  // ${ID}.${suffix}.bam files that the subsequent phases operate on.
  // Auto-detects per-thread files by glob; no-ops if none found.
  runMergePhase(o.id, o.threads, o.verbose);

  // Pick active suffixes — only the ones whose BAM actually exists. No point
  // reserving a thread budget slice for a missing file.
  std::vector<std::string> active;
  for (const auto& s : kAllSuffixes) {
    const std::string bam = o.id + "." + s + ".bam";
    if (fileExists(bam)) active.push_back(s);
  }

  if (active.empty()) {
    std::cerr << "No BAMs found for ID '" << o.id << "' (checked suffixes: ";
    for (size_t i = 0; i < kAllSuffixes.size(); ++i) {
      if (i) std::cerr << ", ";
      std::cerr << kAllSuffixes[i];
    }
    std::cerr << ")" << std::endl;
    return;
  }

  // Split the thread budget across concurrent suffix jobs, matching the
  // shell script's behavior: floor(THREADS / n_active), min 1. If threads <
  // n_active, cap concurrency to threads so we don't oversubscribe.
  const int n_active = static_cast<int>(active.size());
  const int max_parallel   = std::min(o.threads, n_active);
  const int per_job_sort_t = std::max(1, o.threads / n_active);

  if (o.verbose >= 1) {
    std::ostringstream act;
    for (size_t i = 0; i < active.size(); ++i) {
      if (i) act << " ";
      act << active[i];
    }
    std::cerr << "svaba postprocess: id=" << o.id
              << " threads=" << o.threads
              << " mem=" << o.mem
              << " parallel=" << max_parallel
              << " per-job-sort-threads=" << per_job_sort_t
              << " sort=" << (o.dedup_only ? "off" : "on")
              << " dedup=" << (o.sort_only  ? "off" : "on")
              << "\n  active: " << act.str() << std::endl;
  }

  // Two-phase execution:
  //
  //   Phase 1 (PARALLEL): run samtools sort across all active suffixes
  //     concurrently. Sort is disk+CPU bound and perfectly parallelizable
  //     across files, so each worker gets `per_job_sort_t = o.threads /
  //     n_active` threads — total thread usage ≈ o.threads.
  //
  //   Phase 2 (SERIAL): dedup + reheader + index, one suffix at a time,
  //     each with the FULL `o.threads` budget for its BGZF read/write
  //     pool. Running the dedup phase in parallel would fragment the
  //     thread budget (4 suffixes × 2 threads/pool = 2× oversubscription
  //     once you count read+write pools per suffix), yielding worse wall
  //     time than serial-with-wide-pools because BGZF parallelism has
  //     diminishing returns — 8 threads get ~5x vs single-thread, but 4
  //     parallel workers with 2 threads each get only 4 × 1.8x / 4 ≈ 1.8x
  //     effective per-BAM.
  //
  // Each phase is itself idempotent: Phase 1 skips sort when the BAM is
  // already @HD SO:coordinate; Phase 2 skips dedup when the header
  // already carries a `svaba_postprocess` @PG line. Rerunning a
  // completed postprocess is essentially instant.

  // --- Phase 1: parallel sort ----------------------------------------
  if (!o.dedup_only) {
    if (o.verbose >= 1)
      std::cerr << "svaba postprocess: phase 1 — parallel sort across "
                << active.size() << " suffixes ("
                << per_job_sort_t << " threads each)" << std::endl;
    std::vector<std::thread> workers;
    workers.reserve(active.size());
    for (const auto& suffix : active) {
      workers.emplace_back(
          processSuffix,
          o.id,
          suffix,
          per_job_sort_t,
          o.mem,
          /*do_sort=*/    true,
          /*do_dedup=*/   false,  // deferred to phase 2
          /*do_finalize=*/false,  // phase 1 exits right after sort
          cl,
          o.verbose);
    }
    for (auto& t : workers) t.join();
  }

  // --- Phase 2: serial dedup + reheader + index ----------------------
  // `o.threads` per invocation — one BAM at a time gets the whole budget.
  if (o.verbose >= 1)
    std::cerr << "svaba postprocess: phase 2 — serial dedup+reheader+index"
              << " (" << o.threads << " threads per BAM)" << std::endl;
  for (const auto& suffix : active) {
    processSuffix(
        o.id,
        suffix,
        o.threads,                   // full BGZF pool budget
        o.mem,
        /*do_sort=*/    false,       // already sorted in phase 1 (or skipped by --dedup-only)
        /*do_dedup=*/   !o.sort_only,
        /*do_finalize=*/true,
        cl,
        o.verbose);
  }

  if (o.verbose >= 1)
    std::cerr << "svaba postprocess: BAM processing done (sorted"
              << (o.sort_only ? "" : ", deduped where applicable")
              << ", PG-stamped, indexed)" << std::endl;

  // --- Phase 3: sort + dedup + filter bps.txt.gz -------------------------
  // Auto-detects bps.txt.gz (or .bps.txt); no-ops if absent.
  // Collects PASS/somatic cname sets for Phase 4 r2c filtering.
  std::unordered_set<std::string> pass_cnames, som_cnames;
  sortDedupFilterBps(o.id, o.verbose, pass_cnames, som_cnames);

  // --- Phase 4: filter r2c.txt.gz to PASS / PASS+somatic contigs ----------
  // Auto-detects r2c.txt.gz; no-ops if absent or no PASS cnames.
  filterR2c(o.id, pass_cnames, som_cnames, o.verbose);

  // --- Phase 5: split-by-source (optional, --split) --------------------------
  // Demux each dedup-eligible BAM by the first 4 chars of QNAME.
  if (o.split) {
    if (o.verbose >= 1)
      std::cerr << "svaba postprocess: phase 5 — split-by-source" << std::endl;
    for (const auto& suffix : kDedupSuffixes) {
      splitBySource(o.id, suffix, o.threads, o.mem, o.verbose);
    }
  }

  if (o.verbose >= 1)
    std::cerr << "svaba postprocess: done" << std::endl;
}
