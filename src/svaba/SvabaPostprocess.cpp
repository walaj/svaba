// SvabaPostprocess.cpp — see SvabaPostprocess.h for the design rationale.
//
// The flow for each suffix runs in its own std::thread:
//
//   [${ID}.${suffix}.bam]
//        |
//        |  (unless --dedup-only)
//        v
//   samtools sort -@ jthreads -m MEM -o sorted.bam ${bam}       # shell out
//   rename sorted.bam -> ${bam}
//        |
//        |  (if suffix in DEDUP set, unless --sort-only)
//        v
//   native streaming dedup+merge  -> deduped.bam                 # htslib
//   rename deduped.bam -> ${bam}
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

#include <getopt.h>
#include <sys/stat.h>

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <vector>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BamWriter.h"

namespace {

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
  bool        sort_only = false;  // skip dedup
  bool        dedup_only = false; // skip sort (assume sorted)
};

enum { OPT_SORT_ONLY = 1000, OPT_DEDUP_ONLY, OPT_MEM };

const char* kShortOpts = "hi:t:m:v:";
const struct option kLongOpts[] = {
  { "help",       no_argument,       nullptr, 'h' },
  { "id",         required_argument, nullptr, 'i' },
  { "threads",    required_argument, nullptr, 't' },
  { "mem",        required_argument, nullptr, 'm' },
  { "verbose",    required_argument, nullptr, 'v' },
  { "sort-only",  no_argument,       nullptr, OPT_SORT_ONLY  },
  { "dedup-only", no_argument,       nullptr, OPT_DEDUP_ONLY },
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
                       const std::string& tag,
                       int verbose) {
  using clock = std::chrono::steady_clock;
  DedupStats st;
  const auto t0 = clock::now();
  auto last = t0;
  const off_t in_size = fileSize(in_bam);

  SeqLib::BamReader r;
  if (!r.Open(in_bam))
    throw std::runtime_error("dedup: cannot open input " + in_bam);

  SeqLib::BamWriter w;
  if (!w.Open(out_bam))
    throw std::runtime_error("dedup: cannot open output " + out_bam);
  w.SetHeader(r.Header());
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
  idx_by_key.reserve(1024);

  int32_t cur_chr = -2;  // -2 is a sentinel that won't match any valid ChrID
  int32_t cur_pos = -2;
  // Cached display name for cur_chr. ChrName() hits the header every call, so
  // we resolve it only on (rare) chromosome change, not per record.
  std::string cur_chr_name = "*";

  // Flush buffered records for the current locus in insertion order.
  auto flushLocus = [&]() {
    for (auto& kept : keep_buf) {
      if (!w.WriteRecord(kept))
        throw std::runtime_error("dedup: WriteRecord failed on " + out_bam);
    }
    keep_buf.clear();
    idx_by_key.clear();
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

    if (verbose >= 1 && (st.in_records % 1'000'000ULL) == 0) {
      const auto now = clock::now();
      const double elapsed = std::chrono::duration<double>(now - t0).count();
      const double dt      = std::chrono::duration<double>(now - last).count();
      const double inst_rps = dt > 0 ? 1'000'000.0 / dt : 0.0;
      const double avg_rps  = elapsed > 0 ? st.in_records / elapsed : 0.0;
      last = now;

      // Show input size as a coarse "how far into the job are we" hint.
      // A true ETA would need bytes-read from htslib, which SeqLib's API
      // doesn't currently surface; skip rather than fabricate.
      std::string size_hint;
      if (in_size > 0) {
        std::ostringstream oss;
        oss << " (input " << (in_size / (1024 * 1024)) << " MiB)";
        size_hint = oss.str();
      }

      // Report the most recently seen locus as a coarse progress mark.
      // cur_pos is 0-based (BAM convention); print 1-based for readability.
      std::string locus_str;
      {
        std::ostringstream oss;
        oss << cur_chr_name << ":"
            << (cur_pos >= 0 ? cur_pos + 1 : 0);
        locus_str = oss.str();
      }

      logLine("[", tag, "] dedup: at ", locus_str, " | ",
              st.in_records / 1'000'000ULL, "M reads | ",
              st.kept, " kept, ", st.merged, " merged (",
              st.tag_edits, " tag-edits) | ",
              std::fixed, std::setprecision(2),
              inst_rps / 1e6, "M/s last-1M, ",
              avg_rps / 1e6,  "M/s avg, ",
              std::setprecision(1), elapsed, "s elapsed",
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

// Run the full sort + dedup pipeline for a single suffix. Safe to invoke
// from its own thread. Never throws out — caller joins; any failure is
// logged and the pipeline aborts for this suffix only.
void processSuffix(const std::string& id,
                   const std::string& suffix,
                   int  per_job_threads,
                   const std::string& mem,
                   bool do_sort,
                   bool do_dedup,
                   int  verbose) {
  const std::string bam     = id + "." + suffix + ".bam";
  const std::string sorted  = id + "." + suffix + ".sorted.bam";
  const std::string deduped = id + "." + suffix + ".deduped.bam";

  if (!fileExists(bam)) {
    logLine("[", suffix, "] skipping: ", bam, " not found");
    return;
  }

  try {
    if (do_sort) {
      const int rc = runSort(bam, sorted, per_job_threads, mem, suffix, verbose);
      if (rc != 0)
        throw std::runtime_error("samtools sort exited " + std::to_string(rc));
      renameOrThrow(sorted, bam);
    }

    if (do_dedup && isDedupSuffix(suffix)) {
      if (verbose >= 1)
        logLine("[", suffix, "] starting dedup on ", bam);
      const DedupStats ds = streamDedup(bam, deduped, suffix, verbose);
      renameOrThrow(deduped, bam);
      logLine("[", suffix, "] dedup complete: ",
              ds.in_records, " in, ",
              ds.kept, " kept, ",
              ds.merged, " merged (",
              (ds.in_records ?
                100.0 * ds.merged / static_cast<double>(ds.in_records) : 0.0),
              "% dup, ",
              ds.tag_edits, " bi/bz tag-edits) in ",
              std::fixed, std::setprecision(1), ds.seconds, "s");
    } else if (do_dedup && !isDedupSuffix(suffix)) {
      if (verbose >= 2)
        logLine("[", suffix, "] dedup skipped (not a dedup suffix)");
    }
  } catch (const std::exception& e) {
    logLine("[", suffix, "] ERROR: ", e.what());
  }
}

}  // namespace

void runPostprocess(int argc, char** argv) {
  const Opts o = parseOpts(argc, argv);

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

  // Launch one thread per active suffix. max_parallel is only used to cap the
  // samtools sort thread count — since we only launch n_active workers total,
  // we don't need an external throttle. The heavy lifting inside each worker
  // is the samtools sort subprocess (which itself respects per_job_sort_t),
  // plus a CPU-light streaming dedup loop.
  std::vector<std::thread> workers;
  workers.reserve(active.size());
  for (const auto& suffix : active) {
    workers.emplace_back(
        processSuffix,
        o.id,
        suffix,
        per_job_sort_t,
        o.mem,
        /*do_sort=*/  !o.dedup_only,
        /*do_dedup=*/ !o.sort_only,
        o.verbose);
  }
  for (auto& t : workers) t.join();

  if (o.verbose >= 1)
    std::cerr << "svaba postprocess: all suffixes done" << std::endl;
}
