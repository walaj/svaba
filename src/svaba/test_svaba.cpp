// test_svaba.cpp
//
// `svaba test` — a "secret" subcommand (not listed in the top-level help)
// that runs a minimal microbenchmark of svaba's five main hot paths:
//
//   1. walk     — read N reads from a BAM
//   2. correct  — BFC error-correct those reads
//   3. realign  — BWA-MEM align each *corrected* read back to the reference.
//                 This mirrors the per-read native alignment that
//                 BreakPoint::splitCoverage needs for its r2c-vs-native
//                 comparative gate (SVABA_R2C_NATIVE_GATE). On a real
//                 svaba run this is currently an underrated cost: N
//                 BWA-MEM calls PER WINDOW, where N is a few thousand
//                 on a typical WGS, times ~120k windows for a full run.
//   4. assemble — fermi-lite assemble the reads
//   5. align    — BWA-MEM align the resulting contigs back to a reference
//
// Each phase is timed independently across a sweep of N values (default
// 100, 500, 1 000, 5 000, 10 000, 50 000) so you can eyeball how each
// phase scales. Results are not checked — the point is wall-clock per-
// phase, not correctness. Repeat each trial -r times; the median is
// reported (first-trial cold-cache effects wash out).
//
// Designed to isolate the bottleneck profile independently of the real
// svaba pipeline's bookkeeping, scoring, and BAM writing — so when
// someone says "svaba is slow" we can point at a phase.
//
// Usage:
//   svaba test -i BAM -G REF [-n N ...] [-k REGION] [-r REPEATS] [--skip-*]
//
//   -i FILE     input BAM (any BAM with reads; aligned status is irrelevant)
//   -G FILE     BWA-indexed reference FASTA (needs .bwt / .pac / etc)
//   -n N        reads-per-trial; may repeat. Default sweep is
//               100 / 500 / 1000 / 5000 / 10000 / 50000.
//   -k REGION   samtools-style region (chr:start-end) to restrict the BAM
//               read. Defaults to the beginning of the file.
//   -r N        repetitions per trial (default 3). Reports median.
//   --skip-walk       (benchmark still does the read-gathering;
//                     just doesn't include it in the phase table)
//   --skip-ec         skip BFC error correction
//   --skip-assemble   skip fermi-lite assembly
//   --skip-align      skip BWA-MEM alignment
//   -v, --verbose
//   -h, --help

#include <getopt.h>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BFC.h"
#include "SeqLib/BWAAligner.h"
#include "SeqLib/BWAIndex.h"
#include "SeqLib/FermiAssembler.h"
#include "SeqLib/UnalignedSequence.h"

namespace {

struct TestOpts {
  std::string         bam;
  std::string         ref;
  std::string         region;
  std::vector<size_t> n_sweep;
  int                 repeats       = 3;
  bool                skip_walk     = false;
  bool                skip_ec       = false;
  bool                skip_realign  = false;
  bool                skip_assemble = false;
  bool                skip_align    = false;
  int                 verbose       = 0;
};

constexpr const char* kUsage =
"svaba test — microbenchmark the four main svaba hot paths.\n\n"
"Usage:\n"
"  svaba test -i BAM -G REF [-n N ...] [-k REGION] [-r REPEATS] [options]\n\n"
"  -i, --bam FILE        input BAM (reads source)\n"
"  -G, --reference FILE  BWA-indexed reference (for alignment phase)\n"
"  -n, --num-reads N     trial size; may be passed multiple times.\n"
"                        Default sweep: 500, 1000, 2000, 3000, 5000, 10000 —\n"
"                        the 2k-3k range spans typical svaba 25kb-window loads.\n"
"  -k, --region REGION   samtools-style region to restrict BAM reading\n"
"  -r, --repeats N       per-trial repetitions (median reported, default 3)\n"
"      --skip-walk       exclude walk phase from the table (reads are still fetched)\n"
"      --skip-ec         skip BFC error correction\n"
"      --skip-realign    skip BWA-MEM realignment of corrected reads\n"
"      --skip-assemble   skip fermi-lite assembly\n"
"      --skip-align      skip BWA-MEM alignment of the produced contigs\n"
"  -v, --verbose         chatty diagnostics\n"
"  -h, --help            this message\n";

constexpr int OPT_SKIP_WALK     = 3001;
constexpr int OPT_SKIP_EC       = 3002;
constexpr int OPT_SKIP_REALIGN  = 3003;
constexpr int OPT_SKIP_ASSEMBLE = 3004;
constexpr int OPT_SKIP_ALIGN    = 3005;

const struct option LONGOPTS[] = {
  { "bam",            required_argument, nullptr, 'i' },
  { "reference",      required_argument, nullptr, 'G' },
  { "num-reads",      required_argument, nullptr, 'n' },
  { "region",         required_argument, nullptr, 'k' },
  { "repeats",        required_argument, nullptr, 'r' },
  { "skip-walk",      no_argument,       nullptr, OPT_SKIP_WALK },
  { "skip-ec",        no_argument,       nullptr, OPT_SKIP_EC },
  { "skip-realign",   no_argument,       nullptr, OPT_SKIP_REALIGN },
  { "skip-assemble",  no_argument,       nullptr, OPT_SKIP_ASSEMBLE },
  { "skip-align",     no_argument,       nullptr, OPT_SKIP_ALIGN },
  { "verbose",        no_argument,       nullptr, 'v' },
  { "help",           no_argument,       nullptr, 'h' },
  { nullptr, 0, nullptr, 0 }
};

TestOpts parseCli(int argc, char** argv) {
  TestOpts o;
  bool die = (argc <= 2);
  optind = 1;
  for (int c; (c = getopt_long(argc, argv, "hvi:G:n:k:r:", LONGOPTS, nullptr)) != -1;) {
    switch (c) {
      case 'h': die = true; break;
      case 'v': o.verbose = 1; break;
      case 'i': o.bam     = optarg ? optarg : ""; break;
      case 'G': o.ref     = optarg ? optarg : ""; break;
      case 'k': o.region  = optarg ? optarg : ""; break;
      case 'r': o.repeats = std::max(1, std::atoi(optarg ? optarg : "3")); break;
      case 'n': {
        try { o.n_sweep.push_back(std::stoull(optarg ? optarg : "")); }
        catch (...) {
          std::cerr << "ERROR: -n expects a positive integer\n";
          std::exit(EXIT_FAILURE);
        }
        break;
      }
      case OPT_SKIP_WALK:     o.skip_walk     = true; break;
      case OPT_SKIP_EC:       o.skip_ec       = true; break;
      case OPT_SKIP_REALIGN:  o.skip_realign  = true; break;
      case OPT_SKIP_ASSEMBLE: o.skip_assemble = true; break;
      case OPT_SKIP_ALIGN:    o.skip_align    = true; break;
      default: die = true; break;
    }
  }
  if (!die) {
    if (o.bam.empty()) { std::cerr << "ERROR: -i BAM required\n"; die = true; }
    // Reference is needed for both the contig-align AND the corrected-read
    // realign phase. Require it unless BOTH are skipped.
    if ((!o.skip_align || !o.skip_realign) && o.ref.empty()) {
      std::cerr << "ERROR: -G REF required unless --skip-align AND --skip-realign\n";
      die = true;
    }
  }
  if (die) { std::cerr << "\n" << kUsage; std::exit(1); }
  if (o.n_sweep.empty())
    // Default sweep targets the real per-window range (~2k-3k reads/window
    // on typical WGS) with brackets below and above to show scaling.
    o.n_sweep = { 500, 1000, 2000, 3000, 5000, 10000 };
  return o;
}

// Read up to N BamRecord from the BAM, optionally restricted to `region`.
// Anything shorter than `min_read_len` (e.g. unmapped stubs) is skipped so
// the benchmark doesn't get polluted by degenerate records.
SeqLib::BamRecordVector readNReads(const std::string& bam,
                                   const std::string& region,
                                   size_t N,
                                   int min_read_len = 30) {
  SeqLib::BamReader r;
  if (!r.Open(bam)) {
    std::cerr << "ERROR: cannot open " << bam << "\n";
    std::exit(EXIT_FAILURE);
  }
  if (!region.empty()) {
    SeqLib::GenomicRegion gr(region, r.Header());
    if (!r.SetRegion(gr)) {
      std::cerr << "ERROR: SetRegion failed for '" << region << "'\n";
      std::exit(EXIT_FAILURE);
    }
  }
  SeqLib::BamRecordVector out;
  out.reserve(N);
  while (auto opt = r.Next()) {
    SeqLib::BamRecord& rec = *opt;
    if (rec.Length() < min_read_len) continue;
    // BamRecord is move-only (holds a shared_ptr<bam1_t>); transfer rather
    // than copy. r.Next() just handed us ownership via the optional, so
    // we move out of it and drop it before the next loop iteration.
    out.push_back(std::move(rec));
    if (out.size() >= N) break;
  }
  return out;
}

// Per-trial timing result. Fields are wall-clock seconds; 0.0 means the
// phase was skipped.
struct Phase {
  double walk     = 0.0;
  double ec       = 0.0;
  double realign  = 0.0;  // BWA-MEM of corrected reads -> reference
  double assem    = 0.0;  // "asm" is GCC-reserved (inline assembly)
  double align    = 0.0;  // BWA-MEM of contigs -> reference

  double total() const { return walk + ec + realign + assem + align; }
};

double since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

Phase oneTrial(const TestOpts& opts,
               size_t N,
               SeqLib::BWAIndexPtr& bwa_idx) {
  Phase p;
  using clock = std::chrono::steady_clock;

  // ---- Phase 1: walk (always runs; --skip-walk just hides it in the table).
  auto t0 = clock::now();
  SeqLib::BamRecordVector reads = readNReads(opts.bam, opts.region, N);
  p.walk = since(t0);

  if (reads.empty()) {
    std::cerr << "WARNING: no reads fetched for N=" << N << "\n";
    return p;
  }

  // ---- Phase 2: BFC error correction ---------------------------------------
  // Mirrors what SvabaRegionProcessor does: add sequences, train, correct.
  // Kept alive across phases so we can drain the corrected sequences for
  // the realign phase below.
  SeqLib::BFC bfc;
  if (!opts.skip_ec) {
    auto t1 = clock::now();
    for (const auto& rec : reads) {
      bfc.AddSequence(rec.Sequence(), rec.Qualities(0), rec.Qname());
    }
    bfc.Train();
    bfc.ErrorCorrect();
    p.ec = since(t1);
  }

  // ---- Phase 3: realign corrected reads -> reference -----------------------
  // Feeds BreakPoint::splitCoverage's r2c-vs-native comparative gate. Real
  // svaba runs ONE BWA-MEM call per corrected read per window, which on a
  // typical WGS is a few thousand calls per window × ~120k windows. This
  // is the phase the user was curious about — the "additional alignment"
  // that the comparative split-coverage gate introduced.
  //
  // We drain the corrected sequences from BFC via GetSequence() rather
  // than re-reading the original BAM sequence, so timing reflects the
  // same payload svaba's gate operates on.
  if (!opts.skip_realign && !opts.skip_ec) {
    auto t_r = clock::now();
    SeqLib::BWAAligner aligner(bwa_idx);
    aligner.allocBuffer(4096);
    bfc.ResetGetSequence();
    std::string s, name;
    size_t nr = 0;
    while (bfc.GetSequence(s, name)) {
      SeqLib::BamRecordPtrVector hits;
      // Only need the best hit for a native-score comparison, so keep
      // maxSecondary low and keepSecFrac at 0. Matches the intent of
      // the SvabaRegionProcessor corrected_native_cig realign.
      aligner.alignSequence(s, name.empty() ? std::to_string(nr) : name,
                            hits,
                            /*hardclip=*/false,
                            /*keepSecFrac=*/0.0,
                            /*maxSecondary=*/0);
      ++nr;
    }
    p.realign = since(t_r);
    if (opts.verbose) {
      std::cerr << "  [N=" << N << "] realigned " << nr << " corrected reads\n";
    }
  } else if (!opts.skip_realign && opts.skip_ec) {
    // --skip-ec pre-empts realign (no corrected sequences to align).
    if (opts.verbose) {
      std::cerr << "  [N=" << N << "] realign skipped (requires --skip-ec unset)\n";
    }
  }

  // ---- Phase 4: fermi-lite assembly ----------------------------------------
  // Feed raw BamRecords; the SeqLib wrapper extracts sequences.
  std::vector<std::string> contigs;
  if (!opts.skip_assemble) {
    auto t2 = clock::now();
    SeqLib::FermiAssembler fml;
    fml.AddReads(reads);
    fml.PerformAssembly();
    contigs = fml.GetContigs();
    p.assem = since(t2);
    if (opts.verbose) {
      std::cerr << "  [N=" << N << "] assembled " << contigs.size()
                << " contigs\n";
    }
  }

  // ---- Phase 4: BWA-MEM align the contigs back to the reference ------------
  // If no contigs were produced (either --skip-assemble was set or fermi
  // produced nothing), fall back to aligning the raw read sequences so the
  // phase still exercises BWA. This is deliberately a coarse measure —
  // svaba's real align phase has a bunch of contig-sanity logic on top.
  if (!opts.skip_align) {
    auto t3 = clock::now();
    SeqLib::BWAAligner aligner(bwa_idx);
    aligner.allocBuffer(4096);
    std::vector<std::string>* targets = &contigs;
    std::vector<std::string> raw_seqs;
    if (targets->empty()) {
      raw_seqs.reserve(reads.size());
      for (const auto& rec : reads) raw_seqs.push_back(rec.Sequence());
      targets = &raw_seqs;
    }
    for (size_t i = 0; i < targets->size(); ++i) {
      SeqLib::BamRecordPtrVector hits;
      aligner.alignSequence((*targets)[i], "t" + std::to_string(i),
                            hits, /*hardclip=*/false,
                            /*keepSecFrac=*/0.9, /*maxSecondary=*/10);
    }
    p.align = since(t3);
  }

  return p;
}

double median(std::vector<double> v) {
  if (v.empty()) return 0.0;
  std::sort(v.begin(), v.end());
  const size_t m = v.size() / 2;
  return (v.size() % 2) ? v[m] : 0.5 * (v[m-1] + v[m]);
}

// Compact fixed-width table printer. Each row = one (N, phase) cell with
// the median timing and a per-read microseconds column so you can eyeball
// scaling linearity. Scaling column is the ratio vs the smallest N.
void printTable(const TestOpts& opts,
                const std::vector<std::pair<size_t, std::vector<Phase>>>& results) {
  std::cout << "\n";
  std::cout << std::left << std::setw(10) << "phase"
            << std::right << std::setw(10) << "N"
            << std::setw(12) << "median(s)"
            << std::setw(14) << "per-read(us)"
            << std::setw(10) << "scale"
            << "\n";
  std::cout << std::string(56, '-') << "\n";

  struct PhaseInfo { const char* label; double Phase::*field; bool skipped; };
  const PhaseInfo phases[] = {
    { "walk",     &Phase::walk,    opts.skip_walk },
    { "correct",  &Phase::ec,      opts.skip_ec },
    { "realign",  &Phase::realign, opts.skip_realign || opts.skip_ec },
    { "assemble", &Phase::assem,   opts.skip_assemble },
    { "align",    &Phase::align,   opts.skip_align },
  };

  for (const auto& ph : phases) {
    if (ph.skipped) continue;
    double base_per_read = -1.0;   // baseline for the scaling column (smallest N)
    for (const auto& [N, trials] : results) {
      std::vector<double> v;
      v.reserve(trials.size());
      for (const auto& t : trials) v.push_back(t.*ph.field);
      const double med = median(v);
      const double us_per_read = (N > 0) ? 1e6 * med / static_cast<double>(N) : 0.0;
      if (base_per_read < 0) base_per_read = us_per_read;
      const double scale =
          (base_per_read > 0) ? us_per_read / base_per_read : 1.0;

      std::cout << std::left  << std::setw(10) << ph.label
                << std::right << std::setw(10) << N
                << std::setw(12) << std::fixed << std::setprecision(4) << med
                << std::setw(14) << std::fixed << std::setprecision(1) << us_per_read
                << std::setw(10) << std::fixed << std::setprecision(2) << scale
                << "\n";
    }
    std::cout << "\n";
  }

  // Totals row: sum of non-skipped phases, median across trials.
  std::cout << std::left << std::setw(10) << "TOTAL"
            << std::right << std::setw(10) << " "
            << std::setw(12) << "median(s)"
            << std::setw(14) << "per-read(us)"
            << "\n";
  std::cout << std::string(46, '-') << "\n";
  for (const auto& [N, trials] : results) {
    std::vector<double> v;
    for (const auto& t : trials) v.push_back(t.total());
    const double med = median(v);
    const double us_per_read = (N > 0) ? 1e6 * med / static_cast<double>(N) : 0.0;
    std::cout << std::left  << std::setw(10) << ""
              << std::right << std::setw(10) << N
              << std::setw(12) << std::fixed << std::setprecision(4) << med
              << std::setw(14) << std::fixed << std::setprecision(1) << us_per_read
              << "\n";
  }
  std::cout << "\nScaling column is (per-read us at this N) / (per-read us at smallest N).\n"
            << "Near 1.00 = linear scaling. >1 means per-read cost grows with N\n"
            << "(= superlinear algorithm or cache effects). <1 means setup cost\n"
            << "amortizes as N grows.\n";
}

}  // namespace

// Forward-declared in svaba.cpp dispatch.
void runTest(int argc, char** argv) {
  const TestOpts opts = parseCli(argc, argv);

  std::cerr << "svaba test:\n"
            << "  bam:     " << opts.bam    << "\n"
            << "  ref:     " << opts.ref    << "\n"
            << "  region:  " << (opts.region.empty() ? "(all)" : opts.region) << "\n"
            << "  repeats: " << opts.repeats << "\n"
            << "  N sweep:";
  for (auto n : opts.n_sweep) std::cerr << " " << n;
  std::cerr << "\n";

  // Load BWA index once up front; it's large and expensive to load
  // repeatedly. Reuse across every trial and every N.
  SeqLib::BWAIndexPtr bwa_idx;
  if (!opts.skip_align) {
    std::cerr << "loading BWA index from " << opts.ref << " ...\n";
    auto t0 = std::chrono::steady_clock::now();
    bwa_idx = std::make_shared<SeqLib::BWAIndex>();
    bwa_idx->LoadIndex(opts.ref);
    std::cerr << "  loaded in " << std::fixed << std::setprecision(1)
              << since(t0) << "s\n";
  }

  std::vector<std::pair<size_t, std::vector<Phase>>> results;
  results.reserve(opts.n_sweep.size());

  for (size_t N : opts.n_sweep) {
    std::vector<Phase> trials;
    trials.reserve(opts.repeats);
    for (int r = 0; r < opts.repeats; ++r) {
      if (opts.verbose)
        std::cerr << "  [N=" << N << " trial " << (r+1) << "/" << opts.repeats << "]\n";
      trials.push_back(oneTrial(opts, N, bwa_idx));
    }
    results.emplace_back(N, std::move(trials));
  }

  printTable(opts, results);
}
