#pragma once

#include <string>
#include <vector>
#include <map>
#include <cstddef>

// version & date  --  OWNED BY CLAUDE (see CLAUDE.md "Versioning"). On every
// meaningful change, bump this per semver and add a matching CHANGELOG.md entry;
// these two constants drive the startup banner / `svaba --version` / @PG lines.
// The build commit is stamped separately (SvabaGitVersion.h).
//
// 2.4.0: recurrent 5' soft-clip "tag" detection + auto-trim. An untrimmed 5'
//   molecular tag (UMI / inline barcode / adapter / library spacer) mismatches
//   the reference at the read's 5' end on nearly every read; that mismatch sits
//   at the leading edge of each read's overlap with its neighbor and breaks the
//   overlap (string-graph) assembly genome-wide (measured: 133->94 calls on a
//   2 Mb region with a synthetic 5 bp tag). During insert-size learning,
//   BamReadGroup::addRead now accumulates a 5'-soft-clip-length histogram and
//   detectTag() flags a sharp short-clip spike (>=20% of reads at the modal
//   length) as a tag, classifying fixed-seq vs random/UMI. It prints a per-RG
//   WARNING, records `tag_trim_5p` in the .learn.tsv.gz and the --bam-params
//   cache (so scatter-gather shards reuse it without re-detecting), and svaba
//   then auto-trims min(5'-softclip, L) bases off each read's 5' end before
//   assembly (svabaRead::TrimTag5p, gated on SvabaSharedConfig::tag_trim_5p).
//   The trim is self-targeting (only clipped/non-genomic bases, up to L; long
//   split-read clips keep everything past L) and assembly-only (original BAM
//   record/CIGAR untouched -> r2c/scoring/coords unaffected). Verified: auto-
//   trim fully recovers the broken assembly (94->133, exact), inert + byte-
//   identical on clean BAMs (no false positives). Override: --tag-trim N forces
//   N bp, --no-tag-trim disables (detection/warning still run).
// 2.3.3: `--max-cov` no longer perturbs the error model of clean reads. The
//   weird-coverage subsample demotes pileup reads out of fermi assembly
//   (to_assemble=false) to cut the super-linear assembly cost, but it had ALSO
//   been excluding them from BFC k-mer training -- so capping a pileup changed
//   the auto-learned k / k-mer spectrum for the WHOLE assembly window, which
//   silently re-error-corrected (and sometimes dropped) clean calls elsewhere
//   in that window. Example: a clean somatic Alu-Alu deletion at chr1:44147919
//   (somlod 2.14, no breakend reads demoted) vanished only because a pileup ~9kb
//   away was capped. Fix: demoted reads now carry `svabaRead::cov_demoted` and
//   stay in the BFC *training* pool (Phase 1) while remaining excluded from the
//   correction pool and fermi. Training is ~linear/cheap, so the speedup is kept
//   (measured chr1 give-back ~+27s CPU = ~1.5%); the call returns byte-identical
//   (scores unchanged, only the contig serial differs). Residual cap-vs-no-cap
//   churn is now only the unavoidable fermi-graph effect of not assembling a
//   pileup, confined to the marginal somlod<=3 tier.
// 2.3.2: DETERMINISM FIX, part 2 -- byte-identical output. With the calls now
//   deterministic (2.3.1), the only remaining run-to-run variation was the
//   `bp_id` column (bps col 52) + the line ORDER of equal-key rows in the sorted
//   bps. (1) bp_id was "bpTTTNNNNNNNN" = worker-thread ID + per-thread counter,
//   so the SAME variant got a different id depending on which thread/order
//   processed it (and the bi:Z BAM tag, r2c.db split_bps/disc_bps, and VCF EVENT
//   join on it). It's now derived from the BP's stable content (cname + both
//   breakends + contig positions + kind) via FNV-1a -> "bp"+16 hex; the same BP
//   yields the same id across runs/threads/-p. (2) `svaba postprocess` bps sort
//   (std::sort, unstable) now has a final total-order tie-break on the raw line
//   bytes, so equal-key rows can't reorder by thread-emission order. Result: a
//   plain `diff` of .bps.sorted[.dedup[.pass[.somatic]]].txt.gz is now 0 across
//   runs. NOTE: bp_id format changed -- tools that treat it opaquely (all of
//   svaba's own) are unaffected; anything parsing the old thread+counter layout
//   must adapt.
// 2.3.1: DETERMINISM FIX (multithreading). svaba was nondeterministic at -p>1:
//   identical inputs/command produced ~0.3% different calls run-to-run (and ~8%
//   in the somatic-PASS subset, since somatic classification is a threshold on
//   the sensitive somlod). Root cause: the BFC error-corrector's auto-learned
//   k-mer size was cached in the reused per-thread BFC object and frozen from
//   whichever region a worker processed FIRST; region->thread assignment varies
//   at -p>1, so the EC k (hence error-correction, assembly, and the call set)
//   varied between runs. Fix: re-learn k per region in BFC::Train() (SeqLib),
//   honoring an explicit SetKmer(). Verified: -p8 now byte-matches -p1. NOTE:
//   output changes slightly vs 2.3.0 for regions whose own k differed from the
//   frozen one -- this is the corrected, reproducible behavior.
// 2.3.0: `--max-normal-weird-cov N` (opt-in, 0 = off). Somatic-safe artifact
//   skip: in any sub-region where the NORMAL weird-read coverage exceeds N, the
//   reads are dropped from assembly + error-correction + r2c + scoring entirely
//   (DP coverage + discordant clustering already done, so unaffected). Such
//   pileups are shared artifacts/repeats — a somatic event needs a CLEAN normal,
//   so this can never drop a somatic call (keyed on NORMAL only; a tumor-only
//   pileup, which could be somatic, is never skipped). Targets the few
//   high-pileup regions that dominate runtime even after --max-cov; demonstrated
//   ~40%+ region-time cut when it fires. Threshold must sit above germline
//   weird-coverage (scales with normal depth); ~200 is safe for deep WGS
//   normals — validate recall on the sim panel for your depth before relying.
// 2.2.1: `--max-cov` actually works again. It was dead code: the
//   subSampleToWeirdCoverage() call was commented out in svabaBamWalker::readBam
//   AND walker->max_cov was never wired to opts.maxCov, so high-weird-coverage
//   pileups flooded BFC+fermi regardless of the (parsed, logged) cap. Fix: wire
//   walker->max_cov = opts.maxCov in SvabaRegionProcessor, restore the call
//   (guarded `if (get_coverage && max_cov > 0)`). Rather than DELETE excess
//   reads, they are DEMOTED to second-class (to_assemble=false): excluded from
//   fermi assembly + BFC (the read-count-scaling cost = the speedup) but kept in
//   `reads` so they're still r2c-aligned, discordant-clustered, and scored. This
//   preserves support for very-high-coverage real events (amplicons) and never
//   samples away a normal read that could refute a somatic call. Changes results
//   only where per-position WEIRD-read coverage exceeds --max-cov (default 100);
//   deterministic hash-by-qname (mate pairs together), symmetric tumor/normal.
// 2.2.0: `--bam-params FILE` caches learned per-RG insert-size params (+
//   `--learn-only` precompute step), so a scatter-gather run learns once and
//   every shard reuses it instead of re-doing the genome-wide insert sweep —
//   removes a large redundant per-shard CPU + shared-disk cost.
// 2.1.2: contig-alignment confidence (contig_conf / WEAKCONTIG) now penalizes
//   alignment non-uniqueness -- BWA XS≈AS (another near-equal placement, i.e. a
//   repeat/paralog/segdup the contig BLATs to many sites) -- via a new Rule E in
//   scoreContigAlignment. The as_xs_gap signal was computed but unused.
// 2.1.1: LOWSPANDSCRD now catches span-0 FR discordant clusters (degenerate
//   zero-length "deletions" from mildly-over-insert FR pairs) that leaked to
//   PASS as false somatic calls -- `getSpan() > 0` -> `>= 0` in score_dscrd.
// 2.1.0 (since 2.0.0): somatic-safety junction-kmer net (KC FORMAT field),
// `--annotation` BED + per-breakend repeat_anno/poly_a (v6 bps schema),
// DUPREADS via unique split-read starts, hasAdapter 3'-clip fix, and assorted
// false-somatic / SV-coordinate correctness fixes. See CHANGELOG.md.
inline constexpr char SVABA_VERSION[] = "2.5.0";
inline constexpr char SVABA_DATE[]    = "08/2026";

// from AlignmentFragment.h
inline constexpr std::size_t MAX_CONTIG_SIZE = 5'000'000;

// from run_svaba.cpp
inline constexpr int MIN_CLIP_FOR_LOCAL           = 40;
inline constexpr int MAX_NM_FOR_LOCAL             = 10;

// from BreakPoint
inline constexpr double MAX_ERROR                   = 0.08; // for very repetitive regions, what is the maximum expected error rate for an indel
inline constexpr double MIN_ERROR                   = 0.0005; 
inline constexpr int    T_SPLIT_BUFF                = 5;
inline constexpr int    N_SPLIT_BUFF                = 5;
inline constexpr int    INSERT_SIZE_TOO_BIG_SPAN_READS = 16;
// (Removed in SvABA2.0) R2C_MAX_NONEDGE_SOFTCLIP_{TUMOR,NORMAL} —
// replaced by the principled per-read score comparison gate in
// BreakPoint::splitCoverage(). See svaba::readAlignmentScore() in
// ContigAlignmentScore.h: a read is now credited as a variant supporter
// only when its r2c alignment scores strictly higher than its native
// read-to-reference alignment, which subsumes both the interior-clip
// rejection and the duplicated-reference equally-clean case.

// SvABA2.0 (v3): r2c-better-than-native gate. A read is credited as a
// split-supporter iff r2c_score > native_score (strict greater-than).
// Ties don't credit the read — if the contig provides no better
// explanation than the native alignment, the read is not evidence for
// the variant.
//
// v3 originally used a per-sample percentage margin (10% for tumor,
// 0% for normal) to filter junction-homology borderline cases on the
// tumor side. v3.1 removed the margin entirely:
//
//   - The 10% margin was read-length-dependent and mathematically
//     impossible to clear for small indels (a 1bp del on a 150bp read
//     gives only 4.9% improvement; on a 250bp read only 2.9%).
//   - Junction-homology cases where r2c barely beats native by 1-2
//     points: normal already credits these reads (margin was always 0
//     for normal). If both samples credit borderline reads equally,
//     the downstream LOD model sees similar split support in both →
//     low somlod → correctly not called somatic. The somatic/germline
//     distinction is the LOD model's job, not the split-coverage gate's.
//   - The margin was belt-and-suspenders with a worse failure mode
//     (killing real small indels) than the problem it prevented.
//
// This replaces the older both_split / one_split / homlen branching
// gate — when junction homology is on the order of the read length,
// r2c and native tie and such reads fall out naturally via the strict >
// comparison, which is the correct conservative behavior.
inline constexpr double T_R2C_MIN_MARGIN          = 0.0;
inline constexpr double N_R2C_MIN_MARGIN          = 0.0;

// ---------------------------------------------------------------------------
// SVABA_R2C_NATIVE_GATE — compile-time kill-switch for the r2c-vs-native
// alignment-score split-coverage gate in BreakPoint::splitCoverage.
//
//   1  (default)  — gate enabled. Correct, recommended. Each candidate
//                   split-supporting read is scored against both its r2c
//                   alignment and its native (read→reference) alignment,
//                   and must win by the per-sample-prefix margin above
//                   (T_R2C_MIN_MARGIN / N_R2C_MIN_MARGIN). Fixes the
//                   long-homology false-positive somatic bug.
//
//   0  (opt-in)   — gate disabled. Falls back to counting any r2c-spanning
//                   read as a split supporter regardless of how it
//                   compares to the native alignment. Reintroduces the
//                   homology-trap bug (long junction homology + clean
//                   normal-side reads → spurious somatic calls), so do
//                   NOT use for real calling. Provided solely to isolate
//                   the CPU cost of native_score computation on dense
//                   contigs — flip to 0, rebuild, time a reference
//                   region. Delta vs default-build tells you what the
//                   gate costs.
//
// Build with:
//   cmake .. -DCMAKE_CXX_FLAGS='-DSVABA_R2C_NATIVE_GATE=0'
// ---------------------------------------------------------------------------
#ifndef SVABA_R2C_NATIVE_GATE
#define SVABA_R2C_NATIVE_GATE 1
#endif
inline constexpr int    HOMOLOGY_FACTOR             = 4;
inline constexpr int    MIN_SOMATIC_RATIO           = 15;
inline constexpr int    COVERAGE_AVG_BUFF           = 10;

// from DiscordantCluster
inline constexpr int DISC_PAD                 = 150;
inline constexpr int MIN_PER_CLUST            = 3;
inline constexpr int DEFAULT_ISIZE_THRESHOLD  = 2000;

// from run_svaba
inline constexpr std::size_t THREAD_READ_LIMIT      = 20'000;
inline constexpr int         THREAD_CONTIG_LIMIT    =   5'000;

// from svabaAssemblerEngine
inline constexpr std::size_t MAX_OVERLAPS_PER_ASSEMBLY = 20'000;
inline constexpr int         MIN_CONTIG_MATCH           =    35;
inline constexpr int         MATE_LOOKUP_MIN            =     3;
inline constexpr int         SECONDARY_CAP              =    10;
inline constexpr int         MAX_MATE_ROUNDS            =     1;
inline constexpr std::size_t MAX_NUM_MATE_WINDOWS      = 50'000'000;
inline constexpr int         GERMLINE_CNV_PAD           =    10;
inline constexpr int         GET_MATES                  =     1;
inline constexpr int         LARGE_INTRA_LOOKUP_LIMIT   = 50'000;
inline constexpr double      SECONDARY_FRAC             =  0.90;

// from svabaBamWalker — mate-region lookup
//
// A discordant read contributes its mate locus as a candidate mate region
// iff its own MAPQ >= MIN_MATE_MAPQ AND the mate locus is on a chromosome
// with ChrID <= MAX_MATE_CHR_ID (set to -1 to disable the chromosome
// gate, e.g. for non-human genomes; default 23 = through chrY). After
// candidate regions are merged,
// a region must have >= MATE_REGION_MIN_COUNT supporting reads to survive
// the initial BamWalker filter, and then >= mateLookupMin (runtime option,
// default MATE_LOOKUP_MIN=3) to trigger the actual somatic mate lookup
// in SvabaRegionProcessor.
//
// Defaults:
//   MIN_MATE_MAPQ          = -1   (no MAPQ gate; any mapped read qualifies)
//   MAX_MATE_CHR_ID        =  23  (skip chrM/alt/decoy in human; -1 = no limit)
//                                  0-indexed: 0=chr1 .. 21=chr22, 22=chrX, 23=chrY
//   MATE_REGION_MIN_COUNT  =   2  (at least 2 reads to form a candidate)
//   MATE_REGION_PAD        = 250  (bp padding on each side of mate locus)
inline constexpr int    MIN_MATE_MAPQ                      =    -1;
inline constexpr int    MAX_MATE_CHR_ID                    =    23;
inline constexpr int    MATE_REGION_MIN_COUNT              =     2;
inline constexpr int    MATE_REGION_PAD                    =   250;

// from svabaBamWalker — other
inline constexpr int MIN_DSCRD_READS_DSCRD_ONLY          = 6;
inline constexpr int MIN_ISIZE_FOR_DISCORDANT_REALIGNMENT = 1'000;
inline constexpr int DISC_REALIGN_MATE_PAD                =   100;
inline constexpr int MAX_SECONDARY_HIT_DISC               =    10;

// coverage buffer
inline constexpr int INFORMATIVE_COVERAGE_BUFFER = 0;

// from vcf
inline constexpr int VCF_SECONDARY_CAP = 200;
inline constexpr int SOMATIC_LOD       =   1;
inline constexpr int DEDUPEPAD         = 200;

class SvabaLogger;

class SvabaOptions {

 public:
  // high-level flags
  bool   help        = false;
  int    verbose     = 0;
  int    numThreads  = 1;
  std::string analysisId;
  bool hp            = false;
  int perRgLearnLimit = 1'000;
  size_t weird_read_limit = 15'000;
  size_t mate_region_lookup_limit = 5'000;

  // dumping
  //
  // weird reads: compile-time-only toggle. Flip this to `true` here and
  // rebuild if you really want the weird-reads BAM. Deliberately NOT
  // exposed on the CLI — it's a large, niche output mostly used for
  // debugging the read-collection phase, not for routine runs.
  static constexpr bool dump_weird_reads = false;

  // All three below are opt-in at runtime via a single --dump-reads flag.
  // Default off so routine runs don't pay the (substantial) I/O and disk
  // cost of emitting this per-read detail on deep samples:
  //
  //   dump_discordant_reads  -> ${ID}.discordant.bam
  //   dump_corrected_reads   -> ${ID}.corrected.bam
  //   dump_alignments        -> ${ID}.alignments.txt.gz
  //                          -> ${ID}.r2c.txt.gz
  //
  // All three flip together under --dump-reads. The fields are kept
  // separate so an individual callsite can still key off its own narrow
  // concern (e.g. SvabaOutputWriter gates the alignments file streams on
  // dump_alignments only), but there is intentionally no way to toggle
  // them individually at runtime — that would bloat the CLI surface
  // without serving a real workflow.
  bool dump_discordant_reads = false;
  bool dump_corrected_reads  = false;
  bool dump_alignments       = false;

  // --bam-params FILE: cache the learned per-RG insert-size params. If the file
  // exists and covers every input BAM, svaba LOADS it and SKIPS the genome-wide
  // insert-size sweep (big win for scatter-gather: learn once, reuse per shard);
  // otherwise svaba learns normally and WRITES the file. --learn-only learns,
  // writes the file, and exits before any region processing (the precompute step).
  std::string bamParamsFile;
  bool learnOnly = false;

  // --r2c-min-somlod: gate on r2c.db size. A contig's r2c reads are written
  // only if the MAX somlod (BreakPoint::LO_s) across its breakpoints is
  // STRICTLY greater than this. Default -1e9 = write everything (preserves
  // behavior). Set to 0 to keep only somlod>0 (~somatic) events — LO_s
  // defaults to 0 for unscored/germline BPs, so strict `>` is required to
  // actually drop them. Only consulted when dump_alignments is on.
  double r2cMinSomlod = -1e9;

    // inputs
  std::vector<std::string> caseBams;
  std::vector<std::string> controlBams;
  std::string refGenome;
  std::string regionFile;

  // make the log verbose, but will invoke a lot of mutex locks
  bool verbose_log = false;

  int windowpad          = 500;

  std::string main_bam;
  
  // mode flags
  bool singleEnd         = false;
  bool allContigs        = false;
  bool discClusterOnly   = false;
  bool overrideRefCheck  = false;

  // numeric thresholds
  double sdDiscCutoff       = 3.92;  // SD multiplier for tumor discordant cutoff
  double sdDiscCutoffNormal = 3.60;  // SD multiplier for normal (lower = more sensitive)
  int    chunkSize          = 25000;
  int32_t maxReadsPerAssem  = -1;

  // this is site-level Log-odds cutoff for PASS that is variant 
  double lod                = 1.0; //8.0;
  // this is site-level Log-odds cuttof for PASS that is variant, if also has supporting DBSNP site
  double lodDb              = 1.0; //6.0;
  // this is log-odds that this is REF in the "worst" normal (the one with the most potential alt reads)
  double lodSomatic         = 0.0; //6.0;
  // same, but be more strict if this somatic variant overlaps a dbsnp site
  double lodSomaticDb       = 2.0; // 10.0;
  
  int    maxCov             = 100;
  // --max-normal-weird-cov N (0 = off): drop reads in sub-regions where the
  // NORMAL weird-read coverage exceeds N from assembly + r2c entirely. Such
  // pileups are shared artifacts (somatic needs a clean normal), so this is
  // somatic-safe and targets the high-pileup regions that dominate runtime.
  int    maxNormalWeirdCov  = 0;

  // 5' soft-clip tag trimming. -1 = AUTO: use the recurrent-tag length detected
  // during learning (svaba trims an untrimmed UMI/adapter/spacer off read 5'
  // ends before assembly). --tag-trim N forces N bp; --no-tag-trim sets 0
  // (disable). The effective value lands in SvabaSharedConfig::tag_trim_5p.
  int    tagTrimOverride    = -1;

  // SGA / assembly
  int    sgaMinOverlap      = 0;
  float  sgaErrorRate       = 0.0f;
  int    sgaNumRounds       = 3;

  // error correction
  std::string ecCorrectType = "f";   //s, f, or 0
  double      ecSubsample   = 0.50;

  // BWA-MEM tuning
  int   bwaGapOpen       = 32;
  int   bwaGapExt        = 1;
  int   bwaMismatch      = 18;
  int   bwaMatchScore    = 2;
  int   bwaZdrop         = 100;
  int   bwaBandwidth     = 1000;
  float bwaReseedTrigger = 1.5f;
  int   bwaClip3         = 5;
  int   bwaClip5         = 5;

  // filtering
  std::string rulesJson;
  size_t      mateLookupMin       = MATE_LOOKUP_MIN;
  size_t      mateRegionLookupLim = 400;
  bool        noBadAvoid          = true;

  // Mate-region lookup tuning (runtime overrides for compile-time defaults)
  int         minMateMAPQ          = MIN_MATE_MAPQ;          // --min-mate-mapq
  int         maxMateChrID         = MAX_MATE_CHR_ID;        // --max-mate-chr (or -1 via --non-human)
  int         mateRegionMinCount   = MATE_REGION_MIN_COUNT;  // --mate-min-count

  // Non-human genome mode: removes all hardcoded human chromosome
  // assumptions (currently just the mate-lookup ChrID > 23 gate and
  // the LearnBamParams chr1-chrY sampling preference).
  // Sets maxMateChrID = -1 (no limit) and logs a banner.
  bool        nonHuman            = false;

  // When true, skip the high-NM salvage path that pulls in reads with
  // NM/len > 0.02 even when they don't pass the normal weird-read rules.
  // Those reads are useful for catching NM-only SVs (no clips/discordance)
  // but substantially increase read count and r2c cost. --no-nm disables
  // the path entirely for maximum efficiency.
  bool noNmSalvage = false;

  // When false (default), reads whose corrected sequence matches the
  // original BAM sequence reuse the input BAM's CIGAR/NM for the
  // native-vs-r2c gate in splitCoverage, skipping the expensive
  // full-reference BWA re-alignment. This is correct when the input
  // BAM was aligned with BWA-MEM (same scoring model as svaba's
  // internal aligner). When the input BAM was aligned with a
  // different aligner (e.g. bowtie2, STAR, minimap2), the native
  // CIGAR/NM may not be comparable — set this to true to force
  // re-alignment of every corrected read against the reference.
  bool alwaysRealignCorrected = false;

  // external DBs
  std::vector<std::string> blacklistFile;
  // optional non-filtering annotation track(s): labeled BED (RepeatMasker /
  // SegDup / custom). Each breakend is overlap-annotated into the bps
  // `repeat_anno` column. Repeatable. Does NOT filter calls.
  std::vector<std::string> annotationFile;
  std::string germlineSvFile;
  std::string dbsnpVcf;

  // BAM map (idpath)
  // e.g. t001, /path/to/
  std::map<std::string,std::string> bams;

  // ----------------------------------------------------
  // Print help/usage
  static void printUsage();

  void printLogger(SvabaLogger& logger) const;  

  // Parse argc/argv into an SvabaOptions; throws on error
  static SvabaOptions parse(int argc, char** argv);

};
