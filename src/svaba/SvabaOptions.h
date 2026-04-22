#pragma once

#include <string>
#include <vector>
#include <map>
#include <cstddef>

// version & date
//
// Bumped to 2.0.0 to mark the SvABA2.0 overhaul: v3 bps.txt schema with
// per-BP bp_id (col 52), r2c.txt.gz structured r2c emission replacing the
// old alignments.txt.gz, comparative split-coverage gate, standalone
// `svaba tovcf` + `svaba postprocess` subcommands, VCFv4.5 output, and
// fermi-lite as the default local-assembly engine. See CLAUDE.md for the
// full set of changes.
inline constexpr char SVABA_VERSION[] = "2.0.0";
inline constexpr char SVABA_DATE[]    = "04/2026";

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

// SvABA2.0 (v3): per-sample-prefix margin on the r2c-better-than-native
// gate. A read is credited as a split-supporter iff
//     r2c_score >  native_score * (1.0 + margin)
// where margin is T_R2C_MIN_MARGIN for tumor-prefix reads (t***) and
// N_R2C_MIN_MARGIN for normal-prefix reads (n***). Tumor is held to
// a stricter standard (10% better) because a somatic call depends on
// clean separation between the samples; normal only needs strictly
// greater so we retain sensitivity for germline/LOH reads that help
// rule out a somatic event. This replaces the older both_split /
// one_split / homlen branching gate — when junction homology is on
// the order of the read length, neither r2c nor native wins
// decisively and such reads fall out here naturally, which is the
// correct conservative behavior.
inline constexpr double T_R2C_MIN_MARGIN          = 0.10;
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

// from svabaBamWalker
inline constexpr int MIN_DSCRD_READS_DSCRD_ONLY          = 6;
inline constexpr int MIN_MAPQ_FOR_MATE_LOOKUP            =     0;
inline constexpr int MIN_ISIZE_FOR_DISCORDANT_REALIGNMENT = 1'000;
inline constexpr int DISC_REALIGN_MATE_PAD                =   100;
inline constexpr int MAX_SECONDARY_HIT_DISC               =    10;
inline constexpr int MATE_REGION_PAD                      =   250;

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
  double sdDiscCutoff       = 3.92;
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
  size_t      mateLookupMin       = 3;
  size_t      mateRegionLookupLim = 400;
  bool        noBadAvoid          = true;

  // external DBs
  std::vector<std::string> blacklistFile;
  std::string germlineSvFile;
  std::string dbsnpVcf;

  // BAM map (idpath)
  // e.g. t001, /path/to/
  std::map<std::string,std::string> bams;

  // ----------------------------------------------------
  // Print help/usage
  static void printUsage();

  void printLogger(SvabaLogger& logger) const;  

  void addFRRule(const std::string &rgName, int N);
  
  // Parse argc/argv into an SvabaOptions; throws on error
  static SvabaOptions parse(int argc, char** argv);

};
