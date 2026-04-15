// ContigAlignmentScore.h
//
// Composite quantitative reliability score for a single contig-to-genome
// alignment. The idea is: once a contig has been aligned back to the
// reference (via BWAWrapper in SvabaRegionProcessor), we want a single
// float in [0, 1] that summarizes "is this a clean alignment I should
// trust" vs "is this a wonky multi-indel pile that BWA forced into a
// repeat region".
//
// The score flows through svaba in three places:
//   1. Written to the contig BAM record as tag `zc:i` (milli-confidence,
//      integer in [0, 1000] so it fits SeqLib's AddIntTag API).
//   2. Mirrored onto the BreakEnd struct (field `contig_conf`) inside
//      BreakEnd::transferContigAlignmentData so every BreakPoint carries
//      both sides' scores forward.
//   3. Emitted as the `contig_conf` column in bps.txt.gz and used as a
//      co-gate in the PASS/filter logic.

#pragma once

#include "SeqLib/BamRecord.h"
#include <string>

namespace svaba {

// Raw components + the final composite confidence.
struct ContigAlignScore {
  // CIGAR-derived
  int n_indel_ops    = 0;   // count of I + D ops
  int max_indel_bp   = 0;   // largest single I or D length
  int total_indel_bp = 0;   // sum of I + D bp
  int clip_bp        = 0;   // total soft-clip bp (S ops)
  int aligned_bp     = 0;   // query bp in M/I/=/X (excl. clips)
  int ref_span       = 0;   // reference bp in M/D/=/X/N

  // Tag-derived (may be missing; -1 = not present)
  int NM             = -1;
  int AS             = -1;
  int XS             = -1;

  // Derived metrics
  int    mm          = 0;   // ~ NM - total_indel_bp (BWA convention)
  double mm_rate     = 0.0; // mm / aligned_bp
  double as_per_bp   = 0.0; // AS / max(aligned_bp, ref_span)
  int    as_xs_gap   = 0;   // AS - XS (0 if XS absent)

  // Composite confidence in [0, 1]. 1.0 = pristine single-event alignment.
  double confidence  = 1.0;

  // First rule that pushed the score down; empty means nothing tripped.
  // Purely for logging/debug; not used by downstream filters.
  std::string reason;
};

// Pure function: compute the score from a contig alignment record.
// Does not mutate the record.
ContigAlignScore scoreContigAlignment(const SeqLib::BamRecord& r);

// Write the score back onto the record as integer tag `zc` in
// milli-confidence units (round(1000 * confidence), clamped to [0, 1000]).
// Integer tag keeps us within SeqLib's existing AddIntTag API.
// Also writes `zi` (n_indel_ops), `zm` (max_indel_bp), `zr` (mm_rate*1000)
// for diagnostic visibility in IGV / samtools view.
void tagContigAlignment(SeqLib::BamRecord& r, const ContigAlignScore& s);

// Inverse of tagContigAlignment's zc: pull the milli-confidence off a
// record and return it as a double in [0, 1]. Returns 1.0 if tag absent
// (so legacy records default to "trusted").
double readContigConfTag(const SeqLib::BamRecord& r);

// Default PASS threshold. Anything with contig_conf below this is
// demoted from PASS. 0.50 is a reasonable default; tune via opts later.
inline constexpr double kContigConfPassThreshold = 0.50;

// ---------------------------------------------------------------------------
// Per-read SV-focused alignment score
// ---------------------------------------------------------------------------
//
// Used by BreakPoint::splitCoverage to decide whether a read's r2c
// (read-to-contig) alignment is genuinely "better than" the read's
// corrected-native (read-to-reference) alignment. Only when r2c is
// strictly better do we credit the read as variant-supporting. Both sides
// must use the SAME corrected read sequence and the SAME BWA parameters
// (svaba's internal aligner) — see svabaRead::corrected_native_cig.
//
// Scoring philosophy: an indel-bearing alignment should NOT be allowed to
// approach the same score as a clean ungapped one. We therefore use a
// BWA-MEM-style gap penalty (gap-open dominant, modest gap-extend) that
// is much heavier than the per-bp mismatch cost, so a single indel op
// definitively pushes the score below a clean alignment of the same
// length. Mismatches are penalized at a smaller per-bp rate. Soft-clips
// are penalized 1:1 with matched bp (a clipped base contributes nothing
// AND costs as much as it would have earned, mirroring the intuition that
// a clip is "wasted" query sequence).
//
// Penalty schedule (units = "effective matched bp"; loosely modeled on
// BWA-MEM defaults match=+1 / mm=-4 / gapO=-6 / gapE=-1):
//
//   score = aligned_bp                                 (M / I / = / X query bp)
//         - kReadClipPenaltyPerBp       * clip_bp      (S query bp)
//         - kReadGapOpenPenalty         * indel_ops    (count of I + D ops)
//         - kReadGapExtendPenaltyPerBp  * indel_bp     (sum of I + D bp)
//         - kReadMismatchPenaltyPerBp   * mismatches   (NM - indel_bp; >=0)
//
// Worked examples (k_clip=1, k_gapO=6, k_gapE=1, k_mm=4):
//   150M   NM=0    -> 150                              clean baseline
//   80M70S NM=0    ->  80 - 70                  =  10  big clip
//   90M60S NM=0    ->  90 - 60                  =  30
//   38M3I109M NM=3 -> 150 - 6 - 3               = 141  one 3bp indel,
//                                                       no mismatches
//   75M5I70M NM=5  -> 150 - 6 - 5               = 139
//   75M5D75M NM=5  -> 150 - 6 - 5               = 139  D doesn't consume
//                                                       query but NM counts
//   150M   NM=4    -> 150 - 4*4                 = 134  4 mismatches, no
//                                                       indels (still well
//                                                       above an indel-
//                                                       bearing alignment
//                                                       of equal length —
//                                                       that's the point)
//
// Decision rule in splitCoverage: accept iff
//   r2cScore > correctedNativeScore (apples-to-apples comparison).
// Strict ">" so ties don't credit the read.
//
// `nm` may be -1 (tag absent); the mismatch term is dropped in that case
// but the CIGAR-derived gap-open / gap-extend / clip terms still apply.
inline constexpr double kReadClipPenaltyPerBp       = 1.0;
inline constexpr double kReadGapOpenPenalty         = 6.0;
inline constexpr double kReadGapExtendPenaltyPerBp  = 1.0;
inline constexpr double kReadMismatchPenaltyPerBp   = 4.0;

inline double readAlignmentScore(const SeqLib::Cigar& cig, int nm) {
  int aligned_bp = 0;  // M / I / = / X (consumes query)
  int clip_bp    = 0;  // S
  int indel_ops  = 0;  // count of I + D ops
  int indel_bp   = 0;  // sum   of I + D bp
  for (const auto& c : cig) {
    const int len = static_cast<int>(c.Length());
    switch (c.Type()) {
      case 'M': case '=': case 'X':
        aligned_bp += len;
        break;
      case 'I':
        aligned_bp += len;       // I consumes query
        indel_ops  += 1;
        indel_bp   += len;
        break;
      case 'D':
        indel_ops  += 1;
        indel_bp   += len;       // D doesn't consume query
        break;
      case 'S':
        clip_bp    += len;
        break;
      // 'N', 'H', 'P' contribute nothing.
      default: break;
    }
  }
  double score = static_cast<double>(aligned_bp)
               - kReadClipPenaltyPerBp      * clip_bp
               - kReadGapOpenPenalty        * indel_ops
               - kReadGapExtendPenaltyPerBp * indel_bp;
  if (nm >= 0) {
    // BWA convention: NM = mismatches + inserted_bp + deleted_bp. Subtract
    // out indel_bp to recover the pure mismatch count; clamp at 0 in case
    // CIGAR and NM disagree (e.g. NM tag absent / stale).
    const int mm = std::max(0, nm - indel_bp);
    score -= kReadMismatchPenaltyPerBp * mm;
  }
  return score;
}

} // namespace svaba
