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

} // namespace svaba
