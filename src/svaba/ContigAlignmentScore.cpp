// ContigAlignmentScore.cpp
//
// See header for what this does and where it plugs into svaba. The scoring
// rules below map directly to the filter analysis:
//
//   Rule A: n_indel_ops >= 3 is almost certainly mis-alignment in a repeat.
//           Each additional indel beyond 1 adds a 0.20 penalty.
//   Rule B: 2+ indels with no single big event (<=20 bp) and noisy flanks
//           (mm_rate > 0.02) is scattered-small-indel pattern. Flat 0.25
//           penalty. A single large indel exempts the contig from this
//           rule, so real big SVs/indels pass.
//   Rule C: non-indel mismatch rate. Threshold 1% (well above per-base
//           sequencing error ~0.1-0.5%), slope 10 so 3% mm = 0.20,
//           5% = 0.40, 10% = 0.90. Extra quadratic kicker above 5% so
//           very noisy alignments crash to zero.
//   Rule D: as_per_bp < 0.70 means BWA paid a lot of penalties for this
//           placement. Smooth penalty scaled 2x.
//
// confidence = clamp(1 - sum(penalties), 0, 1).

#include "ContigAlignmentScore.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>

namespace svaba {

namespace {

template <class T>
inline T clip_range(T v, T lo, T hi) {
  return std::max(lo, std::min(hi, v));
}

} // anon

ContigAlignScore scoreContigAlignment(const SeqLib::BamRecord& r) {
  ContigAlignScore s;

  if (!r.MappedFlag()) return s;
  
  // --- parse CIGAR --------------------------------------------------------
  const SeqLib::Cigar cig = r.GetCigar();
  for (const auto& c : cig) {
    const int len = static_cast<int>(c.Length());
    const char op = c.Type();
    switch (op) {
      case 'M': case '=': case 'X':
        s.aligned_bp += len;
        s.ref_span   += len;
        break;
      case 'I':
        s.n_indel_ops++;
        s.total_indel_bp += len;
        if (len > s.max_indel_bp) s.max_indel_bp = len;
        s.aligned_bp += len;
        break;
      case 'D':
        s.n_indel_ops++;
        s.total_indel_bp += len;
        if (len > s.max_indel_bp) s.max_indel_bp = len;
        s.ref_span += len;
        break;
      case 'N':
        s.ref_span += len;
        break;
      case 'S':
        s.clip_bp += len;
        break;
      // 'H' and 'P' contribute nothing
      default: break;
    }
  }

  // --- pull NM / AS / XS --------------------------------------------------
  int32_t v = 0;
  if (r.GetIntTag("NM", v)) s.NM = v;
  if (r.GetIntTag("AS", v)) s.AS = v;
  if (r.GetIntTag("XS", v)) s.XS = v;

  if (s.NM >= 0) {
    // NM in BWA includes mismatches + inserted + deleted bp
    s.mm = std::max(0, s.NM - s.total_indel_bp);
  }
  if (s.aligned_bp > 0 && s.NM >= 0) {
    s.mm_rate = static_cast<double>(s.mm) / s.aligned_bp;
  }
  if (s.AS >= 0) {
    const int denom = std::max(s.aligned_bp, s.ref_span);
    if (denom > 0) s.as_per_bp = static_cast<double>(s.AS) / denom;
  }
  if (s.AS >= 0 && s.XS >= 0) {
    s.as_xs_gap = s.AS - s.XS;
  }

  // --- composite penalty --------------------------------------------------
  double pen = 0.0;

  // Rule A: multi-indel
  if (s.n_indel_ops >= 2) {
    pen += 0.20 * (s.n_indel_ops - 1);
    if (s.reason.empty())
      s.reason = "multi_indel(" + std::to_string(s.n_indel_ops) + ")";
  }

  // Rule B: scattered small indels with noisy flanks. Threshold lowered to
  // 0.02 to match the new mm-rate sensitivity in Rule C.
  if (s.n_indel_ops >= 2 && s.max_indel_bp <= 20 && s.mm_rate > 0.02) {
    pen += 0.25;
    if (s.reason.empty()) s.reason = "small_indel_cluster_noisy";
  }

  // Rule C: non-indel mismatch rate. A clean sequencing + assembly pipeline
  // should have ~0.1-0.5% per-base error; 1% is already suspicious, 3% is
  // almost always mis-assembly or mis-alignment into a repeat/paralog.
  //
  //   linear slope 10 above 0.01:  1% -> 0, 2% -> 0.10, 3% -> 0.20,
  //                                5% -> 0.40, 10% -> 0.90
  //   plus quadratic above 0.05:   5% -> +0, 10% -> +0.075, 15% -> +0.30,
  //                                20% -> +0.675
  //
  // Real indels and soft clips are excluded from mm (NM minus indel bp), so
  // this rule only fires on genuine base mismatches.
  if (s.mm_rate > 0.01) {
    pen += 10.0 * (s.mm_rate - 0.01);
    if (s.mm_rate > 0.05) {
      pen += 30.0 * (s.mm_rate - 0.05) * (s.mm_rate - 0.05);
    }
    if (s.reason.empty() && s.mm_rate > 0.02) s.reason = "high_mm_rate";
  }

  // Rule D: low AS-density (only if AS was actually present)
  if (s.AS >= 0 && s.as_per_bp < 0.70) {
    pen += 2.0 * (0.70 - s.as_per_bp);
    if (s.reason.empty() && s.as_per_bp < 0.60) s.reason = "low_AS_density";
  }

  s.confidence = clip_range(1.0 - pen, 0.0, 1.0);
  return s;
}

void tagContigAlignment(SeqLib::BamRecord& r, const ContigAlignScore& s) {
  // zc:i  = round(1000 * confidence), clamped to [0, 1000]
  // zi:i  = n_indel_ops
  // zm:i  = max single-indel length
  // zr:i  = round(1000 * mm_rate), clamped to [0, 1000]

  if (!r.MappedFlag()) return;   // don't touch aux on unmapped/empty records
   
  const int32_t zc_val =
    clip_range(static_cast<int32_t>(std::lround(1000.0 * s.confidence)),
               int32_t{0}, int32_t{1000});
  const int32_t zr_val =
    clip_range(static_cast<int32_t>(std::lround(1000.0 * s.mm_rate)),
               int32_t{0}, int32_t{1000});

  r.AddIntTag("zc", zc_val);
  r.AddIntTag("zi", s.n_indel_ops);
  r.AddIntTag("zm", s.max_indel_bp);
  r.AddIntTag("zr", zr_val);
}

double readContigConfTag(const SeqLib::BamRecord& r) {
  int32_t zc = -1;
  if (!r.GetIntTag("zc", zc)) return -1.0;   // tag truly absent
  return std::clamp(zc / 1000.0, 0.0, 1.0);
}
 
} // namespace svaba
