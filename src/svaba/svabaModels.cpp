#include "SvabaModels.h"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstring>

namespace SvabaModels {
  

  double LogLikelihood(double ref,
		       double alt,
		       double f,
		       double e_fwd,
		       double e_back) {
    
    // less negative log-likelihoods means more likely
  // eg for low error rate, odds that you see 5 ALT and 5 REF
  // if you are testing for AF = 0 is going to be very low (eg -40)
  // To test if something is true, we want to test the log-likelihood that
  // it's AF is != 0, so we test LL(ref, alt, AF=0, er). If this is 
  // a large negative number, it means that AF = 0 is very unlikely.
  // If we use a larger error rate, then it is more likely that we will
  // see ALT reads even if true AF = 0, so as ER goes up, then 
  // LL(ref, alt, AF=0, er) becomes less negative.

  // Unnormalized log10-likelihood of seeing 'ref' ref-reads and 'alt' alt-reads
  // under a simple two-state error/mutation model.
    double ll = 0.0;
    
  // Clamp negative counts to zero
  ref = std::max(0.0, ref);
  alt = std::max(0.0, alt);

  // P(read = REF):
  //   (1-f)*(1-e_fwd)    true-ref, no forward error
  // +  f*e_back          true-alt, mis-read back to ref
  double p_ref = (1.0 - f) * (1.0 - e_fwd) // prior arg1
    + f         *       e_back;
  
  // std::cerr << " pref " << p_ref << " ref " << ref << std::endl;
  if (p_ref > 0.0) {
    ll += ref * std::log10(p_ref);
  } else if (ref > 0) {
    return -1e12;  // zero probability but nonzero observations -> -Inf 
  }
  
  // P(read = ALT):
  //   f*(1-e_back)      true-alt, no back mutation
  // + (1-f)*e_fwd       true-ref, forward error to alt
  
  double p_alt = f         * (1.0 - e_back)
    + (1.0 - f) *       e_fwd;
  // std::cerr << " palt " << p_alt << " alt " << alt << std::endl;    
  if (p_alt > 0.0) {
    ll += alt * std::log10(p_alt);
  } else if (alt > 0) {
    return -1e12;
  }  
  
  return ll;
  
  }

// Forward: keep your public signature e_fwd is the high artifact rate.
static inline double SomaticLOD_withSplitErrors(double aN, double dN, double aT, double dT,
                                                double e_art_fwd, double e_rev,
                                                double eN_fwd,   double eT_fwd);

double SomaticLOD(double aN, double dN, double aT, double dT,
                  double e_fwd, double e_rev)
{
    // Map old params to split errors:
    // - e_fwd (caller) is your high artifact rate -> e_art_fwd
    // - Choose conservative, lower empirical forward error caps for normal/tumor
    const double e_art_fwd = e_fwd;

    // You can tune these caps globally if you like
    const double eN_fwd = std::min(e_art_fwd, 0.005);          // 0.5% in NORMAL
    const double eT_fwd = std::min(std::max(e_art_fwd, 1e-4), 0.02); // 2% in TUMOR

    return SomaticLOD_withSplitErrors(aN, dN, aT, dT, e_art_fwd, e_rev, eN_fwd, eT_fwd);
}

// Implementation with split errors
static inline double SomaticLOD_withSplitErrors(double aN, double dN, double aT, double dT,
                                                double e_art_fwd, double e_rev,
                                                double eN_fwd,   double eT_fwd)
{
    const double epsilon = 1e-6;
    const double dN_safe = std::max(1.0, dN);
    const double dT_safe = std::max(1.0, dT);

    // MLE tumor AF under H1
    double fT_hat_raw = aT / dT_safe;
    double fT_hat     = std::clamp(fT_hat_raw, epsilon, 1.0 - epsilon);

    // Likelihood helpers (use *different* forward error by sample)
    auto LL_N = [&](double d, double a, double f) {
        return LogLikelihood(d - a, a, f, eN_fwd, e_rev); // log10
    };
    auto LL_T = [&](double d, double a, double f) {
        return LogLikelihood(d - a, a, f, eT_fwd, e_rev); // log10
    };
    // Artifact-branch helper: use the generous artifact forward rate
    auto LL_art = [&](double d, double a, double f) {
        return LogLikelihood(d - a, a, f, e_art_fwd, e_rev); // log10
    };

    // =========================
    // Artifact allele fraction
    // =========================
    // Keep artifact AF small; rate of errors (e_art_fwd) can be high, but true AF under artifact is low.
    const double f_art = std::clamp(e_art_fwd, 1e-5, 1e-2);

    // =========================
    // SOMATIC (H1): max of
    //   (A) True somatic: N~0, T~fT_hat
    //   (B) Artifact in both: N~f_art, T~f_art
    // =========================
    double llN_som_true = LL_N(dN_safe, aN, 0.0);
    double llT_som_true = LL_T(dT_safe, aT, fT_hat);
    double ll_som_true  = llN_som_true + llT_som_true;

    double llN_som_art  = LL_art(dN_safe, aN, f_art);
    double llT_som_art  = LL_art(dT_safe, aT, f_art);
    double ll_som_art   = llN_som_art + llT_som_art;

    double ll_som = std::max(ll_som_true, ll_som_art);

    // =========================
    // GERMLINE (H0): max of
    //   Het (0.5), Hom-alt (1-eps), Shared MLE f, Artifact (low f)
    // =========================
    const double f_het = std::clamp(0.5, epsilon, 1.0 - epsilon);
    const double f_hom = std::clamp(1.0, epsilon, 1.0 - epsilon);
    // Upper edge of "looks like a low-level artifact" VAF band. Anything
    // above this is either heterozygous (~0.5), homozygous (~1.0), or
    // an LOH-shifted germline event -- not a low-VAF shared artifact.
    const double germ_plausible_floor = 0.15;
    
// =========================================================================
    // GERM_shared: pooled-MLE allele fraction across N+T.
    //
    // Design goals:
    //   - aN = 0 (clean normal)  -> GERM_shared must be held back so that a
    //                               low-VAF tumor signal can win as somatic.
    //   - aN = 1-2 in a clean normal  -> GERM_shared should compete strongly,
    //                                    because a shared low-VAF artifact is
    //                                    exactly what that data looks like.
    //   - aN >= ~3  -> GERM_shared dominates normally; shared artifact or
    //                  real (possibly LOH) germline event explains the data.
    //
    // We use TWO shaping knobs on top of the raw pooled-MLE log-likelihood:
    //
    //   (1) "Cleanliness" penalty:  a decaying function of aN relative to
    //       expected errors in the normal (dN * eN_fwd). Near-maximum penalty
    //       when the normal looks clean (aN at or below the error floor),
    //       drops off steeply as excess normal alt reads accumulate.
    //
    //   (2) VAF similarity bonus:  when the normal does have real alt
    //       evidence, *matching* VAFs between normal and tumor are a strong
    //       "shared" fingerprint. Gated by (1 - cleanliness) so this bonus
    //       only fires once the normal actually has signal -- this prevents
    //       accidentally penalizing low-VAF-tumor somatic calls where the
    //       VAF-ratio would be noisy and meaningless.
    //
    // f_shared_hat floats to the true pooled MLE (no 0.15 floor). The old
    // clamp turned the "free-MLE" branch into a "germline-only-MLE" branch,
    // which is exactly what we don't want when explaining low-VAF shared
    // signal.
    // =========================================================================

    // Pooled MLE across both samples, free to land anywhere in (0,1).
    const double f_shared_hat_raw = (aN + aT) / (dN_safe + dT_safe);
    const double f_shared_hat     = std::clamp(f_shared_hat_raw, epsilon, 1.0 - epsilon);

    // Start with the raw pooled-MLE log-likelihood -- never -inf, always let
    // GERM_shared compete; shaping below decides how strongly.
    double ll_germ_shared = LL_N(dN_safe, aN, f_shared_hat)
                          + LL_T(dT_safe, aT, f_shared_hat);

    // ----- (1) cleanliness penalty on GERM_shared -----------------------------
    // Expected normal alt count under pure forward error at the (capped)
    // normal error rate. Anything at or below this is indistinguishable from
    // a clean normal.
    const double expected_errors_N = dN_safe * eN_fwd;
    const double excess_alt_N      = std::max(0.0, aN - expected_errors_N);
    // Steep exponential decay: coefficient=3 means a single excess alt read
    // already drops cleanliness by ~20x. At dN=79, eN_fwd=0.005:
    //   aN=0 -> cleanliness=1.00  (full hold-back, clean somatic wins)
    //   aN=1 -> cleanliness~0.16  (penalty already mostly gone)
    //   aN=2 -> cleanliness~0.008 (essentially no hold-back)
    // This is the "even 1-2 normal alt reads should steeply consider shared"
    // behavior from the design discussion.
    const double kCleanDecay   = 3.0;
    const double cleanliness   = std::exp(-kCleanDecay * excess_alt_N);

    // Max hold-back magnitude at aN=0.
    const double kCleanPenalty = 2.0;
    ll_germ_shared -= kCleanPenalty * cleanliness;

    // ----- (2) VAF similarity bonus --------------------------------------------
    // Matching tumor/normal VAFs are the fingerprint of a shared artifact or
    // a shared (possibly LOH-shifted) germline event. Gated by (1-cleanliness)
    // so this bonus only really fires once the normal has alt evidence --
    // otherwise the ratio is dominated by numerical epsilons.
    const double f_n_hat = std::clamp(aN / dN_safe, epsilon, 1.0 - epsilon);
    const double f_t_hat = fT_hat;
    const double ratio   = (std::max(f_n_hat, f_t_hat) > 0.0)
                         ? std::min(f_n_hat, f_t_hat) / std::max(f_n_hat, f_t_hat)
                         : 0.0;
    const double sim_weight = 1.0 - cleanliness;     // ~0 at aN=0, ~1 at aN>=2
    const double kSimBonus  = 3.0;                   // log10
    ll_germ_shared += sim_weight * ratio * kSimBonus;

    // ----- (3) "both samples low VAF" shared-artifact bonus --------------------
    // If both normal and tumor carry nonzero alt evidence AND both VAFs sit
    // below the germline-plausible band, that is the canonical signature of a
    // low-level shared artifact (repeat tract, low-complexity region, etc).
    // We add a bonus proportional to how similar the two VAFs are. This is
    // independent of (2) -- (2) rewards *matching*, (3) rewards *both-low*.
    const double low_vaf_cutoff = germ_plausible_floor; // 0.15
    const bool both_nonzero_alt = (aN > 0.0) && (aT > 0.0);
    const bool both_low_vaf     = (f_n_hat < low_vaf_cutoff) &&
                                  (f_t_hat < low_vaf_cutoff);
    if (both_nonzero_alt && both_low_vaf) {
        const double kSharedLowVafBonus = 2.0; // log10
        ll_germ_shared += kSharedLowVafBonus * ratio;
    }
    
    // other model options (germline het, germline homozygous, germline artifact)
    double ll_germ_het    = LL_N(dN_safe, aN, f_het)        + LL_T(dT_safe, aT, f_het);
    double ll_germ_hom    = LL_N(dN_safe, aN, f_hom)        + LL_T(dT_safe, aT, f_hom);
    double ll_germ_art    = LL_art(dN_safe, aN, f_art)      + LL_art(dT_safe, aT, f_art);

    double ll_germ = std::max({ ll_germ_shared, ll_germ_het, ll_germ_hom, ll_germ_art });

    double lod_somatic_vs_germ = ll_som - ll_germ;

    // --------- DEBUG DUMP ----------
    if (false) {
        // SOM/H0 winners (for readability)
        const bool som_true_better = (ll_som_true >= ll_som_art);
        const char* som_branch = som_true_better ? "SOM:true" : "SOM:art";

        double ll_germ_best = ll_germ_het;
        const char* germ_branch = "GERM:het";
        if (ll_germ_hom    > ll_germ_best) { ll_germ_best = ll_germ_hom;    germ_branch = "GERM:hom"; }
        if (ll_germ_art    > ll_germ_best) { ll_germ_best = ll_germ_art;    germ_branch = "GERM:art"; }
        if (ll_germ_shared > ll_germ_best) { ll_germ_best = ll_germ_shared; germ_branch = "GERM:shared"; }

        std::cerr
            << "SomaticLOD"
            << " | aN=" << aN << " dN=" << dN
            << " aT=" << aT << " dT=" << dT
            << " e_art_fwd=" << e_art_fwd
            << " eN_fwd=" << eN_fwd
            << " eT_fwd=" << eT_fwd
            << " e_rev=" << e_rev
            << " eps=" << epsilon
            << " fT_hat_raw=" << fT_hat_raw << " fT_hat=" << fT_hat
            << " f_art=" << f_art
            << " f_shared_raw=" << f_shared_hat_raw << " f_shared=" << f_shared_hat
            << " expN_err=" << expected_errors_N << " excessN=" << excess_alt_N
            << " cleanliness=" << cleanliness
            << " clean_penalty=" << (kCleanPenalty * cleanliness)
            << " f_n_hat=" << f_n_hat << " f_t_hat=" << f_t_hat
            << " ratio=" << ratio << " sim_weight=" << sim_weight
            << " sim_bonus=" << (sim_weight * ratio * kSimBonus)
            << " both_low=" << (both_nonzero_alt && both_low_vaf ? 1 : 0)
            << " low_vaf_bonus=" << ((both_nonzero_alt && both_low_vaf) ? (2.0 * ratio) : 0.0)
            << " || SOM_TRUE(N,T,sum)=(" << llN_som_true << "," << llT_som_true << "," << ll_som_true << ")"
            << " SOM_ART(N,T,sum)=("     << llN_som_art  << "," << llT_som_art  << "," << ll_som_art  << ")"
            << " SOM=" << ll_som << " [" << som_branch << "]"
            << " || GERM_HET(N,T,sum)=("    << (LL_N(dN_safe,aN,f_het))    << "," << (LL_T(dT_safe,aT,f_het))    << "," << ll_germ_het    << ")"
            << " GERM_HOM(N,T,sum)=("       << (LL_N(dN_safe,aN,f_hom))    << "," << (LL_T(dT_safe,aT,f_hom))    << "," << ll_germ_hom    << ")"
            << " GERM_ART(N,T,sum)=("       << (LL_art(dN_safe,aN,f_art))  << "," << (LL_art(dT_safe,aT,f_art))  << "," << ll_germ_art    << ")"
            << " GERM_SHARED(N,T,sum)=("    << (LL_N(dN_safe,aN,f_shared_hat)) << "," << (LL_T(dT_safe,aT,f_shared_hat)) << "," << ll_germ_shared << ")"
            << " GERM=" << ll_germ << " [" << germ_branch << "]"
            << " || LOD10=" << lod_somatic_vs_germ
            << std::endl;
    }
    // -------------------------------

    return lod_somatic_vs_germ; // log10
}
  

  int GenotypeQuality(const std::vector<int>& PLs) {
    int best = std::numeric_limits<int>::max();
    int second = std::numeric_limits<int>::max();
    for (int p : PLs) {
      if (p < best) {
	second = best;
	best   = p;
      }
      else if (p < second) {
	second = p;
      }
    }
    return std::min(second, 99);
  }

  // g is the number of reference alleles (e.g. g = 2 is homozygous reference)
  // assumes biallelic model
  // http://bioinformatics.oxfordjournals.org/content/early/2011/09/08/bioinformatics.btr509.full.pdf+html
  double GenotypeLikelihoods(int g,
			     double er,
			     int alt,
			     int cov)
  {
    // g = number of reference alleles: 2 = hom_ref, 1 = het, 0 = hom_alt
    // er = per-read error rate (0 <= er <= 1)
    // alt = count of alt reads
    // cov = total read depth
    // Returns log10-likelihood under a biallelic model.
    
    // Clamp inputs
    g   = std::min(2, std::max(0, g));
    alt = std::max(0, alt);
    cov = std::max(0, cov);
    er  = std::clamp(er, 0.0, 1.0);
    if (alt > cov) return -1e12;
    
    int ref = cov - alt;
    double norm = -cov * std::log10(2.0);  // accounts for dividing by 2
    
    // Compute numerators (before /2)
    double num_ref = (2 - g) * er + g * (1.0 - er);
    double num_alt = (2 - g) * (1.0 - er) + g * er;
    
    // True probabilities (divide by 2)
    double p_ref = num_ref / 2.0;
    double p_alt = num_alt / 2.0;
    
    double ll = norm;
    
    // Add ref-read contribution
    if (p_ref > 0.0) {
      ll += ref * std::log10(p_ref);
    } else if (ref > 0) {
      return -1e12;
    }
    
    // Add alt-read contribution
    if (p_alt > 0.0) {
      ll += alt * std::log10(p_alt);
    } else if (alt > 0) {
      return -1e12;
    }
    
    return ll;
  }
  
}
