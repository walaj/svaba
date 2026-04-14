#include "svabaModels.h"

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
    
// =========================================================================
    // GERM_shared: pooled-MLE allele fraction across N+T.
    //
    // We do NOT gate this branch on raw counts (e.g. `aN >= 2`). In a high-
    // artifact region, 2-3 normal alt reads can be perfectly genuine artifact,
    // and shutting GERM_shared off in that regime is the bug that inflates
    // somlod for things like  aN=2/dN=68  vs  aT=16/dT=151  in repeat tracts.
    //
    // Instead, suppress GERM_shared only when *both* of the following hold:
    //   (a) the pooled MLE f_shared sits below the "germline-plausible" band,
    //       AND
    //   (b) the normal contributes essentially no evidence for any nonzero f,
    //       i.e. fitting the normal's own MLE barely improves on f=0
    //       (normal_evidence < ~1 log10  ==  less than ~10x likelihood gain).
    //
    // If either (a) is false (shared looks germline-y -> could be LOH), or (b)
    // is false (normal does have real alt-supporting evidence -> shared could
    // be artifact-in-both), then GERM_shared is allowed to compete and rescue
    // the call from a spurious somlod.
    //
    // We also let f_shared_hat float to its true MLE (no 0.15 floor). The old
    // `clamp(..., 0.15, ...)` turned the "free-MLE" branch into a "germline-
    // only-MLE" branch, which is exactly what we don't want when explaining
    // low-VAF shared signal.
    // =========================================================================

    // Pooled MLE across both samples, free to land anywhere in (0,1).
    const double f_shared_hat_raw = (aN + aT) / (dN_safe + dT_safe);
    const double f_shared_hat     = std::clamp(f_shared_hat_raw, epsilon, 1.0 - epsilon);

    // Threshold below which an AF is not plausibly germline (het ~0.5,
    // hom ~1.0, with some slack for LOH-driven shifts in tumor). This is a
    // statement about *what germline AFs look like*, not about counts.
    const double germ_plausible_floor = 0.15;
    const bool   shared_in_germline_band = (f_shared_hat >= germ_plausible_floor);

    // Normal-side MLE and "evidence for any nonzero f" in log10. This is the
    // error-rate-aware replacement for the old `aN >= 2` count gate -- in a
    // high-artifact region (large eN_fwd), normal_evidence stays small even
    // for aN=2 or 3, because those reads are well explained as errors.
    const double f_n_mle = std::clamp(aN / dN_safe, epsilon, 1.0 - epsilon);
    const double normal_evidence =
        LL_N(dN_safe, aN, f_n_mle) - LL_N(dN_safe, aN, 0.0);
    const double normal_evidence_threshold = 1.0;  // log10  (~10x likelihood)
    const bool   normal_has_real_signal = (normal_evidence >= normal_evidence_threshold);

    double ll_germ_shared = -std::numeric_limits<double>::infinity();
    if (shared_in_germline_band || normal_has_real_signal) {
        ll_germ_shared = LL_N(dN_safe, aN, f_shared_hat)
                       + LL_T(dT_safe, aT, f_shared_hat);
    }
    // Else: pure low-AF somatic with no real normal signal -> GERM_shared
    // stays at -inf and the somatic branch wins on its own merits, just as
    // intended for clean somatic calls.
    
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
