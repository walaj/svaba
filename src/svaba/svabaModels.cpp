#include "svabaModels.h"

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

  /*  Somatic Log-Odds
   *
   *  epsilon : small floor on allele fraction to keep f in (0,1) and
   *            avoid log(0).  Typical 1e4 to 1e5.
   */
  double SomaticLOD(double aN,
		    double dN,
		    double aT,
		    double dT,
		    double e_fwd,
		    double e_rev)
  {
    
    /* --- MLE for tumour AF under somatic hypothesis ---------------- *
     *     Clamp to [eps,1-eps] to keep log(f_T) and log(1-f_T) finite.   */
    const double epsilon = 1e-6;
    double fT_hat_raw = aT / std::max(1.0, dT);
    double fT_hat     = std::clamp(fT_hat_raw, epsilon, 1.0 - epsilon);
    
    /* MLE for shared AF under H0 */
    double f_shared_raw = (aN + aT) / std::max(1.0, dN + dT);
    double f_shared     = std::clamp(f_shared_raw, epsilon, 1.0 - epsilon);
    
    /* --- Log10-likelihoods (two-error model) ----------------------- */
    // H1 is that this is true tumor variant
    // Normal under H1: true normal AF = 0, so alt reads arise only from forward error.
    double llN_som_log10 = LogLikelihood(dN - aN, aN, 0.0, e_fwd, e_rev);
    
    // Tumour under H1: true tumor AF = f_T.
    double llT_som_log10 = LogLikelihood(dT - aT, aT, fT_hat, e_fwd, e_rev);
    
    // Normal under H0: shared AF = f
    // H0 is that this is shared (germline / artefact) hypothesis (H0)
    double llN_shared_log10 = LogLikelihood(dN - aN, aN, f_shared, e_fwd, e_rev);
    
    // Tumour under H0: shared AF = f
    double llT_shared_log10 = LogLikelihood(dT - aT, aT, f_shared, e_fwd, e_rev);
    
    /* Convert all log10 values to natural log for the ratio */
    double ll_som    = (llN_som_log10    + llT_som_log10   );
    double ll_shared = (llN_shared_log10 + llT_shared_log10);
    
    /* LOD_som > 0 -->  evidence for somatic */
    return ll_som - ll_shared;
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
