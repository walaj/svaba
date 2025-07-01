#pragma once

#include <vector>

namespace SvabaModels {
  
  double LogLikelihood(double ref,
		       double alt,
		       double f,
		       double e_fwd,
		       double e_back);
  
  double SomaticLOD(double aN,
		    double dN,
		    double aT,
		    double dT,
		    double e_fwd,
		    double e_rev);
  
  double GenotypeLikelihoods(int g, double er, int alt, int cov);
  
  int GenotypeQuality(const std::vector<int>& PLs);
}

