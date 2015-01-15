//-----------------------------------------------
// Copyright 2009-2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGAStats - Common statistics functions used 
// throughout the program
//

#ifndef SGASTATS_H
#define SGASTATS_H

namespace SGAStats
{
//
// Probability
//
double logPoisson(unsigned int k, double m);
double logFactorial(unsigned int k);
double logBinomial(unsigned int k, unsigned int n, double p);

// The logarithm of the Beta distribution pdf for the special case
// that the shape parameters are positive integers
double logIntegerBetaDistribution(double x, unsigned int a, unsigned int b);

// The logarithm of the integer beta function
double logIntegerBetaFunction(unsigned int a, unsigned int b);

// The logarithm of the integer beta binomail PMF
double logIntegerBetaBinomialDistribution(unsigned int k, unsigned int n, 
                                          unsigned int a, unsigned int b);

}
#endif
