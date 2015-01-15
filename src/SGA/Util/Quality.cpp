//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Quality - functions for manipulating quality values
//
#include "Quality.h"

// Return a uniform log-scaled quality vector of the given size
DoubleVector Quality::uniformLogProbVector(double p_error, size_t n)
{
    double lp = log(p_error);
    DoubleVector dv;
    dv.reserve(n);
    for(size_t i = 0; i < n; ++i)
        dv.push_back(lp);
    return dv;
}

// Return a string which encodes the log-scaled double vector as characters
std::string Quality::encodeLogProbVector(const DoubleVector& dv)
{
    std::string out;
    out.reserve(dv.size());
    for(size_t i = 0; i < dv.size(); ++i)
        out.push_back(phred2char(lnprob2phred(dv[i])));
    return out;
}

