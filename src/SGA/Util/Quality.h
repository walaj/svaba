//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Quality - functions for manipulating quality values
//
#ifndef QUALITY_H
#define QUALITY_H

#include <inttypes.h>
#include <vector>
#include <math.h>
#include <assert.h>
#include <string>
#include <iostream>

static const int DEFAULT_QUAL_SCORE = 15;
static const int PHRED64_DIFF = 31;

typedef std::vector<double> DoubleVector;
namespace Quality
{
    // Convert the quality character from phred64 to phred33 encoding
    inline char phred64toPhred33(char c)
    {
        return (int)c - PHRED64_DIFF;
    }

    // Returns true if the character c is a valid phred 33 encoding
    // of a quality value in the range [0, 60]
    inline bool isValidPhred33(char c)
    {
        int p = (int)c - 33;
        return p >= 0 && p <= 60;
    }

    // Phred score transformations
    inline int char2phred(char b)
    {
        uint8_t v = b;
        assert(v >= 33);
        return v - 33;
    }

    inline char phred2char(int p)
    {
        uint8_t v = (p <= 93) ? p : 93;
        char c = v + 33;
        return c;
    }

    inline int prob2phred(double p)
    {
        return static_cast<int>(round(-10.0f * log10(p)));
    }

    inline int lnprob2phred(double lp)
    {
        // change base
        static double transform = log(10);
        lp /= transform;
        return static_cast<int>(round(-10.0f * lp));
    }

    // Return a uniform log-scaled quality vector of the given size
    DoubleVector uniformLogProbVector(double p_error, size_t n);

    // Return a string which encodes the log-scaled double vector as characters
    std::string encodeLogProbVector(const DoubleVector& dv);
};


#endif
