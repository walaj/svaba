//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Pileup - An array of strings containing all the seen
// bases for a given region/read
//
#ifndef PILEUP_H
#define PILEUP_H
#include "Util.h"
#include "Alphabet.h"
#include "DNADouble.h"

#define DEFAULT_PROB 0.01
#define DEFAULT_LOG_PROB -4.605170186f

struct PUElem
{
    PUElem() : base('A'), lp(0) {}
    PUElem(char b, double l) : base(b), lp(l) {}

    char base;

    // log probability this base is correct
    double lp; 
};

typedef std::vector<PUElem> PUElemVector;

class Pileup
{
    public:
        Pileup() {}
        Pileup(size_t n) { m_data.reserve(n); }

        char calculateSimpleConsensus() const;
        AlphaCount64 getAlphaCount() const;
        DNADouble calculateSimpleAlphaProb() const;
        DNADouble calculateLikelihoodNoQuality(double p_error) const;

        inline void add(char b) { m_data.push_back(PUElem(b, DEFAULT_LOG_PROB)); }
        inline void add(char b, double lp) { m_data.push_back(PUElem(b, lp)); }

        char getBase(size_t idx) const;
        char getCount(char base) const;
        size_t getDepth() const;

        std::string toStr() const;

    private:
        PUElemVector m_data;
};

typedef std::vector<Pileup> PileupVector;

#endif
