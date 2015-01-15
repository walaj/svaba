///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// StatsProcess - Compute statistics about the reads
//
#ifndef STATSPROCESS_H
#define STATSPROCESS_H

#include "Util.h"
#include "BWT.h"
#include "SequenceProcessFramework.h"
#include "SequenceWorkItem.h"
#include "OverlapAlgorithm.h"
#include "KmerDistribution.h"

class StatsResult
{
    public:
        StatsResult() : bases_counted(0), bases_wrong(0), mean_depth(0.0f) {}

        std::vector<int> kmerCoverage;
        int bases_counted;
        int bases_wrong;
        double mean_depth;
};

//
class StatsProcess
{
    public:
        StatsProcess(const BWT* pBWT, const BWT* pRBWT, int kmerLength, int minOverlap, int branchCutoff, bool bNoOverlap);
        ~StatsProcess();
        StatsResult process(const SequenceWorkItem& item);

    private:
        
        const BWT* m_pBWT;
        const BWT* m_pRBWT;
        const int m_kmerLength;
        const int m_minOverlap;
        const int m_branchCutoff;
        const bool m_bNoOverlap;
        static const int m_errorThreshold = 3;

        const OverlapAlgorithm* m_pAllOverlapper;
};

// Write the results from the overlap step to an ASQG file
class StatsPostProcess
{
    public:
        StatsPostProcess(bool bPrintKmer);
        ~StatsPostProcess();

        void process(const SequenceWorkItem& item, const StatsResult& result);

    private:

        KmerDistribution m_kmerDist;
        bool m_bPrintKmer;
        int m_basesCounted;
        int m_basesWrong;
        double m_depthSum;
        int m_numReads;
        int m_numPerfect;
};

#endif
