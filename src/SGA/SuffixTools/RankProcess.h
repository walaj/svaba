///-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// RankProcess - Compute a vector of BWT ranks for
// SequenceWorkItems
//
#ifndef RANKPROCESS_H
#define RANKPROCESS_H

#include "Util.h"
#include "BWT.h"
#include "SequenceWorkItem.h"
#include "GapArray.h"

typedef std::vector<int64_t> RankVector;
struct RankResult
{
    RankResult() : numRanksProcessed(0) {}

    RankVector overflowVec;
    size_t numRanksProcessed;
};

// Compute the overlap blocks for reads
class RankProcess
{
    public:
        RankProcess(const BWT* pBWT, GapArray* pSharedGapArray, bool doReverse, bool removeMode);
        ~RankProcess();

        RankResult process(const SequenceWorkItem& item);
    
    private:

        int64_t parseRankFromID(const std::string& id);

        const BWT* m_pBWT;
        GapArray* m_pSharedGapArray;

        bool m_doReverse;
        bool m_removeMode;
};

// Update the gap array with 
class RankPostProcess
{
    public:
        RankPostProcess(GapArray* pGapArray);
        ~RankPostProcess();

        void process(const SequenceWorkItem& item, const RankResult& result);
        size_t getNumStringsProcessed() const { return num_strings; }
        size_t getNumSymbolsProcessed() const { return num_symbols; }

    private:
        GapArray* m_pGapArray;
        size_t num_strings;
        size_t num_symbols;
        size_t num_serial_updates;
};

#endif
