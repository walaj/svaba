///-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// RankProcess - Compute a vector of BWT ranks for
// SequenceWorkItems
//
#include "RankProcess.h"

//
//
//
RankProcess::RankProcess(const BWT* pBWT, 
                         GapArray* pSharedGapArray, 
                         bool doReverse, 
                         bool removeMode) : m_pBWT(pBWT), 
                                            m_pSharedGapArray(pSharedGapArray),
                                            m_doReverse(doReverse), 
                                            m_removeMode(removeMode)
{

}

//
RankProcess::~RankProcess()
{

}

// Calculate the ranks of the given sequence.
// We attempt to update the count of the rank in the gap array 
// using an atomic compare and swap. The update will fail if the
// small-storage maximum count is exceeded. In this case
// we push the value to the overflow array for a serial
// update in the post-processing step.
RankResult RankProcess::process(const SequenceWorkItem& workItem)
{
    RankResult out;
    DNAString w = workItem.read.seq;
    if(m_doReverse)
        w.reverse();

    size_t l = w.length();
    int i = l - 1;

    // In add mode, the initial rank is zero and we calculate the rank
    // for the last base of the sequence using just C(a). In remove
    // mode we use the index of the read (in the original read table) as
    // the rank so that ranks calculate correspond to the correct
    // entries in the BWT for the read to remove.
    int64_t rank = 0; // add mode
    if(m_removeMode)
    {
        // Parse the read index from the read id
        rank = parseRankFromID(workItem.read.id);
    }

    out.numRanksProcessed += 1;
    if(!m_pSharedGapArray->attemptBaseIncrement(rank))
        out.overflowVec.push_back(rank);

    // Compute the starting rank for the last symbol of w
    char c = w.get(i);

    // In the case that the starting rank is zero (default
    // in add mode, or if we are removing the first read)
    // there can no occurrence of any characters before this
    // suffix so we just calculate the rank from C(a)
    if(rank == 0)
        rank = m_pBWT->getPC(c);
    else
        rank = m_pBWT->getPC(c) + m_pBWT->getOcc(c, rank - 1);
    
    out.numRanksProcessed += 1;
    if(!m_pSharedGapArray->attemptBaseIncrement(rank))
        out.overflowVec.push_back(rank);
    --i;

    // Iteratively compute the remaining ranks
    while(i >= 0)
    {
        char c = w.get(i);
        rank = m_pBWT->getPC(c) + m_pBWT->getOcc(c, rank - 1);
        //std::cout << "c: " << c << " rank: " << rank << "\n";
        out.numRanksProcessed += 1;
        if(!m_pSharedGapArray->attemptBaseIncrement(rank))
            out.overflowVec.push_back(rank);
        --i;
    }
    return out;
}

// Parse the rank of a read from its ID. This must be set by the process
// which discards the read.
int64_t RankProcess::parseRankFromID(const std::string& id)
{
    static const std::string rank_str = "seqrank=";
    // Find the position of the rank expression in the string
    size_t rank_pos = id.rfind(rank_str);
    if(rank_pos == std::string::npos)
    {
        std::cout << "Error: rank token not found in the discarded read with id: " << id << "\n";
        exit(EXIT_FAILURE);
    }

    assert(rank_pos + rank_str.length() < id.length());
    std::stringstream rp(id.substr(rank_pos + rank_str.length()));
    int64_t rank;
    rp >> rank;
    return rank;
}
//
//
//
RankPostProcess::RankPostProcess(GapArray* pGapArray) : m_pGapArray(pGapArray), num_strings(0), num_symbols(0), num_serial_updates(0)
{

}

RankPostProcess::~RankPostProcess()
{
    //printf("Debug: num serial updates: %zu\n", num_serial_updates); 
}

//
void RankPostProcess::process(const SequenceWorkItem& /*item*/, const RankResult& result)
{
    ++num_strings;
    num_symbols += result.numRanksProcessed;

    // We update any overflowed ranks here. This call is serial and only updates
    // the Overflow table in the gap array.
    for(RankVector::const_iterator iter = result.overflowVec.begin(); iter != result.overflowVec.end(); ++iter)
        m_pGapArray->incrementOverflowSerial(*iter);
    num_serial_updates += result.overflowVec.size();
}
