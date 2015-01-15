///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// FMMergeProcess - Merge unambiguously overlapping sequences
//
#ifndef FMMERGEPROCESS_H
#define FMMERGEPROCESS_H

#include "Util.h"
#include "OverlapAlgorithm.h"
#include "SequenceProcessFramework.h"
#include "BitVector.h"
#include "Bigraph.h"
#include "SGUtil.h"

struct FMMergeResult
{
    std::vector<std::string> mergedSequences;
    std::vector<BWTInterval> usedIntervals;
    bool isMerged;
};

// A merge candidate is a read that is a unique extension
// from an existing read. If it has a single edge
// in the direction pointing back to the existing read
// it can be merged with that read. if it has multiple
// edges back to the read it cannot and will be marked as
// invalid.
// 
struct FMMergeCandidate
{
    Vertex* pVertex;
    Edge* pEdge; // Edge from the existing vertex to pVertex
    BWTInterval interval; // interval containing this candidate
};
typedef std::queue<FMMergeCandidate> FMMergeQueue;

// Compute the overlap blocks for reads
class FMMergeProcess
{
    public:
        FMMergeProcess(const OverlapAlgorithm* pOverlapper, 
                       int minOverlap, BitVector* pMarkedReads);

        ~FMMergeProcess();

        FMMergeResult process(const SequenceWorkItem& item);
    
    private:

        // Add the edges of X as candidates to the graph
        void addCandidates(StringGraph* pGraph, 
                           const Vertex* pX, 
                           const Edge* pEdgeToX, 
                           const OverlapBlockList* pBlockList, 
                           FMMergeQueue& candidateQueue);

        // Check whether the candidate can be merged into the current graph
        bool checkCandidate(const FMMergeCandidate& candidate, 
                            const OverlapBlockList* pBlockList) const;

        //
        std::string makeVertexID(BWTInterval interval);

        const OverlapAlgorithm* m_pOverlapper;
        const int m_minOverlap;
        BitVector* m_pMarkedReads;
};

// Write the results from the overlap step to an ASQG file
class FMMergePostProcess
{
    public:
        FMMergePostProcess(std::ostream* pWriter, BitVector* pMarkedReads);
        ~FMMergePostProcess();
        
        void process(const SequenceWorkItem& item, const FMMergeResult& result);

    private:
        size_t m_numMerged;
        size_t m_numTotal;
        size_t m_totalLength;

        std::ostream* m_pWriter;
        BitVector* m_pMarkedReads;
};

#endif
