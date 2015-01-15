///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ConnectProcess - Algorithm to resolve the sequence
// connecting two ends of a paired end read
//
#include "ConnectProcess.h"
#include "StringGraphGenerator.h"

//
//
//
ConnectProcess::ConnectProcess(const OverlapAlgorithm* pOverlapper, 
                               int minOverlap,
                               int maxDistance) : 
                                                 m_pOverlapper(pOverlapper), 
                                                 m_minOverlap(minOverlap),
                                                 m_maxDistance(maxDistance)
{

}

//
ConnectProcess::~ConnectProcess()
{

}

//
ConnectResult ConnectProcess::process(const SequenceWorkItemPair& workItemPair)
{
    assert(getPairID(workItemPair.first.read.id) == workItemPair.second.read.id);
    ConnectResult result;

    StringGraphGenerator localGraph(m_pOverlapper, workItemPair.first.read, workItemPair.second.read, m_minOverlap, ED_SENSE, m_maxDistance);
    SGWalkVector walks = localGraph.searchWalks();
    //std::cout << "Found " << walks.size() << " walk between " << workItemPair.first.read.id << " and " << workItemPair.second.read.id << "\n";
    if(walks.size() == 1)
    {
        SGWalk& solution = walks.front();
        result.resolvedSequence = solution.getString(SGWT_START_TO_END);
        //std::cout << "Solution has length " << solution.getStartToEndDistance() << "\n";
    }
    return result;
}

//
//
//
ConnectPostProcess::ConnectPostProcess(std::ostream* pWriter,
                                       std::ostream* pUnconnectedWriter) : 
                                                                m_pWriter(pWriter),
                                                                m_pUnconnectedWriter(pUnconnectedWriter),
                                                                m_numPairsTotal(0),
                                                                m_numPairsResolved(0)
{

}

//
ConnectPostProcess::~ConnectPostProcess()
{
    printf("connect: Resolved %d out of %d pairs (%lf)\n", m_numPairsResolved, m_numPairsTotal, (double)m_numPairsResolved / m_numPairsTotal);
}

//
void ConnectPostProcess::process(const SequenceWorkItemPair& workItemPair, const ConnectResult& result)
{
    ++m_numPairsTotal;

    if(result.resolvedSequence.length() > 0)
    {
        ++m_numPairsResolved;

        SeqRecord record;
        record.id = getPairBasename(workItemPair.first.read.id);
        record.seq = result.resolvedSequence;
        record.write(*m_pWriter);
    }
    else
    {
        workItemPair.first.read.write(*m_pUnconnectedWriter);
        workItemPair.second.read.write(*m_pUnconnectedWriter);
    }
}
