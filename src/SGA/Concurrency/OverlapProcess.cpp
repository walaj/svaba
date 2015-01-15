///-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapProcess - Wrapper for the overlap computation
//
#include "OverlapProcess.h"

//
//
//
OverlapProcess::OverlapProcess(const std::string& outFile, 
                               const OverlapAlgorithm* pOverlapper, 
                               int minOverlap) : m_pOverlapper(pOverlapper), 
                                                 m_minOverlap(minOverlap)
{
    m_pWriter = createWriter(outFile);
}

//JEREMIAH
OverlapProcess::OverlapProcess(std::stringstream& stream, 
                               const OverlapAlgorithm* pOverlapper, 
                               int minOverlap) : m_pOverlapper(pOverlapper), 
                                                 m_minOverlap(minOverlap)
{
  stream << "TRUE\n";
  //sstream = stream;
}


//
OverlapProcess::~OverlapProcess()
{
    delete m_pWriter;
}

//
OverlapResult OverlapProcess::process(const SequenceWorkItem& workItem)
{
    OverlapResult result = m_pOverlapper->overlapRead(workItem.read, m_minOverlap, &m_blockList);
    m_pOverlapper->writeOverlapBlocks(*m_pWriter, workItem.idx, result.isSubstring, &m_blockList);
    //m_pOverlapper->writeOverlapBlocks(*sstream, workItem.idx, result.isSubstring, &m_blockList);
    m_blockList.clear();
    return result;
}

//
//
//
OverlapPostProcess::OverlapPostProcess(std::ostream* pASQGWriter, 
                                       const OverlapAlgorithm* pOverlapper) : m_pASQGWriter(pASQGWriter),
                                                                              m_pOverlapper(pOverlapper)
{

}

//
void OverlapPostProcess::process(const SequenceWorkItem& item, const OverlapResult& result)
{
    m_pOverlapper->writeResultASQG(*m_pASQGWriter, item.read, result);
}
