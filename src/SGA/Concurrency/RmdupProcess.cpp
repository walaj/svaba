///-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// RmdupProcess - Implements the logic for computing hits
// for the Rmdup program
//
#include "RmdupProcess.h"

//
//
//
RmdupProcess::RmdupProcess(const std::string& outFile, 
                           const OverlapAlgorithm* pOverlapper) : m_pOverlapper(pOverlapper)
{
    m_pWriter = createWriter(outFile);
}

//
RmdupProcess::~RmdupProcess()
{
    delete m_pWriter;
}

//
OverlapResult RmdupProcess::process(const SequenceWorkItem& workItem)
{
    OverlapResult result = m_pOverlapper->alignReadDuplicate(workItem.read, &m_blockList);
    // Write the read sequence and the overlap blocks to the file
    *m_pWriter << workItem.read.id << "\t" << workItem.read.seq.toString() << "\t";
    m_pOverlapper->writeOverlapBlocks(*m_pWriter, workItem.idx, result.isSubstring, &m_blockList);
    m_blockList.clear();
    return result;
}
