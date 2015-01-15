///-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// RmdupProcess - Compute full-length hits for 
// SequenceWorkItems
//
#ifndef RMDUPPROCESS_H
#define RMDUPPROCESS_H

#include "Util.h"
#include "OverlapAlgorithm.h"
#include "SequenceProcessFramework.h"

// Compute the overlap blocks for reads
class RmdupProcess
{
    public:
        RmdupProcess(const std::string& outFile, 
                       const OverlapAlgorithm* pOverlapper);

        ~RmdupProcess();

        OverlapResult process(const SequenceWorkItem& item);
    
    private:
        std::ostream* m_pWriter;
        OverlapBlockList m_blockList;
        const OverlapAlgorithm* m_pOverlapper;
};

// The rmdup process does not have a post-processing step, this just passes-through the data
class RmdupPostProcess
{
    public:
        RmdupPostProcess() {}
        void process(const SequenceWorkItem& /*item*/, const OverlapResult& /*result*/) {}
};

#endif
