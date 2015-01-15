///-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapProcess - Compute overlap hits for SequenceWorkItems
//
#ifndef OVERLAPPROCESS_H
#define OVERLAPPROCESS_H

#include "Util.h"
#include "OverlapAlgorithm.h"
#include "SequenceProcessFramework.h"

// Compute the overlap blocks for reads
class OverlapProcess
{
    public:
        OverlapProcess(const std::string& outFile, 
                       const OverlapAlgorithm* pOverlapper, 
                       int minOverlap);
        OverlapProcess(std::stringstream& stream,
                       const OverlapAlgorithm* pOverlapper, 
                       int minOverlap);


        ~OverlapProcess();

        OverlapResult process(const SequenceWorkItem& item);
    
    private:
        std::ostream* m_pWriter;
        std::stringstream* sstream;
        OverlapBlockList m_blockList;
        const OverlapAlgorithm* m_pOverlapper;
        const int m_minOverlap;
};

// Write the results from the overlap step to an ASQG file
class OverlapPostProcess
{
    public:
        OverlapPostProcess(std::ostream* pASQGWriter, const OverlapAlgorithm* pOverlapper);
        void process(const SequenceWorkItem& item, const OverlapResult& result);

    private:
        std::ostream* m_pASQGWriter;
        const OverlapAlgorithm* m_pOverlapper;
};

#endif
