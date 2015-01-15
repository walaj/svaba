///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// VCFTester - Run dindel on each variant in a vcf file
//
#ifndef VCFTESTER_H
#define VCFTESTER_H

#include <list>
#include <stack>
#include <queue>
#include "GraphCompare.h"
#include "DindelRealignWindow.h"
#include "DindelUtil.h"

//
//
//
class VCFTester
{
    public:

        //
        // Functions
        //
        VCFTester(const GraphCompareParameters& params);
        ~VCFTester();
        
        // Process a read and all its kmers
        void process(const VCFFile::VCFEntry& record);

    private:
 
        // Functions
        std::string applyVariant(const std::string& in, int pos,
                                 const std::string& ref, const std::string& alt);
        
        //
        // Data
        //
        GraphCompareParameters m_parameters;
        GraphCompare m_graphComparer;

        int m_returnCodes[DRC_NUM_CODES];

        // Output file
        VCFFile m_baseVCFFile;
        VCFFile m_variantVCFFile;
        VCFFile m_callsVCFFile;
};

#endif
