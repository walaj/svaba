///----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// PairedDeBruijnHaplotypeBuilder - Build haplotypes
// using a de Bruijn graph and read pair constraints
//
#ifndef PAIRED_DEBRUIJN_HAPLOTYPE_BUILDER_H
#define PAIRED_DEBRUIJN_HAPLOTYPE_BUILDER_H

#include "BWT.h"
#include "SampledSuffixArray.h"
#include "BWTInterval.h"
#include "BWTIntervalCache.h"
#include "SGUtil.h"
#include "SGWalk.h"
#include "VariationBuilderCommon.h"
#include "HaplotypeBuilder.h"
#include "multiple_alignment.h"
#include "GraphCompare.h"
#include "ErrorCorrectProcess.h"
#include "SGWalk.h"
#include "DBGPathGuide.h"
#include <queue>

// Build haplotypes starting from a given sequence.
class PairedDeBruijnHaplotypeBuilder
{
    public:

        PairedDeBruijnHaplotypeBuilder(const GraphCompareParameters& parameters);
        ~PairedDeBruijnHaplotypeBuilder();
        
        // Set the string to start from
        void setInitialHaplotype(const std::string& str);

        // Run the bubble construction process
        // Returns true if the graph was successfully built between the two sequences
        HaplotypeBuilderReturnCode run(StringVector& out_haplotypes);
        
    private:
 
        // Get the IDs of reads containing this kmer       
        std::vector<size_t> getReadIDs(const std::string& kmer) const;

        // 
        void selectGuideAndTargetKmers(const std::string& kmer, bool is_kmer_reversed, DBGPathGuide& guide, std::set<std::string>& target_set) const;

        //
        // Data
        //
        GraphCompareParameters m_parameters;
        std::string m_startingKmer;
};

#endif
