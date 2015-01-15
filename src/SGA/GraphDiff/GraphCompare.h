///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// GraphCompare - Compare two (abstractly represented)
// graphs against each other to find strings
// only present in one.
//
// The graphs are abstractly represented as
// an FM-index.
//
#ifndef GRAPH_COMPARE_H
#define GRAPH_COMPARE_H

#include <list>
#include <stack>
#include <queue>
#include "BWTIndexSet.h"
#include "BWTInterval.h"
#include "SGUtil.h"
#include "SGWalk.h"
#include "BitVector.h"
#include "VariationBuilderCommon.h"
#include "SequenceProcessFramework.h"
#include "BWTIntervalCache.h"
#include "SampledSuffixArray.h"
#include "DindelRealignWindow.h"
#include "BloomFilter.h"
#include "api/BamWriter.h"

enum GraphCompareAlgorithm
{
    GCA_STRING_GRAPH,
    GCA_DEBRUIJN_GRAPH,
    GCA_PAIRED_DEBRUIJN_GRAPH
};

// Parameters structure
class GraphCompareAggregateResults;
struct GraphCompareParameters
{
    //
    // Data structures
    //

    // Base reads index
    BWTIndexSet baseIndex;
    
    // Variant reads index
    BWTIndexSet variantIndex;
    
    // Reference genome index
    BWTIndexSet referenceIndex;
    const ReadTable* pRefTable;

    // Bloom filter to mark user kmers
    BloomFilter* pBloomFilter;

    //
    // Parameters
    //

    // Compare against reference genome?
    bool bReferenceMode;

    // Haplotype assembly
    GraphCompareAlgorithm algorithm;
    size_t kmer;
    size_t minDiscoveryCount;
    size_t maxDiscoveryCount;
    size_t minDBGCount;
    int minOverlap;

    // Read extraction parameters
    size_t maxExtractionIntervalSize;
    size_t maxReads;

    // Filtering
    size_t maxHaplotypes;
    
    // Dindel params
    DindelRealignParameters dindelRealignParameters;
};

// Intermediate result
struct GraphBuildResult
{
    StringVector variant_haplotypes;
    StringVector base_haplotypes;
};

// Final result
struct GraphCompareResult
{
    StringVector varStrings;
    StringVector baseStrings;
    DoubleVector varCoverages;
    DoubleVector baseCoverages;

    StringVector baseVCFStrings;
    StringVector variantVCFStrings;
    StringVector calledVCFStrings;

    DindelReadReferenceAlignmentVector projectedReadAlignments;
};

//
// Statistics tracking object
//
struct GraphCompareStats
{
    // functions
    GraphCompareStats() { clear(); }
    void clear();
    void add(const GraphCompareStats& other);
    void print() const;

    // data
    int numBubbles;
    int numAttempted;
    int numTargetBranched;
    int numSourceBranched;
    int numTargetBroken;
    int numSourceBroken;
    int numWalkFailed;
    int numHBFailed;
    int numNoSolution;

    int numInsertions;
    int numDeletions;
    int numSubs;
};

//
// GraphCompare
//
class GraphCompare
{
    public:

        //
        // Functions
        //
        GraphCompare(const GraphCompareParameters& params);
        ~GraphCompare();
        
        // Process a read and all its kmers
        GraphCompareResult process(const SequenceWorkItem& item) const;
        void debug(const std::string& debugFilename);
        void testKmersFromFile(const std::string& kmerFilename);
        void testKmer(const std::string& kmer);

        //
        void updateSharedStats(GraphCompareAggregateResults* pSharedStats);

    private:
        
        //
        // Functions
        //

        // Attempt to assemble a variant kmer into haplotypes
        GraphBuildResult processVariantKmer(const std::string& str, int count) const;
        
        // Perform quality checks on the variant haplotypes
        void qcVariantHaplotypes(bool bReferenceMode, StringVector& variant_haplotypes) const;

        // Build haplotypes in the base sequence that are parallel to the variant haplotypes
        void buildParallelBaseHaplotypes(const StringVector& variant_haplotypes, StringVector& base_haplotypes) const;

        // Mark all the kmers in str as being visited
        void markVariantSequenceKmers(const std::string& str) const;
        
        // Calculate the number of high coverage branches off a haplotype path through the de Bruijn graph
        size_t calculateHaplotypeBranches(const std::string& sequence, size_t k, size_t min_branch_depth, const BWTIndexSet& indices);

        // Update statistics 
        void updateVariationCount(const BubbleResult& result);
        
        // Debug functions
        bool transformVariantString(const std::string& inStr, std::string& outStr);
        void showMappingLocations(const std::string& str);

        //
        // Data
        //
        GraphCompareParameters m_parameters;

        // Results stats
        GraphCompareStats m_stats;

};

// Shared result object that the threaded
// GraphCompare instances write to. Protected
// by a mutex
class GraphCompareAggregateResults
{

    public:
        GraphCompareAggregateResults(const std::string& fileprefix, 
                                     const StringVector& samples, 
                                     const ReadTable& refTable,
                                     const std::string& referenceFilename);

        ~GraphCompareAggregateResults();

        void process(const SequenceWorkItem& item, const GraphCompareResult& result);

        void updateShared(const GraphCompareStats stats);
        void printStats() const;

    private:
        pthread_mutex_t m_mutex;
        GraphCompareStats m_stats;
        
        // Fasta output file
        std::ostream* m_pWriter;
        
        // VCF output files
        VCFFile m_baseVCFFile;
        VCFFile m_variantVCFFile;
        VCFFile m_callsVCFFile;

        // Bam file output
        BamTools::BamWriter m_evidenceBamFile;
        std::map<std::string, size_t> m_refNameToIndexMap;

        size_t m_numVariants;
};

#endif
