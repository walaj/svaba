///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// GapFillProcess - Fill in gaps in a scaffold
//
#ifndef GAPFILL_PROCESS_H
#define GAPFILL_PROCESS_H

#include <list>
#include <stack>
#include <queue>
#include "BWT.h"
#include "BWTInterval.h"
#include "SGUtil.h"
#include "SGWalk.h"
#include "BitVector.h"
#include "HaplotypeBuilder.h"
#include "SequenceProcessFramework.h"
#include "BWTIntervalCache.h"
#include "SampledSuffixArray.h"

// Structures and typedefs

// Parameters structure
struct GapFillParameters
{
    // BWTS
    const BWT* pBWT; 
    const BWT* pRevBWT;

    const BWTIntervalCache* pBWTCache;
    const BWTIntervalCache* pRevBWTCache;
    
    size_t startKmer;
    size_t endKmer;
    size_t stride;
    size_t kmerThreshold;

    int verbose;
};

// 
struct GapFillResult
{
    std::string scaffold;
};

enum GapFillReturnCode
{
    GFRC_UNKNOWN,
    GFRC_OK,
    GFRC_NO_HAPLOTYPE_PATH,
    GFRC_NO_HAPLOTYPE_ABORTED,
    GFRC_NO_HAPLOTYPE_WALK_FAIL,
    GFRC_NO_ANCHOR,
    GFRC_AMBIGUOUS,
    GFRC_BAD_SIZE,
    GFRC_BAD_TRIM,
    GFRC_NUM_CODES
};

// Statistics tracking object
struct GapFillStats
{
    GapFillStats(); 
    size_t numGapsAttempted;
    size_t numGapsFilled;

    // Failure stats
    size_t numFails[GFRC_NUM_CODES];

    void print() const;
};

//
//
//
class GapFillProcess
{
    public:

        //
        // Functions
        //
        GapFillProcess(const GapFillParameters& params);
        ~GapFillProcess();
        
        // Generate haplotypes from chromosome refName, position [start, end]
        GapFillResult processScaffold(const std::string& scaffold) const;

    private:
        
        //
        // Functions
        //
        
        // Attempt to fill in the sequence between the two anchors
        GapFillReturnCode processGap(size_t k, 
                                     int estimatedSize,
                                     const AnchorSequence& leftAnchor, 
                                     const AnchorSequence& rightAnchor, 
                                     std::string& outSequence) const;

        // Find an anchor sequence to start the process of building the gap sequence
        AnchorSequence findAnchor(size_t k, const std::string& scaffold, int64_t position, bool upstream) const;

        // Attempt to select one of the passed in strings as the gap sequence. If none fit the constraints,
        // this sets gapSequence to the empty string and returns an error code
        GapFillReturnCode selectGapSequence(int estiamtedSize, const StringVector& sequences, std::string& gapSequence) const;

        //
        // Data
        //
        GapFillParameters m_parameters;
        mutable GapFillStats m_stats;
};

#if 0
// 
class GapFillPostProcess
{

    public:
        GapFillPostProcess(const std::string& filename);
        ~GapFillPostProcess();

        //void process(const SequenceWorkItem& item, const GraphCompareResult& result);


    private:
        std::ostream* m_pWriter;
};
#endif

#endif
