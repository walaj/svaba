///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// HapgenUtil - Utility functions used for 
// working with the haplotype generation
// modules and output
//
#ifndef HAPGENUTIL_H
#define HAPGENUTIL_H

#include "BWTIndexSet.h"
#include "StdAlnTools.h"
#include "overlapper.h"

// Simple alignment object representing
// the placement of a string onto a reference genome
struct HapgenAlignment
{
    //
    HapgenAlignment() : referenceID(-1), position(-1), length(0), score(-1), numEvents(0), isRC(false), isSNP(false) {}
    HapgenAlignment(const int& id, int p, int l, int s, int ne, bool rc, bool snp) : referenceID(id), 
                                                                                     position(p), 
                                                                                     length(l), 
                                                                                     score(s),
                                                                                     numEvents(ne), 
                                                                                     isRC(rc), 
                                                                                     isSNP(snp) {}

    // Sort
    // Kees: added sorting on alignment score
    friend bool operator<(const HapgenAlignment& a, const HapgenAlignment& b)
    {
        if(a.referenceID != b.referenceID)
            return a.referenceID < b.referenceID;
        else if(a.position != b.position)
            return a.position < b.position;
        else
            return a.score < b.score;
    }

    // Output
    friend std::ostream& operator<<(std::ostream& o, const HapgenAlignment& align) 
    { 
        o << "R_ID: " << align.referenceID << 
             " P: [" << align.position << ", " << 
             align.position + align.length << "]" << 
             " length: " << align.length << 
             " score: " << align.score;
        return o; 
    }

    //
    int referenceID;
    int position;
    int length; // alignment length
    int score;
    int numEvents;
    bool isRC;
    bool isSNP;
};
typedef std::vector<HapgenAlignment> HapgenAlignmentVector;

//
namespace HapgenUtil
{
    //
    void alignHaplotypeToReferenceKmer(size_t k,
                                       const std::string& haplotype,
                                       const BWTIndexSet& referenceIndex,
                                       const ReadTable* pReferenceTable,
                                       HapgenAlignmentVector& outAlignments);

    


    // Coalesce a set of alignments into distinct locations
    void coalesceAlignments(HapgenAlignmentVector& alignments);

    // Extract a substring for a reference defined by the alignment
    std::string extractReference(const HapgenAlignment& aln, const ReadTable* pRefTable, int flanking);

    // Extract substrings of a reference described by an alignment
    // This writes to 3 out parameters, one for the defined position
    // of the reference, one for the upstream flanking bases and one for the downstream flanking
    void extractReferenceSubstrings(const HapgenAlignment& aln, 
                                    const ReadTable* pRefTable, 
                                    int flanking,
                                    std::string& outUpstream, 
                                    std::string& outDefined, 
                                    std::string& outDownstream);
    
    // Make haplotypes for the given reference alignment, including flanking sequences
    // A haplotype for the reference is also generated
    bool makeFlankingHaplotypes(const HapgenAlignment& aln, 
                                const ReadTable* pRefTable, 
                                int flanking,
                                const StringVector& inHaplotypes,
                                StringVector& outFlankingHaplotypes,
                                StringVector& outHaplotypes);


    // Extract reads from an FM-index that have a k-mer match to any given haplotypes
    // If the number of reads to extract exceeds maxReads, false is returned
    bool extractHaplotypeReads(const StringVector& haplotypes, 
                               const BWTIndexSet& indices,
                               int k,
                               bool doReverse,
                               size_t maxReads,
                               int64_t maxIntervalSize,
                               SeqRecordVector* pOutReads, 
                               SeqRecordVector* pOutMates);

    // Extract reads from an FM-index that have a k-mer match to AT MOST one haplotype
    // If the number of reads to extract exceeds maxReads, false is returned
    bool extractHaplotypeSpecificReads(const StringVector& haplotypes, 
                                       const BWTIndexSet& indices,
                                       int k,
                                       bool doReverse,
                                       size_t maxReads,
                                       int64_t maxIntervalSize,
                                       SeqRecordVector* pOutReads, 
                                       SeqRecordVector* pOutMates);

    
    // Perform an alignment between a haplotype sequence and a reference sequence
    SequenceOverlap alignHaplotypeToReference(const std::string& reference,
                                              const std::string& haplotype);

    // Compute the best local alignment for each read in the array to the given sequence
    LocalAlignmentResultVector alignReadsLocally(const std::string& target, const SeqItemVector& reads);

    // Check that all the strings in the vector align to the same coordinates
    // of the passed in sequence
    bool checkAlignmentsAreConsistent(const std::string& refString, const StringVector& queries);

    // Calculate the dust score around a reference position
    double calculateDustScoreAtPosition(const std::string& name, 
                                        int position, 
                                        const ReadTable* pRefTable,
                                        int window_size = 64);

    // Return the count of the most frequent string within one-edit distance of the given string 
    size_t getMaximumOneEdit(const std::string& str, const BWTIndexSet& indices);
    
    // Calculate the largest k such that every k-mer in the sequence is present at least min_depth times in the BWT
    size_t calculateMaxCoveringK(const std::string& sequence, int min_depth, const BWTIndexSet& indices);

    // Count the number of times each k-mer of str is in the BWT
    std::vector<int> makeCountProfile(const std::string& str, size_t k, int max, const BWTIndexSet& indices);

    // Print an alignment to a reference
    void printAlignment(const std::string& query, const HapgenAlignment& aln, const ReadTable* pRefTable);

};

#endif
