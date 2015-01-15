//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// BWTAlgorithms.h - Algorithms for aligning to a bwt structure
//
#ifndef BWT_ALGORITHMS_H
#define BWT_ALGORITHMS_H

#include "STCommon.h"
#include "BWTIndexSet.h"
#include "BWTInterval.h"
#include "GraphCommon.h"

#include <queue>
#include <list>

#define LEFT_INT_IDX 0
#define RIGHT_INT_IDX 1

// structures

// A (partial) prefix of a string contained in the BWT
// and its lexicographic rank
struct RankedPrefix
{
    size_t rank;
    std::string prefix;
};
typedef std::vector<RankedPrefix> RankedPrefixVector;

// functions
namespace BWTAlgorithms
{

// get the interval(s) in pBWT/pRevBWT that corresponds to the string w using a backward search algorithm
BWTInterval findInterval(const BWT* pBWT, const std::string& w);
BWTInterval findIntervalWithCache(const BWT* pBWT, const BWTIntervalCache* pIntervalCache, const std::string& w);
BWTInterval findInterval(const BWTIndexSet& indices, const std::string& w);

BWTIntervalPair findIntervalPair(const BWT* pBWT, const BWT* pRevBWT, const std::string& w);
BWTIntervalPair findIntervalPairWithCache(const BWT* pBWT, 
                                          const BWT* pRevBWT, 
                                          const BWTIntervalCache* pFwdCache, 
                                          const BWTIntervalCache* pRevCache,
                                          const std::string& w);

// Count the number of times the sequence w appears in the collection, including
// its reverse complement
size_t countSequenceOccurrences(const std::string& w, const BWT* pBWT);
size_t countSequenceOccurrencesWithCache(const std::string& w, const BWT* pBWT, const BWTIntervalCache* pIntervalCache);
size_t countSequenceOccurrences(const std::string& w, const BWTIndexSet& indices);

// Count the occurrences of w, not including the reverse complement
size_t countSequenceOccurrencesSingleStrand(const std::string& w, const BWTIndexSet& indices);

// Update the given interval using backwards search
// If the interval corrsponds to string S, it will be updated 
// for string bS
inline void updateInterval(BWTInterval& interval, char b, const BWT* pBWT)
{
    size_t pb = pBWT->getPC(b);
    interval.lower = pb + pBWT->getOcc(b, interval.lower - 1);
    interval.upper = pb + pBWT->getOcc(b, interval.upper) - 1;
}

// Update the interval pair for the right extension to symbol b.
// In this version the AlphaCounts for the upper and lower intervals
// have been calculated
inline void updateBothR(BWTIntervalPair& pair, char b, const BWT* pRevBWT,
                        AlphaCount64& l, AlphaCount64& u)
{
    AlphaCount64 diff = u - l;

    pair.interval[0].lower = pair.interval[0].lower + diff.getLessThan(b);
    pair.interval[0].upper = pair.interval[0].lower + diff.get(b) - 1;

    // Update the right index directly
    size_t pb = pRevBWT->getPC(b);
    pair.interval[1].lower = pb + l.get(b);
    pair.interval[1].upper = pb + u.get(b) - 1;
}

//
// Update the interval pair for the right extension to symbol b.
// 
inline void updateBothR(BWTIntervalPair& pair, char b, const BWT* pRevBWT)
{
    // Update the left index using the difference between the AlphaCounts in the reverse table
    AlphaCount64 l = pRevBWT->getFullOcc(pair.interval[1].lower - 1);
    AlphaCount64 u = pRevBWT->getFullOcc(pair.interval[1].upper);
    updateBothR(pair, b, pRevBWT, l, u);
}

// Update the interval pair for the left extension to symbol b.
// In this version the AlphaCounts for the upper and lower intervals
// have been calculated.
inline void updateBothL(BWTIntervalPair& pair, char b, const BWT* pBWT, 
                        AlphaCount64& l, AlphaCount64& u)
{
    AlphaCount64 diff = u - l;
    // Update the left index using the difference between the AlphaCounts in the reverse table
    pair.interval[1].lower = pair.interval[1].lower + diff.getLessThan(b);
    pair.interval[1].upper = pair.interval[1].lower + diff.get(b) - 1;

    // Update the left index directly
    size_t pb = pBWT->getPC(b);
    pair.interval[0].lower = pb + l.get(b);
    pair.interval[0].upper = pb + u.get(b) - 1;
}

//
// Update the interval pair for the left extension to symbol b.
//
inline void updateBothL(BWTIntervalPair& pair, char b, const BWT* pBWT)
{
    // Update the left index using the difference between the AlphaCounts in the reverse table
    AlphaCount64 l = pBWT->getFullOcc(pair.interval[0].lower - 1);
    AlphaCount64 u = pBWT->getFullOcc(pair.interval[0].upper);
    updateBothL(pair, b, pBWT, l, u);
}


// Initialize the interval of index idx to be the range containining all the b suffixes
inline void initInterval(BWTInterval& interval, char b, const BWT* pB)
{
    interval.lower = pB->getPC(b);
    interval.upper = interval.lower + pB->getOcc(b, pB->getBWLen() - 1) - 1;
}

// Initialize the interval of index idx to be the range containining all the b suffixes
inline void initIntervalPair(BWTIntervalPair& pair, char b, const BWT* pBWT, const BWT* pRevBWT)
{
    initInterval(pair.interval[LEFT_INT_IDX], b, pBWT);
    initInterval(pair.interval[RIGHT_INT_IDX], b, pRevBWT);
}

// Return the counts of the bases between the lower and upper interval in pBWT
inline AlphaCount64 getExtCount(const BWTInterval& interval, const BWT* pBWT)
{
    return pBWT->getOccDiff(interval.lower - 1, interval.upper);
}

// Return the count of all the possible one base extensions of the string w.
// This returns the number of times the suffix w[i, l]A, w[i, l]C, etc 
// appears in the FM-index for all i s.t. length(w[i, l]) >= minOverlap.
AlphaCount64 calculateExactExtensions(const unsigned int overlapLen, const std::string& w, const BWT* pBWT, const BWT* pRevBWT);

// Calculate de Bruijn graph extensions of the given sequence using an index pair
// Returns an AlphaCount64 with the count of each extension base
// This function optionally takes in an interval cache to speed up the computation
AlphaCount64 calculateDeBruijnExtensions(const std::string str, 
                                         const BWT* pBWT, 
                                         const BWT* pRevBWT, 
                                         EdgeDir direction,
                                         const BWTIntervalCache* pFwdCache = NULL,
                                         const BWTIntervalCache* pRevCache = NULL);

// Calculate de Bruijn graph extensions of the given sequence using a single index
// This version is more computationally expensive than above but allows
// only one index to be held in memory. 
// Returns an AlphaCount64 with the count of each extension base
// This function optionally takes in an interval cache to speed up the computation
AlphaCount64 calculateDeBruijnExtensionsSingleIndex(const std::string str, 
                                                    const BWT* pBWT, 
                                                    EdgeDir direction,
                                                    const BWTIntervalCache* pFwdCache = NULL);

// Extract the complete string starting at idx in the BWT
std::string extractString(const BWT* pBWT, size_t idx);

// Extract the next len bases of the string starting at idx
std::string extractString(const BWT* pBWT, size_t idx, size_t len);

// Extract the substring from start, start+length of the sequence starting at position idx
std::string extractSubstring(const BWT* pBWT, uint64_t idx, size_t start, size_t length = std::string::npos);

// Extract all prefixes of the suffixes for the given interval, along
// with their lexicographic rank.
RankedPrefixVector extractRankedPrefixes(const BWT* pBWT, BWTInterval interval);

// Extract symbols from the starting index of the BWT until the index lies
// within the given interval. If the extraction hits the start of a string
// without finding a prefix, the empty string is returned.
std::string extractUntilInterval(const BWT* pBWT, int64_t start, const BWTInterval& interval);

// Returns a randomly chosen string from the BWT
std::string sampleRandomString(const BWT* pBWT);

// Returns a randomly chosen substring from the BWT 
std::string sampleRandomSubstring(const BWT* pBWT, size_t len);

};

#endif
