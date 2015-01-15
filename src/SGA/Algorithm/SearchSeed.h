//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SearchSeed.h - Data structure holding a partial 
// alignment to a FM-index. 
//
#ifndef SEARCHSEED_H
#define SEARCHSEED_H

#include <queue>
#include <list>
#include "BWTInterval.h"
#include "SearchHistory.h"

// types
enum ExtendDirection
{
    ED_LEFT,
    ED_RIGHT
};

// A search seed describes a partial alignment to a BWT
struct SearchSeed
{
    //
    // Functions
    //
    inline int length() const { return right_index - left_index + 1; }
    inline bool isSeed() const { return length() < seed_len; }
    inline bool isIntervalValid(int idx) { return ranges.interval[idx].isValid(); }
    inline bool allowMismatches() const { return z < maxDiff; }
    inline double calcErrorRate() const { return static_cast<double>(z) / static_cast<double>(length()); }
     
    // Sort the alignments based on their r_lower/r_upper
    static inline bool compareLeftRange(const SearchSeed& a, const SearchSeed& b)
    {
        return BWTInterval::compare(a.ranges.interval[0], b.ranges.interval[0]);
    }

    // Compare for equality based on the left range
    // If the length of the alignment is equal, then if the left ranges
    // are a match, the two alignment objects are redundant and one can be removed
    static inline bool equalLeftRange(const SearchSeed& a, const SearchSeed& b)
    {
#ifdef VALIDATE
        if(BWTInterval::equal(a.ranges.inteval[0], b.ranges.interval[0]))
            assert(a.length() == b.length() && a.z == b.z);
#endif
        return BWTInterval::equal(a.ranges.interval[0], b.ranges.interval[0]);
    }
    
    friend bool operator==(const SearchSeed& a, const SearchSeed& b)
    {
        return a.ranges == b.ranges && a.left_index == b.left_index &&
               a.right_index == b.right_index && a.z == b.z && a.dir == b.dir; 
    }

    //
    void print() const;
    void print(const std::string& w) const;
     
    //
    // Data
    //

    // BWT interval coordinates, the first element of the pair is the left range
    BWTIntervalPair ranges;
    SearchHistoryLink historyLink;

    // Index range is inclusive on both ends
    int left_index;
    int right_index;
    int seed_len;

    // The direction the alignment is currently being extended
    ExtendDirection dir; 
    
    // The current number of mismatches in the block
    int z;

    // The total number of allowed mismatches
    int maxDiff;

};

// Collections
typedef std::vector<SearchSeed> SearchSeedVector;
typedef std::queue<SearchSeed> SearchSeedQueue;

#endif
