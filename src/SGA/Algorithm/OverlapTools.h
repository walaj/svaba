//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapTools - Wrapper for the overlap machinery 
// to perform an overlap computation for two strings
//
#ifndef OVERLAPTOOLS_H
#define OVERLAPTOOLS_H

#include "Match.h"
#include "DPAlignment.h"
#include <list>

namespace OverlapTools
{
    // Datatypes
    typedef std::vector<int> IntVector;
    typedef std::vector<IntVector> DPMatrix;

    //typedef EditDistanceScoring DPScoring;
    typedef SimilarityScoring DPScoring;

    // Find a match between s1 and s2 using a dynamic programming alignment
    // This is something of a heavy function and the longest possible match must be provided
    // If there is a match between minOverlap and maxOverlap with edit distance less than maxErrorRate, it will be returned.
    // If there are multiple valid matches like this, the highest scoring match will be returned under the DPSS_SIMILARITY
    // scoring scheme. If true is returned, a valid match has been found and outMatch is filled in appropriately.
    // This function should not be used to find overlaps greater than 500bp
    bool boundedOverlapDP(const std::string& s1, const std::string& s2, int minOverlap, int maxOverlap, double maxErrorRate, Match& outMatch);

    // Find the best overlap using the dpAlign structure. Returns true
    // if a valid overlap has been found and fills in the match structure
    bool findBestOverlapByScore(const std::string& s1, const std::string& s2, 
                                int minOverlap, double maxErrorRate, const DPAlignment& dpAlign,
                                Match& outMatch);

};

#endif
