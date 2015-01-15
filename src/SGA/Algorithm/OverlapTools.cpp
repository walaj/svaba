//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapTools - Wrapper for the overlap machinery 
// to perform an overlap computation for two strings
//
#include "OverlapTools.h"
#include "SuffixArray.h"
//#define DPOVERLAPPRINT 1

// Find a suffix of s1 that matches a prefix of s2 with minimal edit distance
bool OverlapTools::boundedOverlapDP(const std::string& s1, const std::string& s2, int minOverlap, int maxOverlap, double maxErrorRate, Match& outMatch)
{
    const size_t MAX_REGION_LENGTH = 500; // do not try to align strings that are longer than this number
    if(s1.empty() || s2.empty())
        return false;
 
    if(maxOverlap > (int)MAX_REGION_LENGTH)
    {
        std::cerr << "Maximum requested overlap " << maxOverlap << " is larger than search limit (" << MAX_REGION_LENGTH << "), not aligning\n";
        return false;
    }

    // The dynamic programming is too heavy to compute matches for very large strings
    // so we only try to align the ends of s1 and s2 up to MAX_REGION_LENGTH bases.
    size_t sublen = min3(maxOverlap, s1.size(), s2.size());
    if(sublen > MAX_REGION_LENGTH)
        return false;

    std::string sub1 = s1.substr(s1.size() - sublen);
    std::string sub2 = s2.substr(0, sublen);

    DPAlignment dpAlign(sub1, sub2, DPM_OVERLAP, DPSS_SIMILARITY);
    Match match;
    bool validMatch = findBestOverlapByScore(sub1, sub2, minOverlap, maxErrorRate, dpAlign, match);
    if(validMatch)
    {
        // If we are matching a substring of both s1 and s2, and the match
        // is a containment (one string was completely matched) we call the
        // match invalid. 
        if(match.isContainment() && sublen < s1.size() && sublen < s2.size())
        {
            return false;
        }
        else
        {
            // A valid match has been found which is maximal
            // Expand the matches to the coordinates of s1 and s2 and return
            int s1_add = s1.size() - sublen;
            match.coord[0].interval.start += s1_add;
            match.coord[0].interval.end += s1_add;
            match.coord[0].seqlen = s1.size();

            match.coord[1].seqlen = s2.size();

#ifdef DPOVERLAPPRINT
            match.printMatch(s1, s2);
#endif
            outMatch = match;
            return true;
        }
    }
    else
    {
        return false;
    }
}    

//
bool OverlapTools::findBestOverlapByScore(const std::string& s1, const std::string& s2, 
                                          int minOverlap, double maxErrorRate, 
                                          const DPAlignment& dpAlign, Match& outMatch)
{
    assert(dpAlign.getScoringScheme() == DPSS_SIMILARITY);

    // The last row of the table holds the edit distance for the prefix-suffix overlaps of S1
    // Iterate over this row and calculate the length and edit % of each overlap
    int overlapRowIdx = dpAlign.getNumRows() - 1;
    int numCols = dpAlign.getNumColumns();
    int bestOverlapLen = -1;
    int bestScore = 0;
    for(int j = minOverlap; j < numCols; ++j)
    {
        int overlapLen = j;

        // Do not allow overlaps to end with insertions
        // For example, this case is a valid alignment but would be disallowed:
        // S1 AGTACATTTACA--
        // S2    ACATTTACATA
        // This is an overlap of length 11 but the "true" overlap is only 9 bases
        // Under a reasonable scoring system this should not happen
        DPOperation lastOp = dpAlign.getPathOperationToCell(s1, s2, overlapRowIdx, j);
        
        if(lastOp != DPO_INSERT)
        {
            int score = dpAlign.getScore(overlapRowIdx,j);
            if(score > bestScore)
            {
                // A higher scoring overlap has been found. Make sure that it is valid
                // by checking that the error rate and overlap length are within the expected
                // counts and that the coordinates are extreme w.r.t both sequences
                int i_match = 0;
                int num_edits = 0;
                DPPath overlapPath = dpAlign.calculatePath(s1, s2, overlapRowIdx, overlapLen);
                for(DPPath::iterator iter = overlapPath.begin(); iter != overlapPath.end(); ++iter)
                {
                    if(iter->j == 0 && iter->i > i_match)
                        i_match = iter->i;

                    // Count the mismatches but don't count deletions in the first column which are "free"
                    if(iter->op != DPO_MATCH && !(iter->j == 0 && iter->op == DPO_DELETE))
                        ++num_edits;
                }
                
                Match currMatch;
                currMatch.isReverse = false;
                currMatch.coord[0].interval.start = i_match;
                currMatch.coord[0].interval.end = overlapRowIdx - 1; // match coords are inclusive
                currMatch.coord[0].seqlen = s1.size();
                currMatch.coord[1].interval.start = 0;
                currMatch.coord[1].interval.end = overlapLen - 1;
                currMatch.coord[1].seqlen = s2.size();
                currMatch.setNumDiffs(num_edits);
                double errorRate = (double)num_edits / currMatch.getMinOverlapLength();
                bool isValid = errorRate < maxErrorRate && currMatch.getMinOverlapLength() >= minOverlap &&
                               currMatch.coord[0].isExtreme() && currMatch.coord[1].isExtreme();           
                
                if(isValid)
                {
                    bestOverlapLen = overlapLen;
                    bestScore = score;
                    outMatch = currMatch;
                }
            }
         
            /*   
            std::cout << "\n";
            dpAlign.printAlignment(s1, s2, overlapRowIdx, j);
            std::cout << "OL: " << overlapLen << "\n";
            std::cout << "LO: " << (int)lastOp << " ED: " << score << "\n";
            */
        }
    }

    if(bestOverlapLen > 0)
    {

#ifdef DPOVERLAPPRINT
        std::cout << "\nBest overlap: " << bestOverlapLen << "\n";
        dpAlign.printAlignment(s1, s2, overlapRowIdx, bestOverlapLen);
        outMatch.printMatch(s1, s2);
#endif
        return true;
    }
    else
    {
        return false;
    }
}
