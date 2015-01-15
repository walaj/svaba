///-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// StdAlnTools - Collection of wrappers around the
// stdaln dynamic programming alignment functions
#ifndef STDALNTOOLS_H
#define STDALNTOOLS_H

#include <string>
#include <vector>
#include <inttypes.h>
#include <iostream>
#include "stdaln.h"

// Parameters object 
struct GlobalAlnParams
{
    GlobalAlnParams() { setDefaults(); }

    void setDefaults()
    {
        gap_open = 5;
        gap_ext = 2;
        mismatch = 3;
        match = 1;
        gap_open_extend = gap_open + gap_ext;
        threshold = 30;
        bandwidth = 50;
    }

    // PacBio specific parameters
    void setPacBio()
    {
        setDefaults();

        gap_open = 2;
        gap_ext = 1;
        mismatch = 5;
        gap_open_extend = gap_open + gap_ext;
    }

    int gap_open;
    int gap_ext;
    int mismatch;
    int match;
    int gap_open_extend;
    int threshold;
    int bandwidth;
};

struct LocalAlignmentResult
{
    // Indices of the start/end base in the aligning
    // substrings. INCLUSIVE coordinates.
    int64_t targetStartIndex;
    int64_t targetEndIndex;
    int64_t queryStartIndex;
    int64_t queryEndIndex;
    std::string cigar;
    int score;

    friend std::ostream& operator<<(std::ostream& out, const LocalAlignmentResult& a)
    {
        out << "S:" << a.score << " [" << a.targetStartIndex << " " << a.targetEndIndex << "] C:" << a.cigar;
        return out;
    }
};
typedef std::vector<LocalAlignmentResult> LocalAlignmentResultVector;

namespace StdAlnTools
{
    // Perform a global alignment between target and query using stdaln
    // If bPrint is true, the padded alignment is printed // to stdout.
    // The alignment score is returned.
    int globalAlignment(const std::string& target, const std::string& query, bool bPrint = false);

    // Perform a global alignment between the two strings and return a CIGAR string
    std::string globalAlignmentCigar(const std::string& target, const std::string& query);

    // Perform global alignment and return the cigar and score as out parameters
    void globalAlignment(const std::string& target, const std::string& query, std::string& cigar, int& score);

    // Perform a local alignment
    LocalAlignmentResult localAlignment(const std::string& target, const std::string& query);

    // Expand a Cigar string so there is one symbol per code
    std::string expandCigar(const std::string& cigar);
    
    // Compact an expanded cigar string
    std::string compactCigar(const std::string& cigar);

    // Convert a std::string into the stdAln required packed format.
    // This function allocates memory which the caller must free.
    uint8_t* createPacked(const std::string& s, size_t start = 0, size_t length = std::string::npos);

    // Create a path array representing the global alignment between target and query
    // Caller must free the path
    void createGlobalAlignmentPath(const std::string& target, const std::string& query,
                                            path_t** path, int* path_len, int* score);

    // Calculate the maximum aligned target length possible for a query of length ql
    // under the scoring scheme given by params
    size_t calculateMaxTargetLength(int ql, const GlobalAlnParams& params);

    // Fill in the stdaln AlnParam data, using GlobalAlnParams
    void setAlnParam(AlnParam& par, int matrix[25], const GlobalAlnParams& params);

    // Make a padded strings from a global alignment between the pair of sequences
    void makePaddedStrings(const std::string& s1, const std::string& s2, 
                           std::string& out1, std::string& out2);

    // Convert a stdaln path to a pair of padded strings representing the alignment
    void makePaddedStringsFromPath(const std::string& s1, const std::string& s2, 
                                   path_t* path, int path_len,
                                   std::string& out1, std::string& out2, std::string& outm);

    // Returns a new copy of str with padding characters removed
    std::string unpad(const std::string& str);
    
    // Make a cigar string from a path
    std::string makeCigar(path_t* path, int path_len);

    // Print the padded aligned strings
    void printPaddedStrings(const std::string& s1, const std::string& s2, const std::string& m, int colSize = 120);
};

#endif
