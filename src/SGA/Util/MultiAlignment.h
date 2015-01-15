//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// MultiAlignment - Class constructing a multiple
// alignment between a root sequence and a set of sequences
//
#ifndef MULTIALIGNMENT_H
#define MULTIALIGNMENT_H

#include <string>
#include <vector>
#include "Util.h"

struct MAlignData
{
    public:
        std::string str;
        std::string expandedCigar; // cigar string with symbol counts expanded
        int position; // start position of the alignment to the root
        std::string name; // named identifier for this sequence in the MA

    private:

        friend class MultiAlignment;
        std::string padded;

        // Comparator
        static bool sortPosition(const MAlignData& a, const MAlignData& b) { return a.position < b.position; }
};
typedef std::vector<MAlignData> MAlignDataVector;

class MultiAlignment
{
    public:
        MultiAlignment(std::string rootStr, const MAlignDataVector& inData, std::string rootName = "root");

        // Experimental function to generate a consensus sequence from the MA
        std::string generateConsensus();

        // Get the index into the m_alignData vector for a named row
        size_t getIdxByName(const std::string& name) const;
        size_t getRootIdx() const { return 0; }
        
        // Get the symbol at a particular column and row
        char getSymbol(size_t rowIdx, size_t colIdx) const;

        // Get the index of the base corresponding to colIdx for sequence at rowIdx
        size_t getBaseIdx(size_t rowIdx, size_t colIdx) const;

        // Get a substring of the padded string for the given row
        std::string getPaddedSubstr(size_t rowIdx, size_t start, size_t length) const;

        // Count the length of the homopolymer run in the substring [from, to] inclusive
        // If to is npos, then it is assumed to be until the end of the string
        // if to < from, then the count proceeds backwards (towards the beginning of the string.
        size_t countHomopolymer(size_t rowIdx, int from, int to = -1) const;

        // Return the total number of columns in the MA
        size_t getNumColumns() const;
        
        // Return the cigar string for one of the alignments
        std::string getCigar(size_t rowIdx) const;

        // Filter the multiple alignment by removing sequences that are more than max_edit_distance from the root sequence
        void filterByEditDistance(int max_edit_distance);

        // Print the multiple alignment, optionally with a consensus sequence
        void print(int col_size = 80, const std::string* pConsensus = NULL, bool sorted = true, bool masked = false) const;

    private:
        
        // Generate a string representing the columns that differ between the strings
        std::string generateMatchString() const;

        // data
        MAlignDataVector m_alignData;
        int m_verbose;
};

namespace MultiAlignmentTools
{
    // Construct a multiple alignment of the input strings
    MultiAlignment alignSequencesGlobal(const SeqItemVector& sequences);
    MultiAlignment alignSequencesLocal(const SeqItemVector& sequences);
}

#endif
