///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
//
// ExtensionDP - Class implementing an iterative
// dynamic programming alignment. The alignment
// starts by globally constructing an alignment
// between two strings in a seed region. This initial
// alignment can then be extended row-by-row as one of the
// strings is extended. This is the core data structure for the 
// StringThreader class
// 
#ifndef EXTENSION_DP_H
#define EXTENSION_DP_H

#include <vector>
#include "StdAlnTools.h"

// Structure holding the score for a single cell
// of the DP matrix
struct DPCell
{
    DPCell() : score(0), ctype('\0') {}
    int score;
    char ctype;
};
typedef std::vector<DPCell> CellVector;

// Result object providing the coordinates
// of the endpoints of the alignment
struct ExtensionDPAlignment
{
    int query_align_length;
    int target_align_length;
};

// Core class providing a single column of the dynamic
// programming matrix
class BandedDPColumn
{
    public:

        //
        // Functions
        //
        BandedDPColumn(int ci, int maxRows, int bandwidth, const BandedDPColumn* prevColumn);

        // Simple getters
        int getColIdx() const { return m_colIdx; }
        int getBandwidth() const { return m_bandwidth; }

        int getRowScore(int row) const;
        char getRowType(int row) const;
        int getCellScore(int colIdx, int rowIdx) const;
        int getMinRow() const { return m_rowStartIdx; }   
        int getMaxRow() const { return m_rowEndIdx; }
        int getQueryRows() const { return m_maxRows; }

        const BandedDPColumn* getPreviousColumn() const;

        // Set the cell information for a given row in the column
        void setRowScore(int row, int score, char ctype);
        
        // Returns the row index of the best scoring cell
        int getBestRowIndex() const;

        // Calculate and set the score for the row
        void fillRowEditDistance(int rowIdx, int matchScore);

    private:
        
        //
        // Functions
        //
        int getVectorIndex(int rowIdx) const;

        //
        // Data
        //
        int m_colIdx;
        int m_maxRows;
        int m_rowStartIdx;
        int m_rowEndIdx;
        int m_bandwidth;
        CellVector m_cells;
        const BandedDPColumn* m_pPrevColumn;
};
typedef std::vector<BandedDPColumn*> BandedDPColumnPtrVector;

// Algorithms to use/extend the extensionDP functionality
namespace ExtensionDP
{
    // Initialize the extension DP by computing a global alignment between extendable and fixed
    // This function allocates memory and stores the created pointers in outPtrVec
    void createInitialAlignment(const std::string& extendable, const std::string& fixed, int bandwidth, BandedDPColumnPtrVector& outPtrVec);

    // Create a new column, representing the extended alignment from pPrevColumn to char b
    BandedDPColumn* createNewColumn(char b, const std::string& fixed, const BandedDPColumn* pPrevColumn);

    // Calculate the best alignment through the matrix. Assumes that
    // path_t* is pre-allocated with maxPathLength entries.
    void backtrack(const BandedDPColumn* pLastColumn, path_t* path, int* pPathLen, const int maxPathLength);

    // Check if the extension has terminated. This is true if there are more
    // than insertionThreshold insertions at the end of the alignment.
    bool isExtensionTerminated(const BandedDPColumn* pLastColumn, int insertionThreshold);

    // Return the best trimmed alignment of the extensionDP object requiring minMatches matches
    // at the end of the alignment
    ExtensionDPAlignment findTrimmedAlignment(const BandedDPColumn* pLastColumn, int minMatches);

    // Return the best alignment that contains the entire fixed string
    ExtensionDPAlignment findGlocalAlignment(const BandedDPColumn* pLastColumn);
    
    // Calculate mismatch/error rates between the query and the string up to pStartColumn
    double calculateLocalEditPercentage(const BandedDPColumn* pStartColumn, int numBases);
    double calculateGlobalEditPercentage(const BandedDPColumn* pStartColumn);
    void countEditsAndAlignLength(const BandedDPColumn* pStartColumn, int& edits, int& alignLength);

    // Print the alignment starting from the given column
    void printAlignment(const std::string& s1, const std::string& s2, const BandedDPColumn* pLastColumn);

    // Print the full scoring matrix
    void printMatrix(const BandedDPColumnPtrVector& columnPtrVec);
};

#endif
