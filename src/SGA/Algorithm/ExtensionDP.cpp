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
#include "ExtensionDP.h"
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <limits>

static const int MINUS_INF = -1073741823;
static const int POSITIVE_INF = 1073741823;

// Return the largest of the 3 values
inline int max3(int a, int b, int c)
{
    return std::max(std::max(a,b), c);
}

// Return the smalleset of the 3 values
inline int min3(int a, int b, int c)
{
    return std::min(std::min(a,b), c);
}

// Constructor
BandedDPColumn::BandedDPColumn(int ci, int maxRows, int bandwidth, const BandedDPColumn* prevColumn)
{
    // Calculate the lower and upper band indices
    m_bandwidth = bandwidth; // the requested bandwidth, may be greater than actual
    m_colIdx = ci;
    m_rowStartIdx = std::max(0, m_colIdx - (bandwidth / 2) - 1);
    m_rowEndIdx = std::min(maxRows - 1, m_colIdx + (bandwidth / 2));
    m_maxRows = maxRows;
    assert(m_rowStartIdx <= m_rowEndIdx);
    int numRows = m_rowEndIdx - m_rowStartIdx + 1;
    m_cells.resize(numRows);

    m_pPrevColumn = prevColumn;
}

// Get the score for row
int BandedDPColumn::getRowScore(int row) const
{
    int vec_idx = getVectorIndex(row);
    if(vec_idx == -1)
        return POSITIVE_INF;
    else
        return m_cells[vec_idx].score;
}

// Get the type of the row
char BandedDPColumn::getRowType(int row) const
{
    int vec_idx = getVectorIndex(row);
    if(vec_idx == -1)
        return FROM_M;
    else
        return m_cells[vec_idx].ctype;
}

// Return the index of the best-scoring row
int BandedDPColumn::getBestRowIndex() const
{
    int bestScore = POSITIVE_INF;
    int bestIdx = -1;
    for(int i = m_rowStartIdx; i <= m_rowEndIdx; ++i)
    {
        int s = getRowScore(i);
        if(s < bestScore)
        {
            bestScore = s;
            bestIdx = i;
        }
    }
    return bestIdx;
}


// Return the pointer to the previous column
const BandedDPColumn* BandedDPColumn::getPreviousColumn() const
{
    return m_pPrevColumn;
}

// Set the score and type for row 
void BandedDPColumn::setRowScore(int row, int score, char ctype)
{
    int vec_idx = getVectorIndex(row);
    if(vec_idx != -1)
    {
        m_cells[vec_idx].score = score;
        m_cells[vec_idx].ctype = ctype;
    }
}

// Return the score for the requested cell
int BandedDPColumn::getCellScore(int colIdx, int rowIdx) const
{
    // If accessing the col/row with negative index, return
    // the initial gap penalties. For now, we allow free gaps
    // for the first row
    if(colIdx == -1 || rowIdx == -1)
        return 0;

    const BandedDPColumn* pColumn;
    if(colIdx < m_colIdx)
        pColumn = m_pPrevColumn;
    else
        pColumn = this;
    assert(pColumn != NULL);
    
    return pColumn->getRowScore(rowIdx);   
}

// Calculate the score for the given row in the column using an edit distance scheme
void BandedDPColumn::fillRowEditDistance(int rowIdx, int matchScore)
{
    int vec_idx = getVectorIndex(rowIdx);
    if(vec_idx == -1)
        return; // out of band, ignore

    // Get the scores for the cell above, to the left and diagonal
    int above = getCellScore(m_colIdx, rowIdx - 1) + 1;
    int diag = getCellScore(m_colIdx - 1, rowIdx - 1) + matchScore;
    
    // Give free insertions if the last row
    // JS Apr/12 do not give free insertions in the last row
    //int insert_score = rowIdx < (m_maxRows - 1) ? 1 : 0;
    int insert_score = 1;
    int left = getCellScore(m_colIdx - 1, rowIdx) + insert_score;
    int score = min3(above, left, diag);

    char ctype;
    if(score == diag)
        ctype = FROM_M;
    else if(score == above)
        ctype = FROM_D;
    else
        ctype = FROM_I;
//    printf("Cell(%d, %d) a: %d d: %d l: %d s: %d t: %d\n", m_colIdx, rowIdx, above, diag, left, score, ctype);
    setRowScore(rowIdx, score, ctype);
}

// Transform a row index into an index in the scores vector
int BandedDPColumn::getVectorIndex(int rowIdx) const
{
    if(rowIdx < m_rowStartIdx || rowIdx > m_rowEndIdx)
        return -1;
    else
        return rowIdx - m_rowStartIdx;
}

// Initialize the extension DP by computing a global alignment between extendable and fixed
// This function allocates memory and stores the created pointers in outPtrVec, which
// should initially be empty
void ExtensionDP::createInitialAlignment(const std::string& extendable, const std::string& fixed, int bandwidth, BandedDPColumnPtrVector& outPtrVec)
{
    assert(outPtrVec.empty());
    // Initialize the zero-th column
    size_t numRows = fixed.size() + 1;
    size_t numCols = extendable.size() + 1;

    // Set the score of column zero 
    BandedDPColumn* pZeroCol = new BandedDPColumn(0, numRows, bandwidth, NULL);
    pZeroCol->setRowScore(0, 0, FROM_M);
    for(size_t i = 1; i < numRows; ++i)
        pZeroCol->setRowScore(i, i, FROM_I);

    outPtrVec.push_back(pZeroCol);
    BandedDPColumn* pPrevCol = outPtrVec.back();
    for(size_t colIdx = 1; colIdx < numCols; ++colIdx)
    {
        BandedDPColumn* pCurrCol = createNewColumn(extendable[colIdx - 1], fixed, pPrevCol);
        outPtrVec.push_back(pCurrCol);
        pPrevCol = pCurrCol;
    }
}

// Create a new column, representing the extension of the alignment to char b
BandedDPColumn* ExtensionDP::createNewColumn(char b, const std::string& fixed, const BandedDPColumn* pPrevColumn)
{
    int numRows = fixed.size() + 1;
    int colIdx = pPrevColumn->getColIdx() + 1;
    BandedDPColumn* pCurrCol = new BandedDPColumn(colIdx, numRows, pPrevColumn->getBandwidth(), pPrevColumn);

    // Set the first row score
    pCurrCol->setRowScore(0, colIdx, FROM_D);
    
    // Fill in the rest of the rows
    int startRow = std::max(1, pCurrCol->getMinRow());
    for(int rowIdx = startRow; rowIdx <= pCurrCol->getMaxRow(); ++rowIdx)
    {
        int match_score = b == fixed[rowIdx - 1] ? 0 : 1;
        pCurrCol->fillRowEditDistance(rowIdx, match_score);
    }
    return pCurrCol;
}

// Calculate the edit percentage over the last numBases bases
double ExtensionDP::calculateLocalEditPercentage(const BandedDPColumn* pColumn, int numBases)
{
    // Get the row with in the start column with the best score
    int rowIdx = pColumn->getBestRowIndex();
    assert(rowIdx != -1);
    int colIdx = pColumn->getColIdx();
    int startScore = pColumn->getRowScore(rowIdx);

    int alignLength = 1;
    while(alignLength < numBases)
    {
        char ctype = pColumn->getRowType(rowIdx);
        switch(ctype)
        {
   			case FROM_M: 
                pColumn = pColumn->getPreviousColumn();
                --colIdx; 
                --rowIdx; 
                break;
			case FROM_D: 
                --rowIdx; 
                break;
			case FROM_I: 
                pColumn = pColumn->getPreviousColumn();
                --colIdx; 
                break;
        }
        alignLength += 1;
    }
    int endScore = pColumn->getRowScore(rowIdx);

    assert(startScore >= endScore);
    int numDiffs = startScore - endScore;
    return static_cast<double>(numDiffs) / alignLength;
}

// Returns the error percentage over the entire alignment
double ExtensionDP::calculateGlobalEditPercentage(const BandedDPColumn* pStartColumn)
{
    int edits, alignLength;
    countEditsAndAlignLength(pStartColumn, edits, alignLength);
    return (double)edits / alignLength;
}

// Returns the number of edits and the total alignment length for the submatrix starting at pColumn
void ExtensionDP::countEditsAndAlignLength(const BandedDPColumn* pColumn, int& edits, int& alignLength)
{
    // Get the row with in the start column with the best score
    int rowIdx = pColumn->getBestRowIndex();
    assert(rowIdx != -1);
    int colIdx = pColumn->getColIdx();
    edits = pColumn->getRowScore(rowIdx);
    alignLength = 1;

    while(rowIdx != 0 && colIdx != 0)
    {
        char ctype = pColumn->getRowType(rowIdx);
        switch(ctype)
        {
   			case FROM_M: 
                pColumn = pColumn->getPreviousColumn();
                --colIdx; 
                --rowIdx; 
                break;
			case FROM_D: 
                --rowIdx; 
                break;
			case FROM_I: 
                pColumn = pColumn->getPreviousColumn();
                --colIdx; 
                break;
        }
        alignLength += 1;
    }
}

// Calculate the best alignment through the matrix. Assumes that
// path_t* is pre-allocated with max_path entries.
void ExtensionDP::backtrack(const BandedDPColumn* pLastColumn, path_t* path, int* path_len, const int maxPathLength)
{
    const BandedDPColumn* pCurrColumn = pLastColumn;

    // Calculate the starting row and column
    int colIdx = pLastColumn->getColIdx();
    int rowIdx = pLastColumn->getBestRowIndex();
    int pathLength = 0;


    path_t* p = &path[pathLength];
    char ctype = pLastColumn->getRowType(rowIdx);

	p->ctype = ctype; 
    p->i = rowIdx; 
    p->j = colIdx;
	pathLength += 1;
    
	do {
        assert(pathLength < maxPathLength);
        assert(pCurrColumn != NULL);

        p = &path[pathLength];

        // Update indices
		switch (ctype) {
			case FROM_M: 
                pCurrColumn = pCurrColumn->getPreviousColumn();
                --colIdx; 
                --rowIdx; 
                break;
			case FROM_D: 
                --rowIdx; 
                break;
			case FROM_I: 
                pCurrColumn = pCurrColumn->getPreviousColumn();
                --colIdx; 
                break;
		}
		ctype = pCurrColumn->getRowType(rowIdx);
		p->ctype = ctype; 
        p->i = rowIdx; 
        p->j = colIdx;
        pathLength += 1;
	} while (rowIdx || colIdx);

    // We ignore the path character from the first column/row, like stdaln
    *path_len = pathLength - 1;
}

// Return true if the extension can proceed no further. 
// This is based on two conditions:
// 1) There are n insertions at the end of the sequence
// 2) The bandwidth in this column is one 
bool ExtensionDP::isExtensionTerminated(const BandedDPColumn* pLastColumn, int insertionThreshold)
{
    if(pLastColumn->getMinRow() == pLastColumn->getMaxRow())
        return true;

    // Backtrack up to n steps checking for insertions
    const BandedDPColumn* pCurrColumn = pLastColumn;
    int rowIdx = pLastColumn->getBestRowIndex();
    int colIdx = pLastColumn->getColIdx();
    char ctype = pCurrColumn->getRowType(rowIdx);
	do {
        assert(pCurrColumn != NULL);

		switch (ctype) 
        {
			case FROM_M: 
            case FROM_D:
                return false;
			case FROM_I: 
                pCurrColumn = pCurrColumn->getPreviousColumn();
                --colIdx; 
                insertionThreshold -= 1;
                if(insertionThreshold == 0)
                    return true;
		}
		ctype = pCurrColumn->getRowType(rowIdx);
	} while (rowIdx || colIdx);
    return false;
}

// Parse the score matrix and return a trimmed alignment between the query and target strings.
// The trimming condition requires the last minMatches bases of the two sequences to be perfectly aligned
ExtensionDPAlignment ExtensionDP::findTrimmedAlignment(const BandedDPColumn* pLastColumn, int minMatches)
{
    const BandedDPColumn* pBestColumn = pLastColumn;
    bool columnFound = false;

    // Iterate backwards over the columns searching for one that 
    // has minMatches matches to the query at the end
    while(!columnFound)
    {
        // Check if the best alignment starting from this column
        // ends with n matches.
        int matches = 0;
        int mismatches = 0;
        const BandedDPColumn* pCurrColumn = pBestColumn;
        int rowIdx = pCurrColumn->getBestRowIndex();
        int colIdx = pCurrColumn->getColIdx();
        int startScore = pCurrColumn->getRowScore(rowIdx);
        char ctype = pCurrColumn->getRowType(rowIdx);

        do {
            if(matches == minMatches)
            {
                // Ensure the edit distance is the same
                columnFound = true;
                break;    
            }

            switch (ctype) 
            {
                case FROM_M: 
                    pCurrColumn = pCurrColumn->getPreviousColumn();
                    --colIdx; 
                    --rowIdx;
                    matches += 1;
                    break;

                case FROM_D:
                case FROM_I:
                    mismatches += 1;
                    columnFound = false;
                    break;
            }
            
            ctype = pCurrColumn->getRowType(rowIdx);
            
            // Update mismatches
            if(pCurrColumn->getRowScore(rowIdx) != startScore)
                ++mismatches;

        } while ((rowIdx || colIdx) && pCurrColumn != NULL && mismatches == 0);

        if(!columnFound)
            pBestColumn = pBestColumn->getPreviousColumn();
    }
    
    // set the result output
    ExtensionDPAlignment result;
    if(pBestColumn != NULL)
    {
        result.target_align_length = pBestColumn->getColIdx();
        result.query_align_length = pBestColumn->getBestRowIndex();
    }
    else
    {
        result.target_align_length = 0;
        result.query_align_length = 0;
    }
    return result;
}

// Parse the score matrix and return a the best alignment that includes the entire query string.
ExtensionDPAlignment ExtensionDP::findGlocalAlignment(const BandedDPColumn* pLastColumn)
{
    const BandedDPColumn* pCurrColumn = pLastColumn;
    int lastRowIdx = pLastColumn->getQueryRows() - 1;
    
    const BandedDPColumn* pBestColumn = NULL;
    int best_score = std::numeric_limits<int>::max();

    // Iterate backwards over the columns searching for one that 
    // has minMatches matches to the query at the end
    while(pCurrColumn != NULL && pCurrColumn->getMaxRow() >= lastRowIdx)
    {
        // Get the score for the last row of the matrix
        int score = pCurrColumn->getRowScore(lastRowIdx);
        if(score < best_score)
        {
            best_score = score;
            pBestColumn = pCurrColumn;
        }

        pCurrColumn = pCurrColumn->getPreviousColumn();
    }
    
    // set the result output
    ExtensionDPAlignment result;
    if(pBestColumn != NULL)
    {
        result.target_align_length = pBestColumn->getColIdx();
        result.query_align_length = lastRowIdx;
    }
    else
    {
        result.target_align_length = 0;
        result.query_align_length = 0;
    }
    return result;
}

// Print the alignment starting from the given column
void ExtensionDP::printAlignment(const std::string& extendable, const std::string& fixed, const BandedDPColumn* pLastColumn)
{
    // allocate space for the path
    int maxPathLength = extendable.size() + fixed.size();
    path_t* path = new path_t[maxPathLength];
    int path_len;

    // Perform the backtrack
    backtrack(pLastColumn, path, &path_len, maxPathLength);

    // Print the strings
    std::string paddedF, paddedE, paddedMatch;
    StdAlnTools::makePaddedStringsFromPath(fixed, extendable, path, path_len, paddedF, paddedE, paddedMatch);
    StdAlnTools::printPaddedStrings(paddedF, paddedE, paddedMatch);

    std::cout << "CIGAR: " << StdAlnTools::makeCigar(path, path_len) << "\n";

    delete [] path;
}

// Print the full scoring matrix
void ExtensionDP::printMatrix(const BandedDPColumnPtrVector& columnPtrVec)
{
    int numCols = columnPtrVec.size();
    int numRows = columnPtrVec[numCols - 1]->getMaxRow() + 1;

    printf("c:\t");
    for(int j = 0; j < numCols; ++j)
    {
        printf("%d\t", j);
    }
    printf("\n");

    for(int i = 0; i < numRows; ++i)
    {
        printf("%d:\t", i);
        for(int j = 0; j < numCols; ++j)
        {
            int score = columnPtrVec[j]->getRowScore(i);
            printf("%d\t", score);
        }
        printf("\n");
    }
}

