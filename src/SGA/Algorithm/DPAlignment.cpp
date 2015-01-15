//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// DPAlignment - Generic class for dynamic programming alignments 
//
#include "DPAlignment.h"

DPAlignment::DPAlignment(const std::string& s1, 
                         const std::string& s2, 
                         DPAlignMode mode, DPScoringScheme scoring) : m_mode(mode), m_scoring(scoring)
{
    if(s1.empty() || s2.empty())
        return;

    m_dimM = s1.size() + 1; // = num rows
    m_dimN = s2.size() + 1; // = num columns
    initializeDPMatrix();
    fillMatrix(s1, s2);
    //printMatrix(s1, s2);
}

//
int DPAlignment::getNumRows() const
{
    return m_dimM;
}

//
int DPAlignment::getNumColumns() const
{
    return m_dimN;
}

// 
int DPAlignment::getScore(int i, int j) const
{
    assert(i < m_dimM && j < m_dimN);
    return m_matrix[i][j];
}

//
DPScoringScheme DPAlignment::getScoringScheme() const
{
    return m_scoring;
}

// Setup the dynamic programming matrix for overlap computation
void DPAlignment::initializeDPMatrix()
{
    m_matrix.resize(m_dimM); // rows
    for(int i = 0; i < m_dimM; ++i)
        m_matrix[i].resize(m_dimN); // columns

    // Set score[i][0] (first column) to be zero for all i if in overlap mode.
    // This is so we don't penalize the portion of S1 that
    // S2 does not overlap
    for(int i = 0; i < m_dimM; ++i)
        m_matrix[i][0] = (m_mode == DPM_OVERLAP) ? 0 : i * gap_penalty();

    // Set score [0][j] (first row) to be j for all j
    // This represents not matching the first j bases of S1
    for(int j = 0; j < m_dimN; ++j)
        m_matrix[0][j] = gap_penalty() * j;
}

// Fill in the DP matrix
void DPAlignment::fillMatrix(const std::string& s1, const std::string& s2)
{
    // Update the value of each cell for all i and j greater than zero
    for(int i = 1; i < m_dimM; ++i)
    {
        for(int j = 1; j < m_dimN; ++j)
        {
            int up = m_matrix[i - 1][j] + gap_penalty(); // skip a base of S1
            int left = m_matrix[i][j - 1] + gap_penalty(); // skip a base of S2

            int idx1 = i - 1;
            int idx2 = j - 1;
            int diag = m_matrix[i - 1][j - 1] + match_score(s1[idx1], s2[idx2]);
            m_matrix[i][j] = compare3(up, left, diag);
        }
    }
}

// Print the best alignment starting from the given cell
void DPAlignment::printAlignment(const std::string& s1, const std::string& s2,
                                 int startI, int startJ) const
{
    if(m_matrix.empty())
        return;
    if(startI == -1 || startJ == -1)
    {
        startI = m_dimM - 1;
        startJ = m_dimN - 1;
    }
    int score = m_matrix[startI][startJ];
    DPPath backtrack = calculatePath(s1, s2, startI, startJ);

    // Print the match strings
    std::string m1;
    std::string m2;
    std::string os;
    for(DPPath::iterator iter = backtrack.begin(); iter != backtrack.end(); ++iter)
    {
        int i = iter->i;
        int j = iter->j;
        //printf("(%d %d)\n", i, j);
        switch(iter->op)
        {
            case DPO_MATCH:
            case DPO_SUB:
                assert(i >= 0 && i < (int)s1.size() && j >= 0 && j < (int)s2.size());
                m1.push_back(s1[i]);
                m2.push_back(s2[j]);
                os.push_back(iter->op == DPO_MATCH ? '.' : 'S');
                break;
            case DPO_INSERT:
                assert(j >= 0 && j < (int)s2.size());
                m1.push_back('-');
                m2.push_back(s2[j]);
                os.push_back('I');
                break;
            case DPO_DELETE:
                assert(i >= 0 && i < (int)s1.size());
                m1.push_back(s1[i]);
                m2.push_back('-');
                os.push_back('D');
                break;
            case DPO_NOOP:
                assert(false);
                break;
        }
    }

    std::cout << "M1: " << m1 << "\n";
    std::cout << "M2: " << m2 << "\n";
    std::cout << "OS: " << os << "\n";
    std::cout << "Score: " << score << "\n";
}

// Get the operation that was performed to arrive at cell i,j
// If there are multiple possible paths, the
// precedence is diagonal, above, left
DPOperation DPAlignment::getPathOperationToCell(const std::string& s1, const std::string& s2, int i, int j) const
{
    assert((i > 0 || j > 0) && i < m_dimM && j < m_dimN);
    // Get the score of cell i,j
    int currScore = getScore(i,j);

    // Check diagonal first
    int a = i - 1;
    int b = j - 1;
    if(a >= 0 && b >= 0)
    {
        if(currScore == m_matrix[a][b] + match_score(s1[a],s2[b]))
            return (s1[a] == s2[b]) ? DPO_MATCH : DPO_SUB;
    }

    // Check above element
    a = i - 1;
    b = j;
    if(a >= 0 && b >= 0)
    {
        // If this is overlap mode, there is no deletion
        // penalty in the first column
        int above_score = (m_mode == DPM_OVERLAP && b == 0) ? 0 : gap_penalty();
        if(currScore == m_matrix[a][b] + above_score)
            return DPO_DELETE;
    }

    // Check left element
    a = i;
    b = j - 1;
    if(a >= 0 && b >= 0)
    {                  
        if(currScore == m_matrix[a][b] + gap_penalty())
            return DPO_INSERT;
    }
    assert(false);
    return DPO_NOOP;
}

// Calculate an optimal path through the alignment starting from the given cell
DPPath DPAlignment::calculatePath(const std::string& s1, const std::string& s2, 
                                  int startI, int startJ) const
{
    assert((int)s1.size() == m_dimM - 1 && (int)s2.size() == m_dimN - 1);
    DPPath out;
    DPPathNode currNode;
    currNode.i = startI;
    currNode.j = startJ;
    currNode.op = DPO_NOOP;

    while(currNode.i != 0 || currNode.j != 0)
    {
        // Set the node to be the cell that preceded the current one
        DPPathNode nextNode;

        // The path node op is the operation leading away from the cell
        nextNode.op = getPathOperationToCell(s1, s2, currNode.i, currNode.j);

        // Set the cell coordinates based on the op
        switch(nextNode.op)
        {
            case DPO_MATCH:
            case DPO_SUB:
                nextNode.i = currNode.i - 1;
                nextNode.j = currNode.j - 1;
                break;
            case DPO_DELETE:
                nextNode.i = currNode.i - 1;
                nextNode.j = currNode.j;
                break;
            case DPO_INSERT:
                nextNode.i = currNode.i;
                nextNode.j = currNode.j - 1;
                break;
            case DPO_NOOP:
                assert(false);
        }

        out.push_back(nextNode);
        currNode = nextNode;
    }
    out.reverse();
    return out;
}

//
void DPAlignment::printMatrix(const std::string& s1, const std::string& s2) const
{
    // Print the DP matrix
    std::cout << " \t-\t";
    for(int j = 1; j < m_dimN; ++j)
    {
        std::cout << s2[j - 1] << "\t";
    }
    std::cout << "\n";

    for(int i = 0; i < m_dimM; ++i)
    {
        if(i == 0)
            std::cout << "-\t";
        else
            std::cout << s1[i - 1] << "\t";

        for(int j = 0; j < m_dimN; ++j)
        {
            std::cout << m_matrix[i][j] << "\t";
        }
        std::cout << "\n";
    }
}

