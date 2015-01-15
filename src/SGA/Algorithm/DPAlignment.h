//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// DPAlignment - Generic class for dynamic programming alignments 
//
#ifndef DPALIGNMENT_H
#define DPALIGNMENT_H
#include "Match.h"
#include <list>

inline int min3(int a, int b, int c)
{
    return std::min(std::min(a,b), c);
}

inline int max3(int a, int b, int c)
{
    return std::max(std::max(a,b), c);
}

// Scoring scheme used to calculate the minimum edit distance
// between two strings
class EditDistanceScoring
{
    public:
        static int match_score(char a, char b) { return (a == b) ? 0 : 1; };
        static int compare3(int a, int b, int c) { return min3(a,b,c); }
        static int gap_penalty() { return 1; }
};

// Scoring scheme used to calculate the similarity of two strings
class SimilarityScoring
{
    public:
        static int match_score(char a, char b) { return (a == b) ? 2 : -1; };
        static int compare3(int a, int b, int c) { return max3(a,b,c); }
        static int gap_penalty() { return -1; }
};

// Path data structures
enum DPOperation
{
    DPO_MATCH,
    DPO_SUB,
    DPO_INSERT,
    DPO_DELETE,
    DPO_NOOP
};

//
struct DPPathNode
{
    DPPathNode() : i(0), j(0), op(DPO_MATCH) {}
    int i;
    int j;

    // The operation leading to the next cell, not to this one
    DPOperation op;
};
typedef std::list<DPPathNode> DPPath;

//
enum DPAlignMode
{
    DPM_OVERLAP,
    DPM_ALIGNMENT  
};

//
enum DPScoringScheme
{
    DPSS_EDIT,
    DPSS_SIMILARITY
};

class DPAlignment
{
    public:
   
        //
        DPAlignment(const std::string& s1, const std::string& s2, DPAlignMode mode, DPScoringScheme scoring);

        //
        void printMatrix(const std::string& s1, const std::string& s2) const; 
        void printAlignment(const std::string& s1, const std::string& s2, int startI = -1, int startJ = -1) const;
        DPPath calculatePath(const std::string& s1, const std::string& s2, int startI, int startJ) const;

        //
        int getNumRows() const;
        int getNumColumns() const;
        int getScore(int i, int j) const;
        DPOperation getPathOperationToCell(const std::string& s1, const std::string& s2, int i, int j) const;

        //
        DPScoringScheme getScoringScheme() const;
    private:

        //
        inline int match_score(char a, char b) const
        {
            switch(m_scoring)
            {
                case DPSS_SIMILARITY:
                    return (a == b) ? 7 : -3; 
                case DPSS_EDIT:
                    return (a == b) ? 0 : 1; 
                default:
                    assert(false);
            }
            return 0;
        }

        //
        inline int gap_penalty() const
        {
            switch(m_scoring)
            {
                case DPSS_SIMILARITY:
                    return -10;
                case DPSS_EDIT:
                    return 1;
                default:
                    assert(false);
            }
            return 0;
        }

        // GCC 4.7 warns on this
#pragma GCC diagnostic ignored "-Wuninitialized"
        inline int compare3(int a, int b, int c) const
        {
            switch(m_scoring)
            {
                case DPSS_SIMILARITY:
                    return max3(a,b,c);
                case DPSS_EDIT:
                    return min3(a,b,c);
                default:
                    assert(false);
                    return a;
            }
            return a;
        }

        void initializeDPMatrix();
        void fillMatrix(const std::string& s1, const std::string& s2);
        
        typedef std::vector<int> IntVector;
        typedef std::vector<IntVector> IntMatrix;

        IntMatrix m_matrix;
        int m_dimM;
        int m_dimN;
        DPAlignMode m_mode;
        DPScoringScheme m_scoring;

};

#endif
