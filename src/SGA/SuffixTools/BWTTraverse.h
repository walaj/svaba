//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// BWTTraverse - bwt traversal algorithms
//
#ifndef BWTTRAVERSE_H
#define BWTTRAVERSE_H

#include "BWT.h"
#include "BWTAlgorithms.h"
#include <stack>

struct TraverseElem
{
    // The current index is set to 0 and we iterate forward until we find the first valid (non-zero) entry
    // This skips elem 0 which is the $ 
    TraverseElem(const BWTInterval& i, AlphaCount ac) : base_range(i), desc(ac), currIdx(0) { }

    // Return true if the current position is valid
    bool isValid() const
    {
        return currIdx < ALPHABET_SIZE;
    }

    // Return the current char
    char getCurrChar() const
    {
        assert(isValid());
        return desc.getBase(currIdx);
    }

    // Return the range
    const BWTInterval& getRange() const { return base_range; }

    // go to the next character with a non-zero count
    // this automatically skips index 0, which we want since its the '$' character
    void goNext()
    {
        do
        {
            ++currIdx;
        }
        while(currIdx < ALPHABET_SIZE && desc.getByIdx(currIdx) == 0);
    }

    // The ranges in the BWT that correspond to the string leading up to this element
    BWTInterval base_range;

    // The counts of each base {A,C,G,T} that are left-extensions of this range
    AlphaCount desc;

    // The index of the current base being processed
    int currIdx;
};

typedef std::stack<TraverseElem> TraverseStack;
typedef std::vector<bool> bool_vec;

namespace BWTTraverse
{

// Extract all strings of length len from the BWT
void extract(const BWT* pBWT, unsigned int len);

// Extract the string graph which minimum component length len from the BWT
void extractSG(const BWT* pBWT, const BWT* pRevBWT, const unsigned int len);
void extendLeft(const unsigned int len, std::string& str, bool_vec& visited, const BWT* pBWT, const BWT* pRevBWT);
void extendRight(const unsigned int len, std::string& str, bool_vec& visited, const BWT* pBWT, const BWT* pRevBWT, bool isReverse);
void markVisited(const std::string& str, bool_vec& visited, const BWT* pBWT);

};

#endif
