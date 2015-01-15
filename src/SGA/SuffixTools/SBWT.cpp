//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SBWT.cpp - Burrows Wheeler transform of a generalized suffix array
//
#include "SBWT.h"
#include "Timer.h"
#include "BWTReaderAscii.h"
#include "BWTReaderBinary.h"
#include "BWTWriterAscii.h"
#include <istream>
#include <queue>
#include <inttypes.h>

// macros
#define OCC(c,i) m_occurrence.get(m_bwStr, (c), (i))
#define PRED(c) m_predCount.get((c))

// Parse a BWT from a file
SBWT::SBWT(const std::string& filename, int sampleRate)
{
    BWTReaderBinary reader(filename);
    reader.read(this);
    initializeFMIndex(sampleRate);
}

// Construct the BWT from a suffix array
SBWT::SBWT(const SuffixArray* pSA, const ReadTable* pRT)
{
    size_t n = pSA->getSize();
    m_numStrings = pSA->getNumStrings();
    m_bwStr.resize(n);

    // Set up the bwt string and suffix array from the cycled strings
    for(size_t i = 0; i < n; ++i)
    {
        SAElem saElem = pSA->get(i);
        const SeqItem& si = pRT->getRead(saElem.getID());

        // Get the position of the start of the suffix
        uint64_t f_pos = saElem.getPos();
        uint64_t l_pos = (f_pos == 0) ? si.seq.length() : f_pos - 1;
        char b = (l_pos == si.seq.length()) ? '$' : si.seq.get(l_pos);
        m_bwStr.set(i, b);
    }

    initializeFMIndex(DEFAULT_SAMPLE_RATE);
}

// Fill in the FM-index data structures
void SBWT::initializeFMIndex(int sampleRate)
{
    // initialize the occurance table
    m_occurrence.initialize(m_bwStr, sampleRate);

    // Calculate the C(a) array
    
    // Calculate the total number of occurances of each character in the BW str
    AlphaCount64 tmp;
    for(size_t i = 0; i < m_bwStr.length(); ++i)
    {
        tmp.increment(m_bwStr.get(i));
    }

    m_predCount.set('$', 0);
    m_predCount.set('A', tmp.get('$')); 
    m_predCount.set('C', m_predCount.get('A') + tmp.get('A'));
    m_predCount.set('G', m_predCount.get('C') + tmp.get('C'));
    m_predCount.set('T', m_predCount.get('G') + tmp.get('G'));
}

// Compute the last to first mapping
size_t SBWT::LF(size_t idx) const
{
    return m_bwStr.get(idx) != '$' ? PRED(m_bwStr.get(idx)) + OCC(m_bwStr.get(idx), idx) : 0;
}

// Perform a exact search for the string w using the backwards algorithm
void SBWT::backwardSearch(std::string w) const
{
    std::cout << "Searching for " << w << "\n";
    int len = w.size();
    int j = len - 1;
    char curr = w[j];
    int r_lower = PRED(curr);
    int r_upper = r_lower + OCC(curr, m_bwStr.length() - 1) - 1;
    --j;
    std::cout << "Starting point: " << r_lower << "," << r_upper << "\n";
    for(;j >= 0; --j)
    {
        curr = w[j];
        printf("RL = C(%c) + O(%c,%d) + %zu\n", curr, curr, r_lower - 1, m_numStrings); 
        printf("RU = C(%c) + O(%c,%d)\n", curr, curr, r_upper); 
        printf("RL = %zu + %zu + %zu\n", (size_t)PRED(curr), (size_t)OCC(curr, r_lower - 1), m_numStrings); 
        printf("RU = %zu + %zu\n", (size_t)PRED(curr), (size_t)OCC(curr, r_upper)); 
        r_lower = PRED(curr) + OCC(curr, r_lower - 1);
        r_upper = PRED(curr) + OCC(curr, r_upper) - 1;
        printf("Curr: %c, Interval now: %d,%d\n", curr, r_lower, r_upper);
    }

    std::cout << "Interval found: " << r_lower << "," << r_upper << "\n";
}

//
void SBWT::validate() const
{
    std::cerr << "Warning BWT validation is turned on\n";
    m_occurrence.validate(m_bwStr);
}

// Print the BWT
void SBWT::print(const ReadTable* pRT, const SuffixArray* pSA) const
{
    std::cout << "i\tL(i)\tF(i)\tO(-,i)\tSUFF\n";
    for(size_t i = 0; i < m_bwStr.length(); ++i)
    {
        assert(getF(i) == pSA->getSuffix(i, pRT)[0]);
        std::cout << i << "\t" << m_bwStr.get(i) << "\t" << getF(i) << "\t" << m_occurrence.get(m_bwStr, i) << pSA->getSuffix(i, pRT) << "\n";
    }
}

// Print information about the BWT
void SBWT::printInfo() const
{
    size_t o_size = m_occurrence.getByteSize();
    size_t p_size = sizeof(m_predCount);

    size_t bwStr_size = m_bwStr.getMemSize();
    size_t offset_size = sizeof(m_numStrings);
    size_t total_size = o_size + p_size + bwStr_size + offset_size;
    double total_mb = ((double)total_size / (double)(1024 * 1024));
    printf("\nSBWT info\n");
    printf("Sample rate: %zu\n", m_occurrence.getSampleRate());
    printf("Memory -- OCC: %zu C: %zu Str: %zu Misc: %zu TOTAL: %zu (%lf MB)\n",
            o_size, p_size, bwStr_size, offset_size, total_size, total_mb);
    printf("N: %zu Bytes per symbol: %lf\n", m_bwStr.length(), (double)total_size / m_bwStr.length());
}
