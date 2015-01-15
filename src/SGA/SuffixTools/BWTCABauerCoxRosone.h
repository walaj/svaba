//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// BWTCABauerCoxRosone - In-memory version of Illumina's
// BWT construction algorithm
#ifndef BWTCA_COX_BAUER_ROSONE_H
#define BWTCA_COX_BAUER_ROSONE_H

#include "SAWriter.h"
#include "ReadTable.h"
#include "Alphabet.h"
#include "EncodedString.h"
#include "BWTWriterBinary.h"

namespace BWTCA
{
    // This structure tracks the state of each read
    // It contains a position field, which stores the relative/absolute
    // position to insert the next symbol of a read into the partial BWT.
    // It also stores the read index and the symbol to insert
    struct BCRElem
    {
        BCRElem() : index(0), sym('\0') {}

        // Comparator. Sort first on the symbol, then break ties using the index
        friend bool operator<(const BCRElem& a, const BCRElem& b)
        {
            return a.position < b.position;
        }

        friend std::ostream& operator<<(std::ostream& out, const BCRElem& elem)
        {
            out << "I: " << elem.index << " S: " << elem.sym << " P: " << elem.position;
            return out;
        }

        // Data
        uint64_t position; // the relative or absolute position into the next partial bwt to insert a symbol
        uint32_t index; // the read index this entry represents
        char sym; // a symbol to be inserted
    };
    typedef std::vector<BCRElem> BCRVector;

    // Construct the burrows-wheeler transform of the set of reads
    void runBauerCoxRosone(const DNAEncodedStringVector* pReadSequences, 
                           const std::string& bwt_out_name, 
                           const std::string& sai_out_name);
    
    // Run the initial special first iteration of the algorithm
    void outputInitialCycle(const DNAEncodedStringVector* pReadSequences, BCRVector& bcrVector, DNAEncodedString& bwt, AlphaCount64& suffixSymbolCounts);

    // Write the final bwt to a file
    size_t outputFinalBWT(BCRVector& bcrVector, 
                          const DNAEncodedString& readBWT, 
                          size_t partial_bwt_symbols, 
                          BWTWriterBinary* pBWTWriter,
                          SAWriter* pSAWriter);

    // Output an intermediate bwt 
    size_t outputPartialCycle(int cycle,
                             const DNAEncodedStringVector* pReadSequences, 
                             BCRVector& bcrVector, 
                             const DNAEncodedString& readBWT, 
                             size_t partial_bwt_symbols,
                             DNAEncodedString& writeBWT, 
                             AlphaCount64& suffixStartCounts);

    // Calculate absolute position for each element of the nvector
    void calculateAbsolutePositions(BCRVector& bcrVector, const AlphaCount64& suffixSymbolCounts);
};

#endif
