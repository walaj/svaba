//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// BWTCABauerCoxRosone - In-memory version of Illumina's
// BWT construction algorithm
#include "BWTCABauerCoxRosone.h"
#include "Timer.h"
#include "BWTWriterBinary.h"
#include <algorithm>
#include <iterator>

void BWTCA::runBauerCoxRosone(const DNAEncodedStringVector* pReadSequences,
                              const std::string& bwt_out_name, 
                              const std::string& sai_out_name)
{
    size_t num_reads = pReadSequences->size();

    size_t num_symbols = 0;
    for(size_t i = 0; i < num_reads; ++i)
        num_symbols += pReadSequences->at(i).length();
    num_symbols += num_reads; // include 1 sentinal per read
    printf("Running BCR on %zu symbols, %zu reads\n", num_symbols, num_reads);

    // Allocate two working BWTs
    DNAEncodedString read_bwt;
    DNAEncodedString write_bwt;

    read_bwt.resize(num_symbols);
    write_bwt.resize(num_symbols);

    // Allocate the bcr vector, which tracks the state of the algorithm
    BCRVector bcrVector(num_reads);
    
    // Track the number of suffixes that start with a given symbol
    AlphaCount64 suffixStartCounts;

    // Iteration 1:
    // Output a BWT with the last symbol of every read, in the order they appear in the read table.
    // This is the ordering of all suffixes that start with the sentinel character
    outputInitialCycle(pReadSequences, bcrVector, write_bwt, suffixStartCounts);
    write_bwt.swap(read_bwt);

    // Iteration 2...n, create new bwts from the bwt of the previous cycle
    size_t partial_bwt_length = num_reads;
    Timer timer("cycles", false);

    size_t maxCycles = pReadSequences->at(0).length();
    for(size_t cycle = 2; cycle <= maxCycles; ++cycle)
    {
        //std::cout << "Starting cycle " << cycle << "\n";
        
        // Convert relative positions to absolute positions
        calculateAbsolutePositions(bcrVector, suffixStartCounts);
       
        // Sort nvector by the absolute position to insert the symbol
        std::sort(bcrVector.begin(), bcrVector.end());
        //std::cout << "  done sorting..." << timer.getElapsedWallTime() << "\n";

        // Output the BWT for this cycle and update the vector and suffix start count
        partial_bwt_length = outputPartialCycle(cycle, pReadSequences, bcrVector, read_bwt, partial_bwt_length, write_bwt, suffixStartCounts);
        //std::cout << "  done writing..." << timer.getElapsedWallTime() << "\n";

        // Swap the in/out bwt
        read_bwt.swap(write_bwt);
    }

    // Write the resulting bwt and suffix array index
    BWTWriterBinary bwtWriter(bwt_out_name);
    bwtWriter.writeHeader(num_reads, num_symbols, BWF_NOFMI);

    SAWriter saWriter(sai_out_name);
    saWriter.writeHeader(num_reads, num_reads);

    // Calculate the positions of the final insertion symbols and write them directly to the file
    calculateAbsolutePositions(bcrVector, suffixStartCounts);
    std::sort(bcrVector.begin(), bcrVector.end());
    size_t num_wrote = outputFinalBWT(bcrVector, read_bwt, partial_bwt_length, &bwtWriter, &saWriter);
    assert(num_wrote == num_symbols);
}

// Update N and the output BWT for the initial cycle, corresponding to the sentinel suffixes
// the symbolCounts vector is updated to hold the number of times each symbol has been inserted
// into the bwt
void BWTCA::outputInitialCycle(const DNAEncodedStringVector* pReadSequences, BCRVector& bcrVector, DNAEncodedString& bwt, AlphaCount64& suffixSymbolCounts)
{
    AlphaCount64 incomingSymbolCounts;

    size_t n = pReadSequences->size();
    size_t first_read_len = pReadSequences->at(0).length();
    for(size_t i = 0; i < n; ++i)
    {
        size_t rl =  pReadSequences->at(i).length();
        
        // Check that all reads are the same length
        if(rl != first_read_len)
        {
            std::cout << "Error: This implementation of BCR requires all reads to be the same length\n";
            exit(EXIT_FAILURE);
        }

        char c = pReadSequences->at(i).get(rl - 1);
        bwt.set(i, c);

        assert(rl > 1);

        // Load the elements of the N vector with the next symbol
        bcrVector[i].sym = c;
        bcrVector[i].index = i;

        // Set the relative position of the symbol that is being inserted
        bcrVector[i].position = incomingSymbolCounts.get(c);

        // Increment the count of the first base of the suffix of the
        // incoming strings. This is $ for the initial cycle
        suffixSymbolCounts.increment('$');

        // Update the inserted symbols
        incomingSymbolCounts.increment(c);
    }

    suffixSymbolCounts += incomingSymbolCounts;
}

// Write out the next BWT for the next cycle. This updates BCRVector
// and suffixSymbolCounts. Returns the number of symbols written to writeBWT
size_t BWTCA::outputPartialCycle(int cycle,
                                 const DNAEncodedStringVector* pReadSequences,
                                 BCRVector& bcrVector, 
                                 const DNAEncodedString& readBWT, 
                                 size_t total_read_symbols,
                                 DNAEncodedString& writeBWT, 
                                 AlphaCount64& suffixStartCounts)
{
    // We track the rank of each symbol as it is copied/inserted
    // into the new bwt
    AlphaCount64 rank;

    // Counters
    size_t num_copied = 0;
    size_t num_inserted = 0;
    size_t num_wrote = 0;

    for(size_t i = 0; i < bcrVector.size(); ++i)
    {
        BCRElem& ne = bcrVector[i];
        
        // Copy elements from the read bwt until we reach the target position
        while(num_copied + num_inserted < ne.position)
        {
            char c = readBWT.get(num_copied++);
            writeBWT.set(num_wrote++, c);
            rank.increment(c);
        }

        // Now insert the incoming symbol
        int rl = pReadSequences->at(ne.index).length();
        char c = '$';

        // If the cycle number is greater than the read length, we are
        // on the final iteration and we just add in the '$' characters
        if(cycle <= rl)
            c = pReadSequences->at(ne.index).get(rl - cycle);
        //std::cout << "Inserting " << c << " at position " << num_copied + num_inserted << "\n";
        writeBWT.set(num_wrote++, c);
        num_inserted += 1;

        // Update the nvector element
        ne.sym = c;

        // Record the rank of the inserted symbol
        ne.position = rank.get(c);

        // Update the rank and the number of suffixes that start with c
        rank.increment(c);
        suffixStartCounts.increment(c);
    }

    // Copy any remaining symbols in the bwt
    while(num_copied < total_read_symbols)
        writeBWT.set(num_wrote++, readBWT.get(num_copied++));

    return num_wrote;
}

// Write the final BWT to a file.
size_t BWTCA::outputFinalBWT(BCRVector& bcrVector, 
                             const DNAEncodedString& readBWT, 
                             size_t partial_bwt_symbols,
                             BWTWriterBinary* pBWTWriter,
                             SAWriter* pSAWriter)
{

    // Counters
    size_t num_copied = 0;
    size_t num_inserted = 0;

    for(size_t i = 0; i < bcrVector.size(); ++i)
    {
        BCRElem& ne = bcrVector[i];
        
        // Copy elements from the read bwt until we reach the target position
        while(num_copied + num_inserted < ne.position)
            pBWTWriter->writeBWChar(readBWT.get(num_copied++));
        
        // Write a single $, terminating this string
        pBWTWriter->writeBWChar('$');
        pSAWriter->writeElem(SAElem(ne.index, 0));
        num_inserted += 1;
    }

    // Copy any remaining symbols in the bwt
    while(num_copied < partial_bwt_symbols)
        pBWTWriter->writeBWChar(readBWT.get(num_copied++));

    pBWTWriter->finalize();
    return num_copied + num_inserted;
}


void BWTCA::calculateAbsolutePositions(BCRVector& bcrVector, const AlphaCount64& suffixSymbolCounts)
{
    // Calculate a predecessor array from the suffix symbol counts
    AlphaCount64 predCounts;
    for(int i = 0; i < BWT_ALPHABET::size; ++i)
    {
        char b = RANK_ALPHABET[i];
        int64_t pc = suffixSymbolCounts.getLessThan(b);
        predCounts.set(b, pc);
    }

    for(size_t i = 0; i < bcrVector.size(); ++i)
        bcrVector[i].position += predCounts.get(bcrVector[i].sym);
}

