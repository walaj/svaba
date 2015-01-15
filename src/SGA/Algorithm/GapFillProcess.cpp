///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// GapFillProcess - Fill in gaps in a scaffold
//
#include <algorithm>
#include "GapFillProcess.h"
#include "BWTAlgorithms.h"
#include "SGAlgorithms.h"
#include "SGSearch.h"
#include "StdAlnTools.h"
#include "HaplotypeBuilder.h"
#include "MultiAlignment.h"

//
GapFillStats::GapFillStats()
{
    numGapsAttempted = numGapsFilled = 0;
    for(size_t i = 0; i < GFRC_NUM_CODES; ++i)
        numFails[i] = 0;
}

//
void GapFillStats::print() const
{
    printf("Num gaps attempted: %zu\n", numGapsAttempted);
    printf("Num gaps filled: %zu\n", numGapsFilled);
    printf("    Failed -- no anchor: %zu\n", numFails[GFRC_NO_ANCHOR]);
    printf("    Failed -- no haplotype path: %zu\n", numFails[GFRC_NO_HAPLOTYPE_PATH]);
    printf("    Failed -- no haplotype, search aborted: %zu\n", numFails[GFRC_NO_HAPLOTYPE_ABORTED]);
    printf("    Failed -- no haplotype, walk failed: %zu\n", numFails[GFRC_NO_HAPLOTYPE_WALK_FAIL]);
    printf("    Failed -- ambiguous solution: %zu\n", numFails[GFRC_AMBIGUOUS]);
    printf("    Failed -- assembled size differs from estimate: %zu\n", numFails[GFRC_BAD_SIZE]);
    printf("    Failed -- anchor sequence mismatch: %zu\n", numFails[GFRC_BAD_TRIM]);
}

//
//
//
GapFillProcess::GapFillProcess(const GapFillParameters& params) : m_parameters(params)
{
}

//
GapFillProcess::~GapFillProcess()
{
    m_stats.print();
}

// Process the given scaffold, filling in any gaps found
GapFillResult GapFillProcess::processScaffold(const std::string& scaffold) const
{
    if(m_parameters.verbose > 0)
        std::cout << "Processing scaffold of length " << scaffold.length() << "\n";
    
    GapFillResult result;

    size_t len = scaffold.length();
    size_t currIdx = 0;
    while(currIdx < len)
    {
        char b = scaffold[currIdx];

        if(b != 'N')
        {
            // Append this non-gap character into the output scaffold
            result.scaffold.append(1,b);
            currIdx += 1;
        }
        else
        {
            // Found the start of a gap
            size_t gapLength = 0;
            while(scaffold[currIdx + gapLength] == 'N')
                gapLength += 1;

            if(m_parameters.verbose >= 1)
                printf("Constructing gap at position %zu GapLength: %zu\n", currIdx, gapLength);
            
            GapFillReturnCode code = GFRC_UNKNOWN;

            // Attempt to fill this gap starting with a long kmer, then relaxing the process
            for(size_t k = m_parameters.startKmer; k >= m_parameters.endKmer; k -= m_parameters.stride)
            {
                // Calculate the left-anchor using the sequence of the scaffold already appended
                AnchorSequence leftAnchor = findAnchor(k, result.scaffold, result.scaffold.length() - k, true);

                // Calculate the right anchor using the sequence of the input scaffold
                AnchorSequence rightAnchor = findAnchor(k, scaffold, currIdx + gapLength, false);

                // Estimate the size of the assembled sequence, including the flanking anchors
                int leftFlanking = result.scaffold.length() - leftAnchor.position;
                int rightFlankingPlusGap = rightAnchor.position + k - currIdx;
                int estimatedSize = leftFlanking + rightFlankingPlusGap;

                // Attempt to build the gap sequence
                std::string gapSequence;
                code = processGap(k, estimatedSize, leftAnchor, rightAnchor, gapSequence);

                if(code == GFRC_OK)
                {
                    // Successfully resolved the gap. Calculate the amount of the anchor sequence
                    // that is already present in the scaffold and must be removed. We remove this amount
                    // of sequence from the current scaffold, not the gapSequence, to account
                    // for the case that the gap sequence is significantly different than the sequence
                    // after the left anchor. 
                    size_t trimLength = result.scaffold.length() - leftAnchor.position;
                    result.scaffold.replace(leftAnchor.position, trimLength, gapSequence);
                        
                    // We need to update currIdx to point to the next base in the 
                    // input scaffold that is not already assembled. This is given
                    // by the position of the rightAnchor, plus a kmer
                    currIdx = rightAnchor.position + k;
                    m_stats.numGapsFilled += 1;
                    break; 
                }
            }

            if(code != GFRC_OK)
            {
                // Failed to resolve the gap. Append the gap into the growing scaffold
                m_stats.numFails[code] += 1;

                while(scaffold[currIdx] == 'N')
                {
                    result.scaffold.append(1, 'N');
                    currIdx += 1;
                }
            }
            m_stats.numGapsAttempted += 1;
        }
    }

    if(m_parameters.verbose >= 2)
        StdAlnTools::globalAlignment(scaffold, result.scaffold, true);
    return result;
}

// Fill in the specified gap
GapFillReturnCode GapFillProcess::processGap(size_t k, int estimatedSize, const AnchorSequence& startAnchor, const AnchorSequence& endAnchor, std::string& outSequence) const
{

    if(m_parameters.verbose > 0)
    {
        std::cout << "\tSTART: " << startAnchor << "\n";
        std::cout << "\tEND: " << endAnchor << "\n";
    }

    if(startAnchor.sequence.empty() || endAnchor.sequence.empty() || startAnchor.sequence == endAnchor.sequence)
        return GFRC_NO_ANCHOR;

    HaplotypeBuilder builder;
    builder.setTerminals(startAnchor, endAnchor);
    builder.setIndex(m_parameters.pBWT, m_parameters.pRevBWT);
    builder.setKmerParameters(k, m_parameters.kmerThreshold);
    HaplotypeBuilderReturnCode code = builder.run();
    
    HaplotypeBuilderResult result;

    // The search was successfull, build strings from the walks
    if(code == HBRC_OK)
        code = builder.parseWalks(result);

    if(code == HBRC_OK)
    {
        assert(!result.haplotypes.empty());
        // Calculate the estimated size of the sequence, including the anchors
        return selectGapSequence(estimatedSize, result.haplotypes, outSequence);
    }
    else 
    {
        // Map the haplotype builder code to a gap fill code
        assert(code != HBRC_OK);
        if(code == HBRC_NO_PATH)
            return GFRC_NO_HAPLOTYPE_PATH;
        else if(code == HBRC_TOO_MANY_VERTICES)
            return GFRC_NO_HAPLOTYPE_ABORTED;
        else if(code == HBRC_WALK_FAILED)
            return GFRC_NO_HAPLOTYPE_WALK_FAIL;
    }

    // should not reach here
    return GFRC_UNKNOWN;
}

// Find an anchor sequence to start the process of building the gap sequence
AnchorSequence GapFillProcess::findAnchor(size_t k, const std::string& scaffold, int64_t position, bool upstream) const
{
    AnchorSequence anchor;
    int64_t stride = upstream ? -1 : 1;
    int MAX_DISTANCE = 50;
    int64_t stop = upstream ? position - MAX_DISTANCE : position + MAX_DISTANCE;

    // Cap the travel distance to avoid out of bounds
    if(stop < 0)
        stop = 0;
    if(stop > (int64_t)(scaffold.length() - k))
        stop = scaffold.length() - k;

    for(; position != stop; position += stride)
    {
        assert(position >= 0);
        std::string testSeq = scaffold.substr(position, k);
        std::transform(testSeq.begin(), testSeq.end(), testSeq.begin(), ::toupper);
        if(testSeq.find_first_not_of("ACGT") != std::string::npos)
            continue;

        size_t count = BWTAlgorithms::countSequenceOccurrencesWithCache(testSeq, m_parameters.pBWT, m_parameters.pBWTCache);
        if(count > m_parameters.kmerThreshold)
        {
            anchor.sequence = testSeq;
            anchor.count = count;
            anchor.position = position;
            return anchor;
        }
    }

    anchor.sequence = "";
    anchor.count = -1;
    return anchor;
}

// Attempt to select one of the passed in strings as the gap sequence. If none fit the constraints,
// this sets gapSequence to the empty string and returns an error code
GapFillReturnCode GapFillProcess::selectGapSequence(int estimatedSize, const StringVector& sequences, std::string& gapSequence) const
{
    assert(!sequences.empty());
    int selectedIdx = -1;
    int selectedSizeDiff = std::numeric_limits<int>::max();

    for(size_t i = 0; i < sequences.size(); ++i)
    {
        int diff = abs(sequences[i].size() - estimatedSize);
        //printf("ES: %d S: %zu D: %d\n", estimatedSize, sequences[i].size(), diff);

        if(diff < selectedSizeDiff)
        {
            selectedSizeDiff = diff;
            selectedIdx = i;
        }
    }

    // Perform checks on the quality of the gap sequences
    int MAX_SIZE_DIFF = 100;
    if(selectedSizeDiff > MAX_SIZE_DIFF)
    {
        gapSequence = "";
        return GFRC_BAD_SIZE;
    }
    
    gapSequence = sequences[selectedIdx];
    return GFRC_OK;
}

