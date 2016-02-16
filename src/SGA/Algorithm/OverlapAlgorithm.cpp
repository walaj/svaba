//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
#include "OverlapAlgorithm.h"
#include "ASQG.h"
#include <math.h>

// Collect the complete set of overlaps in pOBOut
static const AlignFlags sufPreAF(false, false, false);
static const AlignFlags prePreAF(false, true, true);
static const AlignFlags sufSufAF(true, false, true);
static const AlignFlags preSufAF(true, true, false);

//#define DEBUGOVERLAP 2

// Perform the overlap
OverlapResult OverlapAlgorithm::overlapRead(const SeqRecord& read, int minOverlap, OverlapBlockList* pOutList) const
{
  
    OverlapResult r;
    if(static_cast<int>(read.seq.length()) < minOverlap)
        return r;

    if(!m_exactModeOverlap)
        r = overlapReadInexact(read, minOverlap, pOutList);
    else
        r = overlapReadExact(read, minOverlap, pOutList);
    return r;
}

//
OverlapResult OverlapAlgorithm::overlapReadInexact(const SeqRecord& read, int minOverlap, OverlapBlockList* pOBOut) const
{
    OverlapResult result;
    OverlapBlockList obWorkingList;
    std::string seq = read.seq.toString();

#ifdef DEBUGOVERLAP
    std::cout << "\n\n***Overlapping read " << read.id << " suffix\n";
    std::cout << "***Overlapping seq " << read.seq.toString() << " suffix\n";
#endif

    // Match the suffix of seq to prefixes

    // findInexact returns false is the maximum search time was exceeded. In this
    // case we dont run any of the subsequent commands and return no overlaps.
    bool valid = true;

    valid = findOverlapBlocksInexact(seq, m_pBWT, m_pRevBWT, sufPreAF, 
                                 minOverlap, &obWorkingList, pOBOut, result);

    if(valid)
       valid = findOverlapBlocksInexact(complement(seq), m_pRevBWT, m_pBWT, prePreAF, 
                                       minOverlap, &obWorkingList, pOBOut, result);


#ifdef DEBUGOVERLAP
    std::cout << "WORKING OverlapBlockList Length SUFFIX: " << obWorkingList.size() << std::endl;
#endif       

    if(valid)
    {
        if(m_bIrreducible)
        {
            computeIrreducibleBlocks(m_pBWT, m_pRevBWT, &obWorkingList, pOBOut);
            obWorkingList.clear();
        }
        else
        {
            pOBOut->splice(pOBOut->end(), obWorkingList);
            assert(obWorkingList.empty());
        }
    }

#ifdef DEBUGOVERLAP
    std::cout << "\n\n***Overlapping read " << read.id << " prefix\n";
#endif

    // Match the prefix of seq to suffixes
    if(valid)
        valid = findOverlapBlocksInexact(reverseComplement(seq), m_pBWT, m_pRevBWT, sufSufAF, minOverlap, &obWorkingList, pOBOut, result);
    
    if(valid)
        valid = findOverlapBlocksInexact(reverse(seq), m_pRevBWT, m_pBWT, preSufAF, minOverlap, &obWorkingList, pOBOut, result);

#ifdef DEBUGOVERLAP
    std::cout << "WORKING OverlapBlockList Length PREFIX: " << obWorkingList.size() << std::endl;
#endif

    if(valid)
    {
        if(m_bIrreducible)
        {
            computeIrreducibleBlocks(m_pBWT, m_pRevBWT, &obWorkingList, pOBOut);
            obWorkingList.clear();
        }
        else
        {
            pOBOut->splice(pOBOut->end(), obWorkingList);
            assert(obWorkingList.empty());
        }
	}
    
    if(!valid)
    {
        pOBOut->clear();
        result.isSubstring = false;
        result.searchAborted = true;
        return result;
    }

#ifdef DEBUGOVERLAP
     unsigned N = 1;
     if (pOBOut->size() > N)
      {
	std::cout << read.id << " OverlapBlockList Length: " << pOBOut->size() << std::endl;
	std::list<OverlapBlock>::iterator it = pOBOut->begin();
	std::advance(it, N);
        //std::cout << "Olen: " << it->overlapLen << " numDiff: " << it->numDiff << std::endl;
        std::cout << " BWT1L: " << it->ranges.interval[0].lower;
        std::cout << " BWT1U: " << it->ranges.interval[0].upper;
        std::cout << " BWT2L: " << it->ranges.interval[1].lower;
        std::cout << " BWT2U: " << it->ranges.interval[1].upper;
        std::cout << " RBWT1L: " << it->rawRanges.interval[0].lower;
        std::cout << " RBWT1U: " << it->rawRanges.interval[0].upper;
        std::cout << " RBWT2L: " << it->rawRanges.interval[1].lower;
        std::cout << " RBWT2U: " << it->rawRanges.interval[1].upper << std::endl;
      }
#endif


    return result;
}

//
OverlapResult OverlapAlgorithm::alignReadDuplicate(const SeqRecord& read, OverlapBlockList* pOBOut) const
{
    OverlapResult result;
    OverlapBlockList obWorkingList;
    std::string seq = read.seq.toString();
    int readLength = seq.length();

    findOverlapBlocksInexact(seq, m_pBWT, m_pRevBWT, sufPreAF, readLength, &obWorkingList, pOBOut, result);
    findOverlapBlocksInexact(complement(seq), m_pRevBWT, m_pBWT, prePreAF, readLength, &obWorkingList, pOBOut, result);
    return result;
}


// Construct the set of blocks describing irreducible overlaps with READ
// and write the blocks to pOBOut
OverlapResult OverlapAlgorithm::overlapReadExact(const SeqRecord& read, int minOverlap, OverlapBlockList* pOBOut) const
{
    OverlapResult result;
    // The complete set of overlap blocks are collected in obWorkingList
    // The filtered set (containing only irreducible overlaps) are placed into pOBOut
    // by calculateIrreducibleHits
    OverlapBlockList obWorkingList;
    std::string seq = read.seq.toString();

    // We store the various overlap blocks using a number of lists, one for the containments
    // in the forward and reverse index and one for each set of overlap blocks
    OverlapBlockList oblFwdContain;
    OverlapBlockList oblRevContain;
    
    OverlapBlockList oblSuffixFwd;
    OverlapBlockList oblSuffixRev;
    OverlapBlockList oblPrefixFwd;
    OverlapBlockList oblPrefixRev;

    // Match the suffix of seq to prefixes
    findOverlapBlocksExact(seq, m_pBWT, m_pRevBWT, sufPreAF, minOverlap, &oblSuffixFwd, &oblFwdContain, result);
    findOverlapBlocksExact(complement(seq), m_pRevBWT, m_pBWT, prePreAF, minOverlap, &oblSuffixRev, &oblRevContain, result);

    // Match the prefix of seq to suffixes
    findOverlapBlocksExact(reverseComplement(seq), m_pBWT, m_pRevBWT, sufSufAF, minOverlap, &oblPrefixFwd, &oblFwdContain, result);
    findOverlapBlocksExact(reverse(seq), m_pRevBWT, m_pBWT, preSufAF, minOverlap, &oblPrefixRev, &oblRevContain, result);

    // Remove submaximal blocks for each block list including fully contained blocks
    // Copy the containment blocks into the prefix/suffix lists
    oblSuffixFwd.insert(oblSuffixFwd.end(), oblFwdContain.begin(), oblFwdContain.end());
    oblPrefixFwd.insert(oblPrefixFwd.end(), oblFwdContain.begin(), oblFwdContain.end());
    oblSuffixRev.insert(oblSuffixRev.end(), oblRevContain.begin(), oblRevContain.end());
    oblPrefixRev.insert(oblPrefixRev.end(), oblRevContain.begin(), oblRevContain.end());
    
    // Perform the submaximal filter
    removeSubMaximalBlocks(&oblSuffixFwd, m_pBWT, m_pRevBWT);
    removeSubMaximalBlocks(&oblPrefixFwd, m_pBWT, m_pRevBWT);
    removeSubMaximalBlocks(&oblSuffixRev, m_pRevBWT, m_pBWT);
    removeSubMaximalBlocks(&oblPrefixRev, m_pRevBWT, m_pBWT);
    
    // Remove the contain blocks from the suffix/prefix lists
    removeContainmentBlocks(seq.length(), &oblSuffixFwd);
    removeContainmentBlocks(seq.length(), &oblPrefixFwd);
    removeContainmentBlocks(seq.length(), &oblSuffixRev);
    removeContainmentBlocks(seq.length(), &oblPrefixRev);

    // Join the suffix and prefix lists
    oblSuffixFwd.splice(oblSuffixFwd.end(), oblSuffixRev);
    oblPrefixFwd.splice(oblPrefixFwd.end(), oblPrefixRev);

    // Move the containments to the output list
    pOBOut->splice(pOBOut->end(), oblFwdContain);
    pOBOut->splice(pOBOut->end(), oblRevContain);

    // Filter out transitive overlap blocks if requested
    if(m_bIrreducible)
    {
        computeIrreducibleBlocks(m_pBWT, m_pRevBWT, &oblSuffixFwd, pOBOut);
        computeIrreducibleBlocks(m_pBWT, m_pRevBWT, &oblPrefixFwd, pOBOut);
    }
    else
    {
        pOBOut->splice(pOBOut->end(), oblSuffixFwd);
        pOBOut->splice(pOBOut->end(), oblPrefixFwd);
    }

    return result;
}

// Write overlap results to an ASQG file
void OverlapAlgorithm::writeResultASQG(std::ostream& writer, const SeqRecord& read, const OverlapResult& result) const
{
    ASQG::VertexRecord record(read.id, read.seq.toString());
    record.setSubstringTag(result.isSubstring);
    record.write(writer);
}

// Write overlap blocks out to a file
void OverlapAlgorithm::writeOverlapBlocks(std::ostream& writer, size_t readIdx, bool isSubstring, const OverlapBlockList* pList) const
{

    // Write the header info
    size_t numBlocks = pList->size();
    writer << readIdx << " " << isSubstring << " " << numBlocks << " ";
    //std::cout << readIdx << " " << isSubstring << " " << numBlocks << " ";
    //std::cout << "<Wrote> idx: " << count << " count: " << numBlocks << "\n";
    for(OverlapBlockList::const_iterator iter = pList->begin(); iter != pList->end(); ++iter)
    {
        writer << *iter << " ";
	//std::cout << *iter << " ";
    }
    writer << "\n";
    //std::cout << "\n";
}

//JEREMIAH
void OverlapAlgorithm::writeOverlapBlocks(std::stringstream& writer, size_t readIdx, bool isSubstring, const OverlapBlockList* pList) const
{
    // Write the header info
    size_t numBlocks = pList->size();
    writer << readIdx << " " << isSubstring << " " << numBlocks << " ";
    //std::cout << readIdx << " " << isSubstring << " " << numBlocks << " ";
    for(OverlapBlockList::const_iterator iter = pList->begin(); iter != pList->end(); ++iter)
    {
        writer << *iter << " ";
	//std::cout << *iter << " ";
    }
    writer << "\n";
    //std::cout << "\n";
}


// Calculate the ranges in pBWT that contain a prefix of at least minOverlap basepairs that
// overlaps with a suffix of w. The ranges are added to the pOBList
void OverlapAlgorithm::findOverlapBlocksExact(const std::string& w, const BWT* pBWT,
                                              const BWT* pRevBWT, const AlignFlags& af, int minOverlap,
                                              OverlapBlockList* pOverlapList, OverlapBlockList* pContainList, 
                                              OverlapResult& result) const
{
    // The algorithm is as follows:
    // We perform a backwards search using the FM-index for the string w.
    // As we perform the search we collect the intervals 
    // of the significant prefixes (len >= minOverlap) that overlap w.
    BWTIntervalPair ranges;
    size_t l = w.length();
    int start = l - 1;
    BWTAlgorithms::initIntervalPair(ranges, w[start], pBWT, pRevBWT);
    
    // Collect the OverlapBlocks
    for(size_t i = start - 1; i >= 1; --i)
    {
        // Compute the range of the suffix w[i, l]
        BWTAlgorithms::updateBothL(ranges, w[i], pBWT);
        int overlapLen = l - i;
        if(overlapLen >= minOverlap)
        {
            // Calculate which of the prefixes that match w[i, l] are terminal
            // These are the proper prefixes (they are the start of a read)
            BWTIntervalPair probe = ranges;
            BWTAlgorithms::updateBothL(probe, '$', pBWT);
            
            // The probe interval contains the range of proper prefixes
            if(probe.interval[1].isValid())
            {
                assert(probe.interval[1].lower > 0);
                pOverlapList->push_back(OverlapBlock(probe, ranges, overlapLen, 0, af));
            }
        }
    }

    // Determine if this sequence is contained and should not be processed further
    BWTAlgorithms::updateBothL(ranges, w[0], pBWT);

    // Ranges now holds the interval for the full-length read
    // To handle containments, we output the overlapBlock to the final overlap block list
    // and it will be processed later
    // Two possible containment cases:
    // 1) This read is a substring of some other read
    // 2) This read is identical to some other read
    
    // Case 1 is indicated by the existance of a non-$ left or right hand extension
    // In this case we return no alignments for the string
    AlphaCount64 left_ext = BWTAlgorithms::getExtCount(ranges.interval[0], pBWT);
    AlphaCount64 right_ext = BWTAlgorithms::getExtCount(ranges.interval[1], pRevBWT);
    if(left_ext.hasDNAChar() || right_ext.hasDNAChar())
    {
        result.isSubstring = true;
    }
    else
    {
        BWTIntervalPair probe = ranges;
        BWTAlgorithms::updateBothL(probe, '$', pBWT);
        if(probe.isValid())
        {
            // terminate the contained block and add it to the contained list
            BWTAlgorithms::updateBothR(probe, '$', pRevBWT);
            assert(probe.isValid());
            pContainList->push_back(OverlapBlock(probe, ranges, w.length(), 0, af));
        }
    }

    //OverlapBlockList containedWorkingList;
    //partitionBlockList(w.length(), &workingList, pOverlapList, &containedWorkingList);
    
    // Terminate the contained blocks
    //terminateContainedBlocks(containedWorkingList);
    
    // Move the contained blocks to the final contained list
    //pContainList->splice(pContainList->end(), containedWorkingList);

    return;
}

// Seeded blockwise BWT alignment of prefix-suffix for reads
// Each alignment is given a seed region and a block region
// The seed region is the terminal portion of w where maxDiff + 1 seeds are created
// at least 1 of these seeds must align exactly for there to be an alignment with 
// at most maxDiff differences between the prefix/suffix. Only alignments within the
// range [block_start, block_end] are output. The block_end coordinate is inclusive.
bool OverlapAlgorithm::findOverlapBlocksInexact(const std::string& w, const BWT* pBWT, 
                                                const BWT* pRevBWT, const AlignFlags& af, int minOverlap,
                                                OverlapBlockList* pOverlapList, OverlapBlockList* pContainList, 
                                                OverlapResult& result) const
{
    int len = w.length();
    int overlap_region_left = len - minOverlap;
    SearchSeedVector* pCurrVector = new SearchSeedVector;
    SearchSeedVector* pNextVector = new SearchSeedVector;
    OverlapBlockList workingList;
    SearchSeedVector::iterator iter;

    // Create and extend the initial seeds
    int actual_seed_length = m_seedLength;
    int actual_seed_stride = m_seedStride;
#ifdef DEBUGOVERLAP
    std::cout << "SeedLen: " << m_seedLength << " SeedStride: " << m_seedStride << std::endl;
#endif

    if(actual_seed_length == 0)
    {
        // Calculate a seed length and stride that will guarantee all overlaps
        // with error rate m_errorRate will be found
        calculateSeedParameters(w, minOverlap, actual_seed_length, actual_seed_stride);
    }

    assert(actual_seed_stride != 0);

    createSearchSeeds(w, pBWT, pRevBWT, actual_seed_length, actual_seed_stride, pCurrVector);
    extendSeedsExactRight(w, pBWT, pRevBWT, ED_RIGHT, pCurrVector, pNextVector);
    pCurrVector->clear();
    pCurrVector->swap(*pNextVector);
    assert(pNextVector->empty());

    int num_steps = 0;

    // Perform the inexact extensions
    bool fail = false;
    while(!pCurrVector->empty())
    {
        if(m_maxSeeds != -1 && (int)pCurrVector->size() > m_maxSeeds)
        {
            fail = true;
            break;
        }

        iter = pCurrVector->begin();
        while(iter != pCurrVector->end())
        {
            SearchSeed& align = *iter;

            // If the current aligned region is right-terminal
            // and the overlap is greater than minOverlap, try to find overlaps
            // or containments
            if(align.right_index == len - 1)
            {
                double align_error = align.calcErrorRate();

                // Check for overlaps
                if(align.left_index <= overlap_region_left && isErrorRateAcceptable(align_error, m_errorRate))
                {
                    int overlapLen = len - align.left_index;
                    BWTIntervalPair probe = align.ranges;
                    BWTAlgorithms::updateBothL(probe, '$', pBWT);
                    
                    // The probe interval contains the range of proper prefixes
                    if(probe.interval[1].isValid())
                    {
                        assert(probe.interval[1].lower > 0);
                        OverlapBlock nBlock(probe, align.ranges, overlapLen, align.z, af, align.historyLink->getHistoryVector());
                        workingList.push_back(nBlock);
                    }
                }

                // Check for containments
                // If the seed is left-terminal and there are [ACGT] left/right extensions of the sequence
                // this read must be a substring of another read
                if(align.left_index == 0)
                {
                    AlphaCount64 left_ext = BWTAlgorithms::getExtCount(align.ranges.interval[0], pBWT);
                    AlphaCount64 right_ext = BWTAlgorithms::getExtCount(align.ranges.interval[1], pRevBWT);
                    if(left_ext.hasDNAChar() || right_ext.hasDNAChar())
                        result.isSubstring = true;
                }
            }

            // Extend the seed to the right/left
            if(align.dir == ED_RIGHT)
                extendSeedInexactRight(align, w, pBWT, pRevBWT, pNextVector);
            else
                extendSeedInexactLeft(align, w, pBWT, pRevBWT, pNextVector);
            ++iter;
            //pCurrVector->erase(iter++);
        }
        pCurrVector->clear();
        assert(pCurrVector->empty());
        pCurrVector->swap(*pNextVector);

        // Remove identical seeds after we have performed seed_len steps
        // as there now might be redundant seeds
        if(num_steps % actual_seed_stride == 0)
        {
            std::sort(pCurrVector->begin(), pCurrVector->end(), SearchSeed::compareLeftRange);
            SearchSeedVector::iterator end_iter = std::unique(pCurrVector->begin(), pCurrVector->end(), 
                                                                   SearchSeed::equalLeftRange);
            pCurrVector->resize(end_iter - pCurrVector->begin());
        }
        ++num_steps;
    }

    if(!fail)
    {
        // parse the working list to remove any submaximal overlap blocks
        // these blocks correspond to reads that have multiple valid overlaps. 
        // we only keep the longest
        removeSubMaximalBlocks(&workingList, pBWT, pRevBWT);

        OverlapBlockList containedWorkingList;
        partitionBlockList(len, &workingList, pOverlapList, &containedWorkingList);
        
        // Terminate the contained blocks
        terminateContainedBlocks(containedWorkingList);
        
        // Move the contained blocks to the final contained list
        pContainList->splice(pContainList->end(), containedWorkingList);
    }

    delete pCurrVector;
    delete pNextVector;
    return !fail;
}

// Build forward history for the blocks
void OverlapAlgorithm::buildForwardHistory(OverlapBlockList* pList) const
{
    OverlapBlockList terminalList;
    OverlapBlockList potentialContainedList;
    while(!pList->empty())
    {
        // The terminalBlock list contains all the blocks that became right-terminal
        // in the current extension round.

        // Perform a single round of extension, any terminal blocks
        // are moved to the terminated list
        extendActiveBlocksRight(m_pBWT, m_pRevBWT, *pList, terminalList, potentialContainedList);
    }
    *pList = terminalList;
    pList->splice(pList->end(), potentialContainedList);
}

// Calculate the single right extension to the '$' for each the contained blocks
// so that the interval ranges are consistent
void OverlapAlgorithm::terminateContainedBlocks(OverlapBlockList& containedBlocks) const
{
    for(OverlapBlockList::iterator iter = containedBlocks.begin(); iter != containedBlocks.end(); ++iter)
        BWTAlgorithms::updateBothR(iter->ranges, '$', iter->getExtensionBWT(m_pBWT, m_pRevBWT));
}

// Calculate the seed length and stride to ensure that we will find all 
// all overlaps with error rate at most m_errorRate.
// To overlap two reads allowing for d errors, we create d+1 seeds so that at least one seed will match
// exactly. d is a function of the overlap length so we define the seed length using the minimum overlap
// parameter. We then tile seeds across the read starting from the end such that for every overlap length
// x, there are at least floor(x * error_rate) + 1 seeds.
void OverlapAlgorithm::calculateSeedParameters(const std::string& w, const int minOverlap, int& seed_length, int& seed_stride) const
{

    int read_len = w.length();
    seed_length = 0;
    
    // The maximum possible number of differences occurs for a fully-aligned read
    int max_diff_high = static_cast<int>(m_errorRate * read_len);

    // Calculate the seed length to use
    // If the error rate is so low that no differences are possible just seed
    // over the entire minOverlap region
    if(max_diff_high > 0)
    {
        // Calculate the maximum number of differences between two sequences that overlap
        // by minOverlap
        int max_diff_low = static_cast<int>(m_errorRate * minOverlap);

         if(max_diff_low == 0)
            max_diff_low = 1;
         
         int seed_region_length = static_cast<int>(ceil(max_diff_low / m_errorRate));
         int num_seeds_low = max_diff_low + 1;
         seed_length = static_cast<int>(seed_region_length / num_seeds_low);
         if(seed_length > static_cast<int>(minOverlap))
            seed_length = minOverlap;
    }
    else
    {
        seed_length = minOverlap;
    }
    seed_stride = seed_length;    
}

// Create and intialize the search seeds
int OverlapAlgorithm::createSearchSeeds(const std::string& w, const BWT* pBWT, 
                                        const BWT* pRevBWT, int seed_length, int seed_stride,
                                        SearchSeedVector* pOutVector) const
{
    // Start a new chain of history links
    SearchHistoryLink rootLink = SearchHistoryNode::createRoot();

    // The maximum possible number of differences occurs for a fully-aligned read
    int read_len = w.length();
    int max_diff_high = static_cast<int>(m_errorRate * read_len);

    // Start the seeds at the end of the read
    int seed_start = read_len - seed_length;

    while(seed_start >= 0)
    {
        SearchSeed seed;
        seed.left_index = seed_start;
        seed.right_index = seed_start;
        seed.dir = ED_RIGHT;
        seed.seed_len = seed_length;
        seed.z = 0;
        seed.maxDiff = max_diff_high;
        seed.historyLink = rootLink;

        // Initialize the left and right suffix array intervals
        char b = w[seed.left_index];
        BWTAlgorithms::initIntervalPair(seed.ranges, b, pBWT, pRevBWT);        
        pOutVector->push_back(seed);

        seed_start -= seed_stride;
        
        // Only create one seed in the exact case
        if(max_diff_high == 0)
            break;
    }
    return seed_length;
}

// Extend all the seeds in pInVector to the right over the entire seed range
void OverlapAlgorithm::extendSeedsExactRightQueue(const std::string& w, const BWT* /*pBWT*/, const BWT* pRevBWT,
                                             ExtendDirection /*dir*/, const SearchSeedVector* pInVector, 
                                             SearchSeedQueue* pOutQueue) const
{
    for(SearchSeedVector::const_iterator iter = pInVector->begin(); iter != pInVector->end(); ++iter)
    {
        SearchSeed align = *iter;
        bool valid = true;
        while(align.isSeed())
        {
            ++align.right_index;
            char b = w[align.right_index];
            BWTAlgorithms::updateBothR(align.ranges, b, pRevBWT);
            if(!align.isIntervalValid(RIGHT_INT_IDX))
            {
                valid = false;
                break;
            }
        }

        if(valid)
            pOutQueue->push(align);
    }
}

// Extend all the seeds in pInVector to the right over the entire seed range
void OverlapAlgorithm::extendSeedsExactRight(const std::string& w, const BWT* /*pBWT*/, const BWT* pRevBWT,
                                             ExtendDirection /*dir*/, const SearchSeedVector* pInVector, 
                                             SearchSeedVector* pOutVector) const
{
    for(SearchSeedVector::const_iterator iter = pInVector->begin(); iter != pInVector->end(); ++iter)
    {
        SearchSeed align = *iter;
        bool valid = true;
        while(align.isSeed())
        {
            ++align.right_index;
            char b = w[align.right_index];
            BWTAlgorithms::updateBothR(align.ranges, b, pRevBWT);
            if(!align.isIntervalValid(RIGHT_INT_IDX))
            {
                valid = false;
                break;
            }
        }

        //std::cout << "Initial seed: ";
        //align.print(w);

        if(valid)
            pOutVector->push_back(align);
    }
}

//
void OverlapAlgorithm::extendSeedInexactRight(SearchSeed& seed, const std::string& w, const BWT* /*pBWT*/, 
                                              const BWT* pRevBWT, SearchSeedVector* pOutVector) const
{
    // If this alignment has run all the way to the end of the sequence
    // switch it to be a left extension sequence
    if(seed.right_index == int(w.length() - 1))
    {
        seed.dir = ED_LEFT;
        pOutVector->push_back(seed);
        return;
    }

    ++seed.right_index;
    
    if(!seed.allowMismatches())
    {
        char b = w[seed.right_index];
        BWTAlgorithms::updateBothR(seed.ranges, b, pRevBWT);
        if(seed.isIntervalValid(RIGHT_INT_IDX))
            pOutVector->push_back(seed);
    }
    else
    {
        // Calculating the AlphaCounts is the heavy part of the computation so we cache
        // the value outside the loop, it is the same for all bases
        AlphaCount64 lower = pRevBWT->getFullOcc(seed.ranges.interval[1].lower - 1);
        AlphaCount64 upper = pRevBWT->getFullOcc(seed.ranges.interval[1].upper);
        for(int i = 0; i < 4; ++i)
        {
            char b = ALPHABET[i];    
            BWTIntervalPair probe = seed.ranges;
            BWTAlgorithms::updateBothR(probe, b, pRevBWT, lower, upper);
            if(probe.interval[RIGHT_INT_IDX].isValid())
            {
                SearchSeed branched = seed;
                branched.ranges = probe;
                if(b != w[seed.right_index])
                {
                    ++branched.z;
                    branched.historyLink = seed.historyLink->createChild(w.length() - seed.right_index, b);
                }
                pOutVector->push_back(branched);
            }
        }
    }
}

//
void OverlapAlgorithm::extendSeedInexactLeft(SearchSeed& seed, const std::string& w, 
                                             const BWT* pBWT, const BWT* /*pRevBWT*/,
                                             SearchSeedVector* pOutVector) const
{
    //printf("ProcessingLEFT: "); align.print(w);
    --seed.left_index;
    if(seed.left_index >= 0)
    {
        if(!seed.allowMismatches())
        {
            // Extend exact
            char b = w[seed.left_index];
            BWTAlgorithms::updateBothL(seed.ranges, b, pBWT);
            if(seed.isIntervalValid(LEFT_INT_IDX))
                pOutVector->push_back(seed);
        }
        else
        {
            // Calculating the AlphaCounts is the heavy part of the computation so we cache
            // the value outside the loop, it is the same for all bases
            AlphaCount64 lower = pBWT->getFullOcc(seed.ranges.interval[0].lower - 1);
            AlphaCount64 upper = pBWT->getFullOcc(seed.ranges.interval[0].upper);
            for(int i = 0; i < 4; ++i)
            {
                char b = ALPHABET[i];
                BWTIntervalPair probe = seed.ranges;
                BWTAlgorithms::updateBothL(probe, b, pBWT, lower, upper);
                if(probe.interval[LEFT_INT_IDX].isValid())
                {
                    SearchSeed branched = seed;
                    branched.ranges = probe;
                    if(b != w[seed.left_index])
                    {
                        ++branched.z;
                        // The history coordinates are wrt the right end of the read
                        // so that each position corresponds to the length of the overlap
                        // including that position        
                        branched.historyLink = seed.historyLink->createChild(w.length() - seed.left_index, b);                
                        //branched.history.add(w.length() - seed.left_index, b);
                    }
                    pOutVector->push_back(branched);
                }
            }
        }
    }
}

// Extend the seed to the right, return true if the seed is still a valid range
bool OverlapAlgorithm::extendSeedExactRight(SearchSeed& seed, const std::string& w, const BWT* /*pBWT*/, const BWT* pRevBWT) const
{
    // If this alignment has run all the way to the end of the sequence
    // switch it to be a left extension sequence
    if(seed.right_index == int(w.length() - 1))
    {
        seed.dir = ED_LEFT;
        return true;
    }

    ++seed.right_index;
    char b = w[seed.right_index];
    BWTAlgorithms::updateBothR(seed.ranges, b, pRevBWT);
    if(seed.isIntervalValid(RIGHT_INT_IDX))
        return true;
    else
        return false;
}

//
bool OverlapAlgorithm::extendSeedExactLeft(SearchSeed& seed, const std::string& w, const BWT* pBWT, const BWT* /*pRevBWT*/) const
{
    --seed.left_index;
    if(seed.left_index >= 0)
    {
        char b = w[seed.left_index];
        BWTAlgorithms::updateBothL(seed.ranges, b, pBWT);
        if(seed.isIntervalValid(LEFT_INT_IDX))
            return true;
        else
            return false;
    }
    else
    {
        return false;
    }
}

// Calculate the irreducible blocks from the vector of OverlapBlocks
void OverlapAlgorithm::computeIrreducibleBlocks(const BWT* pBWT, const BWT* pRevBWT, 
                                                OverlapBlockList* pOBList, 
                                                OverlapBlockList* pOBFinal) const
{
    // processIrreducibleBlocks requires the pOBList to be sorted in descending order
    pOBList->sort(OverlapBlock::sortSizeDescending);

#ifdef DEBUGOVERLAP
	std::cout << " m_exactModeIrreducible " << m_exactModeIrreducible << "\n";
#endif       

    if(m_exactModeIrreducible)
        _processIrreducibleBlocksExactIterative(pBWT, pRevBWT, *pOBList, pOBFinal);
    else
        _processIrreducibleBlocksInexact(pBWT, pRevBWT, *pOBList, pOBFinal);
    pOBList->clear();
}

// iterate through obList and determine the overlaps that are irreducible. This function is recursive.
// The final overlap blocks corresponding to irreducible overlaps are written to pOBFinal.
// Invariant: the blocks are ordered in descending order of the overlap size so that the longest overlap is first.
// Invariant: each block corresponds to the same extension of the root sequence w.
void OverlapAlgorithm::_processIrreducibleBlocksExactIterative(const BWT* pBWT, const BWT* pRevBWT, 
                                                               OverlapBlockList& inList, 
                                                               OverlapBlockList* pOBFinal) const
{
    if(inList.empty())
        return;
    
    // We store the overlap blocks in groups of blocks that have the same right-extension.
    // When a branch is found, the groups are split based on the extension
    typedef std::list<OverlapBlockList> BlockGroups;

    BlockGroups blockGroups;
    blockGroups.push_back(inList);
    int numExtensions = 0;
    int numBranches = 0;
    while(!blockGroups.empty())
    {
        // Perform one extenion round for each group.
        // If the top-level block has ended, push the result
        // to the final list and remove the group from processing
        BlockGroups::iterator groupIter = blockGroups.begin();
        BlockGroups incomingGroups; // Branched blocks are placed here

        while(groupIter != blockGroups.end())
        {
            OverlapBlockList& currList = *groupIter;
            bool bEraseGroup = false;

            // Count the extensions in the top level (longest) blocks first
            int topLen = currList.front().overlapLen;
            AlphaCount64 ext_count;
            OBLIter blockIter = currList.begin();
            while(blockIter != currList.end() && blockIter->overlapLen == topLen)
            {
                ext_count += blockIter->getCanonicalExtCount(pBWT, pRevBWT);
                ++blockIter;
            }
            
            // Three cases:
            // 1) The top level block has ended as it contains the extension $. Output TLB and end.
            // 2) There is a singular unique extension base for all the blocks. Update the blocks and continue.
            // 3) There are multiple extension bases, split the block group and continue.
            // If some block other than the TLB ended, it must be contained within the TLB and it is not output
            // or considered further. 
            // Likewise if multiple distinct strings in the TLB ended, we only output the top one. The rest
            // must have the same sequence as the top one and are hence considered to be contained with the top element.

#ifdef DEBUGOVERLAP
	    std::cout << "ext_gcount.get($) " << ext_count.get('$') << "\n";
#endif             

            if(ext_count.get('$') > 0)
            {
                // An irreducible overlap has been found. It is possible that there are two top level blocks
                // (one in the forward and reverse direction). Since we can't decide which one
                // contains the other at this point, we output hits to both. Under a fixed 
                // length string assumption one will be contained within the other and removed later.
                OBLIter tlbIter = currList.begin();
                while(tlbIter != currList.end() && tlbIter->overlapLen == topLen)
                {
                    // Ensure the tlb is actually terminal and not a substring block
                    AlphaCount64 test_count = tlbIter->getCanonicalExtCount(pBWT, pRevBWT);
                    if(test_count.get('$') == 0)
                    {
                        std::cerr << "Error: substring read found during overlap computation.\n";
                        std::cerr << "Please run sga rmdup before sga overlap\n";
                        exit(EXIT_FAILURE);
                    }
                    
                    // Perform the final right-update to make the block terminal
                    OverlapBlock branched = *tlbIter;
                    BWTAlgorithms::updateBothR(branched.ranges, '$', branched.getExtensionBWT(pBWT, pRevBWT));
                    pOBFinal->push_back(branched);
#ifdef DEBUGOVERLAP
                    std::cout << "[IE] TLB of length " << branched.overlapLen << " has ended\n";
                    std::cout << "[IE]\tBlock data: " << branched << "\n";
#endif             
                    ++tlbIter;
                } 

                // Set the flag to erase this group, it is finished
                bEraseGroup = true;
            }
            else
            {
                // Count the extension for the rest of the blocks
                while(blockIter != currList.end())
                {
                    ext_count += blockIter->getCanonicalExtCount(pBWT, pRevBWT);
                    ++blockIter;
                }

                if(ext_count.hasUniqueDNAChar())
                {
                    // Update all the blocks using the unique extension character
                    // This character is in the canonical representation wrt to the query
                    char b = ext_count.getUniqueDNAChar();
                    updateOverlapBlockRangesRight(pBWT, pRevBWT, currList, b);
                    numExtensions++;
                    bEraseGroup = false;
                }
                else
                {
                    for(size_t idx = 0; idx < DNA_ALPHABET_SIZE; ++idx)
                    {
                        char b = ALPHABET[idx];
                        if(ext_count.get(b) > 0)
                        {
                            numBranches++;
                            OverlapBlockList branched = currList;
                            updateOverlapBlockRangesRight(pBWT, pRevBWT, branched, b);
                            incomingGroups.push_back(branched);
                            bEraseGroup = true;
                        }
                    }
                }
            }

            if(bEraseGroup)
                groupIter = blockGroups.erase(groupIter);
            else
                ++groupIter;
        }

        // Splice in the newly branched blocks, if any
        blockGroups.splice(blockGroups.end(), incomingGroups);
    }
}

// Classify the blocks in obList as irreducible, transitive or substrings. The irreducible blocks are
// put into pOBFinal. The remaining are discarded.
// Invariant: the blocks are ordered in descending order of the overlap size so that the longest overlap is first.
void OverlapAlgorithm::_processIrreducibleBlocksInexact(const BWT* pBWT, const BWT* pRevBWT, 
                                                        OverlapBlockList& activeList, 
                                                        OverlapBlockList* pOBFinal) const
{
    if(activeList.empty())
        return;

#ifdef DEBUGOVERLAP
	std::cout << " activeList.size() " << activeList.size() << "\n";
#endif       
    
    // The activeList contains all the blocks that are not yet right terminal
    // Count the extensions in the top level (longest) blocks first
    bool all_eliminated = false;
    while(!activeList.empty() && !all_eliminated)
    {
        // The terminalBlock list contains all the blocks that became right-terminal
        // in the current extension round.
        OverlapBlockList terminalList;
        OverlapBlockList potentialContainedList;

        // Perform a single round of extension, any terminal blocks
        // are moved to the terminated list
        extendActiveBlocksRight(pBWT, pRevBWT, activeList, terminalList, potentialContainedList);

        // Compare the blocks in the contained list against the other terminal and active blocks
        // If they are a substring match to any of these, discard them
        OverlapBlockList::iterator containedIter = potentialContainedList.begin();
        for(; containedIter != potentialContainedList.end(); ++containedIter)
        {
           if(!isBlockSubstring(*containedIter, terminalList, m_errorRate) && 
              !isBlockSubstring(*containedIter, activeList, m_errorRate))
           {
                // Not a substring, move to terminal list
                terminalList.push_back(*containedIter);
                //std::cout << "Contained block kept: " << containedIter->overlapLen << "\n";
           }
           else
           {
                //std::cout << "Contained block found and removed: " << containedIter->overlapLen << "\n";
           }
        }

        // Using the terminated blocks, mark as eliminated any active blocks
        // that form a valid overlap to the terminal block. These are transitive edges
        // We do not compare two terminal blocks, we don't consider these overlaps to be
        // transitive
#ifdef DEBUGOVERLAP
	std::cout << " terminalList.sizee() " << terminalList.size() << "\n";
#endif       

        OverlapBlockList::iterator terminalIter = terminalList.begin();
        for(; terminalIter != terminalList.end(); ++terminalIter)
        {
#ifdef DEBUGOVERLAP
            std::cout << "[II] ***TLB of length " << terminalIter->overlapLen << " has ended\n";
#endif       
            all_eliminated = true;
            OverlapBlockList::iterator activeIter = activeList.begin();
            for(; activeIter != activeList.end(); ++activeIter)
            {
                if(activeIter->isEliminated)
                    continue; // skip previously marked blocks
                
                // Two conditions must be met for a block to be transitive wrt terminal:
                // 1) It must have a strictly shorter overlap than the terminal block
                // 2) The error rate between the block and terminal must be less than the threshold
                double inferredErrorRate = calculateBlockErrorRate(*terminalIter, *activeIter);
                if(activeIter->overlapLen < terminalIter->overlapLen && 
                   isErrorRateAcceptable(inferredErrorRate, m_errorRate))
                {
#ifdef DEBUGOVERLAP_2                            
                    std::cout << "Marking block of length " << activeIter->overlapLen << " as eliminated\n";
#endif
                    activeIter->isEliminated = true;
                }
                else
                {
                    all_eliminated = false;
                }
            } 
            
            // Move this block to the final list if it has not been previously marked eliminated
            if(!terminalIter->isEliminated)
            {
#ifdef DEBUGOVERLAP
                std::cout << "[II] Adding block " << *terminalIter << " to final list\n";
                //std::cout << "  extension: " << terminalIter->forwardHistory << "\n";
#endif                
                pOBFinal->push_back(*terminalIter);
            }
        }
    }

    activeList.clear();
}

// Extend all the blocks in activeList by one base to the right
// Move all right-terminal blocks to the termainl list. If a block 
// is terminal and potentially contained by another block, add it to 
// containedList
void OverlapAlgorithm::extendActiveBlocksRight(const BWT* pBWT, const BWT* pRevBWT, 
                                               OverlapBlockList& activeList, 
                                               OverlapBlockList& terminalList,
                                               OverlapBlockList& /*containedList*/) const
{
    OverlapBlockList::iterator iter = activeList.begin();
    OverlapBlockList::iterator next;
    while(iter != activeList.end())
    {
        next = iter;
        ++next;

        // Check if block is terminal
        AlphaCount64 ext_count = iter->getCanonicalExtCount(pBWT, pRevBWT);
        if(ext_count.get('$') > 0)
        {
            // Only consider this block to be terminal irreducible if it has at least one extension
            // or else it is a substring block
            if(iter->forwardHistory.size() > 0)
            {
                OverlapBlock branched = *iter;
                BWTAlgorithms::updateBothR(branched.ranges, '$', branched.getExtensionBWT(pBWT, pRevBWT));
                terminalList.push_back(branched);
#ifdef DEBUGOVERLAP_2            
                std::cout << "Block of length " << iter->overlapLen << " moved to terminal\n";
#endif
            }
        }

        int curr_extension = iter->forwardHistory.size();

        // Perform the right extensions
        
        // Best case, there is only a single extension character
        // Handle this case specially so we don't need to copy the potentially
        // large OverlapBlock structure and its full history
        if(ext_count.hasUniqueDNAChar())
        {
            // Get the extension character with respect to the queried sequence
            char canonical_base = ext_count.getUniqueDNAChar();

            // Flip the base into the frame of reference for the block
            char block_base = iter->flags.isQueryComp() ? complement(canonical_base) : canonical_base;

            // Update the block using the base in its frame of reference
            BWTAlgorithms::updateBothR(iter->ranges, block_base, iter->getExtensionBWT(pBWT, pRevBWT));

            // Add the base to the history in the frame of reference of the query read
            // This is so the history is consistent when comparing between blocks from different strands
            iter->forwardHistory.add(curr_extension, canonical_base);
        }
        else
        {
            for(size_t idx = 0; idx < DNA_ALPHABET_SIZE; ++idx)
            {
                char canonical_base = ALPHABET[idx];
                char block_base = iter->flags.isQueryComp() ? complement(canonical_base) : canonical_base;
                if(ext_count.get(canonical_base) == 0)
                    continue;

                // Branch the sequence. This involves copying the entire history which can be large
                // if the input sequences are very long. This could be avoided by using the SearchHistoyNode/Link
                // structure but branches are infrequent enough to not have a large impact
                OverlapBlock branched = *iter;
                BWTAlgorithms::updateBothR(branched.ranges, block_base, branched.getExtensionBWT(pBWT, pRevBWT));
                assert(branched.ranges.isValid());

                // Add the base in the canonical frame
                branched.forwardHistory.add(curr_extension, canonical_base);

                // Insert the new block after the iterator
                activeList.insert(iter, branched);
            }

            // Remove the original block, which has been superceded by the branches
            activeList.erase(iter);
        }

        iter = next; // this skips the newly-inserted blocks
    }
} 

// Return true if the terminalBlock is a substring of any member of blockList
bool OverlapAlgorithm::isBlockSubstring(OverlapBlock& terminalBlock, const OverlapBlockList& blockList, double maxER) const
{
    OverlapBlockList::const_iterator iter = blockList.begin();
    size_t right_extension_length = terminalBlock.forwardHistory.size();
    for(; iter != blockList.end(); ++iter)
    {
        if(terminalBlock.overlapLen == iter->overlapLen && 
           right_extension_length == iter->forwardHistory.size())
        {
            continue; // same length, cannot be a substring
        }
        
        // Calculate error rate between blocks
        double er = calculateBlockErrorRate(terminalBlock, *iter);
        if(isErrorRateAcceptable(er, maxER))
            return true;
    }
    return false;
}

// Calculate the error rate between two overlap blocks using their history
double OverlapAlgorithm::calculateBlockErrorRate(const OverlapBlock& terminalBlock, const OverlapBlock& otherBlock) const
{
    int back_max = std::min(terminalBlock.overlapLen, otherBlock.overlapLen);
    int backwards_diff = SearchHistoryVector::countDifferences(terminalBlock.backHistory, otherBlock.backHistory, back_max);

    // We compare the forward (right) extension only up to the last position of the terminated block's extension
    int forward_len = terminalBlock.forwardHistory.size();
    int forward_max = forward_len - 1;
    int forward_diff = SearchHistoryVector::countDifferences(terminalBlock.forwardHistory, otherBlock.forwardHistory, forward_max);

    // Calculate the length of the inferred overlap
    int trans_overlap_length = back_max + forward_len;
    double er = static_cast<double>(backwards_diff + forward_diff) / trans_overlap_length;
            
#ifdef DEBUGOVERLAP_2
    std::cout << "OL: " << terminalBlock.overlapLen << "\n";
    std::cout << "TLB BH: " << terminalBlock.backHistory << "\n";
    std::cout << "TB  BH: " << otherBlock.backHistory << "\n";
    std::cout << "TLB FH: " << terminalBlock.forwardHistory << "\n";
    std::cout << "TB  FH: " << otherBlock.forwardHistory << "\n";
    std::cout << "BM: " << back_max << " FM: " << forward_max << "\n";
    std::cout << "IOL: " << trans_overlap_length << " TD: " << (backwards_diff + forward_diff) << "\n";
    std::cout << "Block of length " << otherBlock.overlapLen << " has ier: " << er << "\n";
#endif
    return er;
}

// Update the overlap block list with a righthand extension to b, removing ranges that become invalid
void OverlapAlgorithm::updateOverlapBlockRangesRight(const BWT* pBWT, const BWT* pRevBWT, 
                                                     OverlapBlockList& obList, char canonical_base) const
{
    OverlapBlockList::iterator iter = obList.begin(); 
    while(iter != obList.end())
    {
        char relative_base = iter->flags.isQueryComp() ? complement(canonical_base) : canonical_base;
        BWTAlgorithms::updateBothR(iter->ranges, relative_base, iter->getExtensionBWT(pBWT, pRevBWT));
        // remove the block from the list if its no longer valid
        if(!iter->ranges.isValid())
        {
            iter = obList.erase(iter);
        }
        else
        {
            // Add the base to the extension history
            int currExtension = iter->forwardHistory.size();
            iter->forwardHistory.add(currExtension, canonical_base);
            ++iter;
        }
    }
}

