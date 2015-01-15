//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapBlock - Data structures holding
// the result of the alignment of a sequence read
// to a BWT
// 
#include "OverlapBlock.h"
#include "BWTAlgorithms.h"

//#define DEBUG_RESOLVE 1

// 
OverlapBlock::OverlapBlock(BWTIntervalPair r,
                           BWTIntervalPair rawI,
                           int ol, 
                           int nd, 
                           const AlignFlags& af, 
                           const SearchHistoryVector& backHist) : ranges(r), 
                                                                  rawRanges(rawI),
                                                                  overlapLen(ol), 
                                                                  numDiff(nd),
                                                                  flags(af),
                                                                  isEliminated(false),
                                                                  backHistory(backHist)
{
    backHistory.normalize(af.isQueryComp());
}

// Return a pointer to the BWT that should be used to extend the block
// this is the opposite BWT that was used in the backwards search
const BWT* OverlapBlock::getExtensionBWT(const BWT* pBWT, const BWT* pRevBWT) const
{
    if(!flags.isTargetRev())
        return pRevBWT;
    else
        return pBWT;
}

// 
AlphaCount64 OverlapBlock::getCanonicalExtCount(const BWT* pBWT, const BWT* pRevBWT) const
{
    AlphaCount64 out = BWTAlgorithms::getExtCount(ranges.interval[1], getExtensionBWT(pBWT, pRevBWT));
    if(flags.isQueryComp())
        out.complement();
    return out;
}

// Returns 0 if the BWT used for the overlap step was the forward BWT
int OverlapBlock::getCanonicalIntervalIndex() const
{
    if(!flags.isTargetRev())
        return 0;
    else
        return 1;
}

//
BWTInterval OverlapBlock::getCanonicalInterval() const
{
    return ranges.interval[getCanonicalIntervalIndex()];
}

// Get the string corresponding to the overlap block. This is the string found
// during the backwards search
std::string OverlapBlock::getOverlapString(const std::string& original) const
{
    std::string transformed = backHistory.transform(original, flags.isQueryRev());
    // If the query was reversed, we take the first overlapLen (the search
    // was from the front of the sequence) otherwise we take the last overlapLen
    if(flags.isQueryRev())
        return transformed.substr(0, overlapLen);
    else
        return transformed.substr(transformed.length() - overlapLen);
}

// Get the full string corresponding to this overlapblock using the forward history
std::string OverlapBlock::getFullString(const std::string& original) const
{
    std::string str = getOverlapString(original);
    std::string history = forwardHistory.getBaseString();

    if(history.empty() && overlapLen != (int)original.size())
    {
        WARN_ONCE("getFullString() called on block with no history")
    }
/*
    std::cout << "OVERLAP: " << str << "\n";
    std::cout << "HIST: " << history << "\n";
    std::cout << "QREV: " << flags.isQueryRev() << "\n";
    std::cout << "RC: " << flags.isReverseComplement() << "\n";
    std::cout << "QC: " << flags.isQueryComp() << "\n";
*/
    if(!flags.isQueryRev())
    {
        str.append(history);
    }
    else
    {
        history = reverse(history);
        history.append(str);
        str.swap(history);
    }

    if(flags.isReverseComplement())
        str = reverseComplement(str);
    return str;
}


EdgeDir OverlapBlock::getEdgeDir() const
{
    if(flags.isQueryRev())
        return ED_ANTISENSE;
    else
        return ED_SENSE;
}

//
Overlap OverlapBlock::toOverlap(const std::string queryID, const std::string targetID, int queryLen, int targetLen) const
{
    // Compute the sequence coordinates
    int s1 = queryLen - overlapLen;
    int e1 = s1 + overlapLen - 1;
    SeqCoord sc1(s1, e1, queryLen);

    int s2 = 0; // The start of the second hit must be zero by definition of a prefix/suffix match
    int e2 = s2 + overlapLen - 1;
    SeqCoord sc2(s2, e2, targetLen);

    // The coordinates are always with respect to the read, so flip them if
    // we aligned to/from the reverse of the read
    if(flags.isQueryRev())
    {
        sc1.flip();
    }
    if(flags.isTargetRev())
    {
        sc2.flip();
    }
    bool isRC = flags.isReverseComplement();

    Overlap o(queryID, sc1, targetID, sc2, isRC, numDiff);
    return o;
}

//
std::string OverlapBlock::toCanonicalID() const
{
    std::stringstream ss;
    int ci = getCanonicalIntervalIndex();
    ss << "IDX-" << ranges.interval[ci].lower;
    return ss.str();
}


//
void printBlockList(const OverlapBlockList* pList)
{
    for(OverlapBlockList::const_iterator i = pList->begin(); i != pList->end(); ++i)
    {
        std::cout << "Block: " << *i << "\n";
    }
}

//
void removeSubMaximalBlocks(OverlapBlockList* pList, const BWT* pBWT, const BWT* pRevBWT)
{
    // This algorithm removes any sub-maximal OverlapBlocks from pList
    // The list is sorted by the left coordinate and iterated through
    // if two adjacent blocks overlap they are split into maximal contiguous regions
    // with resolveOverlap. The resulting list is merged back into pList. This process
    // is repeated until each block in pList is a unique range
    // The bookkeeping in the intersecting case could be more efficient 
    // but the vast vast majority of the cases will not have overlapping 
    // blocks.
    pList->sort(OverlapBlock::sortIntervalLeft);
    OverlapBlockList::iterator iter = pList->begin();
    OverlapBlockList::iterator last = pList->end();
    --last;

    while(iter != pList->end())
    {
        OverlapBlockList::iterator next = iter;
        ++next;

        if(next == pList->end())
            break;

        // Check if iter and next overlap
        if(Interval::isIntersecting(iter->ranges.interval[0].lower, iter->ranges.interval[0].upper, 
                                    next->ranges.interval[0].lower, next->ranges.interval[0].upper))
        {
            OverlapBlockList resolvedList = resolveOverlap(*iter, *next, pBWT, pRevBWT);
            
            // Merge the new elements in and start back from the beginning of the list
            pList->erase(iter);
            pList->erase(next);
            pList->merge(resolvedList, OverlapBlock::sortIntervalLeft);
            iter = pList->begin();

            //std::cout << "After splice: \n";
            //printList(pList);
        }
        else
        {
            ++iter;
        }
    }
}

// A tracing interval maps an interval of an overlap block representing
// a smaller overlap to an overlap block with a larger overlap. This is
// used to determine which intervals are redundant and can be removed.
struct TracingInterval
{
    int64_t foundPosForward;
    int64_t sourcePosReverse;
    BWTInterval tracing;
    BWTIntervalPair updateRanges;
};

typedef std::list<TracingInterval> TracingIntervalList;

// In rare cases, the overlap blocks may represent sub-maximal overlaps between reads
// we need to distinguish these cases and remove the sub-optimal hits. This
// function splits two overlapping OverlapBlocks into up to three distinct
// blocks, keeping the maximal (longest) overlap at each stage.
OverlapBlockList resolveOverlap(const OverlapBlock& A, const OverlapBlock& B, const BWT* pBWT, const BWT* pRevBWT)
{
    OverlapBlockList outList;

    // Check if A and B have the same overlap length, if so they must be 
    // identical blocks (resulting from different seeds) and we can remove one
    if(A.overlapLen == B.overlapLen)
    {
        if(A.ranges.interval[0].lower == B.ranges.interval[0].lower &&
           A.ranges.interval[0].upper == B.ranges.interval[0].upper)
        {
            outList.push_back(A);
            return outList;
        }
        else
        {
            std::cerr << "Error in resolveOverlap: Overlap blocks with same length do not "
            "the have same coordinates\n";
            assert(false);
        }    
    }
    // A and B must have different overlap lengths
    assert(A.overlapLen != B.overlapLen);

    // Determine which of A and B have a higher overlap
    const OverlapBlock* pHigher;
    const OverlapBlock* pLower;
    if(A.overlapLen > B.overlapLen)
    {
        pHigher = &A;
        pLower = &B;
    }
    else
    {
        pHigher = &B;
        pLower = &A;
    }

    // Complicated logic follows
    // We always want the entirity of the block with the longer
    // overlap so it is added to outList unmodified
    outList.push_back(*pHigher);

    // The lower block can be split into up to two pieces:
    // Case 1:
    //     Lower  ------ 
    //     Higher    ------
    //     Result ---
    //
    // Case 2:
    //     Lower  -----------
    //     Higher    ------
    //     Result ---      --
    //
    // Case 3:
    //     Lower  ------
    //     Higher ------
    //     Result (empty set)

    // It is unclear whether case 2 can happen in reality but we handle it 
    // here anyway. Further complicating matters is that the BWTIntervalPair
    // keeps track of both the BWT coordinates for the backwards search
    // and forward search and we must take care to ensure that both intervals
    // are updated and the mapping between them is correct. We do this
    // by calculating the new forward interval using interval intersections
    // and directly recalculating the coordinate of the reverse interval
    //
    OverlapBlock split = *pLower;

    // Left-hand split
    if( (pLower->ranges.interval[0].lower < pHigher->ranges.interval[0].lower) ||
        (pLower->ranges.interval[0].upper > pHigher->ranges.interval[0].upper) )
    {
        // The intervals do not perfectly overlap and must be recalculated. 
        // We start from the raw intervals in the lower block (the intervals representing
        // overlaps that are not capped by '$' symbols) and search backwards through the
        // bwt until the start of the sequence has been found. This maps the source reverse
        // index position to the forward index position. We can then decide which intervals
        // are redundant and can be removed.
        //
        // If the index has duplicates, it is possible that a given source reverse position
        // will map to multiple forward positions. To handle this case, we record the used
        // forward positions in a std::map so we can lookup the next lowest index that is available.
        //
        // A better algorithm (that doesn't required so many interval calculations) probably exists
        // but this case is very rare so simplicity wins here.
        //

#ifdef DEBUG_RESOLVE
        std::cout << "LOWER -- capped: " << pLower->ranges << "\n";
        std::cout << "LOWER -- raw: " << pLower->rawRanges << "\n";
        std::cout << "HIGHER -- capped: " << pHigher->ranges << "\n";
        std::cout << "HIGHER -- raw: " << pHigher->rawRanges << "\n";
#endif

        std::map<int64_t, int64_t> usedMap;

        // Remap every reverse position to a forward position
        TracingIntervalList tracingList;
        int64_t j = pLower->ranges.interval[1].lower;
        while(j <= pLower->ranges.interval[1].upper)
        {
            TracingInterval ti;
            ti.sourcePosReverse = j;

            ti.tracing.lower = j;
            ti.tracing.upper = j;

            ti.updateRanges = pLower->rawRanges;

            bool done = false;
            while(!done)
            {
                char trace_base = pRevBWT->getChar(ti.tracing.lower);
                if(trace_base == '$')
                {
                    BWTAlgorithms::updateBothL(ti.updateRanges, '$', pBWT);
                    done = true;
                }
                BWTAlgorithms::updateInterval(ti.tracing, trace_base, pRevBWT);
                BWTAlgorithms::updateBothR(ti.updateRanges, trace_base, pRevBWT);
            }
            
            if(ti.updateRanges.interval[0].lower == ti.updateRanges.interval[0].upper)
            {
                // This read is not duplicated
                ti.foundPosForward = ti.updateRanges.interval[0].lower;
            }
            else
            {
                // This read is duplicated, look up its value in the map
                int64_t basePos = ti.updateRanges.interval[0].lower;
                if(usedMap.find(basePos) != usedMap.end())
                {
                    // Use the value in the map and update it
                    ti.foundPosForward = usedMap[basePos];
                    assert(ti.foundPosForward > ti.updateRanges.interval[0].lower && 
                           ti.foundPosForward <= ti.updateRanges.interval[0].upper);
                    usedMap[basePos]++;
                }
                else
                {
                    // Use the base value and initialize the map
                    ti.foundPosForward = basePos;
                    usedMap[basePos] = basePos + 1;
                }
            }
            ++j;
            tracingList.push_back(ti);
        }

        // Reset the mapping between blocks
        std::list<BWTIntervalPair> retainedIntervals;
        TracingIntervalList::iterator tracingIter = tracingList.begin();
        while(tracingIter != tracingList.end())
        {
            // Check if the forward position intersects the higher block, if so this block
            // is redundant and can be removed.
            if(!Interval::isIntersecting(tracingIter->foundPosForward, tracingIter->foundPosForward,
                                         pHigher->ranges.interval[0].lower, pHigher->ranges.interval[0].upper))
            {
                BWTIntervalPair retained;
                
                retained.interval[0].lower = tracingIter->foundPosForward;
                retained.interval[0].upper = tracingIter->foundPosForward;
                retained.interval[1].lower = tracingIter->sourcePosReverse;
                retained.interval[1].upper = tracingIter->sourcePosReverse;

#ifdef DEBUG_RESOLVE
                std::cout << "Retained coords: " << retained << "\n";
#endif
                retainedIntervals.push_back(retained);
            }
            ++tracingIter;
        }
        
        // Write out the final blocks
        std::list<BWTIntervalPair>::iterator iter = retainedIntervals.begin();
        while(iter != retainedIntervals.end())
        {
#ifdef DEBUG_RESOLVE
            std::cout << "OUTPUT: " << *iter << "\n";
#endif
            split.ranges = *iter;

            // Sanity check
            assert(split.ranges.interval[0].size() == split.ranges.interval[1].size());
            assert(split.ranges.interval[0].isValid());
            assert(split.ranges.interval[1].isValid());

            outList.push_back(split);
            ++iter;
        }
    }

    // Sort the outlist by left coordinate
    outList.sort(OverlapBlock::sortIntervalLeft);
    return outList;
}

// Partition the overlap block list into two lists, 
// one for the containment overlaps and one for the proper overlaps
void partitionBlockList(int readLen, OverlapBlockList* pCompleteList, 
                        OverlapBlockList* pOverlapList, 
                        OverlapBlockList* pContainList)
{
    OverlapBlockList::iterator iter = pCompleteList->begin();
    while(iter != pCompleteList->end())
    {
        if(iter->overlapLen == readLen)
            pContainList->splice(pContainList->end(), *pCompleteList, iter++);
        else
            pOverlapList->splice(pOverlapList->end(), *pCompleteList, iter++);
    }
}

// Filter out full-length (containment) overlaps from the block list
void removeContainmentBlocks(int readLen, OverlapBlockList* pList)
{
    OverlapBlockList::iterator iter = pList->begin();
    while(iter != pList->end())
    {
        if(iter->overlapLen == readLen)
            iter = pList->erase(iter);
        else
            ++iter;
    }    
}

// 
MultiOverlap blockListToMultiOverlap(const SeqRecord& record, OverlapBlockList& blockList)
{
    std::string read_idx = record.id;
    std::string read_seq = record.seq.toString();
    std::string qual = record.qual;
    MultiOverlap out(read_idx, read_seq, qual);

    for(OverlapBlockList::iterator iter = blockList.begin(); iter != blockList.end(); ++iter)
    {
        std::string overlap_string = iter->getOverlapString(read_seq);

        // Compute the endpoints of the overlap
        int s1 = read_seq.length() - iter->overlapLen;
        int e1 = s1 + iter->overlapLen - 1;
        SeqCoord sc1(s1, e1, read_seq.length());

        int s2 = 0; // The start of the second hit must be zero by definition of a prefix/suffix match
        int e2 = s2 + iter->overlapLen - 1;
        SeqCoord sc2(s2, e2, overlap_string.length());

        // The coordinates are always with respect to the read, so flip them if
        // we aligned to/from the reverse of the read
        if(iter->flags.isQueryRev())
            sc1.flip();
        if(iter->flags.isTargetRev())
            sc2.flip();

        bool isRC = false; // since we transformed the original sequence, they are never RC
        if(sc1.isContained())
            continue; // skip containments

        // Add an overlap for each member of the block
        for(int64_t i = iter->ranges.interval[0].lower; i <= iter->ranges.interval[0].upper; ++i)
        {
            Overlap o(read_idx, sc1, makeIdxString(i), sc2, isRC, -1);
            out.add(overlap_string, o);
        }
    }
    return out;
}

// make an id string from a read index
std::string makeIdxString(int64_t idx)
{
    std::stringstream ss;
    ss << idx;
    return ss.str();
}

