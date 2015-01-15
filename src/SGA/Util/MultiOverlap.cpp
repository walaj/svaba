//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// MultiOverlap.h - Data structure containing a set
// of overlaps for a given read
//
#include <algorithm>
#include <iostream>
#include "CorrectionThresholds.h"
#include "MultiOverlap.h"
#include "Alphabet.h"

MultiOverlap::MultiOverlap(const std::string& rootID, 
                           const std::string& rootSeq,
                           const std::string rootQual) : m_rootID(rootID), 
                                                         m_rootSeq(rootSeq),
                                                         m_rootQual(rootQual)
{
    
}

//
void MultiOverlap::add(const std::string& seq, const Overlap& ovr)
{
    MOData mod(seq, ovr);

    // Swap root read into first position if necessary
    if(mod.ovr.id[0] != m_rootID)
        mod.ovr.swap();
    assert(mod.ovr.id[0] == m_rootID);

    // RC the sequence if it is different orientation than the root
    if(mod.ovr.match.isRC())
    {
        mod.seq = reverseComplement(mod.seq);
        mod.ovr.match.canonize();
    }

    // Initialize the offset value, the amount that a coordinate 
    // for the non-root sequence must be shifted so that
    // the sequences are aligned
    mod.offset = mod.ovr.match.inverseTranslate(0);
    m_overlaps.push_back(mod);
}

//
void MultiOverlap::add(const MOData& mod)
{
    assert(mod.ovr.id[0] == m_rootID);
    m_overlaps.push_back(mod);
}

//
void MultiOverlap::updateRootSeq(const std::string& newSeq)
{
    m_rootSeq = newSeq;
}

//
Overlap MultiOverlap::getOverlap(size_t idx) const
{
    assert(idx < m_overlaps.size());
    return m_overlaps[idx].ovr;
}

//
std::string MultiOverlap::simpleConsensus() const
{
    std::string out;
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup p = getPileup(i);
        AlphaCount64 ac = p.getAlphaCount();
        char maxBase;
        BaseCount maxCount;
        ac.getMax(maxBase, maxCount);
        BaseCount rootCount = ac.get(m_rootSeq[i]);

        if(rootCount == maxCount)
            out.push_back(m_rootSeq[i]);
        else
            out.push_back(maxBase);
    }
    return out;
}

//
int MultiOverlap::countPotentialIncorrect(size_t cutoff) const
{
    int count = 0;
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup p = getPileup(i);
        AlphaCount64 ac = p.getAlphaCount();
        char maxBase;
        BaseCount maxCount;
        ac.getMax(maxBase, maxCount);
        BaseCount rootCount = ac.get(m_rootSeq[i]);
        if(maxBase != m_rootSeq[i] && rootCount < cutoff)
            ++count;
    }
    return count;
}

//
int MultiOverlap::countBasesCovered() const
{
    int count = 0;
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup p = getPileup(i);
        if(p.getDepth() > 1)
            ++count;
    }
    return count;
}

// Determine if this multi-overlap has a conflicted position,
// which is a position in the multioverlap where the second-most 
// prevalent base has a frequency greater than cutoff
bool MultiOverlap::isConflicted(size_t cutoff) const
{
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup p = getPileup(i);
        AlphaCount64 ac = p.getAlphaCount();

        char order[5];
        ac.getSorted(order, 5);
    
        // If the second-most prevalent base is greater than the cutoff
        // we consider the MO to be conflicted
        if(ac.get(order[1]) > cutoff)
            return true;
    }
    return false;
}

// Returns true if the overlap at idx has the include flag set
int MultiOverlap::getPartition(size_t idx) const
{
    assert(idx < m_overlaps.size());
    return m_overlaps[idx].partitionID;
}

// Returns true if the overlap at idx has the include flag set
void MultiOverlap::setPartition(size_t idx, int p)
{
    assert(idx < m_overlaps.size());
    m_overlaps[idx].partitionID = p;
}

// Return the total number of bases in the multioverlap
size_t MultiOverlap::getNumBases() const
{
    size_t count = 0;
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup p = getPileup(i);
        count += p.getDepth();
    }
    return count;
}

// Conflict-aware consensus algorithm
// For each position of the read, it calculates whether that position is conflicted (there is 
// sufficient support for an alternative base at that position). All reads that differ from the
// root read at a conflicted position are filtered out. The remaining set of reads match the root read
// at all the conflicted positions and the consensus is called from this set. Bases that are not conflicted
// are called from the entire set of overlaps.
std::string MultiOverlap::consensusConflict(double /*p_error*/, int conflictCutoff)
{
    // Calculate the frequency vector for each base of the read
    std::vector<AlphaCount64> acVec;
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        AlphaCount64 ac = getAlphaCount(i);
        acVec.push_back(ac);
    }

    // Sort the alphacounts by frequency
    StringVector sortedVec;
    sortedVec.reserve(acVec.size());
    for(size_t i = 0; i < acVec.size(); ++i)
    {
        char sorted[ALPHABET_SIZE];
        acVec[i].getSorted(sorted, ALPHABET_SIZE);
        sortedVec.push_back(std::string(sorted, ALPHABET_SIZE));
    }

    // Filter out overlaps that do not match the reference
    // at conflicted positions
    for(size_t j = 0; j < m_overlaps.size(); ++j)
    {
        int numMatch = 0;
        int numMismatch = 0;
        int numConflicted = 0;
        for(size_t i = 0; i < acVec.size(); ++i)
        {
            char rootBase = m_rootSeq[i];
            std::string& sorted = sortedVec[i];
            int second = acVec[i].get(sorted[1]);

            // If the second-most prevelent base is above the conflict cutoff,
            // call this position conflicted
            bool isConflict = second > conflictCutoff;

            // Get the basecall for read j at position i
            char b = getMODBase(m_overlaps[j], i);

            if(isConflict)
            {
                // Check if the read matches the root if:
                // a) the read contains a basecall at this position (it overlaps
                //  the read at this position)
                // b) the root read base is one of the two most frequent bases 
                //  (to filter out sequencing errors at this position in the root).
                int rootCount = acVec[i].get(rootBase);
                if(b != '\0' && rootCount > conflictCutoff)
                {
                    if(b == rootBase)
                        ++numMatch;
                    else
                        ++numMismatch;
                    ++numConflicted;
                }
            }
        }

        // Set the overlap score to be the fraction of conflict bases that this read
        // matches the root read at. 
        double frac;
        if(numConflicted == 0)
            frac = 1.0f;
        else
            frac = (double)numMatch/(double)numConflicted;
        m_overlaps[j].score = frac;

        // Filter out the read if there are any conflicted bases
        // that this read does not match the root read at
        if(numConflicted > 0 && numMismatch > 0)
            m_overlaps[j].partitionID = 1;
        else
            m_overlaps[j].partitionID = 0;
    }

    // Calculate the consensus sequence using all the reads
    // in partition 1

    std::string consensus;
    std::vector<int> supportVector;

    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup p0;
        Pileup p1;
        getPartitionedPileup(i, p0, p1);
        AlphaCount64 ac = p0.getAlphaCount();
        
        size_t minSupport = CorrectionThresholds::Instance().getMinSupportLowQuality();
        if(!m_rootQual.empty())
        {
            int phredScore = Quality::char2phred(m_rootQual[i]);
            minSupport = CorrectionThresholds::Instance().getRequiredSupport(phredScore);
        }

        size_t callSupport = ac.get(m_rootSeq[i]);
        if(callSupport >= minSupport)
        {
            // This base does not require correction
            consensus.push_back(m_rootSeq[i]);
        }
        else
        {
            // Attempt to correct the base with the most frequent base
            // in the partitioned pileup if it has been seen more often
            // than the root base
            char sorted[ALPHABET_SIZE];
            ac.getSorted(sorted, ALPHABET_SIZE);
            size_t bestSupport = ac.get(sorted[0]);
            bool corrected = false;
            if(bestSupport > callSupport)
            {
                consensus.push_back(sorted[0]);
                callSupport = bestSupport;
                corrected = true;
            }
            else
            {
                // A correction could not be made with the 
                // partitioned pileup, try to use the full pileup
                // if this base is not conflicted
                acVec[i].getSorted(sorted, ALPHABET_SIZE);
                int second = acVec[i].get(sorted[1]);
                bool isConflict = second > conflictCutoff;
                if(!isConflict)
                {
                    bestSupport = acVec[i].get(sorted[0]);
                    if(bestSupport > callSupport)
                    {
                        consensus.push_back(sorted[0]);
                        corrected = true;
                    }
                }
            }

            if(!corrected)
            {
                consensus.push_back(m_rootSeq[i]);
            }
        }
    }

    return consensus;
}

//
bool MultiOverlap::qcCheck() const
{
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        AlphaCount64 ac = getAlphaCount(i);
        size_t callSupport = ac.get(m_rootSeq[i]);
        if(callSupport < 2)
            return false;
    }
    return true;
}

//
void MultiOverlap::countOverlaps(size_t& prefix_count, size_t& suffix_count) const
{
    prefix_count = 0;
    suffix_count = 0;
    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        if(!m_overlaps[i].ovr.match.coord[0].isContained())
        {
            if(m_overlaps[i].ovr.match.coord[0].isLeftExtreme())
                ++prefix_count;
            if(m_overlaps[i].ovr.match.coord[0].isRightExtreme())
                ++suffix_count;
        }
    }
}

//
double MultiOverlap::getMeanDepth() const
{
    double depth = 0.0f;
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup p = getPileup(i);
        depth += p.getDepth();
    }
    return depth / m_rootSeq.size();
}

//
int MultiOverlap::calculateCoverageOverlap()
{
    int rightmost_prefix = -1;
    int leftmost_suffix = m_rootSeq.length();

    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        if(!m_overlaps[i].ovr.match.coord[0].isContained())
        {
            if(m_overlaps[i].ovr.match.coord[0].isLeftExtreme())
            {
                // Prefix overlap
                int end = m_overlaps[i].ovr.match.coord[0].interval.end;
                if(end > rightmost_prefix)
                    rightmost_prefix = end;
            }
            else
            {
                // Suffix overlap
                int start = m_overlaps[i].ovr.match.coord[0].interval.start;
                if(start < leftmost_suffix)
                    leftmost_suffix = start;
            }
        }
    }

    if(leftmost_suffix > rightmost_prefix)
        return 0;
    else
        return rightmost_prefix - leftmost_suffix + 1;
}

std::string MultiOverlap::calculateConsensusFromPartition(double p_error)
{
    std::string out;
    out.reserve(m_rootSeq.size());

    // require the best base call to be above this above to correct it
    double epsilon = 0.01;
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup p0;
        Pileup p1;
        getPartitionedPileup(i, p0, p1);
        DNADouble ap = p0.calculateLikelihoodNoQuality(p_error);
        char best_c;
        char curr_c = m_rootSeq[i];
        double max;
        ap.getMax(best_c, max);

        //printf("%zu best: %c curr: %c max: %lf curr_l: %lf\n", i, best_c, curr_c, max, ap.get(curr_c));
        // Require the called value to be substantially better than the
        // current base
        if(best_c == curr_c || (best_c != curr_c && max - ap.get(curr_c) < epsilon))
        {
            //std::cout << "using current\n";
            out.push_back(curr_c);
        }
        else
        {
            //std::cout << "corrected\n";
            out.push_back(best_c);
        }
    }
    return out;
}

// Return the base in mod that matches the base at
// idx in the root seq. If mod does not overlap 
// the root at this position, returns '\0'
char MultiOverlap::getMODBase(const MOData& mod, int idx) const
{
    int trans_idx = idx - mod.offset;
    if(trans_idx >= 0 && size_t(trans_idx) < mod.seq.size())
        return mod.seq[trans_idx];
    else
        return '\0';
}

// Get an AlphaCount64 representing the nucleotides
// observed at the given column
AlphaCount64 MultiOverlap::getAlphaCount(int idx) const
{
    AlphaCount64 ac;

    ac.increment(m_rootSeq[idx]);
    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        const MOData& curr = m_overlaps[i];
        // translate idx into the frame of the current sequence
        int trans_idx = idx - curr.offset;
        if(trans_idx >= 0 && size_t(trans_idx) < curr.seq.size())
        {
            ac.increment(curr.seq[trans_idx]);
        }
    }
    return ac;
}


// Get the "stack" of bases that aligns to
// a single position of the root seq, including
// the root base
Pileup MultiOverlap::getPileup(int idx) const
{
    size_t max_depth = m_overlaps.size() + 1;
    Pileup p(max_depth);
    p.add(m_rootSeq[idx]);

    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        const MOData& curr = m_overlaps[i];
        // translate idx into the frame of the current sequence
        int trans_idx = idx - curr.offset;
        if(trans_idx >= 0 && size_t(trans_idx) < curr.seq.size())
        {
            p.add(curr.seq[trans_idx]);
        }
    }
    return p;
}

// Get the "stack" of bases that aligns to
// a single position of the root seq, including
// the root base only including the first numElems items
Pileup MultiOverlap::getPileup(int idx, int numElems) const
{
    assert((size_t)numElems < m_overlaps.size());
    Pileup p;
    p.add(m_rootSeq[idx]);

    for(int i = 0; i < numElems; ++i)
    {
        const MOData& curr = m_overlaps[i];
        // translate idx into the frame of the current sequence
        int trans_idx = idx - curr.offset;
        if(trans_idx >= 0 && (size_t)trans_idx < curr.seq.size())
        {
            p.add(curr.seq[trans_idx]);
        }
    }
    return p;
}


// Get the pileup between the root sequence and one of the sequences in the MO
Pileup MultiOverlap::getSingletonPileup(int base_idx, int ovr_idx) const
{
    assert(ovr_idx < (int)m_overlaps.size());

    Pileup p;
    p.add(m_rootSeq[base_idx]);

    const MOData& curr = m_overlaps[ovr_idx];
    // translate idx into the frame of the current sequence
    int trans_idx = base_idx - curr.offset;
    if(trans_idx >= 0 && size_t(trans_idx) < curr.seq.size())
    {
        p.add(curr.seq[trans_idx]);
    }
    return p;
}

// Fill in the pileup groups g0 and g1
void MultiOverlap::getPartitionedPileup(int idx, Pileup& g0, Pileup& g1) const
{
    g0.add(m_rootSeq[idx]);

    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        const MOData& curr = m_overlaps[i];
        // translate idx into the frame of the current sequence
        int trans_idx = idx - curr.offset;
        if(trans_idx >= 0 && size_t(trans_idx) < curr.seq.size())
        {
            if(curr.partitionID == 0)
                g0.add(curr.seq[trans_idx]);
            else
                g1.add(curr.seq[trans_idx]);
        }
    }
}


// Fill in the pileup groups g0 and g1
PileupVector MultiOverlap::getPartitionedPileup(int idx, int num_parts) const
{
    assert(false && "untested");
    PileupVector pv(num_parts);
    pv.push_back(Pileup());
    
    // the root always goes in partition 0
    pv.back().add(m_rootSeq[idx]);

    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        const MOData& curr = m_overlaps[i];
        // translate idx into the frame of the current sequence
        int trans_idx = idx - curr.offset;
        if(trans_idx >= 0 && size_t(trans_idx) < curr.seq.size())
        {
            assert(curr.partitionID < num_parts);
            pv[curr.partitionID].add(curr.seq[trans_idx]);
        }
    }
    return pv;
}

// Print the MultiOverlap groups specified by the IntVec to stdout
void MultiOverlap::printGroups()
{
    for(int i = 0; i <= 1; ++i)
    {
        MultiOverlap mo_group(m_rootID, m_rootSeq);
        for(size_t j = 0; j < m_overlaps.size(); ++j)
        {
            const MOData& curr = m_overlaps[j];
            if(curr.partitionID == i)
                mo_group.add(curr);
        }
        std::cout << "MO GROUP " << i << "\n";
        mo_group.print();
    }
}

// Print the MultiOverlap to stdout
void MultiOverlap::print(int default_padding, int max_overhang)
{
    std::sort(m_overlaps.begin(), m_overlaps.end(), MOData::sortOffset);
    std::cout << "\nDrawing overlaps for read " << m_rootID << "\n";
    int root_len = int(m_rootSeq.size());
    
    // Print the root row at the bottom
    printRow(default_padding, max_overhang, root_len, 0, root_len, 0, 0.0f, m_rootSeq, m_rootID);

    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        const MOData& curr = m_overlaps[i];
        int overlap_len = curr.ovr.match.getMaxOverlapLength();
        int nd = curr.ovr.match.countDifferences(m_rootSeq, curr.seq);
        double score = static_cast<double>(nd) / overlap_len;
        printRow(default_padding, max_overhang, root_len, 
                 curr.offset, overlap_len, nd, 
                 score, curr.seq, curr.ovr.id[1].c_str());
    }    
}

// Print a single row of a multi-overlap to stdout
void MultiOverlap::printRow(int default_padding, int max_overhang, 
                            int root_len, int offset, int overlap_len, int pid, 
                            double score, const std::string& seq, const std::string& id)
{
    int c_len = seq.length();

    // This string runs from c_offset to c_offset + len
    // Clip the string at -max_overhang to root_len + max_overhang
    int left_clip = std::max(offset, -max_overhang);
    int right_clip = std::min(offset + c_len, root_len + max_overhang);
    
    // translate the clipping coordinates to the string coords
    int t_left_clip = left_clip - offset;
    int t_right_clip = right_clip - offset;
    
    // Calculate the length of the left padding
    int padding = default_padding + left_clip;
    std::string leader = (t_left_clip > 0) ? "..." : "";
    std::string trailer = (t_right_clip < c_len) ? "..." : ""; 
    std::string clipped = seq.substr(t_left_clip, t_right_clip - t_left_clip);
    padding -= leader.size();

    assert(padding >= 0);
    std::string padding_str(padding, ' ');
    std::string outstr = padding_str + leader + clipped + trailer;
    printf("%s\t%d\t%d\t%lf\tID:%s\n", outstr.c_str(), overlap_len, pid, score, id.c_str());
}

void MultiOverlap::printMasked()
{
    printf("ROOT\t%s\t%s\n", m_rootSeq.c_str(), m_rootID.c_str());

    for(size_t j = 0; j < m_overlaps.size(); ++j)
    {
        std::string out;
        for(size_t i = 0; i < m_rootSeq.length(); ++i)
        {
            char b = getMODBase(m_overlaps[j], i);
            if(b == '\0')
                out.push_back('.');
            else if(b == m_rootSeq[i])
                out.push_back('=');
            else
                out.push_back(b);
        }
        printf("OVRL\t%s\t%s\n", out.c_str(), m_overlaps[j].ovr.id[1].c_str());
    }
}

// Print the MultiOverlap horizontally, in a pileup format
void MultiOverlap::printPileup()
{
    std::cout << "\nDrawing overlap pileup for read " << m_rootID << "\n";
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        Pileup p = getPileup(i);
        printf("%zu\t%s\n", i, p.toStr().c_str());
    }
}

//
bool MultiOverlap::MOData::sortOffset(const MOData& a, const MOData& b)
{
    return a.offset < b.offset;
}

// Sort by the ID of the non-root sequence, which is 
// guarenteed to bethe second id in the overlap structure
bool MultiOverlap::MOData::sortID(const MOData& a, const MOData& b)
{
    assert(a.ovr.id[0] == b.ovr.id[0]);
    return a.ovr.id[1] < b.ovr.id[1];
}

