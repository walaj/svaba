///-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ReadCluster - Generate a cluster of overlapping
// reads using the FM-index
//
#include "ReadCluster.h"

//
ReadCluster::ReadCluster(const OverlapAlgorithm* pOverlapper, int minOverlap) : m_pOverlapper(pOverlapper), 
                                                                                m_minOverlap(minOverlap),
                                                                                m_pLimitKmers(NULL),
                                                                                m_limitK(0)
{

}

//
void ReadCluster::setLimitKmers(const std::set<std::string>* pLimitKmers, size_t limitK)
{
    m_pLimitKmers = pLimitKmers;
    m_limitK = limitK;
}

// Add a seed read to the cluster. Overlaps will be found for 
// each seed read to grow the cluster
ClusterNode ReadCluster::addSeed(const std::string& sequence, bool bCheckInIndex)
{
    // Check if this read is a substring
    SeqRecord tempRecord;
    tempRecord.id = "cluster-seed";
    tempRecord.seq = sequence;

    OverlapBlockList tempBlockList;
    OverlapResult overlapResult = m_pOverlapper->alignReadDuplicate(tempRecord, &tempBlockList);
    if(overlapResult.isSubstring)
    {
        // If bCheckInIndex is true, then we are extending clusters from reads that are in the FM-index
        // In this case, seeding a substring read is an error and we abort. If the seed is NOT in the index
        // we just emit a warning and ignore the seed
        if(bCheckInIndex)
        {
            std::cerr << "Error: The seed used for sga-cluster-extend is a substring of some read.\n";
            std::cerr << "Please run sga rmdup before cluster\n";
            std::cerr << "Sequence: " << sequence << "\n";
            exit(1);
        }
        else
        {
            std::cerr << "Warning: The cluster sequence to extend is a substring of some read. This seed is being skipped\n";
            std::cerr << "Sequence: " << sequence << "\n";
            ClusterNode emptyNode;

            // Set an invalid interval
            emptyNode.interval.lower = 0;
            emptyNode.interval.upper = -1;
            return emptyNode;
        }
    }

    // Find the interval in the fm-index containing the read
    const BWT* pBWT = m_pOverlapper->getBWT();
    BWTInterval readInterval = BWTAlgorithms::findInterval(pBWT, sequence);
    BWTAlgorithms::updateInterval(readInterval, '$', pBWT);

    // When building primary clusters, we require each read to be in the index.
    if(bCheckInIndex)
    {
        if(!readInterval.isValid())
        {
            std::cerr << "sga cluster error: The seed read is not part of the FM-index.\n";
            exit(EXIT_FAILURE);
        }
    }

    ClusterNode node;
    node.sequence = sequence;
    node.interval = readInterval;
    node.isReverseInterval = false;
    m_usedIndex.insert(readInterval.lower);
    m_queue.push(node);
    return node;
}

// Run the cluster process. If the number of total nodes
// exceeds max, abort the search.
void ReadCluster::run(size_t max_size, int max_iterations)
{
    int iteration = 0;
    do
    {
        // Perform one iteration of extending each sequence in the queue
        ClusterNodeQueue next_queue;
        bool aborted = false;
            
        // Have we performed enough iterations?
        if(max_iterations > 0 && iteration++ > max_iterations)
            aborted = true;

        while(!aborted && !m_queue.empty())
        {
            if(m_queue.size() + m_outCluster.size() > max_size)
            {
                // Abort the search
                m_outCluster.clear();
                aborted = true;
                break;
            }
            
            ClusterNode node = m_queue.front();
            m_queue.pop();

            // Add this node to the output
            m_outCluster.push_back(node);
            
            // If we are using limit kmers, only try to extend this sequence if it does not contain a limit kmer
            bool extend_read = canExtendRead(node);
            if(!extend_read)
                continue; 

            // Find overlaps for the current node
            SeqRecord tempRecord;
            tempRecord.id = "cluster";
            tempRecord.seq = node.sequence;
            OverlapBlockList blockList;
            m_pOverlapper->overlapRead(tempRecord, m_minOverlap, &blockList);
            
            // Parse each member of the block list and potentially expand the cluster
            for(OverlapBlockList::const_iterator iter = blockList.begin(); iter != blockList.end(); ++iter)
            {
                // Check if the reads in this block are part of the cluster already
                BWTInterval canonicalInterval = iter->getCanonicalInterval();
                int64_t canonicalIndex = canonicalInterval.lower;
                if(m_usedIndex.count(canonicalIndex) == 0)
                {
                    // This is a new node that isn't in the cluster. Add it.
                    m_usedIndex.insert(canonicalIndex);

                    ClusterNode newNode;
                    newNode.sequence = iter->getFullString(node.sequence);
                    newNode.interval = canonicalInterval;
                    newNode.isReverseInterval = iter->flags.isTargetRev();
                    next_queue.push(newNode);
                }
            }
        }

            
        // Clear the queue
        if(aborted) 
        {
            while(!m_queue.empty())
                m_queue.pop();
        } 
        else
        {
            // Swap the reads to extend in the next iteration
            m_queue = next_queue;
        }

    } while(!m_queue.empty());
}

// Check if the given node contains a kmer that we stop extension at.
// Returns true if we can extend the cluster using the given node
bool ReadCluster::canExtendRead(const ClusterNode& node) const
{
    // No limiting kmers
    if(m_pLimitKmers == NULL)
        return true;
    // Node too short
    if(node.sequence.size() < m_limitK)
        return false;

    size_t nk = node.sequence.size() - m_limitK + 1;
    for(size_t i = 0; i < nk; ++i)
    {
        std::string kmer = node.sequence.substr(i, m_limitK);
        if(m_pLimitKmers->find(kmer) != m_pLimitKmers->end())
            return false; // this read contains a limiting kmer
        
        // check reverse complement too
        kmer = reverseComplement(kmer);
        if(m_pLimitKmers->find(kmer) != m_pLimitKmers->end())
            return false; // this read contains a limiting kmer
    }

    return true;
}

// 
ClusterNodeVector ReadCluster::getOutput() const
{
    // Sort the intervals into ascending order and remove any duplicate intervals (which can occur
    // if the subgraph has a simple cycle)
    ClusterNodeVector retVector = m_outCluster;
    std::sort(retVector.begin(), retVector.end(), ClusterNode::compare);

    std::vector<ClusterNode>::iterator newEnd = std::unique(retVector.begin(),
                                                            retVector.end(),
                                                            ClusterNode::equal);

    retVector.erase(newEnd, retVector.end());
    return retVector;
}
