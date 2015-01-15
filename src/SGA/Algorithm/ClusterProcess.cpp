///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ClusterProcess - Merge unambiguously overlapping sequences
//
#include "ClusterProcess.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"

//
ClusterProcess::ClusterProcess(ClusterParameters params) : m_parameters(params)
{
    // Check parameters are set
    assert(m_parameters.pOverlapper != NULL);
    assert(m_parameters.pMarkedReads != NULL);
    assert(m_parameters.minOverlap >= 12);
}

//
ClusterProcess::~ClusterProcess()
{


}

//
ClusterResult ClusterProcess::process(const SequenceWorkItem& item)
{
    ReadCluster cluster(m_parameters.pOverlapper, m_parameters.minOverlap);
    ClusterNode node = cluster.addSeed(item.read.seq.toString(), true);
    assert(node.interval.isValid());
    // Check if this read is already part of a cluster. If so, return an empty result
    for(int64_t i = node.interval.lower; i <= node.interval.upper; ++i)
    {
        if(m_parameters.pMarkedReads->test(i))
        {
            ClusterResult result;
            return result; // already part of a cluster, return nothing
        }
    }
    
    // Add sequences to be used to stop extension, if requested
    if(m_parameters.pLimitKmers != NULL)
        cluster.setLimitKmers(m_parameters.pLimitKmers, m_parameters.limitK);

    // Run the clustering process
    cluster.run(m_parameters.maxClusterSize, 0);

    // If some work was performed, update the bitvector so other threads do not try to merge the same set of reads.
    // As a given set of reads should all be merged together, we only need to make sure we atomically update
    // the bit for the read with the lowest index in the set.
    ClusterResult result;
    result.clusterNodes = cluster.getOutput();

    // Check if the bit in the vector has already been set for the lowest read index
    // If it has some other thread has already output this set so we do nothing
    int64_t lowestIndex = result.clusterNodes.front().interval.lower;
    bool currentValue = m_parameters.pMarkedReads->test(lowestIndex);
    bool updateSuccess = false;

    if(currentValue == false)
    {
        // Attempt to update the bit vector with an atomic CAS. If this returns false
        // the bit was set by some other thread
        updateSuccess = m_parameters.pMarkedReads->updateCAS(lowestIndex, currentValue, true);
    }

    if(updateSuccess)
    {
        // We successfully atomically set the bit for the first read in this set
        // to true. We can safely update the rest of the bits and keep the merged sequences
        // for output.
        std::vector<ClusterNode>::const_iterator iter = result.clusterNodes.begin();
        for(; iter != result.clusterNodes.end(); ++iter)
        {
            for(int64_t i = iter->interval.lower; i <= iter->interval.upper; ++i)
            {
                if(i == lowestIndex) //already set
                    continue;
                currentValue = m_parameters.pMarkedReads->test(i);
                if(currentValue)
                {
                    // This value should not be true, emit a warning
                    std::cout << "Warning: Bit " << i << " was unexpectedly set by a different thread\n";
                }
                else
                {
                    m_parameters.pMarkedReads->updateCAS(i, currentValue, true);
                }
            }
        }
    }
    else
    {
        // Some other thread merged these reads already, discard the intermediate
        // data and set the result to false
        result.clusterNodes.clear();
    }
    return result;
}

// Generate a new cluster from a previously build cluster
ClusterResult ClusterProcess::process(const ClusterVector& inSequences)
{
    ReadCluster cluster(m_parameters.pOverlapper, m_parameters.minOverlap);
    
    // Add seeds to the cluster generator
    for(size_t i = 0; i < inSequences.size(); ++i)
        cluster.addSeed(inSequences[i].sequence, false);

    // Add sequences to be used to stop extension, if requested
    if(m_parameters.pLimitKmers != NULL)
        cluster.setLimitKmers(m_parameters.pLimitKmers, m_parameters.limitK);

    cluster.run(m_parameters.maxClusterSize, m_parameters.maxIterations);

    ClusterResult result;
    result.clusterNodes = cluster.getOutput();
    return result;
}

//
ClusterPostProcess::ClusterPostProcess(std::ostream* pWriter, 
                                       size_t minClusterSize, 
                                       BitVector* pMarkedReads) : m_minClusterSize(minClusterSize),
                                                                  m_numClusters(0), 
                                                                  m_numTotalReads(0), 
                                                                  m_numTotalReadsClustered(0), 
                                                                  m_pWriter(pWriter), 
                                                                  m_pMarkedReads(pMarkedReads)
{

}

//
ClusterPostProcess::~ClusterPostProcess()
{
    printf("[sga cluster] Clustered %zu reads into %zu clusters (%zu total reads input)\n", m_numTotalReadsClustered, m_numClusters, m_numTotalReads);
}

// This just dispatches to main post-process function
void ClusterPostProcess::process(const SequenceWorkItem& /*item*/, const ClusterResult& result)
{
    process(result);
}

//
void ClusterPostProcess::process(const ClusterVector& /*in*/, const ClusterResult& result)
{
    process(result);
}

void ClusterPostProcess::process(const ClusterResult& result)
{
    m_numTotalReads += 1;

    if(result.clusterNodes.size() > 0)
    {
        // Compute the cluster size
        size_t clusterSize = 0;
        for(size_t i = 0; i < result.clusterNodes.size(); ++i)
        {
            if(result.clusterNodes[i].interval.isValid())
                clusterSize += result.clusterNodes[i].interval.size();
            else
                clusterSize += 1;
        }
               
        if(clusterSize >= m_minClusterSize)
        {
            for(size_t i = 0; i < result.clusterNodes.size(); ++i)
                *m_pWriter << "cluster-" << m_numClusters << "\t" << clusterSize << "\t" << result.clusterNodes[i].sequence << "\t" << result.clusterNodes[i].interval << "\n";
            m_numClusters += 1;
            m_numTotalReadsClustered += clusterSize;
        }
    }
}
