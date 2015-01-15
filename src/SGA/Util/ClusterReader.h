//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ClusterReader - Read in a file of read clusters.
// This class conforms to the interface needed
// for the SequenceProcessFramework concurrency lib
//
#ifndef CLUSTER_READER_H
#define CLUSTER_READER_H

#include <string>
#include <iostream>
#include <vector>

struct ClusterRecord
{
    std::string clusterID;
    int numElements;
    std::string readID;
    std::string sequence;
};
typedef std::vector<ClusterRecord> ClusterVector;

class ClusterReader
{
    public:
        ClusterReader(const std::string& filename);
        ~ClusterReader();

        // Read in a cluster from the file and write the records
        // to out. Returns true if a cluster was successfully read.
        bool generate(ClusterVector& out);
        bool readCluster(ClusterRecord& record);

        // Functions to determine the progress of the read
        inline size_t getConsumedLast() const { return m_numConsumedLast; }
        inline size_t getNumConsumed() const { return m_numConsumedTotal; }
        
    private:

        std::istream* m_pReader;
        size_t m_numConsumedLast;
        size_t m_numConsumedTotal;
};

#endif
