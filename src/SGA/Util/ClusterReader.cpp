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
#include "ClusterReader.h"
#include "Util.h"

//
ClusterReader::ClusterReader(const std::string& filename) : m_numConsumedLast(0), m_numConsumedTotal(0)
{
    m_pReader = createReader(filename);
}

ClusterReader::~ClusterReader()
{
    delete m_pReader;
}

// Read a cluster from the file and write its members
// to out. Returns false if the read failed
bool ClusterReader::generate(ClusterVector& out)
{
    out.clear();

    // Read a single cluster from the file.
    ClusterRecord record;
    bool good = readCluster(record);

    // Read failed, return false
    if(!good)
        return false;

    out.push_back(record);

    // Read the remaining records for this cluster
    int remaining = record.numElements - 1;
    for(int i = 0; i < remaining; ++i)
    {
        bool good = readCluster(record);
        if(!good)
        {
            std::cerr << "Error: expected " << remaining + 1 << " elements in the cluster but only read " << i+1 << "\n";
            exit(1);
        }

        if(record.clusterID != out.front().clusterID)
        {
            std::cerr << "Error: cluster names do not match! " << record.clusterID << " != " << out.front().clusterID << "\n";
            exit(1);
        }
        out.push_back(record);
    }

    m_numConsumedLast = 1;
    m_numConsumedTotal += 1;

    return true;
}

//
bool ClusterReader::readCluster(ClusterRecord& record)
{
    std::string line;
    bool good = getline(*m_pReader, line).good();
    if(!good || line.empty())
        return false;
    std::stringstream parser(line);
    parser >> record.clusterID;
    parser >> record.numElements;
    parser >> record.readID;
    parser >> record.sequence;

    if(record.clusterID.empty() || record.numElements == 0 || record.readID.empty() || record.sequence.empty())
    {
        std::cerr << "Could not parse cluster record from line: " << line << "\n";
        exit(1);
    }
    return true;
}

