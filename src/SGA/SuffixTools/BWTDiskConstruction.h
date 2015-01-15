//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// BWTDiskConstruction - Implementation of the
// BWT disk construction algorithm of Ferragina,
// Gagie, Manzini. See Lightweight Data Indexing
// and Compression in External Memory
//
// This implementation does not strictly follow 
// their method but is modified to handle
// large numbers of short strings
#ifndef BWTDISKCONSTRUCTION_H
#define BWTDISKCONSTRUCTION_H

#include "SuffixArray.h"
#include "BWT.h"

struct BWTDiskParameters
{
    std::string inFile;
    std::string outPrefix;
    std::string bwtExtension;
    std::string saiExtension;
    size_t numReadsPerBatch;
    int numThreads;
    int storageLevel;
    bool bBuildReverse;
    bool bUseBCR;
};

// Construct the burrows-wheeler transform of reads in in_filename
// using the disk storage algorithm
void buildBWTDisk(const BWTDiskParameters& parameters);

// Merge the indices for the readsFile1 and readsFile2
void mergeIndependentIndices(const std::string& readsFile1, const std::string& readsFile2, 
                             const std::string& outPrefix, const std::string& bwt_extension, 
                             const std::string& sai_extension, bool doReverse, int numThreads, int storageLevel);

// Compute new indices from allReadsFile without the reads in readsToRemove
void removeReadsFromIndices(const std::string& allReadsFile, const std::string& readsToRemove,
                             const std::string& outPrefix, const std::string& bwt_extension, 
                             const std::string& sai_extension, bool doReverse, int numThreads);

//
void mergeReadFiles(const std::string& readsFile1, const std::string& readsFile2, const std::string& outPrefix);
#endif
