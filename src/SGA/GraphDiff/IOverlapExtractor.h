///----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// IOverlapExtractor - Abstract interface for
// methods to extract all the reads overlapping
// a given sequence.
//
#ifndef IOVERLAP_EXTRACTOR_H
#define IOVERLAP_EXTRACTOR_H

#include "KmerOverlaps.h"

class IOverlapExtractor
{
    public:

        virtual ~IOverlapExtractor() {};

        // Return a vector of reads that exactly overlap a given query
        virtual SequenceOverlapPairVector getExactOverlaps(const std::string& query) = 0;

};

#endif
