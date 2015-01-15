//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// InverseSuffixArray - The inverse of a generalized suffix array
// Implemented as a Vector of Vectors
//
#ifndef INVERSESUFFIXARRAY_H
#define INVERSESUFFIXARRAH_H
#include "STCommon.h"
#include "SuffixArray.h"

typedef std::vector<uint64_t> RankVector;
typedef std::map<uint64_t, uint32_t> RankCountMap;
typedef std::map<uint64_t, RankVector> RankVectorMap;

class InverseSuffixArray
{
    public:
        
        InverseSuffixArray(const SuffixArray& sa);
        
        //
        uint64_t getRank(size_t id, size_t pos) const;
        void validate() const;
        void print() const;

    private:
        RankVectorMap m_data;
};

#endif
