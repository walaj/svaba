//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// InverseSuffixArray - The inverse of a generalized suffix array
// Implemented as a Vector of Vectors
//
#include "InverseSuffixArray.h"

// Constructor which knows the number of strings that will be in the ISA
InverseSuffixArray::InverseSuffixArray(const SuffixArray& sa)
{
    // Make a temporary vector to count the number of ranks for each string
    RankCountMap rcm;

    // Fill in the max ranks
    size_t n_sa = sa.getSize();

    for(size_t i = 0; i < n_sa; ++i)
    {
        size_t id = sa.get(i).getID();
        rcm[id]++;
    }

    for(RankCountMap::iterator iter = rcm.begin(); iter != rcm.end(); ++iter)
    {
        m_data[iter->first].resize(iter->second);
    }

    // Fill in the data
    for(size_t i = 0; i < n_sa; ++i)
    {
        size_t id = sa.get(i).getID();
        size_t pos = sa.get(i).getPos();
        m_data[id][pos] = i;
    }
}

// get the rank
uint64_t InverseSuffixArray::getRank(size_t id, size_t pos) const
{
    return m_data.find(id)->second[pos];
    //return m_data[id][pos];
}

// validate
void InverseSuffixArray::validate() const
{
    RankVector allRanks;
    for(RankVectorMap::const_iterator iter = m_data.begin(); iter != m_data.end(); ++iter)
    {        
        const RankVector& rv = iter->second;
        for(size_t j = 0; j < rv.size(); ++j)
        {
            allRanks.push_back(rv[j]);
        }
    }

    std::sort(allRanks.begin(), allRanks.end());
    for(size_t i = 0; i < allRanks.size(); ++i)
    {
        assert(i == allRanks[i]);
    }
}

// Print the ISA
void InverseSuffixArray::print() const
{
    for(RankVectorMap::const_iterator iter = m_data.begin(); iter != m_data.end(); ++iter)
    {        
        const RankVector& rv = iter->second;
        for(size_t j = 0; j < rv.size(); ++j)
        {
            printf("ISA(%zu, %zu): %zu\n", (size_t)iter->first, j, (size_t)getRank(iter->first, j));
        }
    }
}

