///-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// BWTIntervalCache - Array of cached bwt intervals for all
// substrings of a fixed length
//
#ifndef BWTINTERVAL_CACHE_H
#define BWTINTERVAL_CACHE_H

#include "BWT.h"
#include "BWTInterval.h"

class BWTIntervalCache
{
    public:

        //
        BWTIntervalCache(size_t k, const BWT* pBWT);
        
        // Look up the bwt interval for the given string
        inline BWTInterval lookup(const char* w) const
        {
            // Convert the string to an integer index in the lookup table
            size_t idx = str2int(w);
            return m_table[idx];
        }

        // 
        size_t getCachedLength() const;

    private:

        // Build the array for the given BWt
        void build(const BWT* pBWT);
        
        // Map a string to an integer
        // Precondition: w must be at least m_kmer symbols long
        inline size_t str2int(const char* w) const
        {
            size_t out = 0;
            for(size_t k = 0; k < m_kmer; ++k) {
                assert(w[k] != '$');
                out |= DNA_ALPHABET::getBaseRank(w[k]) << 2*(m_kmer - k - 1);
            }
            return out;
        }

        // Map an integer to a string
        std::string int2string(size_t i) const;

        size_t m_kmer;
        std::vector<BWTInterval> m_table;
};

#endif
