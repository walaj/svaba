///-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// BWTIntervalCache - Array of cached bwt intervals for all
// substrings of a fixed length
//
#include "BWTIntervalCache.h"
#include "BWTAlgorithms.h"

BWTIntervalCache::BWTIntervalCache(size_t k, const BWT* pBWT) : m_kmer(k)
{
    build(pBWT);
}

// Build the table for the given bwt
void BWTIntervalCache::build(const BWT* pBWT)
{
    // Restrict the kmer parameter to something reasonable
    // so we don't try to allocate an absurdly large array
    assert(m_kmer <= 12);

    size_t num_entries = 1 << 2*m_kmer;
    m_table.resize(num_entries);

    // Construct the table
    for(size_t i = 0; i < num_entries; ++i)
    {
        std::string w = int2string(i);
        //printf("i: %zu w: %s o: %zu\n", i, w.c_str(), string2int(w));
        BWTInterval interval = BWTAlgorithms::findInterval(pBWT, w);
        m_table[i] = interval;
    }
}

// Construct the corresponding string for integer i
std::string BWTIntervalCache::int2string(size_t i) const
{
    std::string out(m_kmer, 'A');
    for(size_t k = 0; k < m_kmer; ++k)
    {
        // Get the character as position k (0 = left-most position)
        size_t code = (i >> 2*(m_kmer - k - 1)) & 3;
        char b = DNA_ALPHABET::getBase(code);
        out[k] = b;
    }
    return out;
}

// Return the length of the cached strings
size_t BWTIntervalCache::getCachedLength() const
{
    return m_kmer;
}
