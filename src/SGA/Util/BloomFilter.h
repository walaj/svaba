//-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
#ifndef BLOOM_FILTER_H
#define BLOOM_FILTER_H

#include <vector>
#include <stdint.h>
#include <stddef.h>
#include <limits>
#include "BitVector.h"

//#define TRACK_OCCUPANCY 1

class BloomFilter
{
    public:

        /**
        * @brief Default constructor.
        */
        BloomFilter();
        BloomFilter(size_t width, size_t num_hashes) { initialize(width, num_hashes); }

        /**
        * @brief Initialize the bloom filter.
        *
        * @param width       The number of bits to use
        * @param num_hashes  The number of hashes to use
        */
        void initialize(size_t width, size_t num_hashes);

        /**
        * @brief Add an object to the collection
        *
        * @param key        A pointer to the key data
        * @param num_bytes  The number of bytes to read from the key
        */
        void add(const void* key, int num_bytes);

        /**
        * @brief Test whether an object is in the collection
        *
        * @param key         A pointer to the key data
        * @param num_bytes   The number of bytes to read from the key
        *
        * @return            true if the object is in the bloom filter
        */
        bool test(const void* key, int num_bytes) const;

        /**
        * @brief Print the amount of memory used to stdout
        */
        void printMemory() const;

        /**
        * @brief Count how many bits are set and print to stdout
        */
        void printOccupancy() const;

    private:
        BitVector m_bitvector;
        std::vector<uint32_t> m_hashes;
        size_t m_width;
        size_t m_occupancy;
#if TRACK_OCCUPANCY
        mutable size_t m_test_counter;
#endif
};

#endif
