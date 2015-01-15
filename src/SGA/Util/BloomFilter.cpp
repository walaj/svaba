//-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
#include "BloomFilter.h"
#include <cstdlib>
#include <assert.h>
#include <stdio.h>
#include "MurmurHash3.h"

//
BloomFilter::BloomFilter() : m_occupancy(0)
{

#if TRACK_OCCUPANCY
    m_test_counter = 0;
#endif

}

//
void BloomFilter::initialize(size_t width, size_t num_hashes)
{
    m_width = width;
    m_bitvector.resize(width);
    m_occupancy = 0;

    // Create seeds for murmer hash
    // TODO: Check how independent this method of selecting
    // hash function is...
    m_hashes.resize(num_hashes);
    for(size_t i = 0; i < m_hashes.size(); ++i)
        m_hashes[i] = rand();
}

//
void BloomFilter::add(const void* key, int num_bytes)
{
    // Set h bits
    for(size_t i = 0; i < m_hashes.size(); ++i) {
        int64_t h[2];
        MurmurHash3_x64_128(key, num_bytes, m_hashes[i], &h);
        size_t idx = h[0] % m_width;

        bool was_set = m_bitvector.updateCAS(idx, false, true);
#if TRACK_OCCUPANCY
        // This call uses an atomic update and returns true if the
        // bit was sucessfully set
        if(was_set)
        {
            // A bit was set, update the occupancy
            bool set = false;
            do {
                size_t old_occ = m_occupancy;
                size_t new_occ = old_occ + 1;
                set = __sync_bool_compare_and_swap(&m_occupancy, old_occ, new_occ);
            } while(!set);
        }
#else
        (void)was_set;
#endif
    }
}

//
bool BloomFilter::test(const void* key, int num_bytes) const
{
    for(size_t i = 0; i < m_hashes.size(); ++i) {
        int64_t h[2];
        MurmurHash3_x64_128(key, num_bytes, m_hashes[i], &h);
        size_t idx = h[0] % m_width;
        if(!m_bitvector.test(idx))
            return false;
    }

#if TRACK_OCCUPANCY
    m_test_counter++;
    if(m_test_counter % 500000 == 0) {
        double p = (double)m_occupancy / m_width;
        printf("WIDTH: %zu OCC: %zu\n", m_width, m_occupancy);
        printf("Access: %zu Occupancy: %zu P: %1.4lf\n", m_test_counter, m_occupancy, p);
    }
#endif
    return true;
}

//
void BloomFilter::printOccupancy() const
{
    size_t set_count = 0;
    for(size_t i = 0; i < m_width; ++i)
        set_count += m_bitvector.test(i);
    printf("%zu out of %zu bits are set\n", set_count, m_width);
}

//
void BloomFilter::printMemory() const
{
    size_t bytes = m_bitvector.capacity();
    double mb = (double)bytes / (1 << 20);
    printf("BloomFilter using %.1lf MB\n", mb);
}
