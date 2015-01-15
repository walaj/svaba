//-------------------------------------------------------------------------------
// 
// CountMinSketch - Implementation of the CountMinSketch probabilistic data
// structure as described here: http://sites.google.com/site/countminsketch/
// 
// Copyright (C) 2012 Jared Simpson (jared.simpson@gmail.com)
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
// ------------------------------------------------------------------------------
#include "count_min_sketch.h"
#include "MurmurHash3.h"
#include <cstdlib>
#include <assert.h>
#include <stdio.h>

//
CountMinSketch::CountMinSketch() : m_max_count(std::numeric_limits<CMSData>::max())
{

}

//
void CountMinSketch::initialize(size_t width, size_t num_hashes)
{
    // Clear old data, if any
    for(size_t i = 0; i < m_vectors.size(); ++i)
        m_vectors[i].clear();
    
    // Resize the vectors
    m_vectors.resize(num_hashes);
    for(size_t i = 0; i < m_vectors.size(); ++i)
        m_vectors[i].resize(width);

    // Create seeds for murmer hash
    // TODO: Check how independent this method of selecting
    // hash function is...
    m_hashes.resize(num_hashes);
    for(size_t i = 0; i < m_hashes.size(); ++i)
        m_hashes[i] = rand();
}

//
void CountMinSketch::increment(const void* key, int num_bytes)
{
    assert(!m_vectors.empty());
    // Increment all vectors
    for(size_t i = 0; i < m_vectors.size(); ++i) {
        int64_t h[2];
        MurmurHash3_x64_128(key, num_bytes, m_hashes[i], &h);
        size_t idx = h[0] % m_vectors[i].size();
        
        // Perform an atomic compare and swap to increment the value
        // If the value has reached saturation, do not update
        while(1) {
            CMSData v = m_vectors[i][idx];
            if(v == m_max_count || __sync_bool_compare_and_swap(&m_vectors[i][idx], v, v + 1))
                break;
        }
    }
}

//
CMSData CountMinSketch::get(const void* key, int num_bytes) const
{
    // Return the min count in all vectors
    CMSData min = std::numeric_limits<uint8_t>::max();
    for(size_t i = 0; i < m_vectors.size(); ++i) {
        int64_t h[2];
        MurmurHash3_x64_128(key, num_bytes, m_hashes[i], &h);
        size_t idx = h[0] % m_vectors[i].size();
        if(m_vectors[i][idx] < min)
            min = m_vectors[i][idx];
    }
    return min;
}

//
void CountMinSketch::printMemory() const
{
    size_t bytes = sizeof(CMSData) * m_vectors.size() * m_vectors.front().size();
    double mb = (double)bytes / (1 << 20);
    printf("CountMinSketch using %.1lf MB\n", mb);
}
