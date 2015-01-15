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
#ifndef COUNT_MIN_SKETCH_H
#define COUNT_MIN_SKETCH_H

#include <vector>
#include <stdint.h>
#include <stddef.h>
#include <limits>

// Typedefs
typedef uint8_t CMSData; //single-byte counts, for now
typedef std::vector<CMSData> CMSDataVector;

class CountMinSketch
{
    public:

        /**
        * @brief Default constructor. The returned sketch is uninitialized.
        */
        CountMinSketch();

        /**
        * @brief Initialize the sketch.
        *
        * @param width       The number of bins to use per hash table
        * @param num_hashes  The number of hash tables to use
        */
        void initialize(size_t width, size_t num_hashes);

        /**
        * @brief Increment the counts for the given key
        *        This using gcc's atomic compare and swap to allow
        *        atomic increments to occur when running with multiple threads
        *
        * @param key        A pointer to the key data
        * @param num_bytes  The number of bytes to read from the key
        */
        void increment(const void* key, int num_bytes);

        /**
        * @brief Get the approximate count for the given key
        *
        * @param key         A pointer to the key data
        * @param num_bytes   The number of bytes to read from the key
        *
        * @return            The approximate count
        */
        CMSData get(const void* key, int num_bytes) const;

        /**
        * @brief Print the amount of memory used to stdout
        */
        void printMemory() const;

    private:
        std::vector<CMSDataVector> m_vectors;
        std::vector<uint32_t> m_hashes;
        const CMSData m_max_count;
};

#endif
