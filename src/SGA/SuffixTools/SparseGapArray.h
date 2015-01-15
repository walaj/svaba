//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SparseGapArray
//
// The sparse gap array uses a small integer type to
// keep the count for each element. A small number
// of positions will have a count that is larger
// than the maximum value the integer can represent.
// These values are automatically overflowed into a 
// hash map with a larger integer type. In the typical
// case this will require far less memory than the
// SimpleGapArray. The storage unit is templated
// to allow a dynamic tradeoff between the amount
// of memory required to store the simple cases 
// and the overflow rate.
//
#ifndef SPARSEGAPARRAY_H
#define SPARSEGAPARRAY_H

#include "HashMap.h"
#include "GapArray.h"
#include "BitVector.h"

// Template base storage for the sparse gap array
template<class IntType>
class SparseBaseStorage
{
    public:
        //
        void resize(size_t n)
        {
            m_data.resize(n);
        }

        //
        inline void set(size_t i, IntType c)
        {
            m_data[i] = c;
        }

        // Attempt to set the value using an atomic update
        // Returns false if the update fails
        inline bool setCAS(size_t i, IntType oldV, IntType newV)
        {
            assert(oldV != getMax() && newV <= getMax() && oldV != newV);
            bool success = __sync_bool_compare_and_swap(&m_data[i], oldV, newV);
            return success;
        }

        //
        inline IntType get(size_t i) const
        {
            return m_data[i];
        }

        inline static size_t getMax()
        {
            static size_t max = std::numeric_limits<IntType>::max();
            return max;
        }

        //
        size_t size() const
        {
            return m_data.size();
        }

    private:
        //
        typedef std::vector<IntType> BaseStorageVector;
        BaseStorageVector m_data;
};

static const uint8_t storage4_offset_mask[]={0xF0,0x0F}; // 11110000, 00001111
static const uint8_t storage4_shift[]={4, 0};

// specialized 4-bit per value base storage for the sparse gap array
class SparseBaseStorage4
{
    public:

        SparseBaseStorage4() : m_numElems(0) {}

        //
        void resize(size_t n)
        {
            // With 4 bits per entry, we need n / 2 bytes of storage
            size_t numBytes = ((n & 1) == 0) ? n / 2 : (n / 2) + 1;
            m_data.resize(numBytes);
            m_numElems = n;
        }

        //
        inline void set(size_t i, uint8_t c)
        {
            assert(i < m_numElems);
            assert(c <= getMax());

            size_t idx = i >> 1; // equivalent to i / 2
            size_t offset = i & 1; // 0 if even, 1 if odd
            uint8_t mask = storage4_offset_mask[offset];
            uint8_t& elem = m_data[idx];

            // Clear current value
            elem &= ~mask;
            
            // Set new value
            elem |= c << storage4_shift[offset];
        }

        // Attempt to set the value using an atomic update
        // Returns false if the update fails
        inline bool setCAS(size_t i, uint8_t oldV, uint8_t newV)
        {
            assert(i < m_numElems);
            assert(oldV < getMax() && newV <= getMax());
            assert(oldV != newV);
            size_t idx = i >> 1; // equivalent to i / 2
            size_t offset = i & 1; // 0 if even, 1 if odd
            uint8_t mask = storage4_offset_mask[offset];

            uint8_t oldElem = m_data[idx];
            
            // If the count in the stored element has changed since this
            // function was called, return false
            uint8_t currV = (oldElem & mask) >> storage4_shift[offset];
            if(currV != oldV)
                return false;
            
            // Calculate the new element to swap in
            uint8_t newElem = oldElem;
            newElem &= ~mask;
            newElem |= newV << storage4_shift[offset];

            // Perform the update
            bool success = __sync_bool_compare_and_swap(&m_data[idx], oldElem, newElem);
            return success;
        }

        //
        inline uint8_t get(size_t i) const
        {
            assert(i < m_numElems);

            size_t idx = i >> 1; // equivalent to i / 2
            size_t offset = i & 1; // 0 if even, 1 if odd
            uint8_t mask = storage4_offset_mask[offset];
            uint8_t elem = m_data[idx];

            elem &= mask;
            elem >>= storage4_shift[offset];
            assert(elem <= getMax());
            return elem;
        }

        inline static size_t getMax()
        {
            return 15;
        }

        //
        size_t size() const
        {
            return m_numElems;
        }

    private:
        
        //

        typedef std::vector<uint8_t> BaseStorageVector;
        BaseStorageVector m_data;
        size_t m_numElems;
};

// specialized 1-bit per value base storage for the sparse gap array
// used for removing entries from an fm-index when the count
// will not exceed 1.
class SparseBaseStorage1
{
    public:

        SparseBaseStorage1() : m_numElems(0) {}

        //
        void resize(size_t n)
        {
            m_numElems = n;
            m_bitVector.resize(n);
        }

        // Set bit i from oldV to newV using a compare and swap
        // Returns false if the update fails
        inline bool setCAS(size_t i, uint8_t oldV, uint8_t newV)
        {
            assert(i < m_numElems);
            assert(oldV == 0);
            assert(newV == 1);
            return m_bitVector.updateCAS(i, oldV, newV);
        }

        //
        inline void set(size_t i, uint8_t c)
        {
            assert(i < m_numElems);
            assert(c <= getMax());

            m_bitVector.set(i, true);
        }

        //
        inline uint8_t get(size_t i) const
        {
            assert(i < m_numElems);
            return m_bitVector.test(i) ? 1 : 0;
        }

        inline static size_t getMax()
        {
            return 1;
        }

        //
        size_t size() const
        {
            return m_numElems;
        }

    private:
        
        //
        BitVector m_bitVector;
        size_t m_numElems;
};

// The SparseGapArray has two levels of storage.
// The first level uses x bits to store counts
// up to 2**x for n elements in the array. 
// If the count for a particular element exceeds 2**x
// then an entry in the overflow hash table is created
// allowing arbitrary values to be stored. This
// class is optimized to allow the base storage to be updated
// concurrently with compare and swap operations. Any overflowed
// updates must be serialized though incrementOverflowSerial
// 
template<class BaseStorage, class OverflowStorage>
class SparseGapArray : public GapArray
{
    public:
        SparseGapArray() : m_rankZeroCount(0) {}
        ~SparseGapArray()
        {
            //printf("SparseGapArray -- n: %zu overflow: %zu (%lf)\n", m_baseStorage.size(), m_overflow.size(), (double)m_overflow.size() / m_baseStorage.size());    
        }

        //
        void resize(size_t n)
        {
            m_baseStorage.resize(n);
        }

        // Attempt to increment a value in the GapArray using a compare and swap function
        // This call can fail to perform the update and return false. In this case
        // the calling code is responsible for serializing access to this data structure
        // and calling incrementSerial
        bool attemptBaseIncrement(size_t i)
        {
            assert(i < m_baseStorage.size());

            bool success = false;
            // Rank zero optimization
            // When merging two indices, all reads
            // start at rank 0. This would cause many serial
            // updates to the overflow array so we optimize 
            // for this case by doing a compare and swap update
            // of a large integer value if the rank is zero.
            if(i == 0)
            {
                do
                {
                    size_t count = m_rankZeroCount;
                    success = __sync_bool_compare_and_swap(&m_rankZeroCount, count, count+1);
                }
                while(!success);
                return true;
            }

            success = false;
            do
            {
                size_t count = m_baseStorage.get(i);
                if(count == getBaseMax())
                    return false;
                success = m_baseStorage.setCAS(i, count, count + 1);
            } while(!success);
            return success;
        }

        // Increment the value for rank i in the overflow array. This call must be serialized
        // between any threads.
        void incrementOverflowSerial(size_t i)
        {
            assert(i != 0);
            assert(i < m_baseStorage.size());    
            size_t count = m_baseStorage.get(i);
            assert(count == getBaseMax());
            
            // Check if the overflow map has a value for this index already
            // if not, enter one
            if(m_overflow.count(i) == 0)
                initOverflow(i, count);

            incrementOverflow(i);
        }

        //
        void initOverflow(size_t i, OverflowStorage c)
        {
            m_overflow.insert(std::make_pair(i, c));
        }

        //
        void incrementOverflow(size_t i)
        {
            // precondition: there must be an entry in the hash for this key
            assert(m_overflow.count(i) != 0);
            ++m_overflow[i];
        }

        //
        size_t get(size_t i) const
        {
            // rank zero optimization
            if(i == 0)
                return m_rankZeroCount;

            size_t count = m_baseStorage.get(i);
            if(count == getBaseMax())
            {
                typename OverflowHash::const_iterator iter = m_overflow.find(i);

                // If there is no entry in the overflow table yet
                // the count is exactly the maximum value representable
                // in the base storage
                if(iter == m_overflow.end())
                    return count;
                else
                    return iter->second;
            }
            else
            {
                return count;
            }
        }

        //
        size_t getBaseMax() const
        {
            return BaseStorage::getMax();
        }

        //
        size_t size() const
        {
            return m_baseStorage.size();
        }

   private:

        typedef SparseHashMap<size_t, OverflowStorage> OverflowHash;
        OverflowHash m_overflow;
        BaseStorage m_baseStorage;
        size_t m_rankZeroCount;
};

typedef SparseGapArray<SparseBaseStorage1, size_t> SparseGapArray1;
typedef SparseGapArray<SparseBaseStorage4, size_t> SparseGapArray4;
typedef SparseGapArray<SparseBaseStorage<uint8_t>, size_t> SparseGapArray8;
typedef SparseGapArray<SparseBaseStorage<uint16_t>, size_t> SparseGapArray16;
typedef SparseGapArray<SparseBaseStorage<uint32_t>, size_t> SparseGapArray32;

#endif
