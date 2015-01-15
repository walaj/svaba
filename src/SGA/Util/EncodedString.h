//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// EncodedString - Templated class to store
// a string from a reduced alphabet. The actual
// encoding is done outside this class by the Codec
// template parameter.
//
// Here the storage is in terms of "units" which is
// defined by the codec. For a simple 2-bit encoded
// DNA string this would be a byte but for more complicated
// encodings like the 3-bit BWT this could be a larger
// structure like uint16_t
//
#ifndef ENCODEDSTRING_H
#define ENCODEDSTRING_H
#include <string.h>
#include <string>
#include <assert.h> 
#include "DNACodec.h"
#include "BWTCodec.h"
#include "BWT4Codec.h"
#include "NoCodec.h"

template<class Codec>
class EncodedString
{
    typedef typename Codec::UNIT_TYPE StorageUnit;

    public:
        
        // Constructors/Destructors
        EncodedString() : m_len(0), m_capacity(0), m_data(0) {}

        //
        EncodedString(const EncodedString& other)
        {
            // deep copy
            _alloc(other.m_len);
            _copy(other);
        }

        //
        EncodedString(const std::string& seq)
        {
            size_t n = seq.length();
            _alloc(n);
            _copy(seq.c_str(), n);
        }

        //
        ~EncodedString()
        {
            _dealloc();
        }

        // Assignment op
        EncodedString& operator=(const EncodedString& other)
        {
            if(&other == this)
                return *this; // self-assign

            _dealloc();
            _alloc(other.m_len);
            _copy(other);
            return *this;
        }

        //
        EncodedString& operator=(const std::string& str)
        {
            size_t n = str.length();
            _dealloc();
            _alloc(n);
            _copy(str.c_str(), n);
            return *this;
        }

        // Resize this string to n symbols, setting the 
        // new entries to the default value
        void resize(size_t n)
        {
            if(n > m_capacity)
                _realloc(n);
            m_len = n;
        }

        // Append a std::string
        void append(const std::string& str)
        {
            size_t n = str.length();
            size_t num_total = m_len + n;
            if(num_total > m_capacity)
                _realloc(num_total);
            _append(str.c_str(), n);
        }

        // Append some other EncodedString
        void append(const EncodedString& other)
        {
            size_t n = other.length();
            size_t num_total = m_len + n;
            if(num_total > m_capacity)
                _realloc(num_total);
            _append(other);
        }

        // Swap the contents with another encoded string
        void swap(EncodedString& other)
        {
            size_t tmp;
            tmp = other.m_len;
            other.m_len = m_len;
            m_len = tmp;

            tmp = other.m_capacity;
            other.m_capacity = m_capacity;
            m_capacity = tmp;

            StorageUnit* pTmp = other.m_data;
            other.m_data = m_data;
            m_data = pTmp;
        }

        //
        size_t length() const
        {
            return m_len;
        }

        //
        size_t capacity() const
        {
            return m_capacity;
        }

        //
        bool empty() const
        {
            return m_len == 0;
        }

        // Get the character at idx
        inline char get(size_t idx) const
        {
            assert(idx < m_len);
            return s_codec.get(m_data, idx);
        }

        // Set the character at idx
        inline void set(size_t idx, char b)
        {
            assert(idx < m_len);
            s_codec.store(m_data, idx, b);
        }

        //
        std::string toString() const
        {
            std::string out(m_len, 'A');
            for(size_t i = 0; i < m_len; ++i)
                out[i] = s_codec.get(m_data, i);
            return out;
        }

        //
        std::string substr(size_t start) const
        {
            assert(start < m_len);
            std::string out(m_len - start, 'A');
            for(size_t i = start; i < m_len; ++i)
                out[i - start] = s_codec.get(m_data, i);
            return out;
        }

        //
        std::string substr(size_t start, size_t len) const
        {
            assert(start < m_len);
            assert(start + len <= m_len);
            std::string out(len, 'C');
            for(size_t i = 0; i < len; ++i)
                out[i] = s_codec.get(m_data, i + start);
            return out;
        }

        // 
        friend bool operator==(const EncodedString& a, const EncodedString& b)
        {
            if(a.m_len != b.m_len)
                return false;

            // Would be faster to compare entire units at a time
            size_t n = a.m_len;
            for(size_t i = 0; i < n; ++i)
            {
                if(a.get(i) != b.get(i))
                    return false;
            }
            return true;
        }

        //
        friend std::ostream& operator<<(std::ostream& out, const EncodedString<Codec>& a)
        {
            out << a.toString();
            return out;
        }

        // Return the amount of space this string is using
        size_t getMemSize() const
        {
            return sizeof(*this) + sizeof(StorageUnit) * s_codec.getRequiredUnits(m_capacity);
        }

    private:

        // functions
        void _copy(const char* pData, size_t n)
        {
            // this assumes that storage for n characters has been 
            // allocated
            assert(m_capacity >= n);
            for(size_t i = 0; i < n; ++i)
                s_codec.store(m_data, i, pData[i]);
            m_len = n;
        }
        
        //
        void _copy(const EncodedString& other)
        {
            // storage should have been allocated already
            assert(m_capacity = other.m_capacity);
            size_t num_units = s_codec.getRequiredUnits(other.m_capacity);
            _copyUnitData(other.m_data, num_units);
            m_len = other.m_len;
        }

        // Copy num_units from pData into the internal storage
        void _copyUnitData(StorageUnit* pData, size_t num_units)
        {
            assert(s_codec.getRequiredUnits(m_capacity) >= num_units);
            for(size_t i = 0; i < num_units; ++i)
                m_data[i] = pData[i];
        }

        // append n symbols from pData into the buffer
        void _append(const char* pData, size_t n)
        {
            assert(m_len + n <= m_capacity);
            for(size_t i = 0; i < n; ++i)
                s_codec.store(m_data, m_len + i, pData[i]);
            m_len += n;
        }

        // append 
        void _append(const EncodedString& other)
        {
            size_t n = other.m_len;
            assert(m_len + n <= m_capacity);
            for(size_t i = 0; i < n; ++i)
                s_codec.store(m_data, m_len + i, other.get(i));
            m_len += n;
        }

        // allocate storage for n symbols
        void _alloc(size_t n)
        {
            // Get the number of units that need to be allocated from the codec
            size_t n_units = s_codec.getRequiredUnits(n);
            m_data = new StorageUnit[n_units](); // This zeros the memory
            m_capacity = s_codec.getCapacity(n_units);
        }

        // reallocate the storage so that the capacity is at least n symbols
        // m_len is not changed
        void _realloc(size_t n)
        {
            StorageUnit* oldData = m_data;
            size_t old_units = s_codec.getRequiredUnits(m_len);
            _alloc(n);
            assert(m_capacity >= n);

            _copyUnitData(oldData, old_units);
            delete [] oldData;
        }

        // deallocate storage
        void _dealloc()
        {
            if(m_data != NULL)
                delete [] m_data;
            m_data = 0;
            m_len = 0;
            m_capacity = 0;
        }

        // data
        static Codec s_codec;
        size_t m_len; // the length of the string
        size_t m_capacity; // the maximum length of the string that can be stored
        StorageUnit* m_data;
};

// Initialize the static member
template<class Codec>
Codec EncodedString<Codec>::s_codec;

typedef EncodedString<DNACodec> DNAEncodedString;
typedef EncodedString<BWT4Codec> BWTString;
typedef EncodedString<NoCodec> NoEncodingString;

typedef std::vector<DNAEncodedString> DNAEncodedStringVector;

#endif
