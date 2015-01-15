//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// DNAString - Wrapper for c-string for DNA sequences
// It is tailored for fast suffix operations
//

#ifndef DNASTRING_H
#define DNASTRING_H
#include <string.h>
#include <string>
#include <assert.h> 

class DNAString
{
    public:
        
        // Constructors/Destructors
        DNAString();
        DNAString(const DNAString& other);
        DNAString(std::string seq);
        ~DNAString();

        // Operators
        DNAString& operator=(const DNAString& dna);
        DNAString& operator=(const std::string& str);
        bool operator==(const DNAString& other);

        size_t length() const
        {
            return m_len;
        }

        bool empty() const
        {
            return m_len == 0 || m_data[0] == '\0';
        }

        inline const char* getSuffix(size_t idx) const
        {
            assert(m_len > 0);
            // Force the suffix to point to the empty string if out of bounds
            if(idx <= m_len)
                return m_data + idx;
            else
                return m_data + m_len;
        }
        
        // Return the length of the suffix (not including $) starting at idx
        size_t getSuffixLength(size_t idx) const
        {
            assert(m_len > 0);
            if(idx <= m_len)
                return m_len - idx;
            else
                return 0;
        }

        // Get the character at the given position
        inline char get(size_t idx) const
        {
            assert(idx < m_len + 1);
            return m_data[idx];
        }

        // Set the character at the given position
        inline void set(size_t idx, char b) const
        {
            assert(idx < m_len + 1);
            m_data[idx] = b;
        }

        // Get the substring of length n starting at position pos
        std::string substr(size_t pos, size_t n) const
        {
            assert(m_data != NULL && pos < m_len && pos + n <= m_len);
            return std::string(m_data + pos, n);
        }

        // randomly change the 'N' bases to one of ACGT so that the string is not ambiguous
        void disambiguate();

        void reverse();
        void reverseComplement();
        std::string getSuffixString(size_t idx) const;
        std::string toString() const;

    private:

        // functions
        void _alloc(const char* pData, size_t l);
        void _dealloc();

        // data
        size_t m_len;
        char* m_data;
};

#endif
