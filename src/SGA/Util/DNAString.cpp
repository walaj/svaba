//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// DNAString 
//
#include <iostream>
#include "DNAString.h"
#include "Util.h"

DNAString::DNAString() : m_len(0), m_data(0) {}

//
DNAString::DNAString(std::string seq)
{
    _alloc(seq.c_str(), seq.length());
}

//
DNAString::DNAString(const DNAString& other)
{
    // Deep copy
    _alloc(other.m_data, other.m_len);
}

//
DNAString& DNAString::operator=(const DNAString& dna)
{
    if(&dna == this)
        return *this; // self-assign

    _dealloc();
    _alloc(dna.m_data, dna.m_len);
    return *this;
}

//
DNAString& DNAString::operator=(const std::string& str)
{
    _dealloc();
    _alloc(str.c_str(), str.length());
    return *this;
}

//
bool DNAString::operator==(const DNAString& other)
{
    if(m_len != other.m_len)
        return false;
    return strncmp(m_data, other.m_data, m_len) == 0;
}

//
void DNAString::_alloc(const char* pData, size_t l)
{
    m_len = l;
    size_t alloc_len = m_len + 1;
    m_data = new char[alloc_len];
    strncpy(m_data, pData, m_len);
    m_data[m_len] = '\0';
}

//
void DNAString::_dealloc()
{
    if(m_data != NULL)
        delete [] m_data;
    m_data = 0;
    m_len = 0;
}

//
DNAString::~DNAString()
{
    _dealloc();
}

// Reverse the sequence
void DNAString::reverse()
{
    size_t half = m_len >> 1;
    for(size_t idx = 0; idx < half; ++idx)
    {
        size_t opp = m_len - idx - 1;
        char tmp = m_data[opp];
        m_data[opp] = m_data[idx];
        m_data[idx] = tmp;
    }
}

void DNAString::reverseComplement()
{
    this->reverse();
    for(size_t idx = 0; idx < m_len; ++idx) m_data[idx] = complement(m_data[idx]);    
}


void DNAString::disambiguate()
{
    if(m_len == 0)
        return;
    
    for(size_t idx = 0; idx < m_len; ++idx)
    {
        if(m_data[idx] == 'N')
        {
            m_data[idx] = randomBase();
        }
    }
}

//
std::string DNAString::getSuffixString(size_t idx) const
{
    std::string str;
    str.reserve(getSuffixLength(idx) + 1);
    str += getSuffix(idx);
    str += "$";
    return str;
}

//
std::string DNAString::toString() const
{
    return std::string(m_data, m_len);
}
