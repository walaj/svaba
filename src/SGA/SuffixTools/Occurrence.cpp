//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Occurrence.cpp - Data structure holding the number of times
// the letter b appears in the string S from S[0..i] (inclusive)
//
#include "Occurrence.h"
#include "SBWT.h"

// Initialize the counts from the bwt string b
void Occurrence::initialize(const BWTString& bwStr, int sampleRate)
{
    m_sampleRate = sampleRate;
    m_shift = calculateShiftValue(m_sampleRate);

    size_t l = bwStr.length();
    int num_samples = (l % m_sampleRate == 0) ? (l / m_sampleRate) : (l / m_sampleRate + 1);
    m_values.resize(num_samples);
    
    AlphaCount64 sum;
    for(size_t i = 0; i < l; ++i)
    {
        char currB = bwStr.get(i);
        sum.increment(currB);
        if(i % m_sampleRate == 0)
            m_values[i / m_sampleRate] = sum;
    }
}

// 
int Occurrence::calculateShiftValue(int divisor)
{
    assert(divisor > 0);
    assert(IS_POWER_OF_2(divisor));

    // m_sampleRate is a power of 2, count what bit is set
    unsigned int v = divisor;
    unsigned int c = 0; // c accumulates the total bits set in v

    while(v != 1)
    {
        v >>= 1;
        ++c;
    }
    assert(1 << c == divisor);
    return c;
}

//
void Occurrence::set(char a, size_t i, BaseCount s)
{
    m_values[i].set(a, s);
}

//
size_t Occurrence::getByteSize() const
{
    return m_values.size() * sizeof(AlphaCount64);
}

// Validate that the sampled occurrence array is correct
void Occurrence::validate(const BWTString& bwStr) const
{
    size_t l = bwStr.length();
    AlphaCount64 sum;
    for(size_t i = 0; i < l; ++i)
    {
        char currB = bwStr.get(i);
        sum.increment(currB);
        AlphaCount64 calculated = get(bwStr, i);
        for(int i = 0; i < ALPHABET_SIZE; ++i)
            assert(calculated.get(ALPHABET[i]) == sum.get(ALPHABET[i]));
    }
}

std::ostream& operator<<(std::ostream& out, const Occurrence& o)
{
    out << o.m_sampleRate << "\n";
    out << o.m_values.size() << "\n";
    for(size_t i = 0; i < o.m_values.size(); ++i)
        out << o.m_values[i] << "\n";
    return out;
}

std::istream& operator>>(std::istream& in, Occurrence& o)
{
    in >> o.m_sampleRate;
    size_t n;
    in >> n;
    o.m_values.resize(n);
    for(size_t i = 0; i < n; ++i)
        in >> o.m_values[i];
    o.m_shift = Occurrence::calculateShiftValue(o.m_sampleRate);
    return in;
}


//
void Occurrence::print() const
{
    for(size_t i = 0; i < m_values.size(); i++)
    {
        std::cout << m_values[i];
    }
}
