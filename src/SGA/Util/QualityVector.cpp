//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// QualityVector - Vector of doubles representing
// log-transformed quality values for a sequence
//
#include <assert.h>
#include <algorithm>
#include <iostream>
#include "QualityVector.h"

//
QualityVector::QualityVector()
{

}

//
QualityVector::QualityVector(const QualityVector& vec, int start, int size)
{
    assert(!vec.empty());
    assert(size >= 0);
    assert(size + start <= (int)vec.size());
    m_data.insert(m_data.end(), vec.m_data.begin() + start, vec.m_data.begin() + start + size);
}

//
void QualityVector::add(DNADouble v)
{
    m_data.push_back(v);
}

//
void QualityVector::set(size_t idx, DNADouble v)
{
    if(m_data.size() < idx)
        m_data.resize(idx + 1);
    m_data[idx] = v;
}

//
DNADouble QualityVector::get(size_t idx) const
{
    assert(idx < m_data.size());
    return m_data[idx];
}

//
size_t QualityVector::size() const
{
    return m_data.size();
}

// 
bool QualityVector::empty() const
{
    return m_data.empty();
}

//
void QualityVector::reverseComplement()
{
    std::reverse(m_data.begin(), m_data.end());
    for(size_t i = 0; i < m_data.size(); ++i)
        m_data[i].complement();
}


