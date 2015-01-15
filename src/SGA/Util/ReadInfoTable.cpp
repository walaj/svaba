//
//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ReadInfoTable - A 0-indexed table of ID, length pairs
// Used to convert suffix array hits to overlaps
//
#include <iostream>
#include <algorithm>
#include "ReadInfoTable.h"
#include "SeqReader.h"

//JEREMIAH
ReadInfoTable::ReadInfoTable(ReadTable *pRT) : m_numericIDs(false) {
  for (size_t i = 0; i < pRT->getCount(); i++) {
    m_lengths.push_back(pRT->getReadLength(i));
    m_ids.push_back(pRT->getReadID(i));
  }
}

// Read the sequences from a file
ReadInfoTable::ReadInfoTable(std::string filename, 
                             size_t num_expected, 
                             ReadInfoOption option) : m_numericIDs(false)
{
    // Do not store actual ids, use a numeric id equal to the table index
    if(option == RIO_NUMERICID)
        m_numericIDs = true;

    if(num_expected > 0)
    {
        m_lengths.reserve(num_expected);
        if(!m_numericIDs)
            m_ids.reserve(num_expected);
    }

    SeqReader reader(filename);
    SeqRecord sr;

    // Load the lengths and ids
    while(reader.get(sr))
    {
        m_lengths.push_back(sr.seq.length());
        if(!m_numericIDs)
        {
            m_ids.push_back(sr.id);
        }
    }
}

// 
ReadInfoTable::~ReadInfoTable()
{

}

//
size_t ReadInfoTable::getReadLength(size_t idx) const
{
    assert(idx < m_lengths.size());
    return m_lengths[idx];
}

//
std::string ReadInfoTable::getReadID(size_t idx) const
{
    if(!m_numericIDs)
    {
        assert(idx < m_ids.size());
        return m_ids[idx];
    }
    else
    {
        // Build an identified from the idx
        std::stringstream idss;
        idss << idx;
        return idss.str();
    }
}

//
const ReadInfo ReadInfoTable::getReadInfo(size_t idx) const
{
    ReadInfo ri(getReadID(idx), getReadLength(idx));
    return ri;
}

//
size_t ReadInfoTable::getCount() const
{
    return m_lengths.size();
}

// 
size_t ReadInfoTable::countSumLengths() const
{
    size_t sum = 0;
    for(size_t i = 0; i < m_lengths.size(); ++i)
        sum += m_lengths[i];
    return sum;
}

//
void ReadInfoTable::clear()
{
    m_lengths.clear();
    m_ids.clear();
}

