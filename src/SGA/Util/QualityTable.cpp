//-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// QualityTable - A 0-indexed table of quality strings
//
#include <iostream>
#include <algorithm>
#include "QualityTable.h"
#include "QualityCodec.h"
#include "SeqReader.h"

// Read the sequences from a file
QualityTable::QualityTable() : m_bytes_used(0)
{
    int missing_phred = 20;
    m_missingQualityChar = Quality::phred2char(missing_phred);
}

 
QualityTable::~QualityTable()
{
    for(size_t i = 0; i < m_table.size(); ++i)
        delete m_table[i].encoded_data;
}

//
void QualityTable::loadQualities(const std::string& filename)
{
    SeqReader reader(filename);
    SeqRecord sr;
    while(reader.get(sr))
        addQualityString(sr.qual);
}

//
void QualityTable::addQualityString(const std::string& qual)
{
    size_t num_bytes_required = m_codec.getRequiredUnits(qual.size());
    QualityString incoming;
    incoming.encoded_data = new QualityStorageUnit[num_bytes_required];
    incoming.num_encoded_symbols = qual.size();
    for(size_t i = 0; i < qual.size(); ++i)
        m_codec.store(incoming.encoded_data, i, qual[i]);
    m_table.push_back(incoming);

    //std::cout << "Stored:  " << qual << "\n";
    //std::cout << "Decoded: " << getQualityString(m_table.size() - 1, qual.size()) << "\n";
    m_bytes_used += (sizeof(incoming) + num_bytes_required);
}

//
std::string QualityTable::getQualityString(size_t idx, size_t n) const
{
    // If there is no quality string for this index, return default qualities
    if(idx >= m_table.size())
        return std::string(n, m_missingQualityChar);

    QualityString encoded = m_table[idx];
    std::string out;
    out.reserve(encoded.num_encoded_symbols);
    for(size_t i = 0; i < encoded.num_encoded_symbols; ++i)
        out.push_back(m_codec.get(encoded.encoded_data, i));
    assert(out.length() == n);
    return out;
}

//
size_t QualityTable::getCount() const
{
    return m_table.size();
}

//
void QualityTable::printSize() const
{
    printf("QualityTable: %.2lfGB\n", (double)m_bytes_used / (1000 * 1000 * 1000));
}

