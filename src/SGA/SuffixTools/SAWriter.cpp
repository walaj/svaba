//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// SAWriter.h - Read a suffix array file from disk
//
#include "SAWriter.h"
#include "SuffixArray.h"

//
SAWriter::SAWriter(const std::string& filename) : m_stage(SAIOS_NONE)
{
    m_pWriter = createWriter(filename);
    m_stage = SAIOS_HEADER;
}

//
SAWriter::~SAWriter()
{
    delete m_pWriter;
}

void SAWriter::write(const SuffixArray* pSA)
{
    writeHeader(pSA->m_numStrings, pSA->m_data.size());
    writeElems(pSA->m_data);
}

//
void SAWriter::writeHeader(const size_t& num_strings, const size_t& num_elems)
{
    assert(m_stage == SAIOS_HEADER);
    *m_pWriter << SA_FILE_MAGIC << "\n";
    *m_pWriter << num_strings << "\n";
    *m_pWriter << num_elems << "\n";
    m_stage = SAIOS_ELEM;    
}

//
void SAWriter::writeElems(const SAElemVector& elemVector)
{
    assert(m_stage == SAIOS_ELEM);
    for(size_t i = 0; i < elemVector.size(); ++i)
        writeElem(elemVector[i]);
    m_stage = SAIOS_DONE;
}

//
void SAWriter::writeElem(const SAElem& elem)
{
    *m_pWriter << elem << "\n";
}

