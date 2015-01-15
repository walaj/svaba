//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTWriterAscii - Write a BWT file to disk
//
#include "BWTWriterAscii.h"
#include "SBWT.h"
#include "RLBWT.h"

//
BWTWriterAscii::BWTWriterAscii(const std::string& filename) : m_stage(IOS_NONE)
{
    m_pWriter = createWriter(filename);
    m_stage = IOS_HEADER;
}

//
BWTWriterAscii::~BWTWriterAscii()
{
    delete m_pWriter;
}

void BWTWriterAscii::write(const SBWT* pBWT)
{
    writeHeader(pBWT->m_numStrings, pBWT->m_bwStr.length(), BWF_HASFMI);
    writeBWStr(pBWT->m_bwStr);
    writePred(pBWT->m_predCount);
    writeOccurrence(pBWT->m_occurrence);
}

void BWTWriterAscii::write(const RLBWT* pRLBWT)
{
    writeHeader(pRLBWT->m_numStrings, pRLBWT->m_numSymbols, BWF_NOFMI);
    size_t numRuns = pRLBWT->getNumRuns();
    for(size_t i = 0; i < numRuns; ++i)
    {
        const RLUnit& unit = pRLBWT->m_rlString[i];
        char symbol = unit.getChar();
        size_t length = unit.getCount();
        for(size_t j = 0; j < length; ++j)
            writeBWChar(symbol);
    }
    // Finalize the string
    writeBWChar('\n');
}

//
void BWTWriterAscii::writeHeader(const size_t& num_strings, const size_t& num_symbols, const BWFlag& flag)
{
    assert(m_stage == IOS_HEADER);
    *m_pWriter << BWT_FILE_MAGIC << "\n";
    *m_pWriter << num_strings << "\n";
    *m_pWriter << num_symbols << "\n";
    int temp = flag;
    *m_pWriter << temp << "\n";
    m_stage = IOS_BWSTR;    
}

//
void BWTWriterAscii::writeBWStr(const BWTString& str)
{
    assert(m_stage == IOS_BWSTR);
    size_t n = str.length();
    for(size_t i = 0; i < n; ++i)
        *m_pWriter << str.get(i);
    finalize();
}

// Write a single character of the BWStr
// If the char is '\n' we are finished
void BWTWriterAscii::writeBWChar(char b)
{
    m_pWriter->put(b);
}

void BWTWriterAscii::finalize()
{
    m_pWriter->put('\n');
    m_stage = IOS_PC;
}
//
void BWTWriterAscii::writePred(const AlphaCount64& pc)
{
    assert(m_stage == IOS_PC);
    *m_pWriter << pc << "\n";
    m_stage = IOS_OCC;
}

//
void BWTWriterAscii::writeOccurrence(const Occurrence& occ)
{
    assert(m_stage == IOS_OCC);
    *m_pWriter << occ; // delibrately no newline
    m_stage = IOS_DONE;
}
