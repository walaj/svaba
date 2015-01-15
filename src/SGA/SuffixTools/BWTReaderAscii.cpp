//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTReaderAsciiAscii - Read a BWT file from disk
//
#include "BWTReaderAscii.h"
#include "SBWT.h"
#include "RLBWT.h"

//
BWTReaderAscii::BWTReaderAscii(const std::string& filename) : m_stage(IOS_NONE)
{
    m_pReader = createReader(filename);
    m_stage = IOS_HEADER;
}

//
BWTReaderAscii::~BWTReaderAscii()
{
    delete m_pReader;
}

void BWTReaderAscii::read(SBWT* pBWT)
{
    size_t n;
    BWFlag flag;
    readHeader(pBWT->m_numStrings, n, flag);

    pBWT->m_bwStr.resize(n);
    readBWStr(pBWT->m_bwStr);
 
    // Reading the occurrence array from disk is deprecated
    // we ignore it if it is present
}

void BWTReaderAscii::read(RLBWT* pRLBWT)
{
    size_t n;
    BWFlag flag;
    readHeader(pRLBWT->m_numStrings, n, flag);

    bool done = false;
    while(!done)
    {
        char b = readBWChar();
        if(b != '\n')
        {
            pRLBWT->append(b);
        }
        else
        {
            done = false;
            break;
        }
    }
}

//
void BWTReaderAscii::readHeader(size_t& num_strings, size_t& num_symbols, BWFlag& flag)
{
    assert(m_stage == IOS_HEADER);
    uint16_t magic_number;

    // Ensure the file format is sane
    *m_pReader >> magic_number;
    if(magic_number != BWT_FILE_MAGIC)
    {
        std::cerr << "BWT file is not properly formatted, aborting\n";
        exit(EXIT_FAILURE);
    }
    *m_pReader >> num_strings;
    *m_pReader >> num_symbols;
    int temp;
    *m_pReader >> temp;
    flag = static_cast<BWFlag>(temp);
    // We must explicitly set the stream
    // to the start of the BWStr portion of
    // the file. The above extraction
    // does not discard the trailing newline 
    // so we do it here
    m_pReader->get();

    m_stage = IOS_BWSTR;    
}

//
void BWTReaderAscii::readBWStr(BWTString& out_str)
{
    assert(m_stage == IOS_BWSTR);
    char b;
    size_t idx = 0;
    while(1)
    {
        m_pReader->get(b);
        if(b != '\n')
            out_str.set(idx++, b);
        else
            break;
    }

    assert(idx == out_str.length());
    m_stage = IOS_PC;
}

// Read a single base from the BWStr
char BWTReaderAscii::readBWChar()
{
    assert(m_stage == IOS_BWSTR);
    char b;
    m_pReader->get(b);
    if(b == '\n')
        m_stage = IOS_PC;
    return b;
}

//
void BWTReaderAscii::readPred(AlphaCount64& out_pc)
{
    assert(m_stage == IOS_PC);
    *m_pReader >> out_pc;
    m_stage = IOS_OCC;
}

//
void BWTReaderAscii::readOccurrence(Occurrence& out_occ)
{
    assert(m_stage == IOS_OCC);
    *m_pReader >> out_occ;
    m_stage = IOS_DONE;
}
