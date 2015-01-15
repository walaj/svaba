//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// RLBWTReaderBinary - Read a run length encoded BWT file from disk
//
#include "BWTReaderBinary.h"
#include "SBWT.h"
#include "RLBWT.h"

//
BWTReaderBinary::BWTReaderBinary(const std::string& filename) : m_stage(IOS_NONE), m_numRunsOnDisk(0), m_numRunsRead(0)
{
    m_pReader = createReader(filename, std::ios::binary);
    m_stage = IOS_HEADER;
}

//
BWTReaderBinary::~BWTReaderBinary()
{
    delete m_pReader;
}

void BWTReaderBinary::read(RLBWT* pRLBWT)
{
    BWFlag flag;
    readHeader(pRLBWT->m_numStrings, pRLBWT->m_numSymbols, flag);

    assert(m_numRunsOnDisk > 0);
    readRuns(pRLBWT->m_rlString, m_numRunsOnDisk);

    //pRLBWT->printInfo();
    //pRLBWT->print();
}

void BWTReaderBinary::read(SBWT* pSBWT)
{
    BWFlag flag;
    size_t numSymbols;
    readHeader(pSBWT->m_numStrings, numSymbols, flag);
    size_t numRead = 0;
    pSBWT->m_bwStr.resize(numSymbols);
    while(numRead < numSymbols)
    {
        char b = readBWChar();
        pSBWT->m_bwStr.set(numRead++, b);
    }

    assert(m_numRunsOnDisk > 0);
}

//
void BWTReaderBinary::readHeader(size_t& num_strings, size_t& num_symbols, BWFlag& flag)
{
    assert(m_stage == IOS_HEADER);
    uint16_t magic_number;
    m_pReader->read(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));
    
    if(magic_number != RLBWT_FILE_MAGIC)
    {
        std::cerr << "BWT file is not properly formatted, aborting\n";
        exit(EXIT_FAILURE);
    }

    m_pReader->read(reinterpret_cast<char*>(&num_strings), sizeof(num_strings));
    m_pReader->read(reinterpret_cast<char*>(&num_symbols), sizeof(num_symbols));
    m_pReader->read(reinterpret_cast<char*>(&m_numRunsOnDisk), sizeof(m_numRunsOnDisk));
    m_pReader->read(reinterpret_cast<char*>(&flag), sizeof(flag));
    
    //std::cout << "Read magic: " << magic_number << "\n";
    //std::cout << "strings:" << num_strings << "\n";
    //std::cout << "symbols: " << num_symbols << "\n";
    //std::cout << "runs: " << m_numRunsOnDisk << "\n";
    
    m_stage = IOS_BWSTR;    
}

void BWTReaderBinary::readRuns(RLVector& out, size_t numRuns)
{
    out.resize(numRuns);
    m_pReader->read(reinterpret_cast<char*>(&out[0]), numRuns*sizeof(RLUnit));
    m_numRunsRead = numRuns;
}

// Read a single base from the BWStr
// The BWT is stored as runs on disk, so this class keeps
// an internal buffer of a single run and emits characters from this buffer
// and performs reads as necessary. If all the runs have been read, emit
// a newline character to signal the end of the BWT
char BWTReaderBinary::readBWChar()
{
    assert(m_stage == IOS_BWSTR);

    if(m_currRun.isEmpty())
    {
        // All runs have been read and emitted, return the end marker
        if(m_numRunsRead == m_numRunsOnDisk)
            return '\n';
 
        // Read one run from disk
        m_pReader->read(reinterpret_cast<char*>(&m_currRun), sizeof(RLUnit));
        ++m_numRunsRead;
    }

    // Decrement the current run and emit its symbol
    m_currRun.decrementCount();
    return m_currRun.getChar();
}
