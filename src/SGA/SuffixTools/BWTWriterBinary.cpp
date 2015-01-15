//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTWriterBinary - Write a run-length encoded BWT to a binary file
//
#include "BWTWriterBinary.h"
#include "SBWT.h"
#include "RLBWT.h"

//
BWTWriterBinary::BWTWriterBinary(const std::string& filename) : m_numRuns(0), m_runFileOffset(0), m_stage(IOS_NONE)
{
    m_pWriter = createWriter(filename, std::ios::out | std::ios::binary);
    m_stage = IOS_HEADER;
}

//
BWTWriterBinary::~BWTWriterBinary()
{
    assert(m_stage == IOS_DONE);
    delete m_pWriter;
}

//
void BWTWriterBinary::writeHeader(const size_t& num_strings, const size_t& num_symbols, const BWFlag& flag)
{
    assert(m_stage == IOS_HEADER);
    m_pWriter->write(reinterpret_cast<const char*>(&RLBWT_FILE_MAGIC), sizeof(RLBWT_FILE_MAGIC));
    m_pWriter->write(reinterpret_cast<const char*>(&num_strings), sizeof(num_strings));
    m_pWriter->write(reinterpret_cast<const char*>(&num_symbols), sizeof(num_symbols));

    // Here we do not know the number of runs that are going to be written to the file
    // so we save the offset in the file and write a dummy value. After the bwt string
    // has been written, we return here and fill in the correct value
    m_runFileOffset = m_pWriter->tellp();
    m_numRuns = 0;
    m_pWriter->write(reinterpret_cast<const char*>(&m_numRuns), sizeof(m_numRuns));

    assert(flag == BWF_NOFMI);
    m_pWriter->write(reinterpret_cast<const char*>(&flag), sizeof(flag));

    m_stage = IOS_BWSTR;    
}

// Write a single character of the BWStr
// If the char is '\n' we are finished
void BWTWriterBinary::writeBWChar(char b)
{
    if(m_currRun.isInitialized())
    {
        if(m_currRun.getChar() == b && !m_currRun.isFull())
        {
            m_currRun.incrementCount();
        }
        else
        {
            // Write out the old run and start a new one
            assert(m_currRun.isInitialized());
            writeRun(m_currRun);
            m_currRun = RLUnit(b);
        }        
    }
    else
    {
        // Start a new run
        m_currRun = RLUnit(b);
    }
}

//
void BWTWriterBinary::writeRun(RLUnit& unit)
{
    //std::cout << "Writing " << unit.getChar() << "," << (int)unit.getCount() << "\n";
    m_pWriter->write(reinterpret_cast<const char*>(&unit.data), sizeof(unit.data));
    ++m_numRuns;
}

// write the final run to the stream and fill in the number of runs
void BWTWriterBinary::finalize()
{
    // JS Bugfix 24/05/12
    // Allow empty BWTs to be written out
    if(m_currRun.isInitialized())
        writeRun(m_currRun);

    m_pWriter->seekp(m_runFileOffset);
    m_pWriter->write(reinterpret_cast<const char*>(&m_numRuns), sizeof(m_numRuns));
    m_pWriter->seekp(std::ios_base::end);
    m_stage = IOS_DONE;
}

