//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTWriterBinary - Write a run-length encoded BWT to a binary file
//
#ifndef BWTWRITERBINARY_H
#define BWTWRITERBINARY_H

#include "Util.h"
#include "STCommon.h"
#include "Occurrence.h"
#include "BWTWriter.h"
#include "BWTReaderBinary.h"
#include "EncodedString.h"
#include "RLBWT.h"

class SBWT;
class RLBWT;

class BWTWriterBinary : public IBWTWriter
{
    public:
        BWTWriterBinary(const std::string& filename);
        virtual ~BWTWriterBinary();

        // Write an RLBWT file directly from a suffix array and read table
        virtual void writeHeader(const size_t& num_strings, const size_t& num_symbols, const BWFlag& flag);
        virtual void writeBWChar(char b);
        virtual void finalize(); // this method must be called after writing the BW string

    private:

        void writeRun(RLUnit& unit);

        std::ostream* m_pWriter;
        size_t m_numRuns;
        std::streampos m_runFileOffset;
        RLUnit m_currRun;
        BWIOStage m_stage;
};

#endif
