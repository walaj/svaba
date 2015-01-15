//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTWriterBinary - Read a run length encoded binary BWT file from disk
//
#ifndef BWTREADERBINARY_H
#define BWTREADERBINARY_H

#include "Util.h"
#include "STCommon.h"
#include "Occurrence.h"
#include "EncodedString.h"
#include "BWTReader.h"
#include "RLBWT.h"

class SBWT;
class RLBWT;

class BWTReaderBinary : public IBWTReader
{
    public:
        BWTReaderBinary(const std::string& filename);
        virtual ~BWTReaderBinary();

        //
        virtual void read(RLBWT* pRLBWT);
        virtual void read(SBWT* pSBWT);

        virtual void readHeader(size_t& num_strings, size_t& num_symbols, BWFlag& flag);
        virtual char readBWChar();
        virtual void readRuns(RLVector& out, size_t numRuns);

    private:
        std::istream* m_pReader;
        BWIOStage m_stage;
        RLUnit m_currRun;
        size_t m_numRunsOnDisk;
        size_t m_numRunsRead;
};

#endif
