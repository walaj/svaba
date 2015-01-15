//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTReaderAsciiAscii.h - Read a plain-text BWT file from disk
//
#ifndef BWTREADERASCII_H
#define BWTREADERASCII_H

#include "Util.h"
#include "STCommon.h"
#include "Occurrence.h"
#include "EncodedString.h"
#include "BWTReader.h"

class SBWT;
class RLBWT;

class BWTReaderAscii : public IBWTReader
{
    public:
        BWTReaderAscii(const std::string& filename);
        ~BWTReaderAscii();

        //
        void read(SBWT* pBWT);
        void read(RLBWT* pRLBWT);
        void readHeader(size_t& num_strings, size_t& num_symbols, BWFlag& flag);
        void readBWStr(BWTString& out_str);
        char readBWChar();
        void readPred(AlphaCount64& out_pc);
        void readOccurrence(Occurrence& out_icc);

    private:
        std::istream* m_pReader;
        BWIOStage m_stage;
};

#endif
