//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTReader - Abstract class for reading a BWT file
//
#ifndef BWTREADER_H
#define BWTREADER_H

#include "Util.h"
#include "STCommon.h"
#include "Occurrence.h"
#include "EncodedString.h"

enum BWIOStage
{
    IOS_NONE,
    IOS_HEADER,
    IOS_BWSTR,
    IOS_PC,
    IOS_OCC,
    IOS_DONE
};

enum BWFlag
{
    BWF_NOFMI = 0,
    BWF_HASFMI
};

const uint16_t RLBWT_FILE_MAGIC = 0xCACA;
const uint16_t BWT_FILE_MAGIC = 0xEFEF;

class RLBWT;

class IBWTReader
{
    public:
        IBWTReader() {}
        IBWTReader(const std::string& /*filename*/) {}
        virtual ~IBWTReader() {}

        //
        virtual void read(RLBWT* pRLBWT) = 0;
        virtual void readHeader(size_t& num_strings, size_t& num_symbols, BWFlag& flag) = 0; 
        virtual char readBWChar() = 0;
};

namespace BWTReader
{
    IBWTReader* createReader(const std::string& filename);
};

#endif
