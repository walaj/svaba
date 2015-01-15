//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTWriter - Abstract class for writing a BWT file to disk
//
#ifndef BWTWRITER_H
#define BWTWRITER_H

#include "Util.h"
#include "STCommon.h"
#include "Occurrence.h"
#include "BWTReader.h"
#include "EncodedString.h"
#include "SuffixArray.h"

class IBWTWriter
{
    public:
        IBWTWriter() {}
        IBWTWriter(const std::string& /*filename*/) {}
        virtual ~IBWTWriter() {}

        //
        void write(const SuffixArray* pSA, const ReadTable* pRT);
        virtual void writeHeader(const size_t& num_strings, const size_t& num_symbols, const BWFlag& flag) = 0;
        virtual void writeBWChar(char b) = 0;
        virtual void finalize() = 0;
};

namespace BWTWriter
{
    IBWTWriter* createWriter(const std::string& filename);
}

#endif
