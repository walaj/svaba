//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// SAReader.h - Read a suffix array file from disk
//
#ifndef SAREADER_H
#define SAREADER_H

#include "Util.h"
#include "STCommon.h"
#include "Occurrence.h"

const uint16_t SA_FILE_MAGIC = 0xCACA;

enum SAIOStage
{
    SAIOS_NONE,
    SAIOS_HEADER,
    SAIOS_ELEM,
    SAIOS_DONE
};

class SuffixArray;

class SAReader
{
    public:
        SAReader(const std::string& filename);
        ~SAReader();

        //
        void read(SuffixArray* pSA);

        // Read the header of the file
        void readHeader(size_t& num_strings, size_t& num_elems);
        
        // Read the file into a vector of SAElems
        void readElems(SAElemVector& elemVector);

        // Read the file into a vector of unsigned ints storing
        // the read indices. This is a more compact representation
        // than storing the full SAElems but only works up to 2**32 values
        void readElems(std::vector<uint32_t>& outVector);

        // Read a single element
        SAElem readElem();

    private:
        std::istream* m_pReader;
        SAIOStage m_stage;
};

#endif
