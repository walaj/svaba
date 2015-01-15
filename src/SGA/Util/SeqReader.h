//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SeqReader - Reads fasta or fastq sequence files
//
#ifndef SEQREADER_H
#define SEQREADER_H

#include <fstream>
#include "Util.h"

enum RecordType
{
    RT_FASTA,
    RT_FASTQ,
    RT_UNKNOWN
};

static const uint32_t SRF_NO_VALIDATION = 1;
static const uint32_t SRF_KEEP_CASE = 2;

//
class SeqReader
{
    public:
        SeqReader(std::string filename, uint32_t flags = 0);
        ~SeqReader();
        bool get(SeqRecord& sr);

    private:
        std::istream* m_pHandle;
        uint32_t m_flags;
};

#endif
