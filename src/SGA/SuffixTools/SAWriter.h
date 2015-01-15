//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// SAWriter.h - Read a suffix array file from disk
//
#ifndef SAWRITER_H
#define SAWRITER_H

#include "Util.h"
#include "STCommon.h"
#include "Occurrence.h"
#include "SAReader.h"

class SuffixArray;

class SAWriter
{
    public:
        SAWriter(const std::string& filename);
        ~SAWriter();

        //
        void write(const SuffixArray* pSA);
        void writeHeader(const size_t& num_strings, const size_t& num_elems);
        void writeElems(const SAElemVector& elemVector);
        void writeElem(const SAElem& elem);

    private:
        std::ostream* m_pWriter;
        SAIOStage m_stage;
};

#endif
