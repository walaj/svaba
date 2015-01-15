//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTWriterAsciiAscii - Write a BWT to a plain text file
//
#ifndef BWTWRITERASCII_H
#define BWTWRITERASCII_H

#include "Util.h"
#include "STCommon.h"
#include "Occurrence.h"
#include "BWTWriter.h"
#include "EncodedString.h"
#include "SuffixArray.h"

class SBWT;
class RLBWT;

class BWTWriterAscii : public IBWTWriter
{
    public:
        BWTWriterAscii(const std::string& filename);
        virtual ~BWTWriterAscii();

        //
        virtual void write(const RLBWT* pRLBWT);
        virtual void writeHeader(const size_t& num_strings, const size_t& num_symbols, const BWFlag& flag);
        virtual void writeBWChar(char b);
        virtual void finalize();

        void write(const SBWT* pBWT);
        void writeBWStr(const BWTString& str); 
        void writePred(const AlphaCount64& pc);
        void writeOccurrence(const Occurrence& icc);

    private:
        std::ostream* m_pWriter;
        BWIOStage m_stage;
};

#endif
