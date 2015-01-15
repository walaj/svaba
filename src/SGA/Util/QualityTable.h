//-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// QualityTable - A 0-indexed table of quality scores
//
#ifndef QUALITYTABLE_H
#define QUALITYTABLE_H
#include "Util.h"
#include "SeqReader.h"
#include "QualityCodec.h"
#include <map>

struct QualityString
{
    uint8_t num_encoded_symbols;
    QualityStorageUnit* encoded_data;
};

typedef std::vector<QualityString> QualityStringVector;

class QualityTable
{
    public:

        //
        QualityTable();
        ~QualityTable();

        //
        void loadQualities(const std::string& filename);
        void addQualityString(const std::string& qual);
        std::string getQualityString(size_t idx, size_t n) const;
        size_t getCount() const;
        void clear();

        //
        void printSize() const;

    private:
        QualityCodec<4> m_codec;
        QualityStringVector m_table;
        size_t m_bytes_used;
        char m_missingQualityChar;
};

#endif
