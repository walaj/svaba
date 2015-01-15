//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ReadInfoTable - A 0-indexed table of ID, length pairs
// Used to convert suffix array hits to overlaps
//
#ifndef READINFOTABLE_H
#define READINFOTABLE_H
#include "Util.h"
#include "SeqReader.h"
#include "ReadTable.h"
#include <map>

enum ReadInfoOption
{
    RIO_NONE,
    RIO_NUMERICID
};

struct ReadInfo
{
    ReadInfo() : length(0) {}
    ReadInfo(const std::string& i, const uint32_t& l) : id(i), length(l) {}
    std::string id;
    uint32_t length;
};

typedef std::vector<ReadInfo> InfoVector;

class ReadInfoTable
{
    public:
        //
        ReadInfoTable() {}
        ReadInfoTable(ReadTable *pRT);

        // Load the table using the read in filename
        // If num_expected > 0, reserve room in the table for num_expected reads
        ReadInfoTable(std::string filename, size_t num_expected = 0, ReadInfoOption options = RIO_NONE);
        ~ReadInfoTable();

        //
        const ReadInfo getReadInfo(size_t idx) const;
        std::string getReadID(size_t idx) const;
        size_t getReadLength(size_t idx) const;
        size_t getCount() const;
        size_t countSumLengths() const;
        void clear();

    private:

        std::vector<int> m_lengths;
        std::vector<std::string> m_ids;

        bool m_numericIDs;
};

#endif
