//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// ReadTable - A 0-indexed table of reads
//
#ifndef SNOWS_READTABLE_H
#define SNOWS_READTABLE_H

#include "Util.h"
#include "SeqReader.h"
#include "SeqLib/UnalignedSequence.h"
#include <map>

typedef std::vector<SeqItem> ReadVector;
typedef std::map<std::string, SeqItem*> ReadIndex;

class ReadTable
{
 public:
  //
  ReadTable() : idx(0), m_pIndex(NULL) {}
    //ReadTable(std::string* seq, std::string * id, int length);
    ReadTable(SeqRecordVector &srv);
    ReadTable(const SeqLib::UnalignedSequenceVector &contigs);
    ReadTable(std::string filename, uint32_t reader_flags = 0);
    ~ReadTable();

    // Initialize this read table as the reverse of the passed in read table
    void initializeReverse(const ReadTable* pRT);

    // Reverse all the reads in this table
    void reverseAll();

    // Build a readid -> read index
    void indexReadsByID();

    //
    void addRead(const SeqItem& r);
    const SeqItem& getRead(size_t idx) const;
    const SeqItem& getRead(const std::string& id) const;
    bool getRead(SeqItem &si); // JEREMIAH
    size_t getReadLength(size_t idx) const;
    size_t getCount() const;
    size_t countSumLengths() const;
    std::string getReadID(size_t idx) const;
    void setZero();
    void clear();

    // Get a particular character for a particular read
    inline char getChar(size_t str_idx, size_t char_idx) const
      {
	assert(str_idx < m_table.size());
	return m_table[str_idx].seq.get(char_idx);
      }

    // I/O
    friend std::ostream& operator<<(std::ostream& out, const ReadTable& rt);

    //Jeremiah
    void setReadSequence(const std::string& seqr);

 private:
    ReadVector m_table;
    unsigned int idx = 0;
    // Index of readid -> SeqItem
    // It is not build be default to save memory
    ReadIndex* m_pIndex; 
};

#endif

