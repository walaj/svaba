//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SequenceWorkItem - Definition of the data structure used in the generic
// functions to process a sequence
//
#ifndef SEQUENCEWORKITEM_H
#define SEQUENCEWORKITEM_H

#include "SeqReader.h"

struct SequenceWorkItem
{
    SequenceWorkItem() : idx(0) {}
    SequenceWorkItem(size_t ri, const SeqRecord& sr) : idx(ri), read(sr) {}
    size_t idx;
    SeqRecord read;
};

struct SequenceWorkItemPair
{
    SequenceWorkItem first;
    SequenceWorkItem second;
};

// Genereic class to generate work items using a seq reader
template<class INPUT>
class WorkItemGenerator
{
    public:
        
        WorkItemGenerator(SeqReader* pReader) : m_pReader(pReader), m_numConsumedLast(0), m_numConsumedTotal(0) {}

        // Template specialization for a SequenceWorkItem
        // Returns false when no more sequences could be consumed from the reader
        bool generate(SequenceWorkItem& out)
        {
            SeqRecord read;
            bool valid = m_pReader->get(read);
            if(valid)
            {
                out.idx = m_numConsumedTotal;
                out.read = read;

                m_numConsumedLast = 1;
                m_numConsumedTotal += 1;
                return true;
            }
            else
            {
                return false;
            }
        }

        // Template specialization for a SequenceWorkItemPair
        bool generate(SequenceWorkItemPair& out)
        {
            SeqRecord read1;
            SeqRecord read2;

            bool valid1 = m_pReader->get(read1);
            if(valid1)
            {
                bool valid2 = m_pReader->get(read2);
                assert(valid2);

                out.first.idx = m_numConsumedTotal;
                out.second.idx = m_numConsumedTotal + 1;
                out.first.read = read1;
                out.second.read = read2;

                m_numConsumedLast = 2;
                m_numConsumedTotal += 2;
                return true;
            }
            else
            {
                return false;
            }
        }

        inline size_t getConsumedLast() const { return m_numConsumedLast; }
        inline size_t getNumConsumed() const { return m_numConsumedTotal; }

    private:

        SeqReader* m_pReader;
        size_t m_numConsumedLast;
        size_t m_numConsumedTotal;
};

//Jeremiah
// Genereic class to generate work items using reads in memory
template<class INPUT>
class WorkItemMemGenerator
{
    public:
        
        WorkItemMemGenerator(ReadTable* pRT) : m_pRT(pRT), m_numConsumedLast(0), m_numConsumedTotal(0) {}

        // Template specialization for a SequenceWorkItem
        // Returns false when no more sequences could be consumed from the reader
        bool generate(SequenceWorkItem& out)
        {
	    //SeqRecord read;
	    SeqItem si;
            bool valid = m_pRT->getRead(si);            
	    SeqRecord read;
            read.id = si.id;
            read.seq = si.seq;

            if(valid)
            {
                out.idx = m_numConsumedTotal;
                out.read = read;

                m_numConsumedLast = 1;
                m_numConsumedTotal += 1;
                return true;
            }
            else
            {
                return false;
            }
        }

        inline size_t getConsumedLast() const { return m_numConsumedLast; }
        inline size_t getNumConsumed() const { return m_numConsumedTotal; }

    private:

        ReadTable* m_pRT;
        size_t m_numConsumedLast;
        size_t m_numConsumedTotal;
};


#endif
