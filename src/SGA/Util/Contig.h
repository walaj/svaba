#ifndef CONTIG_H
#define CONTIG_H

#include "Util.h"
#include <iostream>
#include <ostream>
#include <sstream>
#include <iterator>


//
// Enumerations
//
enum UniqueFlag
{
    UF_UNKNOWN = 0,
    UF_UNIQUE,
    UF_REPEAT,
    UF_NOCALL
};

class Contig
{
    public:

        //
        // Constructors
        //
        Contig();

        //
        // Getters
        //
        ContigID getID() const { return m_id; }
        double getCoverage() const { return m_coverage; }
        UniqueFlag getUniqueFlag() const { return m_uniqueFlag; }
        size_t getLength() const { return m_length; }
        Sequence getSequence() const { return m_seq; }
        double getNormalizedCoverage() const { return m_coverage / m_length; }
        bool isUnique() const { return m_uniqueFlag == UF_UNIQUE; }
        unsigned int getAnnotation() { return m_annotation; }
        SequenceVector getVariants() const;

        //
        // Setters
        //
        void setUniqueFlag(UniqueFlag uf) { m_uniqueFlag = uf; }
        void setFromKeyValue(std::string& key, std::string& value);
        void addVariants(const SequenceVector& seqs);

        //
        // Operators
        //
        bool operator()(Contig& c1, Contig& c2) const
        {
            return c1.m_seq < c2.m_seq;
        }

        //
        // Readers
        //
        friend std::istream& readFasta(std::istream& in, Contig& c);
        friend std::istream& readCAF(std::istream& in, Contig& c);
        
        //
        // Writers
        //
        friend std::ostream& writeCAF(std::ostream& out, Contig& c);

    private:
        ContigID m_id;
        size_t m_length;
        double m_coverage;
        UniqueFlag m_uniqueFlag;
        Sequence m_seq;
        unsigned int m_annotation;
        SequenceVector m_variants;
};

#endif
