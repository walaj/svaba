//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// QualityVector - Vector of doubles representing
// log-transformed quality values for a sequence
//
#ifndef QUALITYVECTOR_H
#define QUALITYVECTOR_H

#include <vector>
#include "Alphabet.h"
#include "DNADouble.h"

typedef std::vector<DNADouble> APVec;

class QualityVector
{
    public:
        
        // Constructors
        QualityVector();
        QualityVector(const QualityVector& vec, int start, int size);

        //
        void add(DNADouble v);
        void set(size_t idx, DNADouble v);
        DNADouble get(size_t idx) const;
        size_t size() const;
        bool empty() const;

        //
        void reverseComplement();

    private:
        APVec m_data;

};

#endif
