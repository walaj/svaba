//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SuffixArray - Generalized suffix array
//
#ifndef SUFFIXARRAY_H
#define SUFFIXARRAY_H
#include "STCommon.h"
#include "ReadTable.h"
#include "Match.h"

class LCPArray;

class SuffixArray
{
    public:
        
        //
        SuffixArray() {}
        SuffixArray(const std::string& filename);
        SuffixArray(const ReadTable* pRT, int numThreads, bool silent = false);

        // Construction/Validation functions
        void initialize(const ReadTable& rt);
        void initialize(size_t num_suffixes, size_t num_strings);
        void validate(const ReadTable* pRT) const;
        void sort(const ReadTable* pRT);

        // Remove all the suffixes from the SA that have an id in idSet
        void removeReads(const NumericIDSet& idSet);

        // Simple accessors
        inline const SAElem& get(size_t idx) const { assert(idx < m_data.size()); return m_data[idx]; }
        inline void set(size_t idx, SAElem e) { m_data[idx] = e; }
        size_t getSize() const { return m_data.size(); }
        size_t getNumStrings() const { return m_numStrings; } 
        std::string getSuffix(size_t idx, const ReadTable* pRT) const;
        size_t getSuffixLength(const ReadTable* pRT, const SAElem elem) const;
        

        // Operators
        friend std::ostream& operator<<(std::ostream& out, const SuffixArray& sa);
        friend std::istream& operator>>(std::istream& in, SuffixArray& sa);

        // Output the entire suffix array
        void write(const std::string& filename);

        // Write the BWT directly to disk
        void writeBWT(const std::string& filename, const ReadTable* pRT);

        // Output the suffix array index
        // The suffix array index are the full-length suffixes (the entire string)
        // of the suffix array, sorted into lexographical order
        // The file format is the same as the suffix array
        // This is used during the BWT indexing and is the only part of the
        // suffix array that we need
        void writeIndex(std::string& filename);
        void writeIndex(); //JEREMIAH

        // Print funcs
        void print() const;
        void print(const ReadTable* pRT) const;

        // friends
        friend void saca_induced_copying(SuffixArray* pSA, const ReadTable* pRT, int numThreads, bool silent);
        friend class SAReader;
        friend class SAWriter;

    private:
        
        // Data members
        SAElemVector m_data;
        size_t m_numStrings;
};

#endif 
