//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// SBWT - Uncompressed burrows Wheeler transform 
//        of a set of sequence reads
//          
//
#ifndef SBWT_H
#define SBWT_H

#include "STCommon.h"
#include "Occurrence.h"
#include "SuffixArray.h"
#include "ReadTable.h"
#include "HitData.h"
#include "BWTReader.h"
#include "EncodedString.h"

//
// BWT
//
class SBWT
{
    public:
    
        // Constructors
        SBWT(const std::string& filename, int sampleRate = DEFAULT_SAMPLE_RATE);
        SBWT(const SuffixArray* pSA, const ReadTable* pRT);
        
        //    
        void initializeFMIndex(int sampleRate);

        // Exact match
        void backwardSearch(std::string w) const;

        // L[i] -> F mapping 
        size_t LF(size_t idx) const;

        inline char getChar(size_t idx) const { return m_bwStr.get(idx); }
        inline BaseCount getPC(char b) const { return m_predCount.get(b); }

        // Return the number of times char b appears in bwt[0, idx]
        inline BaseCount getOcc(char b, size_t idx) const { return m_occurrence.get(m_bwStr, b, idx); }

        // Return the number of times each symbol in the alphabet appears in bwt[0, idx]
        inline AlphaCount64 getFullOcc(size_t idx) const { return m_occurrence.get(m_bwStr, idx); }

        // Return the number of times each symbol in the alphabet appears ins bwt[idx0, idx1]
        inline AlphaCount64 getOccDiff(size_t idx0, size_t idx1) const { return m_occurrence.getDiff(m_bwStr, idx0, idx1); }

        inline size_t getNumStrings() const { return m_numStrings; } 
        inline size_t getBWLen() const { return m_bwStr.length(); }

        // Return the first letter of the suffix starting at idx
        inline char getF(size_t idx) const
        {
            size_t ci = 0;
            while(ci < ALPHABET_SIZE && m_predCount.getByIdx(ci) <= idx)
                ci++;
            assert(ci != 0);
            return RANK_ALPHABET[ci - 1];
        }

        // Print the size of the BWT
        void printInfo() const;
        void print(const ReadTable* pRT, const SuffixArray* pSA) const;
        void printRunLengths() const { std::cout << "Using SimpleBWT - No run lengths\n"; }
        void validate() const;

        // IO
        friend class BWTReaderBinary;
        friend class BWTWriterBinary;
        friend class BWTReaderAscii;
        friend class BWTWriterAscii;

        static const int DEFAULT_SAMPLE_RATE = 128;
        static const int DEFAULT_SAMPLE_RATE_SMALL = DEFAULT_SAMPLE_RATE;

    private:

        
        // Default constructor is not allowed
        SBWT() {}

        // The O(a,i) array
        Occurrence m_occurrence;

        // The C(a) array
        AlphaCount64 m_predCount;
        
        // The bw string
        BWTString m_bwStr;

        // The number of strings in the collection
        size_t m_numStrings;
};
#endif
