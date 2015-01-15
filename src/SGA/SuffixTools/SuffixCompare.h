//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SuffixCompare - Suffix Comparator classes
// Used for histogram/bucket sort, mkqs and std::sort
//
#ifndef SUFFIXCOMPARE_H
#define SUFFIXCOMPARE_H
#include "STCommon.h"
#include "ReadTable.h"

// Suffix comparator object for radix-like sorts (bucket/histogram and MKQS)
class SuffixCompareRadix
{
    public:
        SuffixCompareRadix(const ReadTable* pRT);
        SuffixCompareRadix(const ReadTable* pRT, int bucket_len);
        ~SuffixCompareRadix();
    
        //
        void initializeNumSuffixLUT();

        // Bucket function
        int getBucket(SAElem x) const;

        // Get the character at position d for the SAElem
        inline char getChar(SAElem& x, int d) const
        {
            return m_pRT->getChar(x.getID(), x.getPos() + d);
        }

        // Calculate the number of possible suffixes
        int calcNumSuffixes(int maxLen) const;

        // Get the number of possible suffixes of length maxLen using the lookup table
        inline int getNumSuffixes(int maxLen) const
        {
            return m_pNumSuffixLUT[maxLen];
        }

        // Get the suffix character string corresponding to this element
        inline const char* getChrPtr(SAElem& x) const
        {
            return m_pRT->getRead(x.getID()).seq.getSuffix(x.getPos());
        }

        // Calculate the number of suffixes that precede the first instance of b for a 
        // given maximum suffix length
        inline int numPredSuffixes(char b, int maxLen) const
        {
            // base case
            int rb = getBaseRank(b);
            if(rb == 0)
                return 0;
            int block_size = getNumSuffixes(maxLen - 1);
            return block_size * (rb - 1) + 1;
        }

        // Return the number of buckets needed for bucket sort
        int getNumBuckets() const;

        // Return the bucket offset
        size_t getBucketOffset() const { return m_bucketOffset; }

        // Get the bucket length
        size_t getBucketLen() const { return m_bucketLen; }

        // Set the offset
        void setBucketDepth(int depth);

        // Return true if a bucket is degenerate ie it cannot be subdivided any further
        bool isBucketDegenerate(int index) const;


        // Print the element
        void printElem(SAElem& x) const;

    private:

        // Disallow default and copy constructor
        SuffixCompareRadix() {}
        SuffixCompareRadix(const SuffixCompareRadix& /*other*/) { assert(false); }

        const ReadTable* m_pRT;
        size_t m_bucketOffset;
        size_t m_bucketLen;
        int* m_pNumSuffixLUT;
};

// Compare two suffixes by their ID
// This is used for the final pass, after suffixes has been compared by sequence
class SuffixCompareID
{
    public:
        SuffixCompareID(const ReadTable* pRT) : m_pRT(pRT) {}

        // Comparator function
        bool operator()(SAElem x, SAElem y) const;

    private:
        
        // default is not accessible
        SuffixCompareID() : m_pRT(0) {}

        const ReadTable* m_pRT;
};

// Compare two suffixes by their index in the read table
// This is used for the final pass, after suffixes has been compared by sequence
class SuffixCompareIndex
{
    public:
        SuffixCompareIndex() {}

        // Comparator function
        inline bool operator()(SAElem x, SAElem y) const
        {
            return x.getID() < y.getID();
        }

    private:
        
};


#endif
