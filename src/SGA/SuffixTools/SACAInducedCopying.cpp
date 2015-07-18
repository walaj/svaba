//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SACAInducedCopying algorithm
//
#include "SACAInducedCopying.h"
#include "SuffixCompare.h"
#include "mkqs.h"
#include "bucketSort.h"
#include "Util.h"

unsigned char mask[]={0x80,0x40,0x20,0x10,0x08,0x04,0x02,0x01};

#define GET_CHAR(i, j) pRT->getChar((i),(j))
#define isLMS(i, j) ((j) > 0 && getBit(type_array, (i), (j)) && !getBit(type_array, (i), (j-1)))
#define GET_BKT(c) getBaseRank((c))

// Implementation of induced copying algorithm by
// Nong, Zhang, Chan
// Follows implementation given as an appendix to their 2008 paper
// '\0' is the sentinenl in this algorithm
void saca_induced_copying(SuffixArray* pSA, const ReadTable* pRT, int numThreads, bool silent)
{

    // In the multiple strings case, we need a 2D bit array
    // to hold the L/S types for the suffixes
    size_t num_strings = pRT->getCount();
    char** type_array = new char*[num_strings];

    for(size_t i = 0; i < num_strings; ++i)
    {
        size_t s_len = pRT->getReadLength(i) + 1;
        size_t num_bytes = (s_len / 8) + 1;
        type_array[i] = new char[num_bytes];
        assert(type_array[i] != 0);
        memset(type_array[i], 0, num_bytes);
    }

    // Classify each suffix as being L or S type
    for(size_t i = 0; i < num_strings; ++i)
    {
        size_t s_len = pRT->getReadLength(i) + 1;

        // The empty suffix ($) for each string is defined to be S type
        // and hence the next suffix must be L type
        setBit(type_array, i, s_len - 1, 1);
        setBit(type_array, i, s_len - 2, 0);
        for(int64_t j = s_len - 3; j >= 0; --j)
        {
            char curr_c = GET_CHAR(i, j);
            char next_c = GET_CHAR(i, j + 1);

            bool s_type = (curr_c < next_c || (curr_c == next_c && getBit(type_array, i, j + 1) == 1));
            setBit(type_array, i, j, s_type);
        }
    }

    // setup buckets
    const int ALPHABET_SIZE = 5;
    int64_t bucket_counts[ALPHABET_SIZE];
    int64_t buckets[ALPHABET_SIZE];

    // find the ends of the buckets
    countBuckets(pRT, bucket_counts, ALPHABET_SIZE);
    getBuckets(bucket_counts, buckets, ALPHABET_SIZE, true); 

    //std::cout << "initializing SA\n";

    // Initialize the suffix array
    size_t num_suffixes = buckets[ALPHABET_SIZE - 1];
    pSA->initialize(num_suffixes, pRT->getCount());

    // Copy all the LMS substrings into the first n1 places in the SA
    size_t n1 = 0;
    for(size_t i = 0; i < num_strings; ++i)
    {
        size_t s_len = pRT->getReadLength(i) + 1;
        for(size_t j = 0; j < s_len; ++j)
        {
            if(isLMS(i,j))
                pSA->set(n1++, SAElem(i, j));
        }
    }

    /*
    //induceSAl(pRT, pSA, type_array, bucket_counts, buckets, num_suffixes, ALPHABET_SIZE, false);
    //induceSAs(pRT, pSA, type_array, bucket_counts, buckets, num_suffixes, ALPHABET_SIZE, true);
    
    // Compact all the sorted substrings into the first portion of the SA
    size_t n1 = 0;
    for(size_t i = 0; i < num_suffixes; ++i)
    {
        SAElem elem = pSA->get(i);
        if(!elem.isEmpty() && isLMS(elem.getID(), elem.getPos()))
        {
            pSA->set(n1++, elem);
        }
    }
    */

    // Call MKQS, first on the sequence and then on the index in the read table
    SuffixCompareRadix radix_compare(pRT, 6);
    SuffixCompareIndex index_compare;
    //SuffixCompareID id_compare(pRT);

    if(numThreads <= 1)
        mkqs2(&pSA->m_data[0], n1, 0, radix_compare, index_compare);
    //else
    //    parallel_mkqs(&pSA->m_data[0], n1, numThreads, radix_compare, index_compare);
    
    //if(!silent)
    //    std::cout << "[saca] mkqs finished\n";

    // Induction sort the remaining suffixes
    for(size_t i = n1; i < num_suffixes; ++i)
        pSA->set(i, SAElem());
    
    // Find the ends of the buckets
    getBuckets(bucket_counts, buckets, ALPHABET_SIZE, true);

    for(int64_t i = n1 - 1; i >= 0; --i)
    {
        SAElem elem_i = pSA->get(i);
        pSA->set(i, SAElem()); // empty
        char c = GET_CHAR(elem_i.getID(), elem_i.getPos());
        pSA->set(--buckets[GET_BKT(c)], elem_i);
    }

    induceSAl(pRT, pSA, type_array, bucket_counts, buckets, num_suffixes, ALPHABET_SIZE, false);
    induceSAs(pRT, pSA, type_array, bucket_counts, buckets, num_suffixes, ALPHABET_SIZE, true);

    // deallocate t array
    for(size_t i = 0; i < num_strings; ++i)
    {
        delete [] type_array[i];
    }
    delete [] type_array;
}

void induceSAl(const ReadTable* pRT, SuffixArray* pSA, char** p_array, int64_t* counts, int64_t* buckets, size_t n, int K, bool end)
{
    getBuckets(counts, buckets, K, end);
    for(size_t i = 0; i < n; ++i)
    {
        const SAElem& elem_i = pSA->get(i);
        if(!elem_i.isEmpty() && elem_i.getPos() > 0)
        {
            //std::cout << "Curr: " << elem_i << "\n";
            SAElem elem_j(elem_i.getID(), elem_i.getPos() - 1);
            if(!getBit(p_array, elem_j.getID(), elem_j.getPos()))
            {
                char c = GET_CHAR(elem_j.getID(),elem_j.getPos());
                //std::cout << "<iSA1>Placing " << elem_j << " at position " << buckets[GET_BKT(c)] << "\n";
                pSA->set(buckets[GET_BKT(c)]++, elem_j);
            }
        }
    }
}

void induceSAs(const ReadTable* pRT, SuffixArray* pSA, char** p_array, int64_t* counts, int64_t* buckets, size_t n, int K, bool end)
{
    getBuckets(counts, buckets, K, end);
    for(int64_t i = n - 1; i >= 0; --i)
    {
        const SAElem& elem_i = pSA->get(i);
        if(!elem_i.isEmpty() && elem_i.getPos() > 0)
        {
            //std::cout << "<isas>Curr: " << elem_i << "\n";
            SAElem elem_j(elem_i.getID(), elem_i.getPos() - 1);
            if(getBit(p_array, elem_j.getID(), elem_j.getPos()))
            {
                char c = GET_CHAR(elem_j.getID(),elem_j.getPos());
                //std::cout << "<iSAs>Placing " << elem_j << " at position " << buckets[GET_BKT(c)] - 1 << "\n";
                pSA->set(--buckets[GET_BKT(c)], elem_j);
            }
        }
    }
}


// Calculate the number of items that should be in each bucket
void countBuckets(const ReadTable* pRT, int64_t* counts, int K)
{
    for(int i = 0; i < K; ++i)
        counts[i] = 0;

    for(size_t i = 0; i < pRT->getCount(); ++i)
    {
        size_t s_len = pRT->getReadLength(i);
        for(size_t j = 0; j < s_len; ++j)
            counts[getBaseRank(GET_CHAR(i,j))]++;

        counts[getBaseRank('\0')]++;
    }
}

// If end is true, calculate the end of the buckets, otherwise 
// calculate the starts
void getBuckets(int64_t* counts, int64_t* buckets, int K, bool end)
{
    for(int i = 0; i < K; ++i)
        buckets[i] = 0;

    size_t sum = 0; 
    for(int i = 0; i < K; ++i)
    {
        sum += counts[i];
        buckets[i] = end ? sum : sum - counts[i];
    }
}

// set the element to b
void setBit(char** p_array, size_t str_idx, size_t bit_idx, bool b)
{
    char* ba = p_array[str_idx];
    ba[bit_idx / 8] = b ? mask[bit_idx % 8] | ba[bit_idx / 8] : ~mask[bit_idx % 8] & ba[bit_idx / 8];
}

bool getBit(char** p_array, size_t str_idx, size_t bit_idx)
{
    return p_array[str_idx][bit_idx / 8] & mask[bit_idx % 8] ? 1 : 0;
}

void printType(const ReadTable* pRT, char** p_array, size_t str_idx)
{
    std::string suf_string = pRT->getRead(str_idx).seq.getSuffixString(0);
    std::cout << suf_string << "\n";
    for(size_t i = 0; i < suf_string.length(); ++i)
    {
        std::cout << (getBit(p_array, str_idx, i) == 1 ? "S" : "L");
    }
    std::cout << "\n";
}

