//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SACAInducedCopying - Implementation of the 
// induced copying algorithm as described 
// by Nong, Zhang, Chan (2008)
// Modified by JTS to handle multiple strings
#ifndef SACA_INDUCED_COPYING_H
#include "SuffixArray.h"
#include "ReadTable.h"

void saca_induced_copying(SuffixArray* pSA, const ReadTable* pRT, int numThreads, bool silent = false);

void induceSAl(const ReadTable* pRT, SuffixArray* pSA, char** p_array, int64_t* counts, int64_t* buckets, size_t n, int K, bool end);
void induceSAs(const ReadTable* pRT, SuffixArray* pSA, char** p_array, int64_t* counts, int64_t* buckets, size_t n, int K, bool end);

void countBuckets(const ReadTable* pRT, int64_t* buckets, int K);
void getBuckets(int64_t* counts, int64_t* buckets, int K, bool end);
inline void setBit(char** p_array, size_t str_idx, size_t bit_idx, bool b);
inline bool getBit(char** p_array, size_t str_idx, size_t bit_idx);
void printType(const ReadTable* pRT, char** p_array, size_t str_idx);


#endif
