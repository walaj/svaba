//-----------------------------------------------------
// Copyright 2014 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
// Released under the GPL
//-----------------------------------------------------
//
// Data structure for performing proximity queries
// against a set of variants
//
#ifndef VARIANTINDEX_H
#define VARIANTINDEX_H
#include "Util.h"
#include "SeqReader.h"
#include "VCFUtil.h"
#include "ReadTable.h"
#include <map>

// To save space we store a vcf-like record
// with only the required fields
struct VariantRecord
{
    std::string reference;
    std::string ref_sequence;
    std::string alt_sequence;
    size_t position;
};

typedef std::vector<VariantRecord> VariantRecordVector;
typedef std::vector<int> IntVector;
typedef std::vector<IntVector> IntVectorVector;

// map from chromosome to a vector of buckets, one every 10kbp
typedef std::map<std::string, IntVectorVector> VariantIndexMap;

class VariantIndex
{
    public:
        //
        VariantIndex(const std::string& filename, const ReadTable& refTable);

        // Return a vector of the variants that are close to the input variant
        VariantRecordVector getNearVariants(const std::string& reference, 
                                            int position, 
                                            int distance = 30) const;

    private:
        
        //
        void buildIndex(const ReadTable& refTable);

        //
        VariantRecordVector m_records;
        VariantIndexMap m_map;

        // 
        static const int BUCKET_SIZE = 10000;
};

#endif
