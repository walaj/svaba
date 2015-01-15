//-----------------------------------------------------
// Copyright 2014 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
// Released under the GPL
//-----------------------------------------------------
//
// Data structure for performing proximity queries
// against a set of variants
//
#include <iostream>
#include <algorithm>
#include <assert.h>
#include "VariantIndex.h"

VariantIndex::VariantIndex(const std::string& filename, const ReadTable& refTable)
{
    std::ifstream input(filename.c_str());
    std::string line;

    // Parse the VCF file
    while(getline(input, line))
    {
        if(line.empty())
            continue;

        if(line[0] == '#')
            continue;
        
        VCFRecord record(line);

        // Do not allow multi-allelic records
        if(record.isMultiAllelic())
        {
            std::cerr << "Error: multi-allelic Variant found, please run vcfbreakmulti\n";
            exit(EXIT_FAILURE);
        }

        // Convert to a minimal representation of the change
        VariantRecord mvf = { record.refName, record.refStr, record.varStr, record.refPosition };
        m_records.push_back(mvf);
    }

    buildIndex(refTable);
}

void VariantIndex::buildIndex(const ReadTable& refTable)
{
    for(size_t i = 0; i < refTable.getCount(); ++i)
    {
        const SeqItem& seq_record = refTable.getRead(i);
        std::string chromosome = seq_record.id;
        size_t length = seq_record.seq.length();
        int buckets = (length / BUCKET_SIZE) + 1;

        m_map[chromosome].resize(buckets);
    }

    for(size_t i = 0; i < m_records.size(); ++i)
    {
        VariantIndexMap::iterator iter = m_map.find(m_records[i].reference);
        assert(iter != m_map.end());

        int bucket_id = m_records[i].position / BUCKET_SIZE;
        assert(bucket_id < (int)iter->second.size());
        IntVector& bucket = iter->second[bucket_id];
        bucket.push_back(i);
    }
}

VariantRecordVector VariantIndex::getNearVariants(const std::string& reference,
                                                  int position,
                                                  int distance) const
{
    VariantRecordVector out;

    VariantIndexMap::const_iterator iter = m_map.find(reference);
    assert(iter != m_map.end());
    const IntVectorVector& chr_buckets = iter->second;

    int lower_index = std::max(position - distance, 0) / BUCKET_SIZE;
    int upper_index = std::min((position + distance) / BUCKET_SIZE, (int)chr_buckets.size() - 1);

    for(int bi = lower_index; bi <= upper_index; ++bi)
    {
        const IntVector& bucket = chr_buckets[bi];
        for(size_t vi = 0; vi < bucket.size(); ++vi)
        {
            const VariantRecord& record = m_records[bucket[vi]];
            if(abs(record.position - position) < distance)
                out.push_back(record);
        }
    }

    return out;
}

