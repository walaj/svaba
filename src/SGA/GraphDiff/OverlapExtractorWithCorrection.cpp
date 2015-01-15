///----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapExtractorWithCorrection - Extract reads from an 
// FM-index that share a kmer with a query sequence.
// Previously used kmers are cached to avoid duplicate
// lookups.
//
#include "OverlapExtractorWithCorrection.h"
#include "HapgenUtil.h"
#include "Profiler.h"

//
OverlapExtractorWithCorrection::OverlapExtractorWithCorrection(size_t k, 
                                                               const BWTIndexSet& index_set, 
                                                               ErrorCorrectProcess* corrector, 
                                                               int minOverlap, 
                                                               double minIdentity) : m_k(k),
                                                                                     m_index_set(index_set),
                                                                                     m_corrector(corrector),
                                                                                     m_minOverlap(minOverlap),
                                                                                     m_minIdentity(minIdentity)
{

}

// 
SequenceOverlapPairVector OverlapExtractorWithCorrection::getExactOverlaps(const std::string& query)
{
    PROFILE_FUNC("OverlapExtractorWithCorrection::queryOverlaps")

    SequenceOverlapPairVector raw_overlaps;

    // 
    // Get inexact overlaps between the query and uncorrected reads by direct FM-index lookups
    //
    getRawOverlapsDirect(query, &raw_overlaps);

    /*
    //
    // Get inexact overlaps between the query and uncorrected reads using the k-mer cache
    //
    // Extract new reads from the FM-index and update the cache
    extractAndUpdate(query);
    extractAndUpdate(reverseComplement(query));

    // Compute overlaps between the query sequence and the raw extracted reads
    getRawOverlapsCached(query, false, &raw_overlaps);
    getRawOverlapsCached(query, true, &raw_overlaps);
    */
    // Convert raw overlaps to corrected, exact overlaps
    SequenceOverlapPairVector out_overlaps;
    getCorrectedExactOverlaps(query, &raw_overlaps, &out_overlaps);
    return out_overlaps;
}

void OverlapExtractorWithCorrection::getRawOverlapsDirect(const std::string& query, SequenceOverlapPairVector* out_vector)
{
    *out_vector = KmerOverlaps::retrieveMatches(query, m_k, m_minOverlap, m_minIdentity, 2, m_index_set);
}

void OverlapExtractorWithCorrection::getRawOverlapsCached(const std::string& query, bool is_reverse, SequenceOverlapPairVector* out_vector)
{
    PROFILE_FUNC("OverlapExtractorWithCorrection::getRawOverlaps")

    size_t nk = query.size() - m_k + 1;
    std::set<std::string> sequences_seen;

    for(size_t i = 0; i < nk; ++i)
    {
        std::string kmer = query.substr(i, m_k);
        std::string q_kmer = is_reverse ? reverseComplement(kmer) : kmer;
        const std::vector<size_t>& indices = m_cache_map.find(q_kmer)->second;

        for(size_t j = 0; j < indices.size(); ++j)
        {
            size_t index = indices[j];
            std::string match_sequence = m_strings[index];
            if(is_reverse)
                match_sequence = reverseComplement(match_sequence);

            if(sequences_seen.find(match_sequence) != sequences_seen.end())
                continue;
            sequences_seen.insert(match_sequence);
        
            // Ignore identical matches
            if(match_sequence.empty() || match_sequence == query)
                continue;

            // Compute the overlap. If the kmer match occurs a single time in each sequence we use
            // the banded extension overlap strategy. Otherwise we use the slow O(M*N) overlapper.
            SequenceOverlap overlap;
            size_t pos_0 = query.find(kmer);
            size_t pos_1 = match_sequence.find(kmer);

            // If there is a single occurrence of the kmer in each read,
            // use that position to seed the overlap calculation
            if(pos_0 == std::string::npos ||
               pos_1 == std::string::npos ||
               query.find(kmer, pos_0 + 1) != std::string::npos || 
               match_sequence.find(kmer, pos_1 + 1) != std::string::npos) 
            {
                // One of the reads has a second occurrence of the kmer. Use
                // the slow overlapper.
                overlap = Overlapper::computeOverlap(query, match_sequence);
            } else {
                overlap = Overlapper::extendMatch(query, match_sequence, pos_0, pos_1, 2);
            }
            bool bPassedOverlap = overlap.getOverlapLength() >= m_minOverlap;
            bool bPassedIdentity = overlap.getPercentIdentity() / 100 >= m_minIdentity;
            if(bPassedOverlap && bPassedIdentity)
            {
                SequenceOverlapPair op;
                op.sequence[0] = query;
                op.sequence[1] = match_sequence;
                op.overlap = overlap;
                op.is_reversed = is_reverse;
                out_vector->push_back(op);
            }
        }
    }
}

void OverlapExtractorWithCorrection::getCorrectedExactOverlaps(const std::string& query, 
                                                               const SequenceOverlapPairVector* raw_vector,
                                                               SequenceOverlapPairVector* out_vector)
{
    PROFILE_FUNC("OverlapExtractorWithCorrection::getCorrectedExactOverlaps")
    for(size_t i = 0; i < raw_vector->size(); ++i)
    {
        // Get the sequence of the raw read, on its original sequencing strand
        const SequenceOverlapPair& op = raw_vector->at(i);
        std::string sequence;
        if(op.is_reversed)
            sequence = reverseComplement(op.sequence[1]);
        else
            sequence = op.sequence[1];
        
        // Check if this read has been corrected in some other iteration
        std::string corrected_sequence;
        HashMap<std::string, std::string>::iterator find_iter = m_correction_cache.find(sequence);

        if(find_iter != m_correction_cache.end())
        {
            corrected_sequence = find_iter->second;
        }
        else
        {
            // Perform correction
            DNAString dna(sequence);
            SeqRecord record = { "null", dna, "" };
            //SeqRecord record("null", dna.toString(), ""); //JEREMIAH
            SequenceWorkItem wi(0, record);
            ErrorCorrectResult r = m_corrector->correct(wi);
            if(r.kmerQC || r.overlapQC)
                corrected_sequence = r.correctSequence.toString();
            
            // Insert the sequence into the cache
            // We allow it to be empty 
            m_correction_cache.insert(std::make_pair(sequence, corrected_sequence));
        }

        // Skip reads that are empty after correction
        if(corrected_sequence.empty())
            continue;

        // Change the strand of the corrected read back to the same as the query sequence
        if(op.is_reversed)
            corrected_sequence = reverseComplement(corrected_sequence);
        
        // Skip identical matches
        if(corrected_sequence == query)
            continue;

        // Recompute the overlap, using the previous overlap as a guide
        const SequenceOverlap& initial_overlap = op.overlap;
        SequenceOverlap overlap = Overlapper::extendMatch(query, corrected_sequence, initial_overlap.match[0].start, initial_overlap.match[1].start, 2);
        if(overlap.edit_distance == 0 && overlap.getOverlapLength() >= m_minOverlap)
        {
            SequenceOverlapPair sop;
            sop.sequence[1] = corrected_sequence;
            sop.overlap = overlap;
            out_vector->push_back(sop);
        }
    }
}

// Update the cache such that it contains all the reads containing any k-mer in the query
void OverlapExtractorWithCorrection::extractAndUpdate(const std::string& query)
{
    PROFILE_FUNC("OverlapExtractorWithCorrection::extractAndUpdate")

    // Check which kmers have not been cached yet
    std::vector<std::string> not_cached_kmers;
    if(query.size() < m_k)
        return;
            
    size_t nk = query.size() - m_k + 1;
    for(size_t i = 0; i < nk; ++i)
    {
        std::string kmer = query.substr(i, m_k);
        if(m_cache_map.find(kmer) == m_cache_map.end())
            not_cached_kmers.push_back(kmer);
    }

    printf("Kmers not in cache: %zu of %zu\n", not_cached_kmers.size(), nk);

    // For the kmers not in the cache, lookup the reads using the FM-index
    SeqRecordVector left_vector;
    HapgenUtil::extractHaplotypeReads(not_cached_kmers, m_index_set,
                                      m_k, false, 1000000, 500, &left_vector, NULL);
    
    // Insert the new kmers into the cache. The updateStrings() function below will update
    // the strings vector and the kmer -> index cache
    for(size_t i = 0; i < not_cached_kmers.size(); ++i)
        m_cache_map.insert(std::make_pair(not_cached_kmers[i], std::vector<size_t>()));

    printf("Found %zu new reads\n", left_vector.size());
    // Insert the new strings to the collection.
    for(size_t i = 0; i < left_vector.size(); ++i)
        updateStrings(left_vector[i].seq.toString());
}

//
size_t OverlapExtractorWithCorrection::updateStrings(const std::string& incoming)
{
    m_strings.push_back(incoming);
    size_t new_index = m_strings.size() - 1;

    updateCache(incoming, new_index);
    return new_index;
}

// Add the index of a read pair into the cache map.
void OverlapExtractorWithCorrection::updateCache(const std::string& str, size_t index)
{
    assert(str.size() >= m_k);
    size_t nk = str.size() - m_k + 1;
    for(size_t i = 0; i < nk; ++i)
    {
        std::string kmer = str.substr(i, m_k);
        std::map<std::string, std::vector<size_t> >::iterator iter = m_cache_map.find(kmer);
        if(iter != m_cache_map.end())
            iter->second.push_back(index);
    }
}
