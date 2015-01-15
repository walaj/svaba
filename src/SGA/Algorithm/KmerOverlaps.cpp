///-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// KmerOverlaps - Overlap computation functions
// seeded by exact kmer matches
//
#include "KmerOverlaps.h"
#include "HashMap.h"
#include "BWTAlgorithms.h"
#include "Profiler.h"
#include "Timer.h"
#include "Verbosity.h"

//
MultipleAlignment KmerOverlaps::buildMultipleAlignment(const std::string& query, 
                                                       size_t k,
                                                       int min_overlap,
                                                       double min_identity,
                                                       int bandwidth,
                                                       const BWTIndexSet& indices)
{
    SequenceOverlapPairVector overlap_vector = retrieveMatches(query, k, min_overlap, min_identity, bandwidth, indices);
    MultipleAlignment multiple_alignment;
    multiple_alignment.addBaseSequence("query", query, "");
    for(size_t i = 0; i < overlap_vector.size(); ++i)
        multiple_alignment.addOverlap("null", overlap_vector[i].sequence[1], "", overlap_vector[i].overlap);
    return multiple_alignment;
}

// Struct to hold a partial match in the FM-index
// The position field is the location in the query sequence of this kmer.
// The index field is an index into the BWT. 
// The is_reverse flag indicates the strand of the partial match
struct KmerMatch
{
    size_t position:16;
    size_t index:47;
    size_t is_reverse:1;

    friend bool operator<(const KmerMatch& a, const KmerMatch& b)
    {
        if(a.index == b.index)
            return a.is_reverse < b.is_reverse;
        else
            return a.index < b.index;
    }

    friend bool operator==(const KmerMatch& a, const KmerMatch& b)
    {
        return a.index == b.index && a.is_reverse == b.is_reverse;
    }
};

// Return a hash key for a KmerMatch
struct KmerMatchKey
{
    size_t operator()(const KmerMatch& a) const { return a.index; }
};

typedef std::set<KmerMatch> KmerMatchSet;
typedef HashMap<KmerMatch, bool, KmerMatchKey> KmerMatchMap;

//
SequenceOverlapPairVector KmerOverlaps::retrieveMatches(const std::string& query, size_t k, 
                                                        int min_overlap, double min_identity,
                                                        int bandwidth, const BWTIndexSet& indices)
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::retrieveMatches")
    assert(indices.pBWT != NULL);
    assert(indices.pSSA != NULL);

    static size_t n_calls = 0;
    static size_t n_candidates = 0;
    static size_t n_output = 0;
    static double t_time = 0;
    Timer timer("test", true);

    n_calls++;

    int64_t max_interval_size = 200;
    SequenceOverlapPairVector overlap_vector;

    // Use the FM-index to look up intervals for each kmer of the read. Each index
    // in the interval is stored individually in the KmerMatchMap. We then
    // backtrack to map these kmer indices to read IDs. As reads can share
    // multiple kmers, we use the map to avoid redundant lookups.
    // There is likely a faster algorithm which performs direct decompression
    // of the read sequences without having to expand the intervals to individual
    // indices. The current algorithm suffices for now.
    KmerMatchMap prematchMap;
    size_t num_kmers = query.size() - k + 1;
    for(size_t i = 0; i < num_kmers; ++i)
    {
        std::string kmer = query.substr(i, k);
        BWTInterval interval = BWTAlgorithms::findInterval(indices, kmer);
        if(interval.isValid() && interval.size() < max_interval_size) 
        {
            for(int64_t j = interval.lower; j <= interval.upper; ++j)
            {
                KmerMatch match = { i, static_cast<size_t>(j), false };
                prematchMap.insert(std::make_pair(match, false));
            }
        }

        kmer = reverseComplement(kmer);
        interval = BWTAlgorithms::findInterval(indices, kmer);
        if(interval.isValid() && interval.size() < max_interval_size) 
        {
            for(int64_t j = interval.lower; j <= interval.upper; ++j)
            {
                KmerMatch match = { i, static_cast<size_t>(j), true };
                prematchMap.insert(std::make_pair(match, false));
            }
        }
    }

    // Backtrack through the kmer indices to turn them into read indices.
    // This mirrors the calcSA function in SampledSuffixArray except we mark each entry
    // as visited once it is processed.
    KmerMatchSet matches;
    for(KmerMatchMap::iterator iter = prematchMap.begin(); iter != prematchMap.end(); ++iter)
    {
        // This index has been visited
        if(iter->second)
            continue;

        // Mark this as visited
        iter->second = true;

        // Backtrack the index until we hit the starting symbol
        KmerMatch out_match = iter->first;
        while(1) 
        {
            char b = indices.pBWT->getChar(out_match.index);
            out_match.index = indices.pBWT->getPC(b) + indices.pBWT->getOcc(b, out_match.index - 1);

            // Check if the hash indicates we have visited this index. If so, stop the backtrack
            KmerMatchMap::iterator find_iter = prematchMap.find(out_match);
            if(find_iter != prematchMap.end())
            {
                // We have processed this index already
                if(find_iter->second)
                    break;
                else
                    find_iter->second = true;
            }

            if(b == '$')
            {
                // We've found the lexicographic index for this read. Turn it into a proper ID
                out_match.index = indices.pSSA->lookupLexoRank(out_match.index);
                matches.insert(out_match);
                break;
            }
        }
    }

    // Refine the matches by computing proper overlaps between the sequences
    // Use the overlaps that meet the thresholds to build a multiple alignment
    for(KmerMatchSet::iterator iter = matches.begin(); iter != matches.end(); ++iter)
    {
        std::string match_sequence = BWTAlgorithms::extractString(indices.pBWT, iter->index);
        if(iter->is_reverse)
            match_sequence = reverseComplement(match_sequence);
        
        // Ignore identical matches
        if(match_sequence == query)
            continue;

        // Compute the overlap. If the kmer match occurs a single time in each sequence we use
        // the banded extension overlap strategy. Otherwise we use the slow O(M*N) overlapper.
        SequenceOverlap overlap;
        std::string match_kmer = query.substr(iter->position, k);
        size_t pos_0 = query.find(match_kmer);
        size_t pos_1 = match_sequence.find(match_kmer);
        assert(pos_0 != std::string::npos && pos_1 != std::string::npos);

        // Check for secondary occurrences
        if(query.find(match_kmer, pos_0 + 1) != std::string::npos || 
           match_sequence.find(match_kmer, pos_1 + 1) != std::string::npos) {
            // One of the reads has a second occurrence of the kmer. Use
            // the slow overlapper.
            overlap = Overlapper::computeOverlap(query, match_sequence);
        } else {
            overlap = Overlapper::extendMatch(query, match_sequence, pos_0, pos_1, bandwidth);
        }

        n_candidates += 1;
        bool bPassedOverlap = overlap.getOverlapLength() >= min_overlap;
        bool bPassedIdentity = overlap.getPercentIdentity() / 100 >= min_identity;

        if(bPassedOverlap && bPassedIdentity)
        {
            SequenceOverlapPair op;
            op.sequence[0] = query;
            op.sequence[1] = match_sequence;
            op.overlap = overlap;
            op.is_reversed = iter->is_reverse;
            overlap_vector.push_back(op);
            n_output += 1;
        }
    }

    t_time += timer.getElapsedCPUTime();

    if(Verbosity::Instance().getPrintLevel() > 6 && n_calls % 100 == 0)
        printf("[kmer overlaps] n: %zu candidates: %zu valid: %zu (%.2lf) time: %.2lfs\n", 
            n_calls, n_candidates, n_output, (double)n_output / n_candidates, t_time);
    return overlap_vector;
}

struct SeedEdit
{
    SeedEdit(int i, char b) : index(i), base(b) {}
    int index;
    char base;
};

typedef std::vector<SeedEdit> SeedEditVector;

struct ApproxSeed
{
    int query_index; // the index of the last query base included in the interval
    BWTInterval interval;
    int length;
    SeedEditVector edits;

    friend std::ostream& operator<<(std::ostream& o, ApproxSeed& a) {
        o << "QI: " << a.query_index << " IV: " << a.interval;
        o << " Edits: ";
        for(size_t i = 0; i < a.edits.size(); i++)
            o << a.edits[i].index << ":" << a.edits[i].base << ",";
        return o;
    }
};

void _approximateSeededMatch(const std::string& in_query,
                             int min_overlap, 
                             double min_identity,
                             int bandwidth, 
                             int max_interval,
                             bool do_reverse,
                             const BWTIndexSet& indices,
                             SequenceOverlapPairVector& out_vector)
{
    Timer timer("test", true);
    assert(indices.pBWT != NULL);
    assert(indices.pSSA != NULL);
    assert(indices.pCache != NULL);

    static size_t n_calls = 0;
    static size_t n_candidates = 0;
    static size_t n_output = 0;
    static double t_time = 0;

    n_calls++;

    int target_seed_length = 41;
    int seed_stride = target_seed_length / 2;
    size_t d = 1;

    SequenceOverlapPairVector overlap_vector;

    std::string strand_query = do_reverse ? reverseComplement(in_query) : in_query;

    // Initialize seeds 
    int seed_end = strand_query.size();
    int q = indices.pCache->getCachedLength();

    std::queue<ApproxSeed> seeds;

    while(seed_end > target_seed_length)
    {
        // For the last q bases of the seed, create all strings within edit distance d
        std::string qmer = strand_query.substr(seed_end - q, q);
        assert((int)qmer.size() == q);
        // 0-distance seed
        ApproxSeed seed;
        seed.query_index = seed_end - q;
        seed.interval = indices.pCache->lookup(qmer.c_str());
        seed.length = q;
        seeds.push(seed);

        for(int i = 0; i < q; ++i)
        {
            // Switch base to other 3 symbols
            char o = qmer[i];
            for(size_t j = 0; j < 4; ++j)
            {
                char b = "ACGT"[j];
                if(b != o)
                {
                    qmer[i] = b;
                    ApproxSeed seed;
                    seed.query_index = seed_end - q;
                    seed.interval = indices.pCache->lookup(qmer.c_str());
                    seed.length = q;
                    seed.edits.push_back(SeedEdit(i + seed.query_index, b));
                    seeds.push(seed);
                }
            }
            qmer[i] = o;
        }
        seed_end -= seed_stride;
    }
    
    // Extend seeds
    std::vector<ApproxSeed> finished_seeds;
    while(!seeds.empty())
    {
        ApproxSeed& seed = seeds.front();

        // query_index is the index of the last base
        // in the seed. get the next base
        char qb = strand_query[seed.query_index - 1];

        // Branch to inexact match
        if(seed.edits.size() < d)
        {
            for(size_t j = 0; j < 4; ++j)
            {
                char b = "ACGT"[j];
                if(b != qb)
                {
                    ApproxSeed new_seed;
                    new_seed.query_index = seed.query_index - 1;
                    new_seed.interval = seed.interval;
                    BWTAlgorithms::updateInterval(new_seed.interval, b, indices.pBWT);
                    if(new_seed.interval.isValid())
                    {
                        new_seed.length = seed.length + 1;
                        new_seed.edits = seed.edits;
                        new_seed.edits.push_back(SeedEdit(new_seed.query_index, b));

                        if(new_seed.length < target_seed_length && new_seed.query_index > 0)
                            seeds.push(new_seed);
                        else
                            finished_seeds.push_back(new_seed);
                    }
                }
            }
        }
            
        // Extend with the actual query base without branching
        seed.query_index = seed.query_index - 1;
        seed.length += 1;
        BWTAlgorithms::updateInterval(seed.interval, qb, indices.pBWT);
        if(!seed.interval.isValid() || seed.length >= target_seed_length || seed.query_index == 0)
        {
            if(seed.interval.isValid())
                finished_seeds.push_back(seed);
            seeds.pop();
        }
    }
    
    std::set<size_t> rank_set;
    
    for(size_t i = 0; i < finished_seeds.size(); ++i)
    {
        if(finished_seeds[i].interval.size() > max_interval)
            continue;

        //std::cout << finished_seeds[i] << "\n";
        std::string query_seed = strand_query.substr(finished_seeds[i].query_index, target_seed_length);
        std::string match_seed = query_seed;

        // Apply edits to the new sequence
        for(size_t j = 0; j < finished_seeds[i].edits.size(); ++j)
            match_seed[finished_seeds[i].edits[j].index - finished_seeds[i].query_index] = finished_seeds[i].edits[j].base;

        // Flip the seeds to match the strand of the query
        if(do_reverse)
        {
            query_seed = reverseComplement(query_seed);
            match_seed = reverseComplement(match_seed);
        }

        // Extract the prefix of every occurrence of this seed
        RankedPrefixVector extensions = BWTAlgorithms::extractRankedPrefixes(indices.pBWT, finished_seeds[i].interval);

        // Extend the seeds to the full-length string
        for(size_t j = 0; j < extensions.size(); ++j)
        {
            size_t rank = extensions[j].rank;

            // The second element of the returned pair is
            // false if the set already contains this rank
            if(!rank_set.insert(rank).second)
                continue;

            // Extract the reminder of the read
            std::string& prefix = extensions[j].prefix;
            int64_t start_index_of_read = indices.pSSA->lookupLexoRank(rank);
            std::string suffix = BWTAlgorithms::extractUntilInterval(indices.pBWT, 
                                                                     start_index_of_read, 
                                                                     finished_seeds[i].interval);

            std::string match_sequence = prefix + suffix;

            // Ignore identical matches
            if(match_sequence == strand_query)
                continue;

            // Change strands
            if(do_reverse)
                match_sequence = reverseComplement(match_sequence);

            // Compute the overlap
            SequenceOverlap overlap;
            size_t pos_0 = in_query.find(query_seed);
            size_t pos_1 = match_sequence.find(match_seed);
            assert(pos_0 != std::string::npos);
            assert(pos_1 != std::string::npos);

            if(in_query.find(query_seed, pos_0 + 1) != std::string::npos || 
               match_sequence.find(match_seed, pos_1 + 1) != std::string::npos) {
                // One of the reads has a second occurrence of the kmer. Use
                // the slow overlapper.
                overlap = Overlapper::computeOverlap(in_query, match_sequence);
            } else {
                overlap = Overlapper::extendMatch(in_query, match_sequence, pos_0, pos_1, bandwidth);
            }

            n_candidates += 1;
            bool bPassedOverlap = overlap.getOverlapLength() >= min_overlap;
            bool bPassedIdentity = overlap.getPercentIdentity() / 100 >= min_identity;

            if(bPassedOverlap && bPassedIdentity)
            {
                //printf("Rank\t%zu\t%zu\t%s\t%.2lf\t%d\n", n_calls,
                //    rank, match_sequence.c_str(), overlap.getPercentIdentity(), overlap.getOverlapLength());
                
                SequenceOverlapPair op;
                op.sequence[0] = in_query;
                op.sequence[1] = match_sequence;
                op.overlap = overlap;
                op.is_reversed = do_reverse;
                out_vector.push_back(op);
                n_output += 1;
            }                
        }
    }
    t_time += timer.getElapsedCPUTime();
    
    if(Verbosity::Instance().getPrintLevel() > 6 && n_calls % 100 == 0)
        printf("[approx seeds] n: %zu candidates: %zu valid: %zu (%.2lf) time: %.2lfs\n", 
            n_calls, n_candidates, n_output, (double)n_output / n_candidates, t_time);
}

//
SequenceOverlapPairVector KmerOverlaps::approximateMatch(const std::string& query,
                                                         int min_overlap, 
                                                         double min_identity,
                                                         int bandwidth, 
                                                         int max_interval,
                                                         const BWTIndexSet& indices)
{
    SequenceOverlapPairVector opv;
    _approximateSeededMatch(query, min_overlap, min_identity, bandwidth, max_interval, false, indices, opv);
    _approximateSeededMatch(query, min_overlap, min_identity, bandwidth, max_interval, true, indices, opv);
    return opv;
}
