///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// HapgenUtil - Utility functions used for 
// working with the haplotype generation
// modules and output
//

#include "HapgenUtil.h"
#include "LRAlignment.h"
#include "Interval.h"
#include "Profiler.h"
#include "Verbosity.h"
#include "overlapper.h"

struct CandidateKmerAlignment
{
    size_t query_index;
    size_t target_index;
    size_t target_extrapolated_start;
    size_t target_extrapolated_end;
    size_t target_sequence_id;

    static bool sortByStart(const CandidateKmerAlignment& a, const CandidateKmerAlignment& b) { return a.target_extrapolated_start < b.target_extrapolated_start; }
    static bool equalByStart(const CandidateKmerAlignment& a, const CandidateKmerAlignment& b) { return a.target_extrapolated_start == b.target_extrapolated_start; }

};
typedef std::vector<CandidateKmerAlignment> CandidateVector;

// Align the haplotype to the reference genome represented by the BWT/SSA pair
void HapgenUtil::alignHaplotypeToReferenceKmer(size_t k,
                                               const std::string& haplotype,
                                               const BWTIndexSet& referenceIndex,
                                               const ReadTable* pReferenceTable,
                                               HapgenAlignmentVector& outAlignments)
{
    PROFILE_FUNC("HapgenUtil::alignHaplotypesToReferenceKmer")
    int64_t max_interval_size = 4;

    if(haplotype.size() < k)
        return;

    std::vector<int> event_count_vector;
    std::vector<HapgenAlignment> tmp_alignments;
    int min_events = std::numeric_limits<int>::max();

    // Align forward and reverse haplotype to reference
    for(size_t i = 0; i <= 1; ++i)
    {
        bool is_reverse = i == 1;
        std::string query = is_reverse ? reverseComplement(haplotype) : haplotype;

        // Find shared kmers between the haplotype and the reference
        CandidateVector candidates;

        size_t nqk = query.size() - k + 1;
        for(size_t j = 0; j < nqk; ++j)
        {
            std::string kmer = query.substr(j, k);

            // Find the interval of this kmer in the reference
            BWTInterval interval = BWTAlgorithms::findInterval(referenceIndex, kmer);
            if(!interval.isValid() || interval.size() >= max_interval_size)
                continue; // not found or too repetitive

            // Extract the reference location of these hits
            for(int64_t k = interval.lower; k  <= interval.upper; ++k)
            {
                SAElem elem = referenceIndex.pSSA->calcSA(k, referenceIndex.pBWT);

                // Make a candidate alignment
                CandidateKmerAlignment candidate;
                candidate.query_index = j;
                candidate.target_index = elem.getPos();
                candidate.target_extrapolated_start = candidate.target_index - candidate.query_index;
                candidate.target_extrapolated_end = candidate.target_extrapolated_start + query.size();
                candidate.target_sequence_id = elem.getID();
                candidates.push_back(candidate);
            }
        }

        // Remove duplicate candidates
        std::sort(candidates.begin(), candidates.end(), CandidateKmerAlignment::sortByStart);
        CandidateVector::iterator new_end = std::unique(candidates.begin(), candidates.end(), CandidateKmerAlignment::equalByStart);
        candidates.resize(new_end - candidates.begin());
        
        for(size_t j = 0; j < candidates.size(); ++j)
        {
            // Extract window around reference
            size_t window_size = 200;
            int ref_start = candidates[j].target_extrapolated_start - window_size;
            int ref_end = candidates[j].target_extrapolated_end + window_size;
            const SeqItem& ref_record = pReferenceTable->getRead(candidates[j].target_sequence_id);
            const DNAString& ref_sequence = ref_record.seq;
            if(ref_start < 0)
                ref_start = 0;

            if(ref_end > (int)ref_sequence.length())
                ref_end = ref_sequence.length();

            std::string ref_substring = ref_sequence.substr(ref_start, ref_end - ref_start);

            // Align haplotype to the reference
            SequenceOverlap overlap = alignHaplotypeToReference(ref_substring, query);
            if(overlap.score < 0 || !overlap.isValid())
                continue;

            int alignment_start = ref_start + overlap.match[0].start;
            int alignment_end = ref_start + overlap.match[0].end; // inclusive
            int alignment_length = alignment_end - alignment_start + 1;

            // Crude count of the number of distinct variation events
            bool has_indel = false;
            int num_events = overlap.edit_distance;
            std::stringstream c_parser(overlap.cigar);
            int len;
            char t;
            while(c_parser >> len >> t) 
            {
                assert(len > 0);

                // Only count one event per insertion/deletion
                if(t == 'D' || t == 'I')
                {
                    num_events -= (len - 1);
                    has_indel = true;
                }
            }

            // Skip poor alignments
            double mismatch_rate = 1.0f - (overlap.getPercentIdentity() / 100.f);
            if(mismatch_rate > 0.05f || overlap.total_columns < 50)
            {
                if(Verbosity::Instance().getPrintLevel() > 4)
                {
                    printf("Haplotype Alignment - Ignoring low quality alignment (%.3lf, %dbp, %d events) to %s:%d\n", 
                        1.0f - mismatch_rate, overlap.total_columns, num_events, ref_record.id.c_str(), ref_start);
                }
                continue;
            }

            bool is_snp = !has_indel && overlap.edit_distance == 1;

            HapgenAlignment aln(candidates[j].target_sequence_id, 
                                alignment_start, 
                                alignment_length, 
                                overlap.score, 
                                num_events,
                                is_reverse, 
                                is_snp);

            tmp_alignments.push_back(aln);
            event_count_vector.push_back(num_events);
            if(Verbosity::Instance().getPrintLevel() > 4)
            {
                printf("Haplotype Alignment - Accepting alignment (%.3lf, %dbp, %d events) to %s:%d\n", 
                    1.0f - mismatch_rate, overlap.total_columns, num_events, ref_record.id.c_str(), ref_start);
            }            
            // Record the best edit distance
            if(num_events < min_events) 
                min_events = num_events;
        }
    }

    // Copy the best alignments into the output
    int MAX_DIFF_TO_BEST = 10;
    int MAX_EVENTS = 8;
    assert(event_count_vector.size() == tmp_alignments.size());
    for(size_t i = 0; i < event_count_vector.size(); ++i)
    {

        if(event_count_vector[i] <= MAX_EVENTS && event_count_vector[i] - min_events <= MAX_DIFF_TO_BEST)
            outAlignments.push_back(tmp_alignments[i]);
        else if(Verbosity::Instance().getPrintLevel() > 3)
            printf("Haplotype Alignment - Ignoring alignment with too many events (%d)\n", event_count_vector[i]);

    }
}

SequenceOverlap HapgenUtil::alignHaplotypeToReference(const std::string& reference,
                                                      const std::string& haplotype)
{
    SequenceOverlap overlap = Overlapper::computeOverlapAffine(reference, haplotype);
    return overlap;
}


// Coalesce a set of alignments into distinct locations
void HapgenUtil::coalesceAlignments(HapgenAlignmentVector& alignments)
{
    if(alignments.empty())
        return;

    // Sort the alignments by reference id, then position
    std::sort(alignments.begin(), alignments.end());

    HapgenAlignmentVector outAlignments;

    // Iterate over the alignments in sorted order
    // If an alignment is distinct (=does not overlap) from the
    // previous alignment, add it to the output collection.
    
    // First alignment is always ok
    outAlignments.push_back(alignments[alignments.size()-1]);

    // Kees: start from back because alignments are sorted in order of increasing score
    for(size_t i = alignments.size()-1; i-- > 0;)
    {
        // Check this alignment against the last alignment added to the output set
        HapgenAlignment& prevAlign = outAlignments.back();
        const HapgenAlignment& currAlign = alignments[i];

        int s1 = prevAlign.position;
        int e1 = s1 + prevAlign.length;

        int s2 = currAlign.position;
        int e2 = s2 + currAlign.length;
        bool intersecting = Interval::isIntersecting(s1, e1, s2, e2);

        if(prevAlign.referenceID != currAlign.referenceID || !intersecting)
        {
            outAlignments.push_back(currAlign);
        }
        else
        {
            // merge the intersecting alignment into a window that covers both
            prevAlign.position = std::min(s1, s2);
            prevAlign.length = std::max(e1, e2) - prevAlign.position;
        }
    }

    alignments = outAlignments;
}

// Extract a substring for a reference defined by the alignment
std::string HapgenUtil::extractReference(const HapgenAlignment& aln, const ReadTable* pRefTable, int flanking)
{
    const SeqItem& refItem = pRefTable->getRead(aln.referenceID);

    int start = aln.position - flanking;
    int end = aln.position + aln.length + flanking;
    if(start < 0)
        start = 0;
    return refItem.seq.substr(start, end - start);
}


// Extract 3 substrings of a reference described by an alignment
// This writes to 3 out parameters, one for the defined position
// of the reference, one for the upstream flanking bases and one for the downstream flanking
// upstream   defined  downstream
// --------============----------
void HapgenUtil::extractReferenceSubstrings(const HapgenAlignment& aln, 
                                            const ReadTable* pRefTable, 
                                            int flanking,
                                            std::string& outUpstream, 
                                            std::string& outDefined, 
                                            std::string& outDownstream)
{
    const SeqItem& refItem = pRefTable->getRead(aln.referenceID);

    int definedStart = aln.position;
    int definedEnd = aln.position + aln.length;
    int upstreamFlankStart = definedStart - flanking;
    if(upstreamFlankStart < 0)
        upstreamFlankStart = 0;

    if((size_t)definedEnd + flanking > refItem.seq.length())
        flanking = refItem.seq.length() - definedEnd - 1;

    outUpstream = refItem.seq.substr(upstreamFlankStart, definedStart - upstreamFlankStart);
    outDefined = refItem.seq.substr(definedStart, definedEnd - definedStart);
    
    // Bugfix, the substr call will assert if the haplotype
    // aligns to the end of the chromosome
    if(definedEnd < (int)refItem.seq.length())
        outDownstream = refItem.seq.substr(definedEnd, flanking);
    else
        outDownstream = "";
}

bool HapgenUtil::makeFlankingHaplotypes(const HapgenAlignment& aln, 
                                        const ReadTable* pRefTable, 
                                        int flanking,
                                        const StringVector& inHaplotypes,
                                        StringVector& outFlankingHaplotypes,
                                        StringVector& outHaplotypes)
{
    std::string upstream;
    std::string referenceHaplotype;
    std::string downstream;

    extractReferenceSubstrings(aln, pRefTable, flanking, upstream, referenceHaplotype, downstream);

    // Flip reference strings to match the strand of the input haplotypes
    if(aln.isRC)
    {
        // reverse complement each string
        upstream = reverseComplement(upstream);
        referenceHaplotype = reverseComplement(referenceHaplotype);
        downstream = reverseComplement(downstream);

        // Swap up and downstream
        upstream.swap(downstream);
    }

    // Make the reference haplotype w/ flanking sequence
    std::string referenceFlanking = upstream + referenceHaplotype + downstream;
    outFlankingHaplotypes.push_back(referenceFlanking);
    outHaplotypes.push_back(referenceHaplotype);

    // Check that all sequences match the reference haplotype properly
    /*
    bool checkOk = checkAlignmentsAreConsistent(referenceFlanking, inHaplotypes);
    if(!checkOk)
    {
        outHaplotypes.clear();
        return false;
    }
    */

    // Make the flanking sequences for each haplotype
    for(size_t i = 0; i < inHaplotypes.size(); ++i)
    {
        // Skip if the input haplotype exactly matches the reference
        if(inHaplotypes[i] != referenceHaplotype)
        {
            outFlankingHaplotypes.push_back(upstream + inHaplotypes[i] + downstream);
            outHaplotypes.push_back(inHaplotypes[i]);
        }
    }

    return true;
}
    
// Check that all the strings in the vector align to the same coordinates
// of the passed in sequence
bool HapgenUtil::checkAlignmentsAreConsistent(const std::string& refString, const StringVector& queries)
{
    if(queries.empty())
        return true;

    // Perform local alignments of each query to the refString
    LocalAlignmentResultVector alignments;
    for(size_t i = 0; i < queries.size(); ++i)
        alignments.push_back(StdAlnTools::localAlignment(refString, queries[i]));

    size_t i = 0;
    for(size_t j = 1; j < alignments.size(); ++j)
    {
        if(alignments[i].targetStartIndex != alignments[j].targetStartIndex ||
           alignments[j].targetEndIndex != alignments[j].targetEndIndex)
        {
            std::cerr << "Warning: inconsistent alignments found for haplotype realignment\n";
            std::cerr << "A[" << i << "]: " << alignments[i] << "\n";
            std::cerr << "A[" << j << "]: " << alignments[j] << "\n";
            return false;
        }
    }

    return true;
}

// Extract reads from an FM-index that have a k-mer match to any given haplotypes
// Returns true if the reads were successfully extracted, false if there are 
// more reads than maxReads
bool HapgenUtil::extractHaplotypeReads(const StringVector& haplotypes, 
                                       const BWTIndexSet& indices,
                                       int k,
                                       bool doReverse,
                                       size_t maxReads,
                                       int64_t maxIntervalSize,
                                       SeqRecordVector* pOutReads, 
                                       SeqRecordVector* pOutMates)
{
    PROFILE_FUNC("HapgenUtil::extractHaplotypeReads")
    // Extract the set of reads that have at least one kmer shared with these haplotypes
    // This is a bit of a lengthy procedure with a few steps:
    // 1) extract all the kmers in the haplotypes
    // 2) find the intervals for the kmers in the fm-index
    // 3) compute the set of read indices of the reads from the intervals (using the sampled suffix array)
    // 4) finally, extract the read sequences from the index
    // Make a set of kmers from the haplotypes
    std::set<std::string> kmerSet;
    for(size_t i = 0; i < haplotypes.size(); ++i)
    {
        const std::string& h = haplotypes[i];
        if((int)h.size() < k)
            continue;

        for(size_t j = 0; j < h.size() - k + 1; ++j)
        {
            std::string ks = h.substr(j, k);
            if(doReverse)
                ks = reverseComplement(ks);
            kmerSet.insert(ks);
        }
    }

    // Compute suffix array intervals for the kmers
    std::vector<BWTInterval> intervals;
    for(std::set<std::string>::const_iterator iter = kmerSet.begin(); iter != kmerSet.end(); ++iter)
    {
        BWTInterval interval = BWTAlgorithms::findInterval(indices, *iter);
        if(interval.size() < maxIntervalSize)
            intervals.push_back(interval);
    }

    // Compute the set of reads ids 
    std::set<int64_t> readIndices;
    for(size_t i = 0; i < intervals.size(); ++i)
    {
        BWTInterval interval = intervals[i];
        for(int64_t j = interval.lower; j <= interval.upper; ++j)
        {
            // Get index from sampled suffix array
            SAElem elem = indices.pSSA->calcSA(j, indices.pBWT);
            readIndices.insert(elem.getID());
        }
    }

    // Check if we have hit the limit of extracting too many reads
    if(readIndices.size() > maxReads)
        return false;

    for(std::set<int64_t>::const_iterator iter = readIndices.begin(); iter != readIndices.end(); ++iter)
    {
        int64_t idx = *iter;
        
        // Extract the read
        std::stringstream namer;
        namer << "idx-" << idx;
        SeqRecord record;
        record.id = namer.str();
        record.seq = BWTAlgorithms::extractString(indices.pBWT, idx);

        assert(indices.pQualityTable != NULL);
        record.qual = indices.pQualityTable->getQualityString(idx, record.seq.length());
        if(!record.seq.empty())
            pOutReads->push_back(record);

        // Optionally extract its mate
        // If the index is constructed properly, 
        // paired reads are in adjacent indices with the
        // first read at even indices
        if(pOutMates != NULL)
        {
            int64_t mateIdx = idx;
            if(idx % 2 == 0)
                mateIdx += 1;
            else
                mateIdx -= 1;
            
            std::stringstream mateName;
            mateName << "idx-" << mateIdx;
            SeqRecord mateRecord;
            mateRecord.id = mateName.str();
            mateRecord.seq = BWTAlgorithms::extractString(indices.pBWT, mateIdx);
            mateRecord.qual = indices.pQualityTable->getQualityString(mateIdx, mateRecord.seq.length());
            if(!record.seq.empty() && !mateRecord.seq.empty())
                pOutMates->push_back(mateRecord);
        }
    }
    return true;
}

// Extract reads from an FM-index that have a k-mer match to AT MOST one haplotype
// Returns true if the reads were successfully extracted, false if there are 
// more reads than maxReads
bool HapgenUtil::extractHaplotypeSpecificReads(const StringVector& haplotypes, 
                                               const BWTIndexSet& indices,
                                               int k,
                                               bool doReverse,
                                               size_t maxReads,
                                               int64_t maxIntervalSize,
                                               SeqRecordVector* pOutReads, 
                                               SeqRecordVector* pOutMates)
{
    std::map<std::string, int> kmer_map;
    for(size_t i = 0; i < haplotypes.size(); ++i)
    {
        const std::string& h = haplotypes[i];
        printf("Haplotype[%zu]: %s\n", i, h.c_str());
        if((int)h.size() < k)
            continue;

        for(size_t j = 0; j < h.size() - k + 1; ++j)
            kmer_map[h.substr(j,k)]++;
    }

    StringVector specific_kmer_vector;
    for(std::map<std::string, int>::iterator iterator = kmer_map.begin(); 
                                             iterator != kmer_map.end(); ++iterator)
    {
        if(iterator->second == 1)
            specific_kmer_vector.push_back(iterator->first);
    }

    printf("%zu of %zu kmers are haplotype-unique\n", specific_kmer_vector.size(), kmer_map.size());
    return extractHaplotypeReads(specific_kmer_vector, indices, k, doReverse, maxReads, maxIntervalSize, pOutReads, pOutMates);
}


// Align a bunch of reads locally to a sequence
LocalAlignmentResultVector HapgenUtil::alignReadsLocally(const std::string& target, const SeqItemVector& reads)
{
    LocalAlignmentResultVector results;
    for(size_t i = 0; i < reads.size(); ++i)
    {
        LocalAlignmentResult fwdAR = StdAlnTools::localAlignment(target, reads[i].seq.toString());
        LocalAlignmentResult rcAR = StdAlnTools::localAlignment(target, reverseComplement(reads[i].seq.toString()));
        results.push_back(fwdAR.score > rcAR.score ? fwdAR : rcAR);
    }
    return results;
}

//
double HapgenUtil::calculateDustScoreAtPosition(const std::string& name, 
                                                int position, 
                                                const ReadTable* pRefTable,
                                                int window_size)
{
    int dw_start = position - window_size / 2;
    dw_start = std::max(0, dw_start); // clamp to 0

    int dw_end = position + window_size / 2;
    const DNAString& chromosome = pRefTable->getRead(name).seq;
    dw_end = std::min(dw_end, (int)chromosome.length());

    double dust_score = 0.0f;
    if(dw_end - dw_start == window_size)
    {
        std::string ref_dust_window = chromosome.substr(dw_start, dw_end - dw_start);
        dust_score = calculateDustScore(ref_dust_window);
    }
    return dust_score;
}

//
size_t HapgenUtil::getMaximumOneEdit(const std::string& str, const BWTIndexSet& indices)
{
    size_t max = 0;
    std::string t = str;

    //
    for(size_t i = 0; i < str.size(); ++i)
    {
        char tmp = t[i];
        for(int bi = 0; bi < DNA_ALPHABET::size; ++bi)
        {
            char b = DNA_ALPHABET::getBase(bi);
            if(b != tmp)
            {
                t[i] = b;
                size_t count = BWTAlgorithms::countSequenceOccurrences(t, indices);
                if(count > max)
                    max = count;
            }
        }
        t[i] = tmp;
    }

    printf("MaxOneEdit %zu\n", max);
    return max;
}

// Calculate the largest k such that every k-mer in the sequence is present at least min_depth times in the BWT
size_t HapgenUtil::calculateMaxCoveringK(const std::string& sequence, int min_depth, const BWTIndexSet& indices)
{
    size_t min_k = 15;
    for(size_t k = 99; k >= min_k; --k)
    {
        if(sequence.size() < k)
            continue;
        bool covered = true;
        size_t nk = sequence.size() - k + 1;
        for(size_t i = 0; i < nk; ++i)
        {
            std::string kmer = sequence.substr(i, k);
            int c = BWTAlgorithms::countSequenceOccurrences(kmer, indices);

            if(c < min_depth)
            {
                covered = false;
                break;
            }
        }

        if(covered)
            return k;
    }

    return 0;
}

//
std::vector<int> HapgenUtil::makeCountProfile(const std::string& str, size_t k, int max, const BWTIndexSet& indices)
{
    std::vector<int> out;
    if(str.size() < k)
        return out;

    for(size_t i = 0; i < str.size() - k + 1; ++i)
    {
        int count = BWTAlgorithms::countSequenceOccurrences(str.substr(i, k), indices);
        out.push_back(count > max ? max : count);
    }
    return out;
}


// Print an alignment to a reference
void HapgenUtil::printAlignment(const std::string& query, const HapgenAlignment& aln, const ReadTable* pRefTable)
{
    const SeqItem& refItem = pRefTable->getRead(aln.referenceID);
    size_t refStart = aln.position;
    size_t refEnd = refStart + aln.length;
    std::string refSubstring = refItem.seq.substr(refStart, refEnd - refStart);

    std::string q = aln.isRC ? reverseComplement(query) : query;
    StdAlnTools::globalAlignment(q, refSubstring, true);
}

