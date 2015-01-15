///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// GraphCompare - Compare two (abstractly represented)
// graphs against each other to find strings
// only present in one.
//
// The graphs are abstractly represented as
// an FM-index.
//
#include "GraphCompare.h"
#include "BWTAlgorithms.h"
#include "SGAlgorithms.h"
#include "SGSearch.h"
#include "StdAlnTools.h"
#include "LRAlignment.h"
#include "HapgenUtil.h"
#include "DindelRealignWindow.h"
#include "DindelUtil.h"
#include "HaplotypeBuilder.h"
#include "PairedDeBruijnHaplotypeBuilder.h"
#include "DeBruijnHaplotypeBuilder.h"
#include "OverlapHaplotypeBuilder.h"
#include "StringHaplotypeBuilder.h"
#include "Profiler.h"
#include "Verbosity.h"
#include "api/SamHeader.h"
#include "api/BamAlignment.h"

// #define GRAPH_DIFF_DEBUG 1

//
// GraphCompareStats
//
void GraphCompareStats::clear()
{
    numBubbles = 0;
    numAttempted = 0;
    numTargetBranched = 0;
    numSourceBranched = 0;
    numTargetBroken = 0;
    numSourceBroken = 0;
    numWalkFailed = 0;
    numHBFailed = 0;
    numNoSolution = 0;

    numInsertions = 0;
    numDeletions = 0;
    numSubs = 0;
}

//
void GraphCompareStats::add(const GraphCompareStats& other)
{
    numBubbles += other.numBubbles;
    numAttempted += other.numAttempted;
    numTargetBranched += other.numTargetBranched;
    numSourceBranched += other.numSourceBranched;
    numTargetBroken += other.numTargetBroken;
    numSourceBroken += other.numSourceBroken;
    numWalkFailed += other.numWalkFailed;
    numHBFailed += other.numHBFailed;
    numNoSolution += other.numNoSolution;

    numInsertions += other.numInsertions;
    numDeletions += other.numDeletions;
    numSubs += other.numSubs;
}

//
void GraphCompareStats::print() const
{
    /*
    printf("Total attempts: %d\n", numAttempted);
    printf("Total bubbles: %d\n", numBubbles);
    printf("Failed - target branched: %d\n", numTargetBranched);
    printf("Failed - source branched: %d\n", numSourceBranched);
    printf("Failed - target broken: %d\n", numTargetBroken);
    printf("Failed - source broken: %d\n", numSourceBroken);
    printf("Failed - haplotype builder failed: %d\n", numHBFailed);
    printf("Failed - no walk: %d\n", numWalkFailed);
    printf("Failed - no solution: %d\n", numNoSolution);
    
    printf("Num subs found: %d\n", numSubs);
    printf("Num insertions found: %d\n", numInsertions);
    printf("Num deletions found: %d\n", numDeletions);
    */
}

//
//
//
GraphCompare::GraphCompare(const GraphCompareParameters& params) : m_parameters(params)
{
    m_stats.clear();

    // Parameter sanity checking
    assert(m_parameters.maxHaplotypes > 0);
    assert(m_parameters.kmer > 0);
    assert(m_parameters.minDiscoveryCount > 0);
}

//
GraphCompare::~GraphCompare()
{
}

GraphCompareResult GraphCompare::process(const SequenceWorkItem& item) const
{
    PROFILE_FUNC("GraphCompare::process")
    GraphCompareResult result;
    SeqRecord currRead = item.read;
    std::string w = item.read.seq.toString();
    int len = w.size();
    int num_kmers = len - m_parameters.kmer + 1;

    // Check if the read is long enough to be used
    if(w.size() < m_parameters.kmer)
        return result;

    // Process the kmers that have not been previously visited
    // We only try to find variants from one k-mer per read. 
    // This heurestic handles the situation where we process a kmer,
    // spend a lot of time trying to assemble it into a string and fail,
    // then try again with the next k-mer in the read. Since all the k-mers
    // in a read are connected in the de Bruijn graph, if the first attempt
    // to assemble a k-mer into a variant fails, it is very unlikely that
    // subsequent attempts will succeed.
    bool variantAttempted = false;
    for(int j = 0; j < num_kmers; ++j)
    {
        std::string kmer = w.substr(j, m_parameters.kmer);
        std::string rc_kmer = reverseComplement(kmer);

        // Use the lexicographically lower of the kmer and its pair as the key in the bloom filter
        std::string& key_kmer = kmer < rc_kmer ? kmer : rc_kmer;

        // Check if this k-mer is marked as used by the bloom filter
        if(m_parameters.pBloomFilter->test(key_kmer.c_str(), key_kmer.size()))
            continue;
        
    	// Get the interval for this kmer
        BWTInterval interval = BWTAlgorithms::findInterval(m_parameters.variantIndex, kmer);
        BWTInterval rc_interval = BWTAlgorithms::findInterval(m_parameters.variantIndex, rc_kmer);
        
        size_t count = interval.size();
        if(rc_interval.isValid())
            count += rc_interval.size();

        bool both_strands = interval.size() > 0 && rc_interval.size() > 0;

        size_t min_base_coverage = 1;
        if(count >= m_parameters.minDiscoveryCount && count < m_parameters.maxDiscoveryCount && !variantAttempted && both_strands)
        {
            // Update the bloom filter to contain this kmer
            m_parameters.pBloomFilter->add(key_kmer.c_str(), key_kmer.size());

            // Check if this k-mer is present in the other base index
            size_t base_count = BWTAlgorithms::countSequenceOccurrences(kmer, m_parameters.baseIndex);
            
            if(Verbosity::Instance().getPrintLevel() > 6)
                std::cout << "Read: " << currRead.id << " k: " << j << " CV: " << interval.size() << "/" << rc_interval.size() << " " << base_count << "\n";
            
            // k-mer present in the base read set, skip it
            if(base_count >= min_base_coverage)
                continue;
            
            if(Verbosity::Instance().getPrintLevel() > 0)
                std::cout << "Variant read: " << w << "\n";

            // variant k-mer, attempt to assemble it into haplotypes
            GraphBuildResult build_result = processVariantKmer(kmer, count);
            variantAttempted = true;

            // Mark the kmers of the variant haplotypes as being visited
            for(size_t vhi = 0; vhi < build_result.variant_haplotypes.size(); ++vhi)
                markVariantSequenceKmers(build_result.variant_haplotypes[vhi]);
            
            // If we assembled anything, run Dindel on the haplotypes
            if(build_result.variant_haplotypes.size() > 0)
            {
                if(Verbosity::Instance().getPrintLevel() > 0)
                    std::cout << "Running dindel\n";

                std::stringstream baseVCFSS;
                std::stringstream variantVCFSS;
                std::stringstream callsVCFSS;
                DindelReadReferenceAlignmentVector alignments;

                DindelReturnCode drc = DindelUtil::runDindelPairMatePair(kmer,
                                                                         build_result.base_haplotypes,
                                                                         build_result.variant_haplotypes,
                                                                         m_parameters,
                                                                         baseVCFSS,
                                                                         variantVCFSS,
                                                                         callsVCFSS,
                                                                         &alignments);
                
                //
                if(Verbosity::Instance().getPrintLevel() > 0)
                {
                    std::cout << "Dindel returned " << drc << "\n";
                    std::cout << "base vcf records:\n" << baseVCFSS.str() << "\n";
                    std::cout << "variant vcf records:\n" << variantVCFSS.str() << "\n";
                }
                
                // DINDEL ran without error, push its results to the output
                if(drc == DRC_OK)
                {                        
                    result.baseVCFStrings.push_back(baseVCFSS.str());
                    result.variantVCFStrings.push_back(variantVCFSS.str());
                    result.calledVCFStrings.push_back(callsVCFSS.str());

                    result.varStrings.insert(result.varStrings.end(), 
                                             build_result.variant_haplotypes.begin(), build_result.variant_haplotypes.end());
                
                    result.projectedReadAlignments = alignments;
                }
            }
        }
    }
    
    return result;
}

void GraphCompare::updateSharedStats(GraphCompareAggregateResults* pSharedStats)
{
    pSharedStats->updateShared(m_stats);
    m_stats.clear();
}

//
GraphBuildResult GraphCompare::processVariantKmer(const std::string& str, int /*count*/) const
{
    PROFILE_FUNC("GraphCompare::processVariantKmer")

    //
    GraphBuildResult result;

    if(m_parameters.algorithm == GCA_DEBRUIJN_GRAPH)
    {
        DeBruijnHaplotypeBuilder dbg_builder(m_parameters);
        dbg_builder.setInitialHaplotype(str);
        dbg_builder.run(result.variant_haplotypes);
    }
    else if(m_parameters.algorithm == GCA_PAIRED_DEBRUIJN_GRAPH)
    {
        PairedDeBruijnHaplotypeBuilder dbg_builder(m_parameters);
        dbg_builder.setInitialHaplotype(str);
        dbg_builder.run(result.variant_haplotypes);
    }
    else if(m_parameters.algorithm == GCA_STRING_GRAPH)
    {
        //OverlapHaplotypeBuilder overlap_builder(m_parameters);
        StringHaplotypeBuilder overlap_builder(m_parameters);
        overlap_builder.setInitialHaplotype(str);
        overlap_builder.run(result.variant_haplotypes);
    }

    //
    // Haplotype QC
    //
    size_t num_assembled = result.variant_haplotypes.size();
    if(num_assembled > m_parameters.maxHaplotypes)
    {
        result.variant_haplotypes.clear();
        return result;
    }

    qcVariantHaplotypes(m_parameters.bReferenceMode, result.variant_haplotypes);
    size_t num_qc = result.variant_haplotypes.size();

    // If any assembled haplotypes failed QC, do not try to call variants
    if(num_qc < num_assembled)
        result.variant_haplotypes.clear();

    if(!result.variant_haplotypes.empty())
        buildParallelBaseHaplotypes(result.variant_haplotypes, result.base_haplotypes);
    
    return result;
}

// Perform quality checks on each haplotype
void GraphCompare::qcVariantHaplotypes(bool bReferenceMode, StringVector& variant_haplotypes) const
{
    PROFILE_FUNC("GraphCompare::qcVariantHaplotypes")
    if(!bReferenceMode) 
    {
        // If in comparative mode, check whether the haplotype has sparse coverage by the base
        for(size_t i = 0; i < variant_haplotypes.size(); ++i) 
        {
            std::string& haplotype = variant_haplotypes[i];
            size_t k = m_parameters.kmer;
            size_t n = haplotype.size();
            size_t nk = haplotype.size() - k + 1;
            size_t num_covered_kmers = 0;
            std::vector<size_t> start_positions;
            BitVector covered_bases(n);
            std::stringstream cov_ss;

            for(size_t j = 0; j < nk; ++j) 
            {
                size_t c = BWTAlgorithms::countSequenceOccurrences(haplotype.substr(j, k), m_parameters.baseIndex);
                if(c > 0) 
                {
                    num_covered_kmers += 1;
                    for(size_t m = j; m < j + k; ++m)
                        covered_bases.set(m, true);
                    start_positions.push_back(j);
                }

                if(c > 9)
                    cov_ss << 9;
                else
                    cov_ss << c;
            }

            size_t num_uncovered_kmers = nk - num_covered_kmers;
            size_t num_uncovered_bases = 0;
            size_t max_d = std::numeric_limits<size_t>::min();
            size_t sum_d = 0;
            if(!start_positions.empty()) 
            {
                size_t prev_start = start_positions[0];
                for(size_t j = 1; j < start_positions.size(); ++j)
                {
                    size_t d = start_positions[j] - prev_start;
                    if(d > max_d)
                        max_d = d;
                    sum_d += d - 1;
                    prev_start = start_positions[j];
                }
            }

            for(size_t j = 0; j < n; j++)
                num_uncovered_bases += covered_bases.test(j) ? 0 : 1;

            //if(max_d < k / 2)
            //    haplotype = "";

            if(Verbosity::Instance().getPrintLevel() > 2)
            {
                printf("HaplotypeQC for hap[%zu] NK %zu NC %zu NUC %zu NUB %zu MAX_D: %zu SUM_D: %zu\n", 
                    i, nk, num_covered_kmers, num_uncovered_kmers, num_uncovered_bases, max_d, sum_d);
                printf("HaplotypeQC for hap[%zu] %s\n", i, cov_ss.str().c_str());
            }
        }

        // Calculate the maximum k such that every kmer is present in the variant and base BWT
        // The difference between these values must be at least MIN_COVER_K_DIFF
        size_t MIN_COVER_K_DIFF = 10;
        size_t MIN_COVERAGE = bReferenceMode ? 1 : 2;
        for(size_t i = 0; i < variant_haplotypes.size(); ++i)
        {
            if(variant_haplotypes[i] == "")
                continue;

            // Calculate the largest k such that the haplotype is walk through a de Bruijn graph of this k
            size_t max_variant_k = HapgenUtil::calculateMaxCoveringK(variant_haplotypes[i], MIN_COVERAGE, m_parameters.variantIndex);
            size_t max_base_k = HapgenUtil::calculateMaxCoveringK(variant_haplotypes[i], MIN_COVERAGE, m_parameters.baseIndex);
            if(Verbosity::Instance().getPrintLevel() > 2)
                printf("HaplotypeQC hap[%zu] MVK: %zu MBK: %zu\n", i, max_variant_k, max_base_k);

            //
            if( ! (max_variant_k > max_base_k && max_variant_k - max_base_k >= MIN_COVER_K_DIFF && max_base_k < 41))
                variant_haplotypes[i] = "";
        }

    }
    // Copy over the remaining haplotypes
    StringVector tmp_haplotypes;
    for(size_t i = 0; i < variant_haplotypes.size(); ++i) {
        if(!variant_haplotypes[i].empty())
            tmp_haplotypes.push_back(variant_haplotypes[i]);

        if(Verbosity::Instance().getPrintLevel() > 2 && variant_haplotypes[i].empty())
            printf("HaplotypeQC hap[%zu] failed QC\n", i);
    }
    variant_haplotypes.swap(tmp_haplotypes);
}

// Update the bloom filter with the kmers that were assembled into str
void GraphCompare::markVariantSequenceKmers(const std::string& str) const
{
    assert(str.size() >= m_parameters.kmer);
    size_t n = str.size() - m_parameters.kmer + 1;

    for(size_t i = 0; i < n; ++i)
    {
        std::string kmer = str.substr(i, m_parameters.kmer);
        std::string rc_kmer = reverseComplement(kmer);
        std::string& key_kmer = kmer < rc_kmer ? kmer : rc_kmer;
        m_parameters.pBloomFilter->add(key_kmer.c_str(), key_kmer.size());
    }
}

//
void GraphCompare::buildParallelBaseHaplotypes(const StringVector& variant_haplotypes,
                                               StringVector& base_haplotypes) const
{
    size_t haplotype_builder_kmer = m_parameters.kmer;
    // Run haplotype builder on the normal graph
    for(size_t i = 0; i < variant_haplotypes.size(); ++i)
    {
        const std::string& current_variant_haplotype = variant_haplotypes[i];

        std::string startAnchorSeq = current_variant_haplotype.substr(0, haplotype_builder_kmer); 
        std::string endAnchorSeq = current_variant_haplotype.substr(current_variant_haplotype.length() - haplotype_builder_kmer);

        if(startAnchorSeq == endAnchorSeq)
            continue; // degenerate sequence, skip

        AnchorSequence startAnchor;
        startAnchor.sequence = startAnchorSeq;
        startAnchor.count = 0;
        startAnchor.position = 0;

        AnchorSequence endAnchor;
        endAnchor.sequence = endAnchorSeq;
        endAnchor.count = 0;
        endAnchor.position = 0;

        HaplotypeBuilder builder;
        builder.setTerminals(startAnchor, endAnchor);
        builder.setIndex(m_parameters.baseIndex.pBWT, NULL);
        builder.setKmerParameters(haplotype_builder_kmer, m_parameters.bReferenceMode ? 1 : 2);

        // Run the builder
        HaplotypeBuilderReturnCode hbCode = builder.run();
        HaplotypeBuilderResult hbResult;

        // The search was successful, build strings from the walks
        if(hbCode == HBRC_OK) 
        {
            hbCode = builder.parseWalks(hbResult);
            if(hbCode == HBRC_OK) 
                base_haplotypes.insert(base_haplotypes.end(), hbResult.haplotypes.begin(), hbResult.haplotypes.end());
        }
    }
}

//
size_t GraphCompare::calculateHaplotypeBranches(const std::string& sequence, 
                                                size_t k, 
                                                size_t min_branch_depth, 
                                                const BWTIndexSet& indices)
{
    if(sequence.size() < k)
        return 0;

    size_t num_branches = 0;
    size_t nk = sequence.size() - k + 1;
    for(size_t i = 0; i < nk; ++i)
    {
        std::string kmer = sequence.substr(i, k);
        AlphaCount64 extensions = BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(kmer, indices.pBWT, ED_SENSE, indices.pCache);

        // Count number of symbols with coverage >= min_branch_depth
        size_t n = 0;
        for(int j = 0; j < DNA_ALPHABET::size; ++j)
        {
            char b = DNA_ALPHABET::getBase(j);
            if(extensions.get(b) >= min_branch_depth)
                ++n;
        }

        if(n > 0)
            num_branches += (n - 1);
    }

    return num_branches;
}

// test kmers from a file
void GraphCompare::testKmersFromFile(const std::string& kmerFilename)
{
    std::cout << "Kmer file: " << kmerFilename << "\n";
    std::istream* pReader = createReader(kmerFilename);

    std::string line;
    while(getline(*pReader, line))
    {
        std::cout << "Processing " << line << "\n";
        StringVector fields = split(line, '\t');
        std::string kmer = fields[0];

        if(kmer.size() != m_parameters.kmer)
        {
            std::cout << "Incorrect kmer size: " << kmer << " kmer size: " << kmer.size() << " kmer length parameter: " << m_parameters.kmer << "\n";
            continue;
        }
        testKmer(kmer);
    }
}

void GraphCompare::testKmer(const std::string& kmer)
{
    int count = 1;
    GraphBuildResult build_result = processVariantKmer(kmer, count);

    if(build_result.variant_haplotypes.size() > 0 && build_result.base_haplotypes.size() > 0)
    {
        std::cout << "Haplotypes successfully built. Aligning first pair.\n";
        StdAlnTools::globalAlignment(build_result.base_haplotypes.front(), build_result.variant_haplotypes.front(), true);

        std::stringstream baseVCFSS;
        std::stringstream variantVCFSS;
        std::stringstream callsVCFSS;

        DindelReturnCode drc = DindelUtil::runDindelPairMatePair(kmer,
                                                                 build_result.base_haplotypes,
                                                                 build_result.variant_haplotypes,
                                                                 m_parameters,
                                                                 baseVCFSS,
                                                                 variantVCFSS,
                                                                 callsVCFSS,
                                                                 NULL);
        
        std::cout << "base:    " << baseVCFSS.str() << "\n";
        std::cout << "variant: " << variantVCFSS.str() << "\n";
        
        if(drc == DRC_OK)
            std::cout << "DINDEL says: OK.\n";
        else
            std::cout << "DINDEL error: " << drc <<"\n";
    } 
    else 
    {
        
        std::cout << "Error. Not enough haplotypes VH:" << build_result.variant_haplotypes.size() << " BH: " << build_result.base_haplotypes.size() << "\n";
    }
}

void GraphCompare::showMappingLocations(const std::string& str)
{
    BWTInterval ref_interval = BWTAlgorithms::findInterval(m_parameters.referenceIndex, str);
    for(int64_t i = ref_interval.lower; i <= ref_interval.upper; ++i)
        std::cout << "Location: " << m_parameters.referenceIndex.pSSA->calcSA(i, m_parameters.referenceIndex.pBWT) << "\n";
}

//
// GraphCompareAggregateResult
//
GraphCompareAggregateResults::GraphCompareAggregateResults(const std::string& fileprefix, 
                                                           const StringVector& samples,
                                                           const ReadTable& refTable,
                                                           const std::string& referenceFilename) : m_baseVCFFile(fileprefix + ".base.vcf","w"),
                                                                                                   m_variantVCFFile(fileprefix + ".variant.vcf","w"),
                                                                                                   m_callsVCFFile(fileprefix + ".calls.vcf","w"),
                                                                                                   m_numVariants(0)
{
    //
    m_pWriter = createWriter(fileprefix + ".strings.fa");
    
    // Initialize samples
    m_baseVCFFile.setSamples(samples);
    m_variantVCFFile.setSamples(samples);
    m_callsVCFFile.setSamples(samples);

    // Initialize VCF
    m_baseVCFFile.outputHeader(referenceFilename, "");
    m_variantVCFFile.outputHeader(referenceFilename, "");
    m_callsVCFFile.outputHeader(referenceFilename, "");

    // Initialize BAM
    BamTools::SamHeader null_header;

    // Build a RefVector for bamtools
    BamTools::RefVector ref_vector;
    for(size_t i = 0; i < refTable.getCount(); ++i)
    {
        const SeqItem& ref_item = refTable.getRead(i);
        BamTools::RefData rd(ref_item.id, ref_item.seq.length());
        ref_vector.push_back(rd);

        m_refNameToIndexMap.insert(std::make_pair(ref_item.id, i));
    }

    m_evidenceBamFile.Open(fileprefix + ".evidence.bam", null_header, ref_vector);

    // Initialize mutex
    int ret = pthread_mutex_init(&m_mutex, NULL);
    if(ret != 0)
    {
        std::cerr << "Mutex initialization failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }
}

//
GraphCompareAggregateResults::~GraphCompareAggregateResults()
{
    delete m_pWriter;

    int ret = pthread_mutex_destroy(&m_mutex);
    if(ret != 0)
    {
        std::cerr << "Mutex destruction failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }

    m_evidenceBamFile.Close();
}

//
void GraphCompareAggregateResults::updateShared(const GraphCompareStats stats)
{
    pthread_mutex_lock(&m_mutex);
    m_stats.add(stats);
    pthread_mutex_unlock(&m_mutex);
}

void GraphCompareAggregateResults::process(const SequenceWorkItem& /*item*/, const GraphCompareResult& result)
{
    for(size_t i = 0; i < result.varStrings.size(); ++i)
    {
        // Write to the variants file
        std::stringstream varIDMaker;
        varIDMaker << "variant-" << m_numVariants;
        SeqItem item2 = { varIDMaker.str(), result.varStrings[i] };
        item2.write(*m_pWriter);
        m_numVariants += 1;
    }

    assert(result.baseVCFStrings.size() == result.variantVCFStrings.size());
    for(size_t i = 0; i < result.baseVCFStrings.size(); ++i)
    {
        m_baseVCFFile.getOutputStream() << result.baseVCFStrings[i];
        m_variantVCFFile.getOutputStream() << result.variantVCFStrings[i];
    }

    // Write out the final calls
    for(size_t i = 0; i < result.calledVCFStrings.size(); ++i)
        m_callsVCFFile.getOutputStream() << result.calledVCFStrings[i];

    // Write the read-to-reference alignment to the bam file
    for(size_t i = 0; i < result.projectedReadAlignments.size(); ++i)
    {
        DindelReadReferenceAlignment drra = result.projectedReadAlignments[i];
        
        // Set core data
        BamTools::BamAlignment bam_align;
        bam_align.Name = drra.read_name;

        // In BAM, we must feed in the reverse complement of the read when we set the IsReverseStrand flag
        bam_align.QueryBases = drra.is_reference_reverse_strand ? reverseComplement(drra.read_sequence) 
                                                                : drra.read_sequence;
        bam_align.Length = drra.read_sequence.size();
        bam_align.Qualities = "*";
        bam_align.MapQuality = 255;

        // Find reference index
        std::map<std::string, size_t>::iterator ref_iter = m_refNameToIndexMap.find(drra.reference_name);
        assert(ref_iter != m_refNameToIndexMap.end());
        bam_align.RefID = ref_iter->second;
        bam_align.Position = drra.reference_start_position - 1;

        // Build CigarOp
        std::stringstream parser(drra.cigar);
        int length;
        char symbol;
        while(parser >> length >> symbol)
            bam_align.CigarData.push_back(BamTools::CigarOp(symbol, length));

        // Set flags
        bam_align.SetIsMapped(true);
        bam_align.SetIsPrimaryAlignment(true);
        bam_align.SetIsReverseStrand(drra.is_reference_reverse_strand);
        bam_align.SetIsProperPair(false);
        bam_align.SetIsFailedQC(false);

        // write
        m_evidenceBamFile.SaveAlignment(bam_align);
    }
}

//
void GraphCompareAggregateResults::printStats() const
{
    m_stats.print();
}
