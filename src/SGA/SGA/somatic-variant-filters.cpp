//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// variant-filtering - apply various filters
// to a somatic VCF file
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "somatic-variant-filters.h"
#include "BWTAlgorithms.h"
#include "BWTIndexSet.h"
#include "SGACommon.h"
#include "Timer.h"
#include "SequenceProcessFramework.h"
#include "HapgenUtil.h"
#include "SGAStats.h"
#include "HashMap.h"
#include "DindelRealignWindow.h"
#include "VariantIndex.h"
#include "HapgenUtil.h"
#include "MultiAlignment.h"
#include "api/BamReader.h"

// Types
typedef HashMap<std::string, std::string> StringStringHash;

//
// Getopt
//
#define SUBPROGRAM "somatic-variant-filters"
static const char *SOMATIC_VARIANT_FILTERS_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2014 Ontario Institute for Cancer Research\n";

static const char *SOMATIC_VARIANT_FILTERS_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... VCF_FILE\n"
"Count read support for variants in a vcf file\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -t, --threads=NUM                use NUM threads to compute the overlaps (default: 1)\n"
"          --reference=STR              load the reference genome from FILE\n"
"          --tumor-bam=STR              load the aligned tumor reads from FILE\n"
"          --normal-bam=STR             load the aligned normal reads from FILE\n"
"          --annotate-only              only annotate with INFO tags, rather than filter\n"
"\n"
"Filtering cutoffs:\n"
"          --min-af=FLOAT               minimum allele frequency (AF tag)\n"
"          --min-var-dp=INT             minimum variant depth (VarDP)\n"
"          --max-strand-bias=FLOAT      maximum strand bias (SB)\n"
"          --max-hp-length=INT          maximum homopolymer context length (HPlen)\n"
"          --max-dust=FLOAT             maximum dust scores (Dust)\n"
"          --max-normal-reads=INT       maximum normal reads showing variant\n"
"          --min-normal-depth=INT       minimum normal depth at the locus\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static std::string vcfFile;
    static std::string referenceFile;
    static std::string tumorBamFile;
    static std::string normalBamFile;

    static int numThreads = 1;
    static double minAF = 0.0f;
    static double errorBound = 0.02;
    static int minVarDP = 0;
    static double maxStrandBias = 2.0f;
    static int maxHPLen = 7;
    static double maxDust = 2.0f;
    static size_t maxNormalReads = 1;
    static size_t minNormalDepth = 5;
    static int minMedianQuality = 15;
    static size_t maximumAlignments = 500;
    static double minHaplotypeLength = 125.0f;
    static bool annotateOnly = false;
}

static const char* shortopts = "t:v";

enum { OPT_HELP = 1, 
       OPT_VERSION, 
       OPT_REFERENCE, 
       OPT_TUMOR_BAM, 
       OPT_NORMAL_BAM,
       OPT_ANNOTATE,
       OPT_MIN_AF,
       OPT_MIN_VARDP,
       OPT_MAX_SB,
       OPT_MAX_HP,
       OPT_MAX_DUST,
       OPT_MAX_NORMAL_READS,
       OPT_MIN_NORMAL_DEPTH };

static const struct option longopts[] = {
    { "verbose",               no_argument,       NULL, 'v' },
    { "annotate-only",         no_argument,       NULL, OPT_ANNOTATE },
    { "threads",               required_argument, NULL, 't' },
    { "reference",             required_argument, NULL, OPT_REFERENCE },
    { "tumor-bam",             required_argument, NULL, OPT_TUMOR_BAM },
    { "normal-bam",            required_argument, NULL, OPT_NORMAL_BAM },
    { "min-af",                required_argument, NULL, OPT_MIN_AF },
    { "min-var-dp",            required_argument, NULL, OPT_MIN_VARDP },
    { "max-strand-bias",       required_argument, NULL, OPT_MAX_SB },
    { "max-hp-length",         required_argument, NULL, OPT_MAX_HP },
    { "max-dust",              required_argument, NULL, OPT_MAX_DUST },
    { "max-normal-reads",      required_argument, NULL, OPT_MAX_NORMAL_READS },
    { "min-normal-depth",      required_argument, NULL, OPT_MIN_NORMAL_DEPTH },
    { "help",                  required_argument, NULL, OPT_HELP },
    { "version",               required_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

// Calculate the median value of the vector
template<typename T>
double median(const std::vector<T> v)
{
    std::vector<T> c = v;
    std::sort(c.begin(), c.end());
    double r = 0.0f;
    int m = c.size() / 2;
    if(c.size() % 2 == 1)
        r = c[m];
    else
        r = (c[m - 1] + c[m]) / 2.0f;
    return r;
}

struct CoverageStats
{
    CoverageStats() : n_total_reads(0), n_evidence_reads(0), too_many_alignments(false) {} 

    // functions
    double calculateVAF() const
    {
        return n_total_reads > 0 ? 
            n_evidence_reads / (double)n_total_reads : 0.0f;
    }

    //
    size_t n_total_reads;
    size_t n_evidence_reads;
    double median_mapping_quality;
    std::vector<int> snv_evidence_quals;
    bool too_many_alignments;
};

// This struct holds segments
// of a read that match a variant
// The read is broken into (up to) 3 pieces
//  -the segment that aligns to the reference before the variant
//  -the segment that aligns to the variant region
//  -the segment that aligns after the variant
struct VariantReadSegments
{
    std::string preSegment;
    std::string variantSegment;
    std::string postSegment;

    std::string preQual;
    std::string variantQual;
    std::string postQual;
};

VariantReadSegments splitReadAtVariant(const BamTools::BamAlignment& alignment,
                                       const VCFRecord& record)
{
    // Set up a vector of aligned reference positions
    // for every base of the read
    std::vector<int> reference_positions(alignment.QueryBases.size(), -1);
    size_t read_idx = 0;
    size_t ref_idx = alignment.Position;

    std::stringstream cigar;
    // 
    for(size_t ci = 0; ci < alignment.CigarData.size(); ++ci)
    {
        BamTools::CigarOp op = alignment.CigarData[ci];
        switch(op.Type) {
            case 'M':
                for(size_t l = 0; l < op.Length; ++l)
                    reference_positions[read_idx++] = ref_idx++;
                break;
            case 'D':
            case 'N':
                ref_idx += op.Length;
                break;
            case 'I':
            case 'S':
                read_idx += op.Length;
                break;
        }

        cigar << op.Length << op.Type;
    }

    // Find the last aligned position before the variant and the first position after
    int rightmost_before_variant = -1;
    int leftmost_after_variant = std::numeric_limits<int>::max();

    int variant_start = record.refPosition - 1; // VCF is 1-based, bamtools is 0-based
    int variant_end = variant_start + record.refStr.size();

    for(int i = 0; i < (int)reference_positions.size(); ++i) {
        if(reference_positions[i] != -1 && reference_positions[i] < variant_start && i > rightmost_before_variant)
            rightmost_before_variant = i;
        if(reference_positions[i] != -1 && reference_positions[i] >= variant_end && i < leftmost_after_variant)
            leftmost_after_variant = i;
    }

    VariantReadSegments ret;
    ret.preSegment = alignment.QueryBases.substr(0, rightmost_before_variant + 1);
    ret.preQual = alignment.Qualities.substr(0, rightmost_before_variant + 1);

    ret.variantSegment = alignment.QueryBases.substr(rightmost_before_variant + 1, leftmost_after_variant - rightmost_before_variant - 1);
    ret.variantQual = alignment.Qualities.substr(rightmost_before_variant + 1, leftmost_after_variant - rightmost_before_variant - 1);

    if(leftmost_after_variant < (int)alignment.QueryBases.size())
    {
        ret.postSegment = alignment.QueryBases.substr(leftmost_after_variant);
        ret.postQual = alignment.Qualities.substr(leftmost_after_variant);
    }

    return ret;
}

//
//
//
CoverageStats getVariantCoverage(BamTools::BamReader* pReader, const VCFRecord& record, const ReadTable* refTable)
{
    CoverageStats stats;
    
    static const int flankingSize = 100;
    static const double minPercentIdentity = 95.0f;

    bool is_snv = record.refStr.size() == 1 && record.varStr.size() == 1;

    // Grab the reference haplotype
    int zeroBasedPos = record.refPosition - 1;
    int start = zeroBasedPos - flankingSize - 1;
    if(start < 0)
        start = 0;

    int end = zeroBasedPos + record.refStr.length() + 2 * flankingSize;
    const SeqItem& chr = refTable->getRead(record.refName);
    if(end > (int)chr.seq.length())
        end = (int)chr.seq.length();

    std::string reference_haplotype = chr.seq.substr(start, end - start);
    int translatedPos = zeroBasedPos - start;

    std::string variant_haplotype = reference_haplotype;
    
    // Ensure that the reference string at the variant matches the expected
    if(variant_haplotype.substr(translatedPos, record.refStr.length()) != record.refStr)
    {
        std::cerr << "Warning, reference base does not match VCF record.\n";
        std::cerr << record << "\n";
    }
    variant_haplotype.replace(translatedPos, record.refStr.length(), record.varStr);

    // Grab all reads in reference region
    int refID = pReader->GetReferenceID(record.refName);
    if(refID < 0)
        return stats;

    int refStart = record.refPosition;
    int refEnd = record.refPosition;
    pReader->SetRegion(refID, refStart, refID, refEnd);
    BamTools::BamAlignment aln;

    std::vector<double> mapping_quality;
    std::vector<BamTools::BamAlignment> alignments;
    while(pReader->GetNextAlignment(aln)) {
        if(aln.MapQuality > 0)
            alignments.push_back(aln);
        mapping_quality.push_back(aln.MapQuality);

        if(alignments.size() > opt::maximumAlignments)
        {
            stats.too_many_alignments = true;
            return stats;
        }
    }

    if(!mapping_quality.empty())
        stats.median_mapping_quality = median(mapping_quality);
    else
        stats.median_mapping_quality = 60;

    // Shuffle and take the first N alignments only
    std::random_shuffle(alignments.begin(), alignments.end());

#if HAVE_OPENMP
    #pragma omp parallel for
#endif    
    for(size_t i = 0; i < alignments.size(); ++i) {
        BamTools::BamAlignment alignment = alignments[i];

        VariantReadSegments segments = splitReadAtVariant(alignment, record);

        if(opt::verbose > 1)
        {
            fprintf(stderr, "var: %zu %s -> %s\n",  record.refPosition, record.refStr.c_str(), record.varStr.c_str());
            fprintf(stderr, "pos: %d\n",  alignment.Position);
            fprintf(stderr, "strand: %s\n", alignment.IsReverseStrand() ? "-" : "+");
            fprintf(stderr, "read: %s\n", alignment.QueryBases.c_str());
            fprintf(stderr, "qual: %s\n", alignment.Qualities.c_str());
            fprintf(stderr, "alnb: %s\n", alignment.AlignedBases.c_str());
            
            fprintf(stderr, "Pre: %s\n",  segments.preSegment.c_str());
            fprintf(stderr, "Var: %s\n",  segments.variantSegment.c_str());
            fprintf(stderr, "Pos: %s\n",  segments.postSegment.c_str());
            
            fprintf(stderr, "PreQual: %s\n",  segments.preQual.c_str());
            fprintf(stderr, "VarQual: %s\n",  segments.variantQual.c_str());
            fprintf(stderr, "PosQual: %s\n",  segments.postQual.c_str());
        }

        bool aligned_at_variant = segments.variantSegment.size() > 0 && 
                                  (segments.preSegment.size() > 0 || segments.postSegment.size() > 0);

        if(!aligned_at_variant)
            continue;
        
        bool is_evidence_read = false;
        if(segments.variantSegment != record.refStr)
        {
            if(segments.variantSegment == record.varStr)
            {
                // Evidence read via the current alignment
                is_evidence_read = true;
            }
            else
            {
                // Check for evidence via realignment to the variant haplotype
                SequenceOverlap ref_overlap = 
                    Overlapper::computeOverlapAffine(alignment.QueryBases, reference_haplotype);
                SequenceOverlap var_overlap = 
                    Overlapper::computeOverlapAffine(alignment.QueryBases, variant_haplotype);
                
                bool quality_alignment = (ref_overlap.getPercentIdentity() >= minPercentIdentity || 
                                          var_overlap.getPercentIdentity() >= minPercentIdentity);

                is_evidence_read = quality_alignment && var_overlap.score > ref_overlap.score;
            }
        }

        #pragma omp critical
        {
            stats.n_total_reads += 1;
            if(is_evidence_read)
            {
                stats.n_evidence_reads += 1;
                if(is_snv && segments.variantQual.size() == 1)
                {
                    char qb = segments.variantQual[0];
                    int q = Quality::char2phred(qb);
                    stats.snv_evidence_quals.push_back(q);
                }
            }
        }
    }

    return stats;
}

void makeTagHash(const VCFRecord& record, StringStringHash& tagHash)
{
    for(size_t i = 0; i < record.comments.size(); i++)
    {
        StringVector kv = split(record.comments[i], '=');
        assert(kv.size() == 2 || kv.size() == 1);
        if(kv.size() == 2)
            tagHash[kv[0]] = kv[1];
        else
            tagHash[kv[0]] = "true";
    }
}


template<typename T>
bool getTagValue(StringStringHash& tagHash, const std::string& key, T& out)
{
    StringStringHash::iterator iter = tagHash.find(key);
    if(iter == tagHash.end())
        return false; // tag not available
    std::stringstream parser(iter->second);
    parser >> out;
    return true;
}

std::string makeVariantKey(const VCFRecord& record)
{
    std::stringstream ss;
    ss << record.refName << ":" << record.refPosition << ":"
       << record.refStr << ":" << record.varStr;
    return ss.str(); 
}

// Calculate Homopolymer lengths at the variant position
int calculateHomopolymerLength(const VCFRecord& record, const ReadTable* refTable)
{
    static const int flankingSize = 100;

    // Grab the reference haplotype
    int zeroBasedPos = record.refPosition - 1;
    int start = zeroBasedPos - flankingSize - 1;
    if(start < 0)
        start = 0;

    int end = zeroBasedPos + record.refStr.length() + 2 * flankingSize;
    const SeqItem& chr = refTable->getRead(record.refName);
    if(end > (int)chr.seq.length())
        end = (int)chr.seq.length();

    std::string reference_haplotype = chr.seq.substr(start, end - start);

    // we use the MA class for the homopolymer counting code
    MAlignDataVector mav;
    MultiAlignment ma(reference_haplotype, mav);


    // Count run lengths of nucleotides
    std::vector<size_t> run_lengths;

    // Count runs
    char curr_char = reference_haplotype[0];
    size_t curr_run = 1;
    for(size_t i = 1; i < reference_haplotype.size(); ++i)
    {
        if(reference_haplotype[i] != curr_char)
        {
            run_lengths.push_back(curr_run);
            curr_run = 1;
            curr_char = reference_haplotype[i];
        }
        else
        {
            curr_run++;
        }
    }
    run_lengths.push_back(curr_run); // last run

    // set up a map from the reference base to the run length its in
    std::vector<size_t> index_to_run_lengths(reference_haplotype.size());
    size_t j = 0;
    for(size_t i = 0; i < run_lengths.size(); ++i)
    {
        size_t stop = j + run_lengths[i];
        while(j < stop)
        {
            index_to_run_lengths[j] = i;
            j++;
        }
    }
    assert(j == reference_haplotype.size());

    // Calculate the maximum homopolymer run as the max over the variant reference region
    size_t eventStart = zeroBasedPos - start;
    size_t eventEnd = eventStart + record.refStr.size();
    size_t maxhp = 0;
    for(size_t i = eventStart; i <= eventEnd; ++i)
    {
        size_t idx = index_to_run_lengths[i];
        if(run_lengths[idx] > maxhp)
            maxhp = run_lengths[idx];
    }
    return maxhp;
}

// Get the 3-mer preceding and following the variant
void getVariantContext(const VCFRecord& record, const ReadTable* refTable,
                       std::string& prefix, std::string& suffix)
{
    static const int flankingSize = 100;

    // Grab the reference haplotype
    int zeroBasedPos = record.refPosition - 1;
    int start = zeroBasedPos - flankingSize - 1;
    if(start < 0)
        start = 0;

    int end = zeroBasedPos + record.refStr.length() + 2 * flankingSize;
    const SeqItem& chr = refTable->getRead(record.refName);
    if(end > (int)chr.seq.length())
        end = (int)chr.seq.length();

    std::string reference_haplotype = chr.seq.substr(start, end - start);

    // Calculate the maximum homopolymer run as the max over the variant reference region
    size_t eventStart = zeroBasedPos - start;
    size_t eventEnd = eventStart + record.refStr.size();
    static const size_t k = 3;
    prefix = reference_haplotype.substr(eventStart - k, k);
    suffix = reference_haplotype.substr(eventEnd, k);
}

//
// Main
//
int somaticVariantFiltersMain(int argc, char** argv)
{
    parseSomaticVariantFiltersOptions(argc, argv);

    Timer* pTimer = new Timer(PROGRAM_IDENT);

#if HAVE_OPENMP
    omp_set_num_threads(opt::numThreads);
#endif 
    // Load Reference
    ReadTable refTable(opt::referenceFile, SRF_NO_VALIDATION);
    refTable.indexReadsByID();

    // Load BAMs
    BamTools::BamReader* pTumorBamReader = new BamTools::BamReader;
    pTumorBamReader->Open(opt::tumorBamFile);
    pTumorBamReader->LocateIndex();
    assert(pTumorBamReader->HasIndex());

    BamTools::BamReader* pNormalBamReader = new BamTools::BamReader;
    pNormalBamReader->Open(opt::normalBamFile);
    pNormalBamReader->LocateIndex();
    assert(pNormalBamReader->HasIndex());

    // Track duplicated variants
    HashSet<std::string> duplicateHash;

    std::istream* pInput = createReader(opt::vcfFile.c_str());
    std::string line;

    while(getline(*pInput, line))
    {
        if(line.empty())
            continue;

        if(line[0] == '#')
        {
            std::cout << line << "\n";
            continue;
        }
        
        // parse record
        VCFRecord record(line);
        if(record.isMultiAllelic())
        {
            std::cerr << "Error: multi-allelic VCF found, please run vcfbreakmulti\n";
            exit(EXIT_FAILURE);
        }

        // Check if we've seen this variant already
        std::string key = makeVariantKey(record);
        if(duplicateHash.find(key) != duplicateHash.end())
            continue;
        else
            duplicateHash.insert(key);

        if(opt::verbose > 0)
        {
            std::stringstream ss;
            ss << "Variant: " << record << "\n";
            fprintf(stderr, "===============================================\n%s", ss.str().c_str());
        }

        StringStringHash tagHash;
        makeTagHash(record, tagHash);

        StringVector fail_reasons;

        int hplen = 0;
        if(!getTagValue(tagHash, "HPLen", hplen))
            hplen = calculateHomopolymerLength(record, &refTable);
        if(hplen > opt::maxHPLen)
            fail_reasons.push_back("Homopolymer");

        double dust = 0.0f;
        if(!getTagValue(tagHash, "Dust", dust))
            dust = HapgenUtil::calculateDustScoreAtPosition(record.refName, 
                                                            record.refPosition, 
                                                            &refTable);

        if(dust > opt::maxDust)
            fail_reasons.push_back("LowComplexity");
        
        double af;
        if(getTagValue(tagHash, "AF", af) && af < opt::minAF)
            fail_reasons.push_back("LowAlleleFrequency");

        int varDP;
        if(getTagValue(tagHash, "VarDP", varDP) && varDP < opt::minVarDP)
            fail_reasons.push_back("LowVarDP");
        
        double avgHapLen;
        if(getTagValue(tagHash, "AvgHapLen", avgHapLen) && avgHapLen < opt::minHaplotypeLength)
            fail_reasons.push_back("ShortHaplotype");

        double strandBias;
        if(getTagValue(tagHash, "SB", strandBias) && strandBias >= opt::maxStrandBias)
            fail_reasons.push_back("StrandBias");

        CoverageStats tumor_stats = getVariantCoverage(pTumorBamReader, record, &refTable);
        CoverageStats normal_stats = getVariantCoverage(pNormalBamReader, record, &refTable);

        if(opt::verbose > 0)
        {
            fprintf(stderr, "Tumor: [%zu %zu]\n",  tumor_stats.n_total_reads, tumor_stats.n_evidence_reads);
            fprintf(stderr, "Normal: [%zu %zu]\n", normal_stats.n_total_reads, normal_stats.n_evidence_reads);
        }

        if(!tumor_stats.too_many_alignments && !normal_stats.too_many_alignments)
        {
            // Check that there is not evidence for the variant in the normal sample
            if(normal_stats.n_evidence_reads > opt::maxNormalReads)
                fail_reasons.push_back("NormalEvidence");
            
            // If the normal is poorly covered in this region, we may call a germline variant as a variant
            if(normal_stats.n_total_reads < opt::minNormalDepth)
                fail_reasons.push_back("LowNormalDepth");

            if(!tumor_stats.snv_evidence_quals.empty())
            {
                // Check that the base scores of SNVs are reasonably high
                double median_quality = median(tumor_stats.snv_evidence_quals);
                if(median_quality < opt::minMedianQuality)
                    fail_reasons.push_back("LowQuality");

                // For very deep regions, errors can be mistaken for SNVs
                // Check if the variant frequency is very low. This check is not performed
                // for indels and MNPs as they might not be aligned to this region.
                double vf_estimate = tumor_stats.n_evidence_reads / (double)tumor_stats.n_total_reads;
                if(vf_estimate < opt::errorBound)
                    fail_reasons.push_back("PossibleError");
            }
            
            // Check that the mapping quality of reads in the region is reasonable
            if(tumor_stats.median_mapping_quality < opt::minMedianQuality)
                fail_reasons.push_back("LowMappingQuality");
        }
        else 
        {
            // If the depth in the region is excessively high, the statistical models may break down
            fail_reasons.push_back("DepthLimitReached");
        }

        if(!fail_reasons.empty() && !opt::annotateOnly)
        {
            if(record.passStr != "PASS" && record.passStr != ".")
                fail_reasons.insert(fail_reasons.begin(), record.passStr);

            std::stringstream strss;
            std::copy(fail_reasons.begin(), fail_reasons.end(), std::ostream_iterator<std::string>(strss, ";"));
            record.passStr = strss.str();
            record.passStr.erase(record.passStr.size() - 1); // erase trailing ;
        }
        
        // Add INFO tags indicating allele coverage
        if(opt::annotateOnly)
        {
            double tumor_vf = tumor_stats.calculateVAF();
            double normal_vf = normal_stats.calculateVAF();

            std::string prefix;
            std::string suffix;
            getVariantContext(record, &refTable, prefix, suffix);

            record.addComment("TumorVAF", tumor_vf);
            record.addComment("NormalVAF", normal_vf);
            record.addComment("TumorVarDepth", (int)tumor_stats.n_evidence_reads);
            record.addComment("TumorTotalDepth", (int)tumor_stats.n_total_reads);
            record.addComment("NormalVarDepth", (int)normal_stats.n_evidence_reads);
            record.addComment("NormalTotalDepth", (int)normal_stats.n_total_reads);
            record.addComment("5pContext", prefix);
            record.addComment("3pContext", suffix);
        }

        std::cout << record << "\n";
    }
    
    // Cleanup
    delete pInput;
    delete pTumorBamReader;
    delete pNormalBamReader;
    delete pTimer;

    return 0;
}

// 
// Handle command line arguments
//
void parseSomaticVariantFiltersOptions(int argc, char** argv)
{
    std::string algo_str;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case OPT_REFERENCE: arg >> opt::referenceFile; break;
            case OPT_TUMOR_BAM: arg >> opt::tumorBamFile; break;
            case OPT_NORMAL_BAM: arg >> opt::normalBamFile; break;
            case OPT_MIN_AF: arg >> opt::minAF; break;
            case OPT_MIN_VARDP: arg >> opt::minVarDP; break;
            case OPT_MAX_SB: arg >> opt::maxStrandBias; break;
            case OPT_MAX_HP: arg >> opt::maxHPLen; break;
            case OPT_MAX_DUST: arg >> opt::maxDust; break;
            case OPT_MAX_NORMAL_READS: arg >> opt::maxNormalReads; break;
            case OPT_MIN_NORMAL_DEPTH: arg >> opt::minNormalDepth; break;
            case OPT_ANNOTATE: opt::annotateOnly = true; break;
            case 't': arg >> opt::numThreads; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << SOMATIC_VARIANT_FILTERS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << SOMATIC_VARIANT_FILTERS_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 1) 
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } 
    else if (argc - optind > 1) 
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if(opt::numThreads <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
        die = true;
    }

    if(opt::tumorBamFile.empty())
    {
        std::cerr << SUBPROGRAM ": a --tumor-bam must be provided\n";
        die = true;
    }
    
    if(opt::normalBamFile.empty())
    {
        std::cerr << SUBPROGRAM ": a --normal-bam must be provided\n";
        die = true;
    }

    if(opt::referenceFile.empty())
    {
        std::cerr << SUBPROGRAM ": a --reference must be provided\n";
        die = true;
    }

    opt::vcfFile = argv[optind++];

    if (die) 
    {
        std::cout << "\n" << SOMATIC_VARIANT_FILTERS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}
