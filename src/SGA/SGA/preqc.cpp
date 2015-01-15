//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// pre-QC - Perform pre-assembly quality checks on a set
//          of reads
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <queue>
#include "Util.h"
#include "preqc.h"
#include "Timer.h"
#include "BWT.h"
#include "BWTAlgorithms.h"
#include "SGACommon.h"
#include "HashMap.h"
#include "KmerDistribution.h"
#include "KmerOverlaps.h"
#include "BloomFilter.h"
#include "SGAStats.h"
#include "DindelRealignWindow.h"
#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/filestream.h"

#if HAVE_OPENMP
#include <omp.h>
#endif

// Enums
enum BranchClassification
{
    BC_NO_CALL,
    BC_ERROR,
    BC_VARIANT,
    BC_REPEAT
};

// Structs

// struct to store some estimated properties of the genome
struct GenomeEstimates
{
    double average_read_length;
    size_t genome_size;
    size_t total_reads;

    // the mean number of reads starting at a given reference position
    double mean_read_starts;
};

// struct to store the parameters used by our classifier
struct ModelParameters
{
    int k;

    // the modal k-mer depth inferred from the count distribution
    double mode;

    // the kmer count model is a mixture of poissons representing
    // the copy number of a k-mer. this vector gives the mixture
    // proportions:
    // [ error_proportion, het_proportion, diploid_prop, 1-repeat, 2-repeat, ...]
    static const size_t NUM_MIXTURE_STATES = 10;
    static const size_t MAX_REPEAT_COPY = NUM_MIXTURE_STATES - 1;

    double mixture_proportions[NUM_MIXTURE_STATES];
    double mixture_means[NUM_MIXTURE_STATES];

    // This vector contains normalized proportions
    // for the repeat states. The indexing is the same
    // as the above arrays - for non-repeat states the
    // value will be 0.
    double repeat_mixture_normalized[NUM_MIXTURE_STATES];

    // prior on the per-base error rate
    static const double error_rate_prior;
};
const double ModelParameters::error_rate_prior = 0.02;

// struct to store posterior probabilities calculated
// by the model
struct ModelPosteriors
{
    ModelPosteriors() : posterior_error(0.0f), 
                        posterior_variant(0.0f), 
                        posterior_repeat(0.0f),
                        classification(BC_NO_CALL) {}

    double posterior_error;
    double posterior_variant;
    double posterior_repeat;
    BranchClassification classification;
};

// struct to store a description of the neighbors of a vertex in the graph
struct KmerNeighbors
{
    // This stores the occurrence count for each neighbor
    // including both strands.
    AlphaCount64 total_count;

    // Same as above but restricted to the kmers that have
    // coverage on both strands
    AlphaCount64 count_both_strands;

    // This stores the bases that have coverage on both strands
    std::string extensions_both_strands;

    // Return an ordering of the 4 possible extensions by their count
    // For instance A:2 C:5 G:8 T:0 will return "GCAT"
    static std::string getExtensionsFromCount(const AlphaCount64& a)
    {
        char sorted_bases[5] = "ACGT";
        sorted_bases[4] = '\0';
        a.getSorted(sorted_bases, 5);
        std::string r(sorted_bases);
        return r;
    }
};


// Typedefs
typedef rapidjson::PrettyWriter<rapidjson::FileStream> JSONWriter;

// Functions

// Compute the distribution of unipath lengths from the k-de Bruijn graph.
// A branch in the graph with coverage c will be ignored when 
//    c / max_c < coverage_ratio_threshold 
//    where max_c is the highest-coverage branch
// 
void unipath_length_distribution(JSONWriter* pWriter,
                                 const BWTIndexSet& index_set,
                                 size_t k,
                                 double coverage_ratio_threshold, 
                                 size_t n_samples);

//
// Getopt
//
#define SUBPROGRAM "preqc"
static const char *PREQC_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2013 Wellcome Trust Sanger Institute\n";

static const char *PREQC_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Perform pre-assembly quality checks\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -t, --threads=NUM                use NUM threads (default: 1)\n"
"          --simple                     only compute the metrics that do not need the FM-index\n"
"          --max-contig-length=N        stop contig extension at N bp (default: 50000)\n"
"          --reference=FILE             use the reference FILE to calculate GC plot\n"
"          --diploid-reference-mode     generate metrics assuming that the input data\n"
"                                       is a reference genome, not a collection of reads\n"
"          --force-EM                   force preqc to proceed even if the coverage model\n"
"                                       does not converge. This allows the rest of the program to continue\n"
"                                       but the branch and genome size estimates may be misleading\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static size_t maxContigLength = 50000;
    static size_t kmerDistributionSamples = 50000;
    static std::string prefix;
    static std::string readsFile;
    static std::string referenceFile;
    static int diploidReferenceMode = 0;
    static bool forceEM = false;
    static bool simple = false;
}

static const char* shortopts = "p:d:t:o:k:n:b:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_REFERENCE, OPT_MAX_CONTIG, OPT_DIPLOID, OPT_FORCE_EM, OPT_SIMPLE };

static const struct option longopts[] = {
    { "verbose",                no_argument,       NULL, 'v' },
    { "threads",                required_argument, NULL, 't' },
    { "prefix",                 required_argument, NULL, 'p' },
    { "sample-rate",            required_argument, NULL, 'd' },
    { "kmer-size",              required_argument, NULL, 'k' },
    { "num-reads",              required_argument, NULL, 'n' },
    { "branch-cutoff",          required_argument, NULL, 'b' },
    { "max-contig-length",      required_argument, NULL, OPT_MAX_CONTIG },
    { "reference",              required_argument, NULL, OPT_REFERENCE },
    { "simple",                 no_argument,       NULL, OPT_SIMPLE },
    { "force-EM",               no_argument,       NULL, OPT_FORCE_EM },
    { "diploid-reference-mode", no_argument,       NULL, OPT_DIPLOID },
    { "help",                   no_argument,       NULL, OPT_HELP },
    { "version",                no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

// Returns the valid SUFFIX neighbors of the kmer
// in the de Bruijn graph. The neighbors are encoded
// as the single base they add. A coverage threshold
// is applied to filter out low-coverage extensions.
std::string get_valid_dbg_neighbors_ratio(const std::string& kmer,
                                          const BWTIndexSet& index_set,
                                          double coverage_ratio_threshold)
{
    std::string out;
    AlphaCount64 counts = 
        BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(kmer, 
                                                              index_set.pBWT,
                                                              ED_SENSE,
                                                              index_set.pCache);
    
    if(!counts.hasDNAChar())
        return out; // no extensions

    char max_b = counts.getMaxDNABase();
    size_t max_c = counts.get(max_b);

    for(size_t j = 0; j < 4; ++j)
    {
        char b = "ACGT"[j];
        size_t c = counts.get(b);
        if(c > 0 && (double)c / max_c >= coverage_ratio_threshold)
            out.push_back(b);
    }
    return out;
}

std::string get_valid_dbg_neighbors_coverage_and_ratio(const std::string& kmer,
                                                       const BWTIndexSet& index_set,
                                                       size_t min_coverage,
                                                       double min_ratio,
                                                       EdgeDir dir)
{
    std::string out;
    AlphaCount64 counts = 
        BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(kmer, 
                                                              index_set.pBWT,
                                                              dir,
                                                              index_set.pCache);
    
    if(!counts.hasDNAChar())
        return out; // no extensions

    char max_b = counts.getMaxDNABase();
    size_t max_c = counts.get(max_b);

    for(size_t j = 0; j < 4; ++j)
    {
        char b = "ACGT"[j];
        size_t c = counts.get(b);
        if(c >= min_coverage && (double)c / max_c >= min_ratio)
            out.push_back(b);
    }
    return out;    
}

void fill_neighbor_count_by_strand(const std::string& kmer, const BWTIndexSet& index_set, int f_counts[5], int r_counts[5])
{
    // Count the extensions from this kmer
    char f_ext[] = "ACGT";
    char r_ext[] = "TGCA";
    size_t k = kmer.size();

    std::string f_mer = kmer.substr(1) + "A";
    std::string r_mer = "A" + reverseComplement(kmer).substr(0, k - 1);
    assert(f_mer.size() == k);
    assert(r_mer.size() == k);

    // fill in the counts for the 4 dna bases
    for(size_t bi = 0; bi < 4; ++bi)
    {
        f_mer[f_mer.size() - 1] = f_ext[bi];
        r_mer[0] = r_ext[bi];

        f_counts[bi] = BWTAlgorithms::countSequenceOccurrencesSingleStrand(f_mer, index_set);
        r_counts[bi] = BWTAlgorithms::countSequenceOccurrencesSingleStrand(r_mer, index_set);
    }

    // fill in the count for the number of reads that terminated with kmer
    // at the end (or start of the reverse complement
    f_mer = kmer + "$";
    r_mer = "$" + reverseComplement(kmer);
    f_counts[4] = BWTAlgorithms::countSequenceOccurrencesSingleStrand(f_mer, index_set);
    r_counts[4] = BWTAlgorithms::countSequenceOccurrencesSingleStrand(r_mer, index_set);
}

// Get a description of the neighbors of the given kmer
// This is used in the classifier
KmerNeighbors calculate_neighbor_data(const std::string& kmer, const BWTIndexSet& index_set)
{
    KmerNeighbors out;

    // these vectors are in order ACGT$ on the forward strand
    int f_counts[5] = { 0, 0, 0, 0, 0 };
    int r_counts[5] = { 0, 0, 0, 0, 0 };

    fill_neighbor_count_by_strand(kmer, index_set, f_counts, r_counts);

    // Make sure both strands are represented for every branch
    for(size_t bi = 0; bi < 4; ++bi)
    {
        out.total_count.set("ACGT"[bi], f_counts[bi] + r_counts[bi]);
        if(f_counts[bi] >= 1 && r_counts[bi] >= 1)
        {
            out.count_both_strands.set("ACGT"[bi], f_counts[bi] + r_counts[bi]);
            out.extensions_both_strands.append(1, "ACGT"[bi]);
        }
    }

    return out;
}

// calculate the delta statistic from the neighboring kmer data
int calculate_delta(const std::string& kmer, 
                    const KmerNeighbors nd, 
                    const BWTIndexSet& index_set)
{
    // Calculate the number of reads that have both the 
    // kmer and the neighbor kmer
    std::vector<int> double_kmer_counts;
    for(size_t bi = 0; bi < nd.extensions_both_strands.size(); ++bi)
    {
        std::string n_kmer = kmer + nd.extensions_both_strands[bi];
        double_kmer_counts.push_back(
                BWTAlgorithms::countSequenceOccurrences(n_kmer, index_set.pBWT));
    }

    std::string sorted = KmerNeighbors::getExtensionsFromCount(nd.count_both_strands);
    assert(sorted.size() >= 2);

    char b_1 = sorted[0];
    char b_2 = sorted[1];

    size_t c_1 = nd.count_both_strands.get(b_1);
    size_t c_2 = nd.count_both_strands.get(b_2);

    int delta = c_1 + c_2 - double_kmer_counts[0] - double_kmer_counts[1];
    return delta;
}

KmerDistribution sample_kmer_counts(size_t k, size_t n, const BWTIndexSet& index_set)
{
    // Learn k-mer occurrence distribution for this value of k
    KmerDistribution distribution;
    for(size_t i = 0; i < n; ++i)
    {
        std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);
        if(s.size() < k)
            continue;
        size_t nk = s.size() - k + 1;    
        for(size_t j = 0; j < nk; ++j)
        {
            std::string kmer = s.substr(j, k);
            int count = BWTAlgorithms::countSequenceOccurrences(kmer, index_set.pBWT);
            distribution.add(count);
        }
    }

    return distribution;
}

// Find the single-copy peak of the kmer count distribution
// Extremely heterozygous diploid genomes might have multiple
// coverage peaks so we need to take care to not call
// the half-coverage peak
size_t find_single_copy_peak(const KmerDistribution& distribution)
{
    size_t MAX_COUNT = 10000;
   
    // Step 1: find first local minima
    size_t local_min = 1;
    size_t local_count = distribution.getNumberWithCount(local_min);
    while(local_min < MAX_COUNT)
    {
        size_t c = distribution.getNumberWithCount(local_min + 1);
        if(c > local_count)
            break;
        local_count = c;
        local_min += 1;
    }

    // Step 2: find the highest peak at count >= local_count
    size_t global_peak_v = 0;
    size_t global_peak_c = 0;
    for(size_t i = local_min; i < MAX_COUNT; ++i)
    {
        size_t c = distribution.getNumberWithCount(i);
        if(c > global_peak_c)
        {
            global_peak_v = i;
            global_peak_c = c;
        }
    }

    if(opt::verbose)
    {
        // print the kmer distribution to stderr
        fprintf(stderr, "findSingleCopyPeak -- KmerDistribution:\n");
        distribution.print(stderr, MAX_COUNT);

        fprintf(stderr, "Local min: %zu\n", local_min);
        fprintf(stderr, "Global peak: %zu\n", global_peak_v);
    }

    // Check whether the 2*n peak is reasonably close to the same height
    // as the global peak
    size_t n2_count = distribution.getNumberWithCount(2*global_peak_v);
    const double HEIGHT_FRAC = 0.5;
    if(HEIGHT_FRAC * global_peak_c < n2_count)
        return 2*global_peak_v;
    else
        return global_peak_v;
}

// Learn the parameters for the mixture count model
void learn_mixture_parameters(const KmerDistribution& distribution, ModelParameters& params)
{
    size_t lambda_hat = params.mode;
    //printf("\nLearning mixture\n");
    //distribution.print(100);

    assert(ModelParameters::NUM_MIXTURE_STATES > 3);

    // this array stores the lambda parameter for each poisson
    params.mixture_means[0] = (double)lambda_hat * params.error_rate_prior; // errors
    params.mixture_means[1] = lambda_hat / 2.f; // hets
    params.mixture_means[2] = lambda_hat; // homs
    for(size_t i = 3; i < ModelParameters::NUM_MIXTURE_STATES; ++i)
        params.mixture_means[i] = (i - 1) * lambda_hat; // repeats

    // initialize mixture proportions
    for(size_t i = 0; i < ModelParameters::NUM_MIXTURE_STATES; ++i)
    {
        //printf("c(%d): %d\n", (int)mixture_means[i], (int)distribution.getNumberWithCount(mixture_means[i]));
        params.mixture_proportions[i] = 1.0f / ModelParameters::NUM_MIXTURE_STATES;
    }

    // this vector stores the number of k-mers that have been softly
    // classified as each state for the current iteration
    double mixture_counts[ModelParameters::NUM_MIXTURE_STATES];

    // Calculate the likelihood using the current parameters

    // The model only fits counts up to a this copy number
    int max = params.mode * ModelParameters::MAX_REPEAT_COPY;
    std::vector<int> count_vector = distribution.toCountVector(max);
    
    double error_sum = 0;
    double error_n = 0;
    double prev_ll = -std::numeric_limits<double>::max();

    int max_iterations = 30;
    int iteration = 0;
    while(iteration++ < max_iterations) 
    {
        for(size_t i = 0; i < ModelParameters::NUM_MIXTURE_STATES; ++i)
            mixture_counts[i] = 0;

        for(size_t c = 1; c < count_vector.size(); ++c)
        {
            int n_c = count_vector[c];

            // Calculate P(c) and P(c|m_i)P(m_i) by summing over all mixture states
            double mixture_log_p[ModelParameters::NUM_MIXTURE_STATES];
            double log_p_c = 0.0f;
            for(size_t i = 0; i < ModelParameters::NUM_MIXTURE_STATES; ++i)
            {
                // Calculate the zero-truncation term
                double log_zero_trunc = log(1 / (1 - exp(-params.mixture_means[i])));

                // P(c| m_i) * P(m_i)
                mixture_log_p[i] = 
                    SGAStats::logPoisson(c, params.mixture_means[i]) + log_zero_trunc + log(params.mixture_proportions[i]);

                if(i == 0)
                    log_p_c = mixture_log_p[i];
                else
                    log_p_c = addLogs(log_p_c, mixture_log_p[i]);
            }

            //printf("c: %2zu n_c: %6d\n", c, n_c);
            // Calculate P(m_i|c)
            for(size_t i = 0; i < ModelParameters::NUM_MIXTURE_STATES; ++i)
            {
                mixture_log_p[i] -= log_p_c;
                mixture_counts[i] += exp(mixture_log_p[i]) * n_c;
                //double zero_trunc = 1 / (1 - exp(-mixture_means[i]));
                //double tp = exp(SGAStats::logPoisson(c, mixture_means[i])) * zero_trunc;
                //printf("\tp(m_%zu|c): %.3lf lambda: %.1lf p(c): %.3lf count: %.3lf\n", i, mixture_probability[i], mixture_means[i], p_c, mixture_counts[i]);
            }

            error_sum += c * n_c * exp(mixture_log_p[0]);
            error_n += n_c * exp(mixture_log_p[0]);
        }

        // Recalculate mixture proportions
        double sum = 0;
        for(size_t i = 0; i < ModelParameters::NUM_MIXTURE_STATES; ++i)
            sum += mixture_counts[i];
        
        for(size_t i = 0; i < ModelParameters::NUM_MIXTURE_STATES; ++i)
            params.mixture_proportions[i] = mixture_counts[i] / sum;

        // Recalculate lambda for the error state
        params.mixture_means[0] = error_sum / error_n;

        // Calculate the log-likelihood of the new parameters to see if we should stop
        double sum_ll = 0.0f;
        for(size_t c = 1; c < count_vector.size(); ++c)
        {
            int n_c = count_vector[c];
            
            double partial_ll = 0.0f;
            for(size_t i = 0; i < ModelParameters::NUM_MIXTURE_STATES; ++i)
            {
                // Calculate the zero-truncation term
                double log_zero_trunc = log(1 / (1 - exp(-params.mixture_means[i])));
                double ll = SGAStats::logPoisson(c, params.mixture_means[i]) + log_zero_trunc + log(params.mixture_proportions[i]);

                if(i == 0)
                    partial_ll = ll;
                else
                    partial_ll = addLogs(partial_ll, ll);
            }
            
            partial_ll *= n_c;
            sum_ll += partial_ll;
        }        
    
        double diff = sum_ll - prev_ll;
        double ll_improvement = diff / abs(sum_ll);

        if(iteration > 1 && opt::verbose > 1)
            fprintf(stderr, "PREV: %lf NOW: %lf DIFF: %lf I: %lf\n", prev_ll, sum_ll, diff, ll_improvement);

        if(ll_improvement < 0.00001)
            break;
        prev_ll = sum_ll;
    }

    if(iteration == max_iterations)
    {
        if(!opt::forceEM)
        {
            std::cerr << "\nError: coverage model failed to converge\n";
            std::cerr << "you can disable this error with the --force-em flag but the\n";
            std::cerr << "branch classification and genome size results may be erroneous\n";
            exit(EXIT_FAILURE);
        }
        else
        {
            std::cerr << "\nwarning: coverage model failed to converge\n";
        }
    }

    // Calculate the repeat-normalized proportions
    double sum_repeat = 0.0;
    for(size_t r = 3; r < ModelParameters::NUM_MIXTURE_STATES; ++r)
        sum_repeat += params.mixture_proportions[r];

    for(size_t r = 3; r < ModelParameters::NUM_MIXTURE_STATES; ++r)
        params.repeat_mixture_normalized[r] = params.mixture_proportions[r] / sum_repeat;


    /*
    printf("\nprop\t(");
    for(size_t i = 0; i < ModelParameters::NUM_MIXTURE_STATES; ++i)
        printf("%10.4lf, ", params.mixture_proportions[i]);
    printf(")\n");

    printf("counts\t(");
    for(size_t i = 0; i < ModelParameters::NUM_MIXTURE_STATES; ++i)
        printf("%10.4lf, ", mixture_counts[i]);
    printf("]\n");

    printf("means\t(");
    for(size_t i = 0; i < ModelParameters::NUM_MIXTURE_STATES; ++i)
        printf("%10.3lf, ", params.mixture_means[i]);
    printf(")\n");

    printf("repeats\t(");
    for(size_t i = 0; i < ModelParameters::NUM_MIXTURE_STATES; ++i)
        printf("%10.3lf, ", params.repeat_mixture_normalized[i]);
    printf(")\n");
    */
}

// Calculate parameters for the k-mer count model for the given distribution
ModelParameters calculate_model_parameters(size_t k, 
                                           KmerDistribution& distribution)
{
    ModelParameters params;
    
    params.k = k;

    // learn k-mer mode
    params.mode = find_single_copy_peak(distribution);

    // learn the mixture proportions
    learn_mixture_parameters(distribution, params);

    return params;
}


// Calculate parameters for the k-mer count model
ModelParameters calculate_model_parameters(size_t k, 
                                           size_t samples, 
                                           const BWTIndexSet& index_set)
{
    // sample kmers
    KmerDistribution distribution = sample_kmer_counts(k, samples, index_set);
    return calculate_model_parameters(k, distribution);
}

//
void generate_unipath_length_data(JSONWriter* pWriter, const BWTIndexSet& index_set)
{
    pWriter->String("UnipathLength");
    pWriter->StartArray();
    for(size_t k = 16; k < 96; k += 5)
        unipath_length_distribution(pWriter, index_set, k, 0.9, 1000);
    pWriter->EndArray();
}

//
void generate_kmer_coverage(JSONWriter* pWriter, const BWTIndexSet& index_set)
{
    pWriter->String("KmerDistribution");
    pWriter->StartObject();
    size_t n_samples = 10000;
    size_t k = 51;
    
    pWriter->String("k");
    pWriter->Int(k);
    
    KmerDistribution kmerDistribution;
    for(size_t i = 0; i < n_samples; ++i)
    {
        std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);
        int n = s.size();
        int nk = n - k + 1;
        for(int j = 0; j < nk; ++j)
        {
            std::string kmer = s.substr(j, k);
            int count = BWTAlgorithms::countSequenceOccurrences(kmer, index_set.pBWT);
            kmerDistribution.add(count);
        }
    }

    pWriter->String("distribution");
    pWriter->StartArray();
    int max = kmerDistribution.getCutoffForProportion(0.95f);
    std::vector<int> count_vector = kmerDistribution.toCountVector(max);
    for(size_t i = 1; i < count_vector.size(); ++i)
    {
        pWriter->StartObject();
        pWriter->String("kmer-depth");
        pWriter->Int(i);
        pWriter->String("count");
        pWriter->Int(count_vector[i]);
        pWriter->EndObject();
    }

    pWriter->EndArray();
    pWriter->EndObject();
}

//
void generate_position_of_first_error(JSONWriter* pWriter, const BWTIndexSet& index_set)
{
    size_t n_samples = 100000;
    size_t k = 41;
    size_t starting_count = 10;
    size_t min_count = 3;

    std::vector<size_t> position_count;
    std::vector<size_t> error_count;
    for(size_t i = 0; i < n_samples; ++i)
    {
        std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);
        if(s.length() < k)
            continue;

        size_t nk = s.size() - k + 1;
        size_t first_kmer_count = 
            BWTAlgorithms::countSequenceOccurrences(s.substr(0, k), index_set.pBWT);

        // Skip reads with a weak starting kmer
        if(first_kmer_count < starting_count)
            continue;

        for(size_t j = 1; j < nk; ++j)
        {
            size_t kmer_count =
                BWTAlgorithms::countSequenceOccurrences(s.substr(j, k), index_set.pBWT);

            if(j >= position_count.size())
            {
                position_count.resize(j+1);
                error_count.resize(j+1);
            }

            position_count[j] += 1;
            if(kmer_count < min_count)
            {
                error_count[j] += 1;
                break;
            }
        }
    }

    pWriter->String("FirstErrorPosition");
    pWriter->StartObject();

    pWriter->String("indices");
    pWriter->StartArray();
    for(size_t i = 1; i < position_count.size(); ++i)
        pWriter->Int(i);
    pWriter->EndArray();

    pWriter->String("base_count");
    pWriter->StartArray();
    for(size_t i = 1; i < position_count.size(); ++i)
        pWriter->Int(position_count[i]);
    pWriter->EndArray();

    pWriter->String("error_count");
    pWriter->StartArray();
    for(size_t i = 1; i < position_count.size(); ++i)
        pWriter->Int(error_count[i]);
    pWriter->EndArray();
    pWriter->EndObject();
}

//
void generate_errors_per_base(JSONWriter* pWriter, const BWTIndexSet& index_set)
{
    int n_samples = 100000;
    size_t k = 31;

    double max_error_rate = 0.95;
    size_t min_overlap = 50;
    
    std::vector<size_t> position_count;
    std::vector<size_t> error_count;

    Timer timer("test", true);
#if HAVE_OPENMP
        omp_set_num_threads(opt::numThreads);
        #pragma omp parallel for
#endif
    for(int i = 0; i < n_samples; ++i)
    {
        std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);
        if(s.length() < k)
            continue;

        KmerOverlaps::retrieveMatches(s, k, min_overlap, max_error_rate, 2, index_set);
        //KmerOverlaps::approximateMatch(s, min_overlap, max_error_rate, 2, 200, index_set);

        MultipleAlignment ma = 
            KmerOverlaps::buildMultipleAlignment(s, k, min_overlap, max_error_rate, 2, index_set);

        // Skip when there is insufficient depth to classify errors
        size_t ma_rows = ma.getNumRows();
        if(ma_rows <= 1)
            continue;

        size_t ma_cols = ma.getNumColumns();
        size_t position = 0;
        for(size_t j = 0; j < ma_cols; ++j)
        {
            char s_symbol = ma.getSymbol(0, j);

            // Skip gaps
            if(s_symbol == '-' || s_symbol == '\0')
                continue;
            
            SymbolCountVector scv = ma.getSymbolCountVector(j);
            int s_symbol_count = 0;
            char max_symbol = 0;
            int max_count = 0;

            for(size_t k = 0; k < scv.size(); ++k)
            {
                if(scv[k].symbol == s_symbol)
                    s_symbol_count = scv[k].count;
                if(scv[k].count > max_count)
                {
                    max_count = scv[k].count;
                    max_symbol = scv[k].symbol;
                }
            }

            // Call an error at this position if the consensus symbol differs from the read
            //    and the support for the read symbol is less than 4 and the consensus symbol
            //    is strongly supported.
            bool is_error = s_symbol != max_symbol && s_symbol_count < 4 && max_count >= 3;

#if HAVE_OPENMP
            #pragma omp critical
#endif
            {
                if(position >= position_count.size())
                {
                    position_count.resize(position+1);
                    error_count.resize(position+1);
                }

                position_count[position]++;
                error_count[position] += is_error;
            }
            position += 1;
        }
    }
    
    pWriter->String("ErrorsPerBase");
    pWriter->StartObject();
    
    pWriter->String("base_count");
    pWriter->StartArray();
    for(size_t i = 0; i < position_count.size(); ++i)
        pWriter->Int(position_count[i]);
    pWriter->EndArray();
    
    pWriter->String("error_count");
    pWriter->StartArray();
    for(size_t i = 0; i < position_count.size(); ++i)
        pWriter->Int(error_count[i]);
    pWriter->EndArray();

    pWriter->EndObject();
}

// Generate local graph complexity measure
void generate_local_graph_complexity(JSONWriter* pWriter, const BWTIndexSet& index_set)
{
    int n_samples = 50000;
    size_t min_coverage_to_test = 5;
    size_t min_coverage_for_branch = 3;
    double min_coverage_ratio = 0.5f;

    pWriter->String("LocalGraphComplexity");
    pWriter->StartArray();
    for(size_t k = 16; k < 86; k += 5)
    {
        size_t num_branches = 0;
        size_t num_kmers = 0;

#if HAVE_OPENMP
        omp_set_num_threads(opt::numThreads);
        #pragma omp parallel for
#endif
        for(int i = 0; i < n_samples; ++i)
        {
            std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);
            if(s.size() < k)
                continue;
            
            for(size_t j = 0; j < s.size() - k + 1; ++j)
            {
                std::string kmer = s.substr(j, k);
                size_t count = BWTAlgorithms::countSequenceOccurrences(kmer, index_set);
                if(count < min_coverage_to_test)
                    break;

                std::string extensions = 
                    get_valid_dbg_neighbors_coverage_and_ratio(kmer, 
                                                               index_set, 
                                                               min_coverage_for_branch, 
                                                               min_coverage_ratio,
                                                               ED_SENSE);
#if HAVE_OPENMP
                #pragma omp critical
#endif
                {    
                    num_branches += extensions.size() > 1;
                    num_kmers += 1;
                }

            }
        }

        pWriter->StartObject();
        pWriter->String("k");
        pWriter->Int(k);
        pWriter->String("num_kmers");
        pWriter->Int(num_kmers);
        pWriter->String("num_branches");
        pWriter->Int(num_branches);
        pWriter->EndObject();
    }
    pWriter->EndArray();
}

// Measure genome repetitiveness using the rate of k-mers
// that branch on both ends
void generate_double_branch(JSONWriter* pWriter, const BWTIndexSet& index_set)
{
    int n_samples = 50000;
    size_t min_coverage_to_test = 5;
    size_t min_coverage_for_branch = 3;
    double min_coverage_ratio = 0.5f;

    pWriter->String("DoubleBranch");
    pWriter->StartArray();
    for(size_t k = 16; k < 86; k += 5)
    {
        size_t num_branches = 0;
        size_t num_kmers = 0;

#if HAVE_OPENMP
        omp_set_num_threads(opt::numThreads);
        #pragma omp parallel for
#endif
        for(int i = 0; i < n_samples; ++i)
        {
            std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);
            if(s.size() < k)
                continue;

            std::string kmer = s.substr(0, k);
            size_t count = BWTAlgorithms::countSequenceOccurrences(kmer, index_set);
            if(count >= min_coverage_to_test)
            {
                std::string right_extensions = 
                    get_valid_dbg_neighbors_coverage_and_ratio(kmer, 
                                                               index_set,
                                                               min_coverage_for_branch, 
                                                               min_coverage_ratio,
                                                               ED_SENSE);

                std::string left_extensions = 
                    get_valid_dbg_neighbors_coverage_and_ratio(kmer, 
                                                               index_set,
                                                               min_coverage_for_branch, 
                                                               min_coverage_ratio,
                                                               ED_ANTISENSE);
#if HAVE_OPENMP
                    #pragma omp critical
#endif
                    {    
                        num_branches += (left_extensions.size() > 1 && right_extensions.size() > 1);
                        num_kmers += 1;
                    }
            }
        }

        pWriter->StartObject();
        pWriter->String("k");
        pWriter->Int(k);
        pWriter->String("num_kmers");
        pWriter->Int(num_kmers);
        pWriter->String("num_branches");
        pWriter->Int(num_branches);
        pWriter->EndObject();
    }
    pWriter->EndArray();
}

// Generate random walk length
void generate_random_walk_length(JSONWriter* pWriter, const BWTIndexSet& index_set)
{
    int n_samples = 1000;
    size_t max_length = 30000;
    size_t bf_overcommit = 10;

    pWriter->String("RandomWalkLength");
    pWriter->StartArray();

    for(size_t k = 16; k < 86; k += 5)
    {
        // Use a bloom filter to skip previously seen kmers
        BloomFilter* bloom_filter = new BloomFilter;;
        bloom_filter->initialize(n_samples * max_length * bf_overcommit, 3);

        pWriter->StartObject();
        pWriter->String("k");
        pWriter->Int(k);
        pWriter->String("walk_lengths");
        pWriter->StartArray();

#if HAVE_OPENMP
        omp_set_num_threads(opt::numThreads);
        #pragma omp parallel for
#endif
        for(int i = 0; i < n_samples; ++i)
        {
            size_t walk_length = 0;
            std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);
            if(s.size() < k)
                continue;

            std::string kmer = s.substr(0, k);
            std::string rc_kmer = reverseComplement(kmer);

            // Only start a walk from this kmer if is not in the bloom filter and has coverage on both strands
            bool in_filter = bloom_filter->test( (kmer < rc_kmer ? kmer.c_str() : rc_kmer.c_str()), k);
            size_t fc = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, index_set);
            size_t rc = BWTAlgorithms::countSequenceOccurrencesSingleStrand(rc_kmer, index_set);
            if(in_filter || fc == 0 || rc == 0)
                continue;

            while(walk_length < max_length) 
            {
                bloom_filter->add( (kmer < rc_kmer ? kmer.c_str() : rc_kmer.c_str()), k);

                // Get the possible extensions of this kmer
                int f_counts[5] = { 0, 0, 0, 0, 0 };
                int r_counts[5] = { 0, 0, 0, 0, 0 };

                fill_neighbor_count_by_strand(kmer, index_set, f_counts, r_counts);

                // Only allow extensions to vertices that have coverage on both strands
                std::string extensions;
                for(size_t bi = 0; bi < 4; ++bi)
                {
                    if(f_counts[bi] >= 1 && r_counts[bi] >= 1)
                        extensions.append(1, "ACGT"[bi]);
                }

                if(!extensions.empty())
                {
                    kmer.erase(0, 1);
                    kmer.append(1, extensions[rand() % extensions.size()]);
                    rc_kmer = reverseComplement(kmer);
                    walk_length += 1;
                }
                else
                {
                    break;
                }
            }
#if HAVE_OPENMP
        #pragma omp critical
#endif
            pWriter->Int(walk_length);
        }

        pWriter->EndArray();
        pWriter->EndObject();
        delete bloom_filter;
    }

    pWriter->EndArray();
}

void generate_gc_distribution(JSONWriter* pJSONWriter, const BWTIndexSet& index_set)
{
    int n_samples = 100000;
    size_t k = 31;
    double gc_bin_size = 0.05;

    std::vector<double> read_gc_sum;
    size_t read_gc_n = 0;
    
    std::vector<double> ref_gc_sum;
    size_t ref_gc_n = 0;

    std::vector<size_t> coverage_vector;
    std::vector<double> gc_vector;

    read_gc_sum.resize((size_t)(1.0f / gc_bin_size) + 1);
    ref_gc_sum.resize((size_t)(1.0f / gc_bin_size) + 1);

    // Calculate the gc content of sampled reads
    for(int i = 0; i < n_samples; ++i)
    {
        std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);
        if(s.size() < k)
            continue;

        size_t cov = BWTAlgorithms::countSequenceOccurrences(s.substr(0, k), index_set.pBWT);

        // Ignore singletons
        if(cov == 1)
            continue;

        double gc = 0.f;
        double at = 0.f;
        for(size_t j = 0; j < s.size(); ++j)
        {
            if(s[j] == 'C' || s[j] == 'G')
                gc += 1;
            else
                at += 1;
        }
        
        double gc_f = gc / (gc + at);
        size_t bin_idx = (size_t)(gc_f / gc_bin_size);
        read_gc_sum[bin_idx] += gc_f;
        read_gc_n += 1;

        gc_vector.push_back(gc_f);
        coverage_vector.push_back(cov);
    }

    // If a reference file is provided, calculate the reference GC for comparison
    if(!opt::referenceFile.empty())
    {
        size_t window_size = 100;
        SeqReader reader(opt::referenceFile, SRF_NO_VALIDATION);
        SeqRecord record;
        while(reader.get(record))
        {
            DNAString& seq = record.seq;
            size_t l = seq.length();
            for(size_t i = 0; i < l; i += window_size)
            {
                double gc = 0.f;
                double at = 0.f;
                for(size_t j = i; j < i + window_size && j < l; ++j)
                {
                    char c = seq.get(j);
                    if(c == 'C' || c == 'G')
                        gc += 1;
                    else if(c == 'A' || c == 'T') // IUPAC symbols may be found
                        at += 1;
                }
                double gc_f = gc / (gc + at);
                size_t bin_idx = (size_t)(gc_f / gc_bin_size);
                ref_gc_sum[bin_idx] += gc_f;
                ref_gc_n += 1;
            }
        }
    }

    pJSONWriter->String("GCDistribution");
    pJSONWriter->StartObject();

    pJSONWriter->String("gc_bins");
    pJSONWriter->StartArray();
    for(size_t i = 0; i < read_gc_sum.size(); ++i)
        pJSONWriter->Double(i * gc_bin_size + gc_bin_size / 2);
    pJSONWriter->EndArray();

    pJSONWriter->String("read_gc_prop");
    pJSONWriter->StartArray();
    for(size_t i = 0; i < read_gc_sum.size(); ++i)
        pJSONWriter->Double(read_gc_sum[i] / read_gc_n);
    pJSONWriter->EndArray();

    if(!opt::referenceFile.empty())
    {
        pJSONWriter->String("ref_gc_prop");
        pJSONWriter->StartArray();
        for(size_t i = 0; i < read_gc_sum.size(); ++i)
            pJSONWriter->Double(ref_gc_sum[i] / ref_gc_n);
        pJSONWriter->EndArray();
    }

    pJSONWriter->String("gc_samples");
    pJSONWriter->StartArray();
    for(size_t i = 0; i < gc_vector.size(); ++i)
        pJSONWriter->Double(gc_vector[i]);
    pJSONWriter->EndArray();

    pJSONWriter->String("cov_samples");
    pJSONWriter->StartArray();
    for(size_t i = 0; i < gc_vector.size(); ++i)
        pJSONWriter->Double(coverage_vector[i]);
    pJSONWriter->EndArray();

    pJSONWriter->EndObject();
}

void generate_duplication_rate(JSONWriter* pJSONWriter, const BWTIndexSet& index_set)
{
    int n_samples = 10000;
    size_t k = 50;

    size_t total_pairs = index_set.pBWT->getNumStrings() / 2;
    size_t num_pairs_checked = 0;
    size_t num_duplicates = 0;
#if HAVE_OPENMP
    omp_set_num_threads(opt::numThreads);
    #pragma omp parallel for
#endif
    for(int i = 0; i < n_samples; ++i)
    {
        // Choose a read pair
        int64_t source_pair_idx = rand() % total_pairs;
        std::string r1 = BWTAlgorithms::extractString(index_set.pBWT, source_pair_idx * 2);
        std::string r2 = BWTAlgorithms::extractString(index_set.pBWT, source_pair_idx * 2 + 1);

        // Get the interval for $k1/$k2 which corresponds to the 
        // lexicographic rank of reads starting with those kmers
        std::string k1 = "$" + r1.substr(0, k);
        std::string k2 = "$" + r2.substr(0, k);

        BWTInterval i1 = BWTAlgorithms::findInterval(index_set.pBWT, k1);
        BWTInterval i2 = BWTAlgorithms::findInterval(index_set.pBWT, k2);

        std::vector<int64_t> pair_ids;
        for(int64_t j = i1.lower; j <= i1.upper; ++j)
        {
            int64_t read_id = index_set.pSSA->lookupLexoRank(j);
            if(read_id % 2 == 1)
                continue;
            
            int64_t pair_id = read_id % 2 == 0 ? read_id / 2 : (read_id - 1) / 2;
            if(pair_id != source_pair_idx)
                pair_ids.push_back(pair_id);
        }

        for(int64_t j = i2.lower; j <= i2.upper; ++j)
        {
            int64_t read_id = index_set.pSSA->lookupLexoRank(j);
            if(read_id % 2 == 0)
                continue;
            int64_t pair_id = read_id % 2 == 0 ? read_id / 2 : (read_id - 1) / 2;
            if(pair_id != source_pair_idx)
                pair_ids.push_back(pair_id);
        }

        std::sort(pair_ids.begin(), pair_ids.end());
        std::vector<int64_t>::iterator iter = 
            std::adjacent_find(pair_ids.begin(), pair_ids.end());
                                           
        bool has_duplicate = iter != pair_ids.end();
#if HAVE_OPENMP
        #pragma omp critical
#endif
        {
            num_pairs_checked += 1;
            num_duplicates += has_duplicate;
        }
    }

    pJSONWriter->String("PCRDuplicates");
    pJSONWriter->StartObject();
    pJSONWriter->String("num_duplicates");
    pJSONWriter->Int(num_duplicates);
    pJSONWriter->String("num_pairs");
    pJSONWriter->Int(num_pairs_checked);
    pJSONWriter->EndObject();
}

// Write a stream of calculated fragments sizes to the JSON file
void generate_pe_fragment_sizes(JSONWriter* pJSONWriter, const BWTIndexSet& index_set)
{
    int n_samples = 100000;
    size_t k = 51;
    size_t MAX_INSERT = 1500;

    std::vector<size_t> fragment_sizes;

    size_t total_pairs = index_set.pBWT->getNumStrings() / 2;
#if HAVE_OPENMP
    omp_set_num_threads(opt::numThreads);
    #pragma omp parallel for
#endif
    for(int i = 0; i < n_samples; ++i)
    {
        // Choose a read pair
        int64_t source_pair_idx = rand() % total_pairs;
        std::string start_kmer = BWTAlgorithms::extractString(index_set.pBWT, source_pair_idx * 2).substr(0, k);
        std::string end_kmer = BWTAlgorithms::extractString(index_set.pBWT, source_pair_idx * 2 + 1).substr(0, k);
        // We assume that the pairs are orientated F/R therefore we reverse-complement k_end
        end_kmer = reverseComplement(end_kmer);

        // Aggressively walk the de Bruijn graph starting from k_start until k_end is found or we give up
        size_t steps = 0;
        bool found = false;
        while(!found && steps < MAX_INSERT)
        {
            // A coverage ratio of 1.0 will force use to only use the highest-coverage branch
            // This may generate erroneous insert sizes in (rare?) cases but will give a good approximation
            // to the real distribution
            std::string extensions = get_valid_dbg_neighbors_ratio(start_kmer, index_set, 1.0f);
            if(extensions.empty())
                break;

            start_kmer.erase(0, 1);
            start_kmer.append(1, extensions[0]);
            found = start_kmer == end_kmer;
            steps += 1;
        }

        if(found)
        {
#if HAVE_OPENMP
            #pragma omp critical
#endif
            fragment_sizes.push_back(steps + k);
        }
    }

    pJSONWriter->String("FragmentSize");
    pJSONWriter->StartObject();
    pJSONWriter->String("sizes");
    pJSONWriter->StartArray();
    for(size_t i = 0; i < fragment_sizes.size(); ++i)
        pJSONWriter->Int(fragment_sizes[i]);
    pJSONWriter->EndArray();
    pJSONWriter->EndObject();
}

// Generate a report of the quality of each base
void generate_quality_stats(JSONWriter* pJSONWriter, const std::string& filename)
{
    size_t max_reads = 10000000;
    double sample_rate = 0.05;
    SeqReader reader(filename, SRF_KEEP_CASE | SRF_NO_VALIDATION);
    SeqRecord record;

    size_t n_reads = 0;
    std::vector<size_t> bases_checked;
    std::vector<size_t> sum_quality;
    std::vector<size_t> num_q30;

    while(reader.get(record) && n_reads++ < max_reads)
    {
        if((double)rand() / RAND_MAX < sample_rate && record.qual.length() == record.seq.length())
        {
            size_t l = record.seq.length();
            if(l > bases_checked.size())
            {
                bases_checked.resize(l);
                sum_quality.resize(l);
                num_q30.resize(l);
            }

            for(size_t i = 0; i < l; ++i)
            {
                bases_checked[i]++;
                size_t q = record.getPhredScore(i);
                sum_quality[i] += q;
                num_q30[i] += (q >= 30);
            }
        }
    }

    pJSONWriter->String("QualityScores");
    pJSONWriter->StartObject();
    
    pJSONWriter->String("mean_quality");
    pJSONWriter->StartArray();
    for(size_t i = 0; i < bases_checked.size(); ++i)
        pJSONWriter->Double((float)sum_quality[i] / bases_checked[i]);
    pJSONWriter->EndArray();

    pJSONWriter->String("fraction_q30");
    pJSONWriter->StartArray();
    for(size_t i = 0; i < bases_checked.size(); ++i)
        pJSONWriter->Double((float)num_q30[i] / bases_checked[i]);
    pJSONWriter->EndArray();
    pJSONWriter->EndObject();
}

//
GenomeEstimates estimate_genome_size_from_k_counts(size_t k, const BWTIndexSet& index_set)
{
    //
    size_t sum_read_length = 0;
    
    KmerDistribution kmerDistribution;
    for(size_t i = 0; i < opt::kmerDistributionSamples; ++i)
    {
        std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);
        int n = s.size();
        int nk = n - k + 1;
        for(int j = 0; j < nk; ++j)
        {
            std::string kmer = s.substr(j, k);
            int count = BWTAlgorithms::countSequenceOccurrences(kmer, index_set.pBWT);
            kmerDistribution.add(count);
        }
        sum_read_length += s.size();
    }

    // calculate the k-mer count model parameters from the distribution
    // this gives us the estimated proportion of kmers that contain errors
    ModelParameters params = calculate_model_parameters(k, kmerDistribution);

    double prop_kmers_with_error = params.mixture_proportions[0];

    size_t avg_rl = sum_read_length / opt::kmerDistributionSamples;
    size_t n = index_set.pBWT->getNumStrings();
    double total_read_kmers = n * (avg_rl - k + 1);
    double corrected_mode = params.mode / (1.0f - prop_kmers_with_error);
    
    //printf("k: %zu mode: %.2lf c_mode: %.2lf g_m: %.2lf g_cm: %.2lf \n", k, mode, corrected_mode, total_read_kmers / mode, total_read_kmers / corrected_mode);

    GenomeEstimates out;
    out.genome_size = static_cast<size_t>(total_read_kmers / corrected_mode);
    out.average_read_length = avg_rl;
    out.total_reads = n;
    out.mean_read_starts = (double)out.total_reads / out.genome_size;
    return out;
}

GenomeEstimates generate_genome_size(JSONWriter* pJSONWriter, const BWTIndexSet& index_set)
{

    /* Debug 
    for(size_t k = 16; k < 81; k += 5)
        estimate_genome_size_from_k_counts(k, index_set);
    */

    size_t estimate_k = 31;
    GenomeEstimates e = estimate_genome_size_from_k_counts(estimate_k, index_set);
    pJSONWriter->String("GenomeSize");
    pJSONWriter->StartObject();
    pJSONWriter->String("k");
    pJSONWriter->Int64(estimate_k);
    pJSONWriter->String("size");
    pJSONWriter->Int64((size_t)e.genome_size);
    pJSONWriter->EndObject();
    return e;
}

// Run the probablisitic classifier on a two-edge branch in the de Bruijn graph
ModelPosteriors classify_2_branch(const ModelParameters& params,
                                  const GenomeEstimates& genome,
                                  size_t higher_count, 
                                  size_t lower_count, 
                                  int delta)
{
    size_t total = higher_count + lower_count;

    // the expected increase in read count for the variant and error models
    double delta_param_unique = genome.mean_read_starts + params.mixture_means[0] * params.mixture_proportions[0];

    // Error model
    double log_p_delta_error = SGAStats::logPoisson(delta, delta_param_unique);
    double log_p_balance_error = SGAStats::logIntegerBetaBinomialDistribution(higher_count, total, 50, 1);

    // Variation Model
    double log_p_delta_variant = log_p_delta_error;
    double log_p_balance_variant = SGAStats::logBinomial(higher_count, total, 0.5);

    // Repeat model
    // This is explicitly diploid where 2 copies is the normal state for the genome
    // The repeat models (1...n) extra genomic copies
    double log_p_delta_repeat = 0.0f;
    int extra_copies_min = 1;
    int extra_copies_max = ModelParameters::MAX_REPEAT_COPY - 2;

    for(int r = extra_copies_min; r < extra_copies_max; ++r)
    {
        // Calculate the probability of having c copies using our mixture model
        // The first repeat state in the below array starts at entry 3, hence the offset
        double log_p_copies = params.repeat_mixture_normalized[r + 2];

        // Calculate the probability of the increase in read count given the number of extra copies
        double log_delta_given_copies = SGAStats::logPoisson(delta, r * params.mode);
        
        if(r == extra_copies_min)
            log_p_delta_repeat = log_p_copies + log_delta_given_copies;
        else
            log_p_delta_repeat = addLogs(log_p_delta_repeat, log_p_copies + log_delta_given_copies);
    }

    //double log_p_balance_repeat = SGAStats::logIntegerBetaDistribution(norm_allele_balance, 3, 1);
    double log_p_balance_repeat = SGAStats::logIntegerBetaBinomialDistribution(higher_count, total, 5, 1);

    // priors
    double log_prior_error = log(1.0/3);
    double log_prior_variant = log(1.0/3);
    double log_prior_repeat = log(1.0/3);

    // Calculate posterior for each state
    double ls1 = log_p_delta_error + log_p_balance_error + log_prior_error;
    double ls2 = log_p_delta_variant + log_p_balance_variant + log_prior_variant;
    double ls3 = log_p_delta_repeat + log_p_balance_repeat + log_prior_repeat;

    double log_sum = ls1;
    log_sum = addLogs(log_sum, ls2);
    log_sum = addLogs(log_sum, ls3);

    ModelPosteriors out;
    out.posterior_error = exp(ls1 - log_sum);
    out.posterior_variant = exp(ls2 - log_sum);
    out.posterior_repeat = exp(ls3 - log_sum);

    double max = 0;
    out.classification = BC_NO_CALL;
    if(out.posterior_error > max)
    {
        max = out.posterior_error;
        out.classification = BC_ERROR;
    } 

    if(out.posterior_variant > max)
    {
        max = out.posterior_variant;
        out.classification = BC_VARIANT;
    }

    if(out.posterior_repeat > max)
    {
        max = out.posterior_repeat;
        out.classification = BC_REPEAT;
    }

    /*
    double allele_balance = (double)higher_count / total;
    fprintf(stderr, "\tlog_sum: %.4lf\n", log_sum);
    fprintf(stderr, "\t\ttc|e: %.4lf y|e: %.4lf p(e): %.4lf\n", exp(log_p_delta_error), exp(log_p_balance_error), out.posterior_error);
    fprintf(stderr, "\t\ttc|v: %.4lf y|v: %.4lf p(v): %.4lf\n", exp(log_p_delta_variant), exp(log_p_balance_variant), out.posterior_variant);
    fprintf(stderr, "\t\ttc|r: %.4lf y|r: %.4lf p(r): %.4lf\n", exp(log_p_delta_repeat), exp(log_p_balance_repeat), out.posterior_repeat);
    fprintf(stderr, "DF %.0lf %zu %zu %zu %.4lf %.6lf %.6lf %.6lf %d %d\n",
        params.mode, higher_count, lower_count, total,
        allele_balance,
        out.posterior_error, out.posterior_variant, out.posterior_repeat,
        out.classification, delta);
    */
    return out;
}

// Calculate the probability that a kmer seen c times
// is a kmer that is contained once on each parental chromsome
// We assume the probability of a kmer appearing twice
// in one parent and zero times in the other is negligable
double probability_diploid_copy(const ModelParameters& params, size_t c)
{
    double mixture_log_p[ModelParameters::NUM_MIXTURE_STATES];
    double log_p_c = 0.0f;
    for(size_t i = 0; i < ModelParameters::NUM_MIXTURE_STATES; ++i)
    {
        // Calculate the zero-truncation term
        double log_zero_trunc = log(1 / (1 - exp(-params.mixture_means[i])));

        // P(c| m_i) * P(m_i)
        mixture_log_p[i] = 
            SGAStats::logPoisson(c, params.mixture_means[i]) + 
            log_zero_trunc + 
            log(params.mixture_proportions[i]);

        if(i == 0)
            log_p_c = mixture_log_p[i];
        else
            log_p_c = addLogs(log_p_c, mixture_log_p[i]);
    }

    return exp(mixture_log_p[2] - log_p_c);
}

//
void generate_branch_classification(JSONWriter* pWriter, 
                                    GenomeEstimates estimates, 
                                    const BWTIndexSet& index_set)
{
    int classification_samples = 1000000;

    pWriter->String("BranchClassification");
    pWriter->StartArray();

    for(size_t k = 21; k <= 71; k += 5)
    {
        // Estimate parameters to the model
        ModelParameters params = calculate_model_parameters(k, 
                                                            opt::kmerDistributionSamples, 
                                                            index_set);

        // Do not attempt classification if k-mer coverage is too low
        if(params.mode < 15)
            continue;

        // Cache the estimates that a kmer is diploid-unique for count [0,3m]
        size_t max_count = static_cast<size_t>(3 * params.mode);
        std::vector<double> p_unique_by_count(max_count, 0.0f);
        
        // We start at c=4 so that we never use low-count kmers in our model
        for(size_t c = 4; c < max_count; ++c)
            p_unique_by_count[c] = probability_diploid_copy(params, c);

        // Our output counts
        double num_error_branches = 0;
        double num_variant_branches = 0;
        double num_repeat_branches = 0;
        double num_kmers = 0;
        double mean_count = 0;
        size_t n_tests = 0;

#if HAVE_OPENMP
        omp_set_num_threads(opt::numThreads);
        #pragma omp parallel for
#endif
        for(int i = 0; i < classification_samples; ++i)
        {
            std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);
            if(s.size() < k)
                continue;
            
            size_t nk = s.size() - k + 1;    
            for(size_t j = 0; j < nk; ++j)
            {
                std::string kmer = s.substr(j, k);
                size_t count = BWTAlgorithms::countSequenceOccurrences(kmer, index_set.pBWT);
                
                // deep kmer, ignore
                if(count >= p_unique_by_count.size())
                    continue;

                double p_single_copy = p_unique_by_count[count];
                if(p_single_copy < 0.90)
                    continue;

                KmerNeighbors neighbors = calculate_neighbor_data(kmer, index_set);

                // if the neighbor data indicates a branch, classify it
                ModelPosteriors ret;
                if(neighbors.extensions_both_strands.size() > 1)
                {
                    std::string sorted = KmerNeighbors::getExtensionsFromCount(neighbors.count_both_strands);
                    char b_1 = sorted[0];
                    char b_2 = sorted[1];

                    size_t c_1 = neighbors.count_both_strands.get(b_1);
                    size_t c_2 = neighbors.count_both_strands.get(b_2);
                    
                    // Calculate delta, the increase in coverage for the neighboring kmers
                    int delta = calculate_delta(kmer, neighbors, index_set);
                    assert(delta >= 0);
                    ret = classify_2_branch(params, estimates, c_1, c_2, delta);
                }

#if HAVE_OPENMP
                #pragma omp critical
#endif
                {
                    num_error_branches += (p_single_copy * ret.posterior_error);
                    num_variant_branches += (p_single_copy * ret.posterior_variant);
                    num_repeat_branches += (p_single_copy * ret.posterior_repeat);
                    num_kmers += p_single_copy;
                    mean_count += count;
                    n_tests += 1;
                }

                if(ret.classification == BC_VARIANT || ret.classification == BC_REPEAT)
                    break;
            }
        }

        pWriter->StartObject();
        pWriter->String("k");
        pWriter->Int(k);
        pWriter->String("mode");
        pWriter->Int((int)params.mode);
        pWriter->String("num_kmers");
        pWriter->Double(num_kmers);
        pWriter->String("num_error_branches");
        pWriter->Double(num_error_branches);
        pWriter->String("num_variant_branches");
        pWriter->Double(num_variant_branches);
        pWriter->String("num_repeat_branches");
        pWriter->Double(num_repeat_branches);

        pWriter->String("variant_rate");
        double vr = num_variant_branches > 0 ? num_kmers / num_variant_branches : num_kmers;
        pWriter->Double(vr);
        
        pWriter->String("repeat_rate");
        double rr = num_repeat_branches > 0 ? num_kmers / num_repeat_branches : 10e-9;
        pWriter->Double(rr);

        pWriter->EndObject();
    }
    pWriter->EndArray();
}

// This is a version of the branch classification code that operates
// on the de bruijn graph of a diploid reference genome.
// The probablistic classifier is replaced by a simple
// counting classifier.
void generate_reference_branch_classification(JSONWriter* pWriter, const BWTIndexSet& index_set)
{
    int classification_samples = 50000;

    //double min_probability_for_classification = 0.25;
    pWriter->String("BranchClassification");
    pWriter->StartArray();

    for(size_t k = 21; k <= 71; k += 5)
    {
        size_t num_error_branches = 0;
        size_t num_variant_branches = 0;
        size_t num_repeat_branches = 0;
        size_t num_kmers = 0;

#if HAVE_OPENMP
        omp_set_num_threads(opt::numThreads);
        #pragma omp parallel for
#endif
        for(int i = 0; i < classification_samples; ++i)
        {
            std::string s = BWTAlgorithms::sampleRandomSubstring(index_set.pBWT, 100);
            if(s.size() < k)
                continue;
            
            size_t nk = s.size() - k + 1;
            for(size_t j = 0; j < nk; ++j)
            {
                std::string kmer = s.substr(j, k);
                int count = BWTAlgorithms::countSequenceOccurrences(kmer, index_set.pBWT);
                
                // explicitly assuming diploid
                if(count != 2)
                    continue;

                // these vectors are in order ACGT$ on the forward strand
                int f_counts[5] = { 0, 0, 0, 0, 0 };
                int r_counts[5] = { 0, 0, 0, 0, 0 };

                fill_neighbor_count_by_strand(kmer, index_set, f_counts, r_counts);
                
                // Make sure both strands are represented for every branch
                AlphaCount64 sum_counts;
                size_t num_extensions = 0;
                for(size_t bi = 0; bi < 4; ++bi)
                {
                    size_t t = f_counts[bi] + r_counts[bi];
                    sum_counts.set("ACGT"[bi], t);
                    if(t > 0)
                        num_extensions += 1;
                }

                // classify the branch
                BranchClassification classification = BC_NO_CALL;
                if(num_extensions == 2)
                {
                    char sorted_bases[5] = "ACGT";
                    sorted_bases[4] = '\0';
                    sum_counts.getSorted(sorted_bases, 5);

                    size_t c_1 = sum_counts.get(sorted_bases[0]);
                    size_t c_2 = sum_counts.get(sorted_bases[1]);

                    if(c_1 == 1 && c_2 == 1)
                        classification = BC_VARIANT;
                    else
                        classification = BC_REPEAT;
                }

#if HAVE_OPENMP
                #pragma omp critical
#endif
                {    
                    num_error_branches += (classification == BC_ERROR);
                    num_variant_branches += (classification == BC_VARIANT);
                    num_repeat_branches += (classification == BC_REPEAT);
                    num_kmers += 1;
                }

                // Do not continue with this read if we hit a repeat or variant
                if(classification == 1 || classification == 2)
                    break;
            }
        }

        pWriter->StartObject();
        pWriter->String("k");
        pWriter->Int(k);
        pWriter->String("mode");
        pWriter->Int(2);
        pWriter->String("num_kmers");
        pWriter->Int(num_kmers);
        pWriter->String("num_error_branches");
        pWriter->Int(num_error_branches);
        pWriter->String("num_variant_branches");
        pWriter->Int(num_variant_branches);
        pWriter->String("num_repeat_branches");
        pWriter->Int(num_repeat_branches);

        pWriter->String("variant_rate");
        double vr = num_variant_branches > 0 ? (double)num_kmers / num_variant_branches : num_kmers;
        pWriter->Double(vr);
        
        pWriter->String("repeat_rate");
        double rr = num_repeat_branches > 0 ? (double)num_kmers / num_repeat_branches : num_kmers;
        pWriter->Double(rr);

        pWriter->EndObject();
    }
    pWriter->EndArray();
}


//
void generate_de_bruijn_simulation(JSONWriter* pWriter,
                                   GenomeEstimates estimates,
                                   const BWTIndexSet& index_set)
{
    int n_samples = 20000;

    pWriter->String("SimulateAssembly");
    pWriter->StartArray();

    for(size_t k = 21; k <= 91; k += 5)
    {
        // Set up a bloom filter to avoid sampling paths multiple times
        size_t max_expected_kmers = opt::maxContigLength * n_samples;

        BloomFilter* bf = new BloomFilter;
        bf->initialize(5 * max_expected_kmers, 3);

        ModelParameters params = 
            calculate_model_parameters(k, opt::kmerDistributionSamples, index_set);

        // Cache the estimates that a kmer is diploid-unique for count [0,3m]
        size_t max_count = static_cast<size_t>(3 * params.mode);
        std::vector<double> p_unique_by_count(max_count, 0.0f);
        
        // We start at c=4 so that we never use low-count kmers in our model
        for(size_t c = 4; c < max_count; ++c)
            p_unique_by_count[c] = probability_diploid_copy(params, c);

        pWriter->StartObject();
        pWriter->String("k");
        pWriter->Int(k);
        pWriter->String("mode");
        pWriter->Int(static_cast<int>(params.mode));
        pWriter->String("walk_lengths");
        pWriter->StartArray();
#if HAVE_OPENMP
        omp_set_num_threads(opt::numThreads);
        #pragma omp parallel for
#endif
        for(int i = 0; i < n_samples; ++i)
        {
            //
            // Find a new starting point for the walk
            //

            // Get a random read from the BWT
            std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);

            // Use the first-kmer of the read to seed the seach
            if(s.size() < k)
                continue;
            
            std::string start_kmer = s.substr(0, k);

            size_t count = BWTAlgorithms::countSequenceOccurrences(start_kmer, index_set);
            
            // Only start walks from paths that are likely to be unique diploid sequence
            if(count >= 3 * params.mode)
                continue;

            double p_single_copy = p_unique_by_count[count];
            if(p_single_copy < 0.50)
                continue;

            // skip if this kmer has been used in a previous walk
            std::string rc_start_kmer = reverseComplement(start_kmer);
            bool in_filter = bf->test( (start_kmer < rc_start_kmer ? start_kmer.c_str() : rc_start_kmer.c_str()), k);
            if(in_filter)
                continue;

            //
            // All checks pass, start a new walk
            //
            std::set<std::string> kmer_set;

            for(size_t dir = 0; dir <= 1; ++dir)
            {
                std::string curr_kmer = dir == 0 ? start_kmer : reverseComplement(start_kmer);
                kmer_set.insert(curr_kmer);

                bool done = false;
                while(!done && kmer_set.size() < opt::maxContigLength)
                {
                    KmerNeighbors neighbors = calculate_neighbor_data(curr_kmer, index_set);

                    char extension_base = '\0';
                    if(neighbors.extensions_both_strands.size() < 2)
                    {
                        // No ambiguity, just pick the highest coverage extension as the next node
                        // In this case we do not require the picked node to have coverage on both strands
                        std::string sorted = KmerNeighbors::getExtensionsFromCount(neighbors.total_count);
                        char best_extension = sorted[0];
                        if(neighbors.total_count.get(best_extension) > 0)
                            extension_base = best_extension;
                    } 
                    else
                    {
                        std::string sorted = KmerNeighbors::getExtensionsFromCount(neighbors.count_both_strands);
                        char b_1 = sorted[0];
                        char b_2 = sorted[1];

                        size_t c_1 = neighbors.count_both_strands.get(b_1);
                        size_t c_2 = neighbors.count_both_strands.get(b_2);

                        // Calculate delta and classify the branch
                        int delta = calculate_delta(curr_kmer, neighbors, index_set);
                        assert(delta >= 0);
                        ModelPosteriors ret = classify_2_branch(params, estimates, c_1, c_2, delta);

                        if(ret.classification == BC_ERROR || ret.classification == BC_VARIANT)
                        {
                            // if this is an error branch, we take the non-error (higher coverage) option
                            // if this is a variant path we also take the higher coverage option to simulate
                            // a successfully popped bubble
                            extension_base = sorted[0];
                        }
                    }

                    if(extension_base != '\0')
                    {
                        curr_kmer.erase(0, 1);
                        curr_kmer.append(1, extension_base);
                        
                        // the insert call returns true in the second
                        // element of the pair if it succeeds
                        if(!kmer_set.insert(curr_kmer).second)
                            done = true;
                    }
                    else
                    {
                        done = true;
                    }
                }
            }

#if HAVE_OPENMP
            #pragma omp critical
#endif
            {
                // For small genomes there is a very real possibility that multiple threads
                // found the same path simultaneously. We do this update section within
                // a lock to detect this case using the bloom filter. If more than p percentage 
                // kmers in the walk are already in the filter, we reject the path
                size_t total_kmers = kmer_set.size();
                size_t kmers_in_filter = 0;
                for(std::set<std::string>::iterator iter = kmer_set.begin();
                        iter != kmer_set.end(); ++iter)
                {
                    std::string rc_curr = reverseComplement(*iter);
                    const char* bf_key = *iter < rc_curr ? iter->c_str() : rc_curr.c_str();
                    if(bf->test(bf_key, k))
                        kmers_in_filter += 1;
                    else
                        bf->add(bf_key, k);
                }
                
                // Put a threshold on the number of kmers that can
                // already be in the filter to reject the path.
                // This threshold should be way above the false positive
                // rate of the filter
                double f_rate = (double)kmers_in_filter / total_kmers;
                if(f_rate < 0.2)
                    pWriter->Int(total_kmers);
            }
        }
        pWriter->EndArray();
        pWriter->EndObject();
        delete bf;
    }
    pWriter->EndArray();
}

// Main
//
int preQCMain(int argc, char** argv)
{
    parsePreQCOptions(argc, argv);
    Timer* pTimer = new Timer(PROGRAM_IDENT);

    rapidjson::FileStream f(stdout);
    JSONWriter writer(f);

    // Top-level document
    writer.StartObject();
    
    // In simple mode we only compute metrics that do not need the FM-index
    generate_quality_stats(&writer, opt::readsFile);

    // Load the FM-index and compute the rest of the metrics if requested
    if(!opt::simple)
    {
        BWTIndexSet index_set;
        fprintf(stderr, "Loading FM-index of %s\n", opt::readsFile.c_str());
        index_set.pBWT = new BWT(opt::prefix + BWT_EXT);
        index_set.pSSA = new SampledSuffixArray(opt::prefix + SAI_EXT, SSA_FT_SAI);
        index_set.pCache = new BWTIntervalCache(10, index_set.pBWT);

        if(!opt::diploidReferenceMode)
        {
            GenomeEstimates estimates = generate_genome_size(&writer, index_set);

            generate_branch_classification(&writer, estimates, index_set);
            generate_de_bruijn_simulation(&writer, estimates, index_set);
            generate_gc_distribution(&writer, index_set);
            generate_pe_fragment_sizes(&writer, index_set);
            generate_kmer_coverage(&writer, index_set);
            generate_position_of_first_error(&writer, index_set);
            generate_errors_per_base(&writer, index_set);
            generate_unipath_length_data(&writer, index_set);
            generate_duplication_rate(&writer, index_set);
        }
        else
        {
            generate_reference_branch_classification(&writer, index_set);
        }

        delete index_set.pBWT;
        delete index_set.pSSA;
        delete index_set.pCache;
    }

    // End document
    writer.EndObject();

    delete pTimer;
    return 0;
}

//
void unipath_length_distribution(JSONWriter* pWriter,
                                 const BWTIndexSet& index_set, 
                                 size_t k,
                                 double coverage_ratio_threshold, 
                                 size_t n_samples)
{
    pWriter->StartObject();
    pWriter->String("k");
    pWriter->Int(k);
    pWriter->String("walk_lengths");
    pWriter->StartArray();

#if HAVE_OPENMP
    omp_set_num_threads(opt::numThreads);
    #pragma omp parallel for
#endif
    for(int i = 0; i < (int)n_samples; ++i)
    {
        // Get a random read from the BWT
        std::string s = BWTAlgorithms::sampleRandomString(index_set.pBWT);

        // Use the first-kmer of the read to seed the seach
        if(s.size() < k)
            continue;
        
        HashMap<std::string, bool> loop_check;
        std::string start_kmer = s.substr(0, k);
        std::string curr_kmer = start_kmer;

        size_t walk_length = 0;
        bool done = false;
        while(!done)
        {
            loop_check[curr_kmer] = true;
            std::string extensions = get_valid_dbg_neighbors_ratio(curr_kmer, index_set, coverage_ratio_threshold);
            if(extensions.size() == 1)
            {
                curr_kmer.erase(0, 1);
                curr_kmer.append(1, extensions[0]);

                if(loop_check.find(curr_kmer) == loop_check.end())
                    walk_length += 1;
                else
                    done = true;
            }
            else
            {
                done = true;
            }
        }
#if HAVE_OPENMP
        #pragma omp critical
#endif
        pWriter->Int(walk_length);
    }
    pWriter->EndArray();
    pWriter->EndObject();
}
// 
// Handle command line arguments
//
void parsePreQCOptions(int argc, char** argv)
{
    std::string algo_str;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'p': arg >> opt::prefix; break;
            case 't': arg >> opt::numThreads; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_SIMPLE: opt::simple = true; break;
            case OPT_MAX_CONTIG: arg >> opt::maxContigLength; break;
            case OPT_DIPLOID: opt::diploidReferenceMode = true; break;
            case OPT_REFERENCE: arg >> opt::referenceFile; break;
            case OPT_FORCE_EM: opt::forceEM = true; break;
            case OPT_HELP:
                std::cout << PREQC_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << PREQC_VERSION_MESSAGE;
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

    if (die) 
    {
        std::cout << "\n" << PREQC_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::readsFile = argv[optind++];

    if(opt::prefix.empty())
        opt::prefix = stripGzippedExtension(opt::readsFile);
}
