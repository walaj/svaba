//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// haplotype-filter - perform quality checks
// on haplotypes and their associated variants
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "haplotype-filter.h"
#include "BWTAlgorithms.h"
#include "BWTIndexSet.h"
#include "SGACommon.h"
#include "Timer.h"
#include "SequenceProcessFramework.h"
#include "HapgenUtil.h"
#include "SGAStats.h"
#include "HashMap.h"

// Functions
// Add the log-scaled values l1 and l2 using a transform to avoid
// precision errors
inline double addLogs(const double l1, const double l2)
{
    if (l1>l2) {
        double diff=l2-l1;
        return l1+log(1.0+exp(diff));
    } else {
        double diff=l1-l2;
        return l2+log(1.0+exp(diff));
    }
}


bool applyFilter(const std::string& kmer, const BWTIndexSet& indices);

// Get the number of times the kmer appears in each samples
std::vector<size_t> getPopulationCoverageCount(const std::string& kmer, const BWTIndexSet& indices);

// Get the mean depth of a random k-mer in each sample
std::vector<double> getSampleMeanKmerDepth(size_t k, const BWTIndexSet& indices);

//
double LMHaploid(double d, const std::vector<size_t>& sample_count);
double LMHaploidNonUniform(const std::vector<double>& depths, const std::vector<size_t>& sample_count);
double LMDiploid(double d, const std::vector<size_t>& sample_count);
double LMDiploidNonUniform(const std::vector<double>& depths, const std::vector<size_t>& sample_count);

double _haploidNonUniform(const std::vector<double>& depths, const std::vector<size_t>& sample_count, double M);
double _diploidNonUniform(const std::vector<double>& depths, const std::vector<size_t>& sample_count, double M);

//
std::vector<double> loadSampleDepthsFromFile(std::string filename);

// Simulation functions
std::vector<size_t> simulateCoverageNull(double d, size_t o, size_t N);
std::vector<size_t> simulateCoverageHaploid(double d, size_t o, size_t N);
std::vector<size_t> simulateCoverageDiploid(double d, size_t o, size_t N);
std::vector<size_t> simulateCoverageNullNonUniform(std::vector<double> sample_coverage, size_t o);
std::vector<size_t> simulateCoverageHaploidNonUniform(std::vector<double> sample_coverage, size_t o, double& out_num_carriers);
std::vector<size_t> simulateCoverageDiploidNonUniform(std::vector<double> sample_coverage, size_t o, double& out_num_carriers);
std::vector<double> simulateSampleKmerDepths(size_t mean_depth, size_t N);

//
void runSimulation();

//
// Getopt
//
#define SUBPROGRAM "haplotype-filter"
static const char *HAPLOTYPE_FILTER_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2013 Wellcome Trust Sanger Institute\n";

static const char *HAPLOTYPE_FILTER_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... HAPLOTYPE_FILE VCF_FILE\n"
"Remove haplotypes and their associated variants from a data set.\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -r, --reads=FILE                 load the FM-index of the reads in FILE\n"
"          --reference=STR              load the reference genome from FILE\n"
"          --haploid                    force use of the haploid model\n"
"      -o, --out-prefix=STR             write the passed haplotypes and variants to STR.vcf and STR.fa\n" 
"      -t, --threads=NUM                use NUM threads to compute the overlaps (default: 1)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static bool bHaploid = false;
    static int numThreads = 1;
    static size_t k = 31;
    static std::string readsFile;
    static std::string outFile = "hapfilter.vcf";
    static std::string haplotypeFile;
    static std::string vcfFile;
    static std::string referenceFile;
}

static const char* shortopts = "o:r:k:t:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_REFERENCE, OPT_HAPLOID };

static const struct option longopts[] = {
    { "verbose",               no_argument,       NULL, 'v' },
    { "threads",               required_argument, NULL, 't' },
    { "outfile",               required_argument, NULL, 'o' },
    { "reads",                 required_argument, NULL, 'r' },
    { "reference",             required_argument, NULL, OPT_REFERENCE },
    { "haploid",               no_argument,       NULL, OPT_HAPLOID },
    { "help",                  no_argument,       NULL, OPT_HELP },
    { "version",               no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int haplotypeFilterMain(int argc, char** argv)
{
    parseHaplotypeFilterOptions(argc, argv);
    
    //runSimulation();
    //exit(0);

    Timer* pTimer = new Timer(PROGRAM_IDENT);

    // Load FM-index of the reads
    std::cout << "Loading FM-index of " << opt::readsFile << "...";
    BWTIndexSet indices;
    std::string prefix = stripGzippedExtension(opt::readsFile);
    indices.pBWT = new BWT(prefix + BWT_EXT);
    indices.pPopIdx = new PopulationIndex(prefix + POPIDX_EXT);
    indices.pSSA = new SampledSuffixArray(prefix + SAI_EXT, SSA_FT_SAI);
    indices.pQualityTable = new QualityTable();
    std::cout << "done\n";

    // Load FM-index of the reference
    std::cout << "Loading reference index... " << std::flush;
    BWTIndexSet referenceIndex;
    std::string refPrefix = stripExtension(opt::referenceFile);
    referenceIndex.pBWT = new BWT(refPrefix + BWT_EXT);
    std::cout << "done" << std::endl;

    //std::vector<double> depths = loadSampleDepthsFromFile("depths.k31.txt");
    std::vector<double> depths = getSampleMeanKmerDepth(opt::k, indices);

    // Read haplotypes into a kmer hash
    size_t assembly_k = 61; // hack hack
    HashMap<std::string, size_t> kmer_to_haplotype;
    StringVector haplotypes;

    std::cout << "Loading haplotypes\n";
    SeqReader haplotypeReader(opt::haplotypeFile);
    SeqRecord hapRecord;
    while(haplotypeReader.get(hapRecord))
    {
        std::string haplotype = hapRecord.seq.toString();
        haplotypes.push_back(haplotype);

        size_t nk = haplotype.size() - assembly_k + 1;
        for(size_t i = 0; i < nk; ++i)
        {
            std::string kmer = haplotype.substr(i, assembly_k);
            kmer_to_haplotype[kmer] = haplotypes.size() - 1;
        }
    }

    std::ofstream outFile(opt::outFile.c_str());
    std::ifstream inFile(opt::vcfFile.c_str());

    std::string line;
    while(getline(inFile, line))
    {
        assert(line.size() > 0);
        if(line[0] == '#')
        {
            if(line[1] != '#')
            {
                outFile << "##INFO=<ID=LM,Number=1,Type=Float,Description=\"Log-likelihood ratio statistic using diploid segregation model\">" << "\n";
                outFile << "##INFO=<ID=O,Number=1,Type=Integer,Description=\"Number of reads used in segregation test\">" << "\n";
            }
            outFile << line << "\n";
            continue;
        }

        StringVector fields = split(line, '\t');
        std::string vcf_kmer = fields[2];
        
        // Load the haplotype with this kmer
        HashMap<std::string, size_t>::iterator iter = kmer_to_haplotype.find(vcf_kmer);
        if(iter == kmer_to_haplotype.end())
            iter = kmer_to_haplotype.find(reverseComplement(vcf_kmer));

        assert(iter != kmer_to_haplotype.end());
        std::string& haplotype = haplotypes[iter->second];
        
        printf("Kmer --- %s\n", vcf_kmer.c_str());
        printf("Haplotype --- %s\n", haplotype.c_str());

        // Find the highest-depth non-reference kmer to use to calculate the segregation stats
        size_t best_index = 0;
        size_t best_count = 0;
        size_t nk = haplotype.size() - opt::k + 1;
        for(size_t i = 0; i < nk; ++i)
        {
            std::string seg_kmer = haplotype.substr(i, opt::k);
            size_t ref_c = BWTAlgorithms::countSequenceOccurrences(seg_kmer, referenceIndex);
            printf("seg_kmer --- %s ref_c? %zu\n", seg_kmer.c_str(), ref_c);
            if(ref_c == 0)
            {
                size_t read_c = BWTAlgorithms::countSequenceOccurrences(seg_kmer, indices);
                if(read_c > best_count)
                {
                    best_count = read_c;
                    best_index = i;
                }
                printf("read_c: %zu\n", read_c);
            }
        }
        
        double LM = 0.f;
        size_t total_coverage = 0;
        if(best_count > 0)
        {
            std::string kmer = haplotype.substr(best_index, opt::k);
            std::vector<size_t> sample_coverage = getPopulationCoverageCount(kmer, indices);
            std::copy(sample_coverage.begin(), sample_coverage.end(), std::ostream_iterator<size_t>(std::cout, " "));
            std::cout << "\n";
          
            for(size_t i = 0; i < sample_coverage.size(); ++i)
                total_coverage += sample_coverage[i];
            
            if(opt::bHaploid)
                LM = LMHaploidNonUniform(depths, sample_coverage);
            else
                LM = LMDiploidNonUniform(depths, sample_coverage);
        }

        std::stringstream lmss;
        lmss << fields[7];
        lmss << ";LM=" << LM << ";";
        lmss << "O=" << total_coverage << ";";
        fields[7] = lmss.str();
        for(size_t i = 0; i < fields.size()-1; ++i)
            outFile << fields[i] << "\t";
        outFile << fields[fields.size()-1];
        outFile << "\n";
    }
    
    // Cleanup
    delete indices.pBWT;
    delete indices.pPopIdx;
    delete indices.pQualityTable;
    delete referenceIndex.pBWT;
    delete pTimer;

    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

//
void runSimulation()
{
    size_t trials = 500;
    size_t N = 1000;
    double d = 3;

    size_t min_o = 5;
    size_t max_o = (size_t)(d * N);
    srand(time(NULL));
    bool use_diploid = true;

    double fpr = 0.5;
    size_t n_tp = 0;
    size_t n_fp = 0;
    size_t n_tn = 0;
    size_t n_fn = 0;

    //std::vector<double> depths(N, d);
    //std::vector<double> depths = simulateSampleKmerDepths(d, N);
    std::vector<double> depths = loadSampleDepthsFromFile("depths.txt");

    if(opt::verbose > 0)
    {
        printf("Per-sample depth:\n");
        std::copy(depths.begin(), depths.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << "\n";
    }

    for(size_t i = 0; i < trials; ++i)
    {
        double p = (double)rand() / RAND_MAX;
        size_t o = (rand() % (max_o - min_o)) + min_o;
        std::vector<size_t> counts;
        bool is_fp = p < fpr;
        double dbg_num_carriers = 0.0f;
        if(is_fp)
            counts = simulateCoverageNullNonUniform(depths, o);
        else
            counts = simulateCoverageDiploidNonUniform(depths, o, dbg_num_carriers);

        if(opt::verbose > 0)
        {
            std::cout << "Read counts:\n";
            std::copy(counts.begin(), counts.end(), std::ostream_iterator<size_t>(std::cout, " "));
            std::cout << "\n";
        }

        double s_h = 0.0f;
        double s_d_nu = LMDiploidNonUniform(depths, counts);
        double s_h_nu = LMHaploidNonUniform(depths, counts);
        double s_d = LMDiploid(d, counts);
        double s = use_diploid ? s_d : s_h;

        if(is_fp && s < 0)
            n_tn += 1;
        else if(is_fp && s > 0)
            n_fp += 1;
        else if(!is_fp && s < 0)
            n_fn += 1;
        else
            n_tp += 1;
        printf("%zu\t%d\t%zu\t%lf\t%lf\t%lf\t%lf\n", i, is_fp, o, s_d, s_h, s_h_nu, s_d_nu);
    }
    fprintf(stderr, "tp: %zu fp: %zu tn: %zu fn: %zu\n", n_tp, n_fp, n_tn, n_fn);
}

//
bool applyFilter(const std::string& kmer, const BWTIndexSet& indices)
{
    printf("Kmer --- %s\n", kmer.c_str());
    std::vector<size_t> sample_coverage = getPopulationCoverageCount(kmer, indices);
    std::copy(sample_coverage.begin(), sample_coverage.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cout << "\n";

    double d = 4;
    double lr = LMHaploid(d, sample_coverage);
    return lr > 0;
}

double LMHaploid(double d, const std::vector<size_t>& sample_count)
{
    size_t N = sample_count.size();
    double o = 0;
    size_t n0 = 0;
    for(size_t i = 0; i < N; ++i)
    {
        if(sample_count[i] == 0)
            n0 += 1;

        o += sample_count[i];
    }

    //double M = std::max(round(o/d), N - n_0);
    double M = round(o/d);
    // Clamp M at the number of samples
    M = std::min((double)N, M);

    double p = o / N;
    double q = o / M;
    
    double sum = 0;
    for(size_t i = 0; i < sample_count.size(); ++i)
    {
        size_t oi = sample_count[i];
        if(oi == 0)
        {
            //sum += p - q;
            double h1_zero = (M/N) * exp(-q) + ((N - M) / N);
            double log_h0_zero = -p;
            double log_pzero = log(h1_zero) - log_h0_zero;
            sum += log_pzero;
        }
        else
        {
            sum += log(M/N) + oi * log(q) - oi * log(p) + p - q;
        }
    }

    if(opt::verbose > 0)
    {
        printf("LMHaploid:\n");
        printf("\to: %lf M: %lf N: %zu p: %lf q: %lf n0: %zu\n", o, M, N, p, q, n0);
    }

    return sum;
}

double LMHaploidNonUniform(const std::vector<double>& depths, const std::vector<size_t>& sample_count)
{
    size_t N = sample_count.size();
    double o = 0;
    size_t n0 = 0;
    double sum_depth = 0.0;
    for(size_t i = 0; i < N; ++i)
    {
        if(sample_count[i] == 0)
            n0 += 1;

        o += sample_count[i];
        sum_depth += depths[i];
    }
    
    double mean_depth = sum_depth / N;
    double M1 = N - n0;
    double M2 = o / mean_depth;
    double M_est = std::max(M1, M2);
    if(opt::verbose)
        printf("HAP ParamEst o: %zu n0: %zu M1: %lf M2: %lf M_est %lf\n", (size_t)o, n0, M1, M2, M_est);

    return _haploidNonUniform(depths, sample_count, M2);
}

// 
double _haploidNonUniform(const std::vector<double>& depths, const std::vector<size_t>& sample_count, double M)
{
    size_t N = sample_count.size();
    assert(depths.size() == N);

    double o = 0;
    size_t n0 = 0;

    std::vector<double> mean_p(N);
    std::vector<double> mean_q(N);
    double sum_depth = 0.0f;
    for(size_t i = 0; i < N; ++i)
    {
        if(sample_count[i] == 0)
            n0 += 1;

        o += sample_count[i];
        sum_depth += depths[i];
    }

    double expected_reads_per_allele_null = o / N;
    double expected_reads_per_allele_alt = o / M;
    double mean_depth = sum_depth / N;

    // Compute the expected mean depth under the NULL and H1
    for(size_t i = 0; i < N; ++i)
    {
        mean_p[i] = (depths[i] /  mean_depth) * expected_reads_per_allele_null;
        mean_q[i] = (depths[i] /  mean_depth) * expected_reads_per_allele_alt;
    }

    if(opt::verbose > 0)
    {
        printf("HAP_NU Parameters -- mean_p0 %lf\n", mean_p[0]);
        printf("HAP_NU Parameters -- mean_q0 %lf\n", mean_q[0]);
        printf("HAP_NU Parameters -- M %lf\n", M);
    }

    //double M = std::max(round(o/d), N - n_0);
    //M = round(o / mean_depth);

    // Clamp M at the number of samples
    //M = std::min((double)N, M);

    double sum = 0;
    double sum_alt = 0.0f;
    double sum_null = 0.0f;

    for(size_t i = 0; i < sample_count.size(); ++i)
    {
        double L_NULL = 0.0f;
        double L_ALT = 0.0f;

        double p = mean_p[i];
        double q = mean_q[i];
        size_t oi = sample_count[i];
        if(oi == 0)
        {
            double h1_zero = (M/N) * exp(-q) + ((N - M) / N);
            double log_h0_zero = -p;
            double log_pzero = log(h1_zero) - log_h0_zero;
            sum += log_pzero;
            
            L_NULL = log_h0_zero;
            L_ALT = log(h1_zero);
        }
        else
        {
            L_NULL = SGAStats::logPoisson(oi, p);
            L_ALT = (log(M/N) + SGAStats::logPoisson(oi, q));
            sum += L_ALT - L_NULL;
        }

        sum_alt += L_ALT;
        sum_null += L_NULL;

        if(opt::verbose > 1)
            printf("HAP_NU o: %zu LA: %lf LN: %lf S: %lf\n", oi, L_ALT, L_NULL, sum);
    }
    
    if(opt::verbose > 0)
        printf("FINAL_HAP_NU M: %lf O: %lf L_A: %lf L_N: %lf S: %lf\n", M, o, sum_alt, sum_null, sum);        
    return sum;
}

double LMDiploid(double d, const std::vector<size_t>& sample_count)
{
    size_t N = sample_count.size();
    double o = 0;
    size_t n0 = 0;
    for(size_t i = 0; i < N; ++i)
    {
        if(sample_count[i] == 0)
            n0 += 1;

        o += sample_count[i];
    }

    double M = round(o/d);
    double f = M / (2 * N);

    // Clamp M at the number of samples
    M = std::min((double)N, M);

    double p = o / N;
    double q = o / M;

    double sum = 0;
    for(size_t i = 0; i < sample_count.size(); ++i)
    {
        size_t oi = sample_count[i];
        if(oi == 0)
        {
            double h1_zero = pow(1 - f, 2) + 2 * f * (1-f) * exp(-q) + pow(f, 2) * exp(-2*q);
            double log_h0_zero = -p;
            double log_pzero = log(h1_zero) - log_h0_zero;
            sum += log_pzero;
        }
        else
        {
            double log_null = SGAStats::logPoisson(oi, p);
            double t1 = 2 * f * (1 - f) * pow(q, oi) * exp(-q) + pow(f, 2) * pow(2 * q, oi) * exp(-2 * q);
            double t2 = SGAStats::logFactorial(oi);
            double t3 = log(t1) - t2;
            sum += t3 - log_null;
        }
    }

    if(opt::verbose > 0)
    {
        printf("LMDiploid:\n");
        printf("\to: %lf M: %lf N: %zu p: %lf q: %lf n0: %zu f: %lf\n", o, M, N, p, q, n0, f);
    }

    return sum;
}

//
double LMDiploidNonUniform(const std::vector<double>& depths, const std::vector<size_t>& sample_count)
{
    size_t N = sample_count.size();
    double o = 0;
    size_t n0 = 0;
    double sum_depth = 0.0;
    double M3 = 0.0f;

    for(size_t i = 0; i < N; ++i)
    {
        if(sample_count[i] == 0)
            n0 += 1;

        o += sample_count[i];
        sum_depth += depths[i];
        M3 += sample_count[i] / depths[i];
    }
    M3 *= 2;
    
    double mean_depth = sum_depth / (2 * N);
    
    double M1 = N - n0;
    double M2 = o / mean_depth;
    double M_est = std::max(M1, M2);
    if(opt::verbose)
        printf("DIP ParamEst o: %zu n0: %zu M1: %lf M2: %lf M3: %lf M_est: %lf\n", (size_t)o, n0, M1, M2, M3, M_est);
    
    return _diploidNonUniform(depths, sample_count, M3);
}

// 
double _diploidNonUniform(const std::vector<double>& depths, const std::vector<size_t>& sample_count, double M)
{
    size_t N = sample_count.size();
    assert(depths.size() == N);

    double f = M / (2 * N);
    double hom_alt_prop = pow(f, 2.0f);
    double het_prop = 2 * f * (1.0f - f);
    double hom_ref_prop = pow(1.0 - f, 2.0f);

    double o = 0;
    double sum_depth = 0.0f;
    std::vector<double> mean_p(N);
    std::vector<double> mean_q(N);
    
    for(size_t i = 0; i < N; ++i)
    {
        o += sample_count[i];
        sum_depth += depths[i];
    }

    double expected_reads_per_allele_null = o / N;
    double expected_reads_per_allele_alt = o / M;
    double mean_depth = sum_depth / N;

    // Compute the expected mean depth under the NULL and H1
    double sum_p = 0.0f;
    double sum_q = 0.0f;

    for(size_t i = 0; i < N; ++i)
    {
        mean_p[i] = (depths[i] /  mean_depth) * expected_reads_per_allele_null;
        //mean_q[i] = o / M;
        mean_q[i] = (depths[i] /  mean_depth) * expected_reads_per_allele_alt;

        sum_p += mean_p[i];
        sum_q += mean_q[i];
    }


    if(opt::verbose)
    {
        printf("DIP_NU Parameters -- M: %lf\n", M);
        printf("DIP_NU parameters -- mean depth: %lf\n", mean_depth);
        printf("DIP_NU parameters -- exp reads null: %lf\n", expected_reads_per_allele_null);
        printf("DIP_NU parameters -- exp reads alt: %lf\n", expected_reads_per_allele_alt);
        printf("DIP_NU parameters -- sum_p: %lf\n", sum_p);
        printf("DIP_NU parameters -- sum_q: %lf\n", sum_q);
        printf("DIP_NU Parameters -- mean_p0: %lf\n", mean_p[0]);
        printf("DIP_NU Parameters -- mean_q0: %lf\n", mean_q[0]);
        printf("DIP_NU Parameters -- f: %lf\n", f);
        printf("DIP_NU Parameters -- hom_alt_prop: %lf\n", hom_alt_prop);
        printf("DIP_NU Parameters -- het_prop: %lf\n", het_prop);
        printf("DIP_NU Parameters -- hom_ref: %lf\n", hom_ref_prop);
    }

    double sum = 0;
    double sum_alt = 0.0f;
    double sum_null = 0.0f;

    for(size_t i = 0; i < sample_count.size(); ++i)
    {
        double p = mean_p[i];
        double q = mean_q[i];
        size_t oi = sample_count[i];
        double L_ALT = 0.0f;
        double L_NULL = 0.0f;
        if(oi == 0)
        {
            double h1_zero = hom_ref_prop + het_prop * exp(-q) + hom_alt_prop * exp(-2*q);
            double log_h0_zero = -p;
            double log_pzero = log(h1_zero) - log_h0_zero;
            sum += log_pzero;
            L_NULL = log_h0_zero;
            L_ALT = log(h1_zero);
        }
        else
        {
            double log_null = SGAStats::logPoisson(oi, p);
            double t1 = het_prop * pow(q, oi) * exp(-q) + hom_alt_prop * pow(2 * q, oi) * exp(-2 * q);
            double t2 = SGAStats::logFactorial(oi);
            double t3 = log(t1) - t2;
            sum += t3 - log_null;
            L_NULL = log_null;
            L_ALT = t3;
        }

        sum_alt += L_ALT;
        sum_null += L_NULL;

        if(opt::verbose > 1)
            printf("DIP_NU i %zu d: %lf p: %lf q: %lf o: %zu LA: %lf LN: %lf S: %lf\n", i, depths[i], p, q, oi, L_ALT, L_NULL, sum);
    }

    if(opt::verbose > 0)
        printf("FINAL_DIP_NU M: %lf O: %lf f: %lf L_A: %lf L_N: %lf S: %lf\n", M, o, f, sum_alt, sum_null, sum);

    return sum;
}

std::vector<size_t> simulateCoverageNull(double /*d*/, size_t o, size_t N)
{
    // Distribute o reads over the N individuals
    std::vector<size_t> counts(N);
    for(size_t i = 0; i < o; ++i)
    {
        // Choose a random index and give it the count
        counts[rand() % (size_t)N]++;
    }
    return counts;
}

std::vector<size_t> simulateCoverageHaploid(double d, size_t o, size_t N)
{
    double M = round(o / d);
    if(opt::verbose > 0)
        printf("Simulate haploid %zu reads across %lf samples\n", o, M);
    assert(M <= N);

    // Select M indices at random to receive the counts
    std::vector<size_t> indices(N);
    for(size_t i = 0; i < N; ++i)
        indices[i] = i;

    std::random_shuffle(indices.begin(), indices.end());

    // Distribute o reads over M individuals
    std::vector<size_t> counts(N);
    for(size_t i = 0; i < o; ++i)
    {
        size_t idx = indices[rand() % (size_t)M];
        counts[idx]++;
    }

    return counts;
}

std::vector<size_t> simulateCoverageDiploid(double d, size_t o, size_t N)
{
    // Here M is the number of non-reference alleles
    double M = round(o / d);
    assert(M <= N);
    
    // Proportion of non-reference alleles in the population
    double q = M / (2 * N);
    size_t num_hom_alt = (size_t)round(pow(q, 2) * N);
    size_t num_het = (size_t)(M - (2*num_hom_alt));
    
    if(opt::verbose > 0)
    {
        printf("Simulate diploid %zu reads across %lf alleles\n", o, M);
        printf("    q: %lf num_hom_alt: %zu num_het: %zu\n", q, num_hom_alt, num_het);
    }

    // Select indices at random to become carries of the variant
    std::vector<size_t> initial_indices(N);
    for(size_t i = 0; i < N; ++i)
        initial_indices[i] = i;
    std::random_shuffle(initial_indices.begin(), initial_indices.end());

    // The first num_hom_alt entries will carry two alleles
    std::vector<size_t> final_indices((size_t)M);
    for(size_t i = 0; i < num_hom_alt; ++i)
    {
        final_indices[2*i] = initial_indices[i];
        final_indices[2*i + 1] = initial_indices[i];
    }

    // The next num_het entries carry one allele
    for(size_t i = 0; i < num_het; ++i)
        final_indices[2*num_hom_alt + i] = initial_indices[num_hom_alt + i];    

    // Distribute o reads over M individuals
    std::vector<size_t> counts(N);
    for(size_t i = 0; i < o; ++i)
    {
        size_t idx = final_indices[rand() % (size_t)M];
        counts[idx]++;
    }

    return counts;
}

//
std::vector<size_t> simulateCoverageNullNonUniform(std::vector<double> sample_coverage, size_t o)
{
    size_t N = sample_coverage.size();
    
    if(opt::verbose > 0)
        printf("SimulateNullNonUniform %zu reads\n", o);

    // Make a discrete distribution using the per-sample coverage depths
    std::vector<double> distribution(N);
    double sum = 0.f;
    for(size_t i = 0; i < N; ++i)
    {
        double p = sample_coverage[i];
        sum += p;
        distribution[i] = sum;
    }

    // Normalize
    for(size_t i = 0; i < N; ++i)
        distribution[i] /= sum;

    // Distribute o reads over the individuals according to the read depth
    std::vector<size_t> counts(N);
    for(size_t i = 0; i < o; ++i)
    {
        double p = (double)rand() / RAND_MAX;
        size_t j = 0;
        while(p > distribution[j])
            j += 1;
        counts[j]++;
    }
    return counts;
}

//
std::vector<size_t> simulateCoverageHaploidNonUniform(std::vector<double> sample_coverage, size_t o, double& out_M)
{
    size_t N = sample_coverage.size();
    double mean_coverage = 0.0f;
    for(size_t i = 0; i < N; ++i)
        mean_coverage += sample_coverage[i];
    mean_coverage /= N;

    // Here M is the number of carriers
    double M = round(o / mean_coverage);
    M = std::min(M, (double)N);
    assert(M <= N);
    
    if(opt::verbose > 0)
        printf("SimulateNonUniform haploid %zu reads across %lf carriers\n", o, M);

    // Select indices at random to become carriers of the variant
    std::vector<size_t> initial_indices(N);
    for(size_t i = 0; i < N; ++i)
        initial_indices[i] = i;
    std::random_shuffle(initial_indices.begin(), initial_indices.end());

    // The first num_hom_alt entries will carry two alleles
    std::vector<size_t> allele_indices((size_t)M);
    for(size_t i = 0; i < M; ++i)
        allele_indices[i] = initial_indices[i];

    // Make a discrete distribution using the per-sample coverage depths
    std::vector<double> distribution(allele_indices.size());
    double sum = 0.f;
    for(size_t i = 0; i < allele_indices.size(); ++i)
    {
        double p = sample_coverage[allele_indices[i]] / 2;
        distribution[i] = sum + p;
        sum = distribution[i];
    }

    // Normalize
    for(size_t i = 0; i < allele_indices.size(); ++i)
        distribution[i] /= sum;

    // Distribute o reads over the individuals according to the read depth
    std::vector<size_t> counts(N);
    for(size_t i = 0; i < o; ++i)
    {
        double p = (double)rand() / RAND_MAX;
        size_t j = 0; 
        while(p > distribution[j])
            j += 1;

        size_t idx = allele_indices[j];
        counts[idx]++;
    }
    out_M = M;
    return counts;
}


//
std::vector<size_t> simulateCoverageDiploidNonUniform(std::vector<double> sample_coverage, size_t o, double& out_M)
{
    size_t N = sample_coverage.size();
    double mean_coverage_per_allele = 0.0f;
    for(size_t i = 0; i < N; ++i)
        mean_coverage_per_allele += sample_coverage[i];
    mean_coverage_per_allele = mean_coverage_per_allele / (2 * N);

    // Here M is the number of non-reference alleles
    double M = round(o / mean_coverage_per_allele);
    M = std::min(M, (double)2*N);
    
    // Proportion of non-reference alleles in the population
    double q = M / (2 * N);
    size_t num_hom_alt = (size_t)round(pow(q, 2) * N);
    size_t num_het = (size_t)(M - (2*num_hom_alt));
    
    if(opt::verbose > 0)
    {
        printf("SimulateNonUniformDiploid %zu reads across %lf alleles\n", o, M);
        printf("    q: %lf num_hom_alt: %zu num_het: %zu\n", q, num_hom_alt, num_het);
    }

    // Select indices at random to become carriers of the variant
    std::vector<size_t> initial_indices(N);
    for(size_t i = 0; i < N; ++i)
        initial_indices[i] = i;
    std::random_shuffle(initial_indices.begin(), initial_indices.end());

    std::vector<size_t> allele_count(N);

    // The first num_hom_alt entries will carry two alleles
    std::vector<size_t> allele_indices((size_t)M);
    for(size_t i = 0; i < num_hom_alt; ++i)
    {
        allele_indices[2*i] = initial_indices[i];
        allele_indices[2*i + 1] = initial_indices[i];

        allele_count[initial_indices[i]] = 2;
    }

    // The next num_het entries carry one allele
    for(size_t i = 0; i < num_het; ++i)
    {
        allele_indices[2*num_hom_alt + i] = initial_indices[num_hom_alt + i];    
        allele_count[initial_indices[num_hom_alt + i]] = 1;
    }

    if(opt::verbose > 0)
    {
        std::cout << "Allele Counts:\n";
        std::copy(allele_count.begin(), allele_count.end(), std::ostream_iterator<size_t>(std::cout, " "));
        std::cout << "\n";
    }

    // Make a discrete distribution using the per-sample coverage depths
    std::vector<double> distribution(allele_indices.size());
    double sum = 0.f;
    for(size_t i = 0; i < allele_indices.size(); ++i)
    {
        double p = sample_coverage[allele_indices[i]] / 2;
        distribution[i] = sum + p;
        sum = distribution[i];
    }

    // Normalize
    for(size_t i = 0; i < allele_indices.size(); ++i)
        distribution[i] /= sum;

    // Distribute o reads over the individuals according to the read depth
    std::vector<size_t> counts(N);
    for(size_t i = 0; i < o; ++i)
    {
        double p = (double)rand() / RAND_MAX;
        size_t j = 0; 
        while(p > distribution[j])
            j += 1;

        size_t idx = allele_indices[j];
        counts[idx]++;
    }
    out_M = num_het + num_hom_alt;
    return counts;
}

std::vector<double> loadSampleDepthsFromFile(std::string filename)
{
    std::ifstream stream(filename.c_str());
    if(!stream.good())
    {
        std::cerr << "Could not load depths from file " << filename << "\n";
        exit(EXIT_FAILURE);
    }

    std::vector<double> depths;
    double d;
    while(stream >> d)
        depths.push_back(d);
    return depths;
}

std::vector<double> simulateSampleKmerDepths(size_t mean_depth, size_t N)
{
    std::vector<double> depths(N);
    for(size_t i = 0; i < N; ++i)
    {
        // uniform
        //depths[i] = ((double)rand() / RAND_MAX) * mean_depth + mean_depth / 2;
        
        // exponential
        double p = (double)rand() / RAND_MAX;
        depths[i] = log(1 - p) * -1 * mean_depth; 
    }
    return depths;
}

//
std::vector<size_t> getPopulationCoverageCount(const std::string& kmer, const BWTIndexSet& indices)
{
    //
    StringVector kmers;
    kmers.push_back(kmer);
    kmers.push_back(reverseComplement(kmer));
    SeqRecordVector records;
    HapgenUtil::extractHaplotypeReads(kmers, indices, opt::k, false, -1, 100000, &records, NULL);

    // Count the number of reads per samples using their index in the BWT
    std::vector<size_t> sample_counts(indices.pPopIdx->getNumSamples());

    for(size_t i = 0; i < records.size(); ++i)
    {
        std::stringstream parser(records[i].id.substr(4));
        size_t index;
        parser >> index;
        size_t s_idx = indices.pPopIdx->getSampleIndex(index);
        sample_counts[s_idx]++;
    }
    return sample_counts;
}

//
std::vector<double> getSampleMeanKmerDepth(size_t k, const BWTIndexSet& indices)
{
    size_t N = indices.pPopIdx->getNumSamples();
    std::vector<double> average_counts(N);
    
    printf("starting sampling\n");
    size_t target_samples = 1000;
    size_t used_samples = 0;
    for(size_t i = 0; i < target_samples; ++i)
    {
        if(i % 100 == 0)
            printf("sampling iteration %zu\n", i);

        std::string r_str = BWTAlgorithms::sampleRandomString(indices.pBWT);

        // We use the first k-mer in the read
        std::string test_kmer = r_str.substr(0, k);

        // Screen out ultra-low depth k-mers
        size_t count = BWTAlgorithms::countSequenceOccurrences(test_kmer, indices);
        if(count < 5)
            continue;

        //
        std::vector<size_t> incoming_counts = getPopulationCoverageCount(test_kmer, indices);

        for(size_t i = 0; i < incoming_counts.size(); ++i)
            average_counts[i] += incoming_counts[i];
        used_samples += 1;
    }

    printf("average counts:\n");
    for(size_t i = 0; i < average_counts.size(); ++i)
    {
        average_counts[i] /= used_samples;
        printf("%.2lf ", average_counts[i]);
    }
    printf("\n");

    return average_counts;
}

// 
// Handle command line arguments
//
void parseHaplotypeFilterOptions(int argc, char** argv)
{
    std::string algo_str;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'o': arg >> opt::outFile; break;
            case 't': arg >> opt::numThreads; break;
            case 'r': arg >> opt::readsFile; break;
            case OPT_REFERENCE: arg >> opt::referenceFile; break;
            case OPT_HAPLOID: opt::bHaploid = true; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << HAPLOTYPE_FILTER_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << HAPLOTYPE_FILTER_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 2) 
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } 
    else if (argc - optind > 2) 
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if(opt::numThreads <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
        die = true;
    }

    if(opt::referenceFile.empty())
    {
        std::cerr << SUBPROGRAM ": a reference file must be given: " << opt::referenceFile << "\n";
        die = true;
    }

    opt::haplotypeFile = argv[optind++];
    opt::vcfFile = argv[optind++];

    if (die) 
    {
        std::cout << "\n" << HAPLOTYPE_FILTER_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}
