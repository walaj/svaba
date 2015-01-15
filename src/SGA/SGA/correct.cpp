//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// correct - Correct sequencing errors in reads using the FM-index
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "correct.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "ASQG.h"
#include "gzstream.h"
#include "SequenceProcessFramework.h"
#include "ErrorCorrectProcess.h"
#include "CorrectionThresholds.h"
#include "KmerDistribution.h"
#include "BWTIntervalCache.h"
#include "LRAlignment.h"

// Functions
int learnKmerParameters(const BWT* pBWT);

//
// Getopt
//
#define SUBPROGRAM "correct"
static const char *CORRECT_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *CORRECT_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Correct sequencing errors in all the reads in READSFILE\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -p, --prefix=PREFIX              use PREFIX for the names of the index files (default: prefix of the input file)\n"
"      -o, --outfile=FILE               write the corrected reads to FILE (default: READSFILE.ec.fa)\n"
"      -t, --threads=NUM                use NUM threads for the computation (default: 1)\n"
"          --discard                    detect and discard low-quality reads\n"
"      -d, --sample-rate=N              use occurrence array sample rate of N in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
"      -a, --algorithm=STR              specify the correction algorithm to use. STR must be one of kmer, hybrid, overlap. (default: kmer)\n"
"          --metrics=FILE               collect error correction metrics (error rate by position in read, etc) and write them to FILE\n"
"\nKmer correction parameters:\n"
"      -k, --kmer-size=N                The length of the kmer to use. (default: 31)\n"
"      -x, --kmer-threshold=N           Attempt to correct kmers that are seen less than N times. (default: 3)\n"
"      -i, --kmer-rounds=N              Perform N rounds of k-mer correction, correcting up to N bases (default: 10)\n"
"          --learn                      Attempt to learn the k-mer correction threshold (experimental). Overrides -x parameter.\n"
"\nOverlap correction parameters:\n"
"      -e, --error-rate                 the maximum error rate allowed between two sequences to consider them overlapped (default: 0.04)\n"
"      -m, --min-overlap=LEN            minimum overlap required between two reads (default: 45)\n"
"      -c, --conflict=INT               use INT as the threshold to detect a conflicted base in the multi-overlap (default: 5)\n"
"      -l, --seed-length=LEN            force the seed length to be LEN. By default, the seed length in the overlap step\n"
"                                       is calculated to guarantee all overlaps with --error-rate differences are found.\n"
"                                       This option removes the guarantee but will be (much) faster. As SGA can tolerate some\n"
"                                       missing edges, this option may be preferable for some data sets.\n"
"      -s, --seed-stride=LEN            force the seed stride to be LEN. This parameter will be ignored unless --seed-length\n"
"                                       is specified (see above). This parameter defaults to the same value as --seed-length\n"
"      -b, --branch-cutoff=N            stop the overlap search at N branches. This parameter is used to control the search time for\n"
"                                       highly-repetitive reads. If the number of branches exceeds N, the search stops and the read\n"
"                                       will not be corrected. This is not enabled by default.\n"
"      -r, --rounds=NUM                 iteratively correct reads up to a maximum of NUM rounds (default: 1)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static int numOverlapRounds = 1;
    static std::string prefix;
    static std::string readsFile;
    static std::string outFile;
    static std::string discardFile;
    static std::string metricsFile;
    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
    
    static double errorRate = 0.04;
    static unsigned int minOverlap = DEFAULT_MIN_OVERLAP;
    static int seedLength = 0;
    static int seedStride = 0;
    static int conflictCutoff = 5;
    static int branchCutoff = -1;

    static int kmerLength = 31;
    static int kmerThreshold = 3;
    static int numKmerRounds = 10;
    static bool bLearnKmerParams = false;
    static int intervalCacheLength = 10;

    static ErrorCorrectAlgorithm algorithm = ECA_KMER;
}

static const char* shortopts = "p:m:d:e:t:l:s:o:r:b:a:c:k:x:i:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_METRICS, OPT_DISCARD, OPT_LEARN };

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "threads",       required_argument, NULL, 't' },
    { "min-overlap",   required_argument, NULL, 'm' },
    { "rounds",        required_argument, NULL, 'r' },
    { "outfile",       required_argument, NULL, 'o' },
    { "prefix",        required_argument, NULL, 'p' },
    { "error-rate",    required_argument, NULL, 'e' },
    { "seed-length",   required_argument, NULL, 'l' },
    { "seed-stride",   required_argument, NULL, 's' },
    { "algorithm",     required_argument, NULL, 'a' },
    { "sample-rate",   required_argument, NULL, 'd' },
    { "conflict",      required_argument, NULL, 'c' },
    { "branch-cutoff", required_argument, NULL, 'b' },
    { "kmer-size",     required_argument, NULL, 'k' },
    { "kmer-threshold",required_argument, NULL, 'x' },
    { "kmer-rounds",   required_argument, NULL, 'i' },
    { "learn",         no_argument,       NULL, OPT_LEARN },
    { "discard",       no_argument,       NULL, OPT_DISCARD },
    { "help",          no_argument,       NULL, OPT_HELP },
    { "version",       no_argument,       NULL, OPT_VERSION },
    { "metrics",       required_argument, NULL, OPT_METRICS },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int correctMain(int argc, char** argv)
{
    parseCorrectOptions(argc, argv);

    std::cout << "Correcting sequencing errors for " << opt::readsFile << "\n";

    // Load indices
    BWT* pBWT = new BWT(opt::prefix + BWT_EXT, opt::sampleRate);
    BWT* pRBWT = NULL;
    SampledSuffixArray* pSSA = NULL;

    if(opt::algorithm == ECA_OVERLAP || opt::algorithm == ECA_HYBRID)
        pSSA = new SampledSuffixArray(opt::prefix + SAI_EXT, SSA_FT_SAI);

    BWTIntervalCache* pIntervalCache = new BWTIntervalCache(opt::intervalCacheLength, pBWT);

    BWTIndexSet indexSet;
    indexSet.pBWT = pBWT;
    indexSet.pRBWT = pRBWT;
    indexSet.pSSA = pSSA;
    indexSet.pCache = pIntervalCache;

    // Learn the parameters of the kmer corrector
    if(opt::bLearnKmerParams)
    {
        int threshold = learnKmerParameters(pBWT);
        if(threshold != -1)
            CorrectionThresholds::Instance().setBaseMinSupport(threshold);
    }

    // Open outfiles and start a timer
    std::ostream* pWriter = createWriter(opt::outFile);
    std::ostream* pDiscardWriter = (!opt::discardFile.empty() ? createWriter(opt::discardFile) : NULL);
    Timer* pTimer = new Timer(PROGRAM_IDENT);
    pBWT->printInfo();

    // Set the error correction parameters
    ErrorCorrectParameters ecParams;
    ecParams.pOverlapper = NULL;
    ecParams.indices = indexSet;
    ecParams.algorithm = opt::algorithm;

    ecParams.minOverlap = opt::minOverlap;
    ecParams.numOverlapRounds = opt::numOverlapRounds;
    ecParams.minIdentity = 1.0f - opt::errorRate;
    ecParams.conflictCutoff = opt::conflictCutoff;

    ecParams.numKmerRounds = opt::numKmerRounds;
    ecParams.kmerLength = opt::kmerLength;
    ecParams.printOverlaps = opt::verbose > 0;

    // Setup post-processor
    bool bCollectMetrics = !opt::metricsFile.empty();
    ErrorCorrectPostProcess postProcessor(pWriter, pDiscardWriter, bCollectMetrics);

    if(opt::numThreads <= 1)
    {
        // Serial mode
        ErrorCorrectProcess processor(ecParams); 
        SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
                                                         ErrorCorrectResult, 
                                                         ErrorCorrectProcess, 
                                                         ErrorCorrectPostProcess>(opt::readsFile, &processor, &postProcessor);
    }
    else
    {
        // Parallel mode
        std::vector<ErrorCorrectProcess*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            ErrorCorrectProcess* pProcessor = new ErrorCorrectProcess(ecParams);
            processorVector.push_back(pProcessor);
        }
        
        SequenceProcessFramework::processSequencesParallel<SequenceWorkItem,
                                                           ErrorCorrectResult, 
                                                           ErrorCorrectProcess, 
                                                           ErrorCorrectPostProcess>(opt::readsFile, processorVector, &postProcessor);

        for(int i = 0; i < opt::numThreads; ++i)
        {
            delete processorVector[i];
        }
    }

    if(bCollectMetrics)
    {
        std::ostream* pMetricsWriter = createWriter(opt::metricsFile);
        postProcessor.writeMetrics(pMetricsWriter);
        delete pMetricsWriter;
    }

    delete pBWT;
    delete pIntervalCache;
    if(pRBWT != NULL)
        delete pRBWT;

    if(pSSA != NULL)
        delete pSSA;

    delete pTimer;
    
    delete pWriter;
    if(pDiscardWriter != NULL)
        delete pDiscardWriter;

    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

// Learn parameters of the kmer corrector
int learnKmerParameters(const BWT* pBWT)
{
    std::cout << "Learning kmer parameters\n";
    srand(time(0));
    size_t n_samples = 10000;

    //
    KmerDistribution kmerDistribution;
    int k = opt::kmerLength;
    for(size_t i = 0; i < n_samples; ++i)
    {
        std::string s = BWTAlgorithms::sampleRandomString(pBWT);
        int n = s.size();
        int nk = n - k + 1;
        for(int j = 0; j < nk; ++j)
        {
            std::string kmer = s.substr(j, k);
            int count = BWTAlgorithms::countSequenceOccurrences(kmer, pBWT);
            kmerDistribution.add(count);
        }
    }

    //
    kmerDistribution.print(75);

    double ratio = 2.0f;
    int chosenThreshold = kmerDistribution.findErrorBoundaryByRatio(ratio);
    double cumulativeLEQ = kmerDistribution.getCumulativeProportionLEQ(chosenThreshold);

    if(chosenThreshold == -1)
    {
        std::cerr << "[sga correct] Error k-mer threshold learning failed\n";
        std::cerr << "[sga correct] This can indicate the k-mer you choose is too high or your data has very low coverage\n";
        exit(EXIT_FAILURE);
    }

    std::cout << "Chosen kmer threshold: " << chosenThreshold << "\n";
    std::cout << "Proportion of kmer density right of threshold: " << 1.0f - cumulativeLEQ << "\n";
    if(cumulativeLEQ > 0.25f)
    {
        std::cerr << "[sga correct] Warning: Proportion of kmers greater than the chosen threshold is less than 0.75 (" << 1.0f - cumulativeLEQ  << ")\n";
        std::cerr << "[sga correct] This can indicate your chosen kmer size is too large or your data is too low coverage to reliably correct\n";
        std::cerr << "[sga correct] It is suggest to lower the kmer size and/or choose the threshold manually\n";
    }
    
    return chosenThreshold;
}

// 
// Handle command line arguments
//
void parseCorrectOptions(int argc, char** argv)
{
    std::string algo_str;
    bool bDiscardReads = false;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'm': arg >> opt::minOverlap; break;
            case 'p': arg >> opt::prefix; break;
            case 'o': arg >> opt::outFile; break;
            case 'e': arg >> opt::errorRate; break;
            case 't': arg >> opt::numThreads; break;
            case 'l': arg >> opt::seedLength; break;
            case 's': arg >> opt::seedStride; break;
            case 'r': arg >> opt::numOverlapRounds; break;
            case 'a': arg >> algo_str; break;
            case 'd': arg >> opt::sampleRate; break;
            case 'c': arg >> opt::conflictCutoff; break;
            case 'k': arg >> opt::kmerLength; break;
            case 'x': arg >> opt::kmerThreshold; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case 'b': arg >> opt::branchCutoff; break;
            case 'i': arg >> opt::numKmerRounds; break;
            case OPT_LEARN: opt::bLearnKmerParams = true; break;
            case OPT_DISCARD: bDiscardReads = true; break;
            case OPT_METRICS: arg >> opt::metricsFile; break;
            case OPT_HELP:
                std::cout << CORRECT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << CORRECT_VERSION_MESSAGE;
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

    if(opt::numOverlapRounds <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of overlap rounds: " << opt::numOverlapRounds << ", must be at least 1\n";
        die = true;
    }
    
    if(opt::numKmerRounds <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of kmer rounds: " << opt::numKmerRounds << ", must be at least 1\n";
        die = true;
    }

    if(opt::kmerLength <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid kmer length: " << opt::kmerLength << ", must be greater than zero\n";
        die = true;
    }

    if(opt::kmerThreshold <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid kmer threshold: " << opt::kmerThreshold << ", must be greater than zero\n";
        die = true;
    }

    // Determine the correction algorithm to use
    if(!algo_str.empty())
    {
        if(algo_str == "hybrid")
            opt::algorithm = ECA_HYBRID;
        else if(algo_str == "kmer")
            opt::algorithm = ECA_KMER;
        else if(algo_str == "overlap")
            opt::algorithm = ECA_OVERLAP;
        else
        {
            std::cerr << SUBPROGRAM << ": unrecognized -a,--algorithm parameter: " << algo_str << "\n";
            die = true;
        }
    }

    if (die) 
    {
        std::cout << "\n" << CORRECT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Validate parameters
    if(opt::errorRate <= 0)
        opt::errorRate = 0.0f;

    if(opt::errorRate > 1.0f)
    {
        std::cerr << "Invalid error-rate parameter: " << opt::errorRate << "\n";
        exit(EXIT_FAILURE);
    }

    if(opt::seedLength < 0)
        opt::seedLength = 0;

    if(opt::seedLength > 0 && opt::seedStride <= 0)
        opt::seedStride = opt::seedLength;

    // Parse the input filenames
    opt::readsFile = argv[optind++];

    if(opt::prefix.empty())
    {
        opt::prefix = stripFilename(opt::readsFile);
    }

    // Set the correction threshold
    if(opt::kmerThreshold <= 0)
    {
        std::cerr << "Invalid kmer support threshold: " << opt::kmerThreshold << "\n";
        exit(EXIT_FAILURE);
    }
    CorrectionThresholds::Instance().setBaseMinSupport(opt::kmerThreshold);

    std::string out_prefix = stripFilename(opt::readsFile);
    if(opt::outFile.empty())
    {
        opt::outFile = out_prefix + ".ec.fa";
    }

    if(bDiscardReads)
    {
        opt::discardFile = out_prefix + ".discard.fa";
    }
    else
    {
        opt::discardFile.clear();
    }
}
