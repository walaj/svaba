//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// stats - Print statistics about the data set
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "stats.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "ASQG.h"
#include "gzstream.h"
#include "SequenceProcessFramework.h"
#include "StatsProcess.h"
#include "BWTDiskConstruction.h"

// Functions

//
// Getopt
//
#define SUBPROGRAM "stats"
static const char *STATS_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *STATS_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Print statistics about the read set. Currently this only prints a histogram\n"
"of the k-mer counts\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -p, --prefix=PREFIX              use PREFIX for the names of the index files (default: prefix of the input file)\n"
"      -t, --threads=NUM                use NUM threads to compute the overlaps (default: 1)\n"
"      -d, --sample-rate=N              use occurrence array sample rate of N in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
"      -k, --kmer-size=N                The length of the kmer to use. (default: 27)\n"
"      -n, --num-reads=N                Only use N reads to compute the statistics\n"
"      -b, --branch-cutoff=N            stop the overlap search at N branches. This lowers the compute time but will bias the statistics\n"
"                                       away from repetitive reads\n"
"      --run-lengths                    Print the run length distribution of the BWT\n"
"      --kmer-distribution              Print the distribution of kmer counts\n"
"      --no-overlap                     Suppress the overlap-based error statistics (faster if you only want the k-mer distribution)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static std::string prefix;
    static std::string readsFile;
    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
    static int kmerLength = 27;
    static int minOverlap = 45;
    static int branchCutoff = -1;
    static size_t numReads = -1;
    static bool bPrintRunLengths = false;
    static bool bPrintKmerDist = false;
    static bool bNoOverlap = false;
}

static const char* shortopts = "p:d:t:o:k:n:b:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_RUNLENGTHS, OPT_KMERDIST, OPT_NOOVERLAP};

static const struct option longopts[] = {
    { "verbose",            no_argument,       NULL, 'v' },
    { "threads",            required_argument, NULL, 't' },
    { "prefix",             required_argument, NULL, 'p' },
    { "sample-rate",        required_argument, NULL, 'd' },
    { "kmer-size",          required_argument, NULL, 'k' },
    { "num-reads",          required_argument, NULL, 'n' },
    { "branch-cutoff",      required_argument, NULL, 'b' },
    { "kmer-distribution",  no_argument,       NULL, OPT_KMERDIST },
    { "no-overlap",         no_argument,       NULL, OPT_NOOVERLAP },
    { "run-lengths",        no_argument,       NULL, OPT_RUNLENGTHS },
    { "help",               no_argument,       NULL, OPT_HELP },
    { "version",            no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int statsMain(int argc, char** argv)
{
    parseStatsOptions(argc, argv);
    Timer* pTimer = new Timer(PROGRAM_IDENT);

    BWT* pBWT = new BWT(opt::prefix + BWT_EXT, opt::sampleRate);
    BWT* pRBWT = NULL;

    // Do not need the reverse BWT if the overlap statistics are disabled
    if (!opt::bNoOverlap)
        pRBWT = new BWT(opt::prefix + RBWT_EXT, opt::sampleRate);

    if(opt::bPrintRunLengths)
    {
        pBWT->printInfo();
        pBWT->printRunLengths();
    }

    SeqReader reader(opt::readsFile);
    
    StatsPostProcess postProcessor(opt::bPrintKmerDist);
    if(opt::numThreads <= 1)
    {
        // Serial mode
        StatsProcess processor(pBWT, pRBWT, opt::kmerLength, opt::minOverlap, opt::branchCutoff, opt::bNoOverlap);

        SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
                                                         StatsResult, 
                                                         StatsProcess, 
                                                         StatsPostProcess>(reader, &processor, &postProcessor, opt::numReads);
    }
    else
    {
        // Parallel mode
        std::vector<StatsProcess*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            StatsProcess* pProcessor = new StatsProcess(pBWT, pRBWT, opt::kmerLength, opt::minOverlap, opt::branchCutoff, opt::bNoOverlap);
            processorVector.push_back(pProcessor);
        }
        
        SequenceProcessFramework::processSequencesParallel<SequenceWorkItem,
                                                           StatsResult, 
                                                           StatsProcess, 
                                                           StatsPostProcess>(reader, processorVector, &postProcessor, opt::numReads);

        for(int i = 0; i < opt::numThreads; ++i)
        {
            delete processorVector[i];
        }
    }

    delete pBWT;
    if (pRBWT)
        delete pRBWT;
    delete pTimer;

    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

// 
// Handle command line arguments
//
void parseStatsOptions(int argc, char** argv)
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
            case 'd': arg >> opt::sampleRate; break;
            case 'k': arg >> opt::kmerLength; break;
            case 'n': arg >> opt::numReads; break;
            case 'b': arg >> opt::branchCutoff; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_KMERDIST: opt::bPrintKmerDist = true; break;
            case OPT_RUNLENGTHS: opt::bPrintRunLengths = true; break;
            case OPT_NOOVERLAP: opt::bNoOverlap = true; break;
            case OPT_HELP:
                std::cout << STATS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << STATS_VERSION_MESSAGE;
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

    if(opt::kmerLength <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid kmer length: " << opt::kmerLength << ", must be greater than zero\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << STATS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::readsFile = argv[optind++];

    if(opt::prefix.empty())
    {
        opt::prefix = stripFilename(opt::readsFile);
    }
}
