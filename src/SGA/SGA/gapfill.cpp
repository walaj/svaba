//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// gapfill - Fill intrascaffold gaps
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "SequenceProcessFramework.h"
#include "SGACommon.h"
#include "GraphCompare.h"
#include "GapFillProcess.h"
#include "gapfill.h"

// Defines to clarify awful template function calls
#define PROCESS_GAPFILL_SERIAL SequenceProcessFramework::processSequencesSerial<SequenceWorkItem, GraphCompareResult, \
                                                                                GraphCompare, GraphCompareAggregateResults>

#define PROCESS_GAPFILL_PARALLEL SequenceProcessFramework::processSequencesParallel<SequenceWorkItem, GraphCompareResult, \
                                                                                    GraphCompare, GraphCompareAggregateResults>

   
//
// Getopt
//
#define SUBPROGRAM "gapfill"
static const char *GAPFILL_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2011 Wellcome Trust Sanger Institute\n";

static const char *GAPFILL_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] SCAFFOLDS.fa\n"
"Fill in scaffold gaps using walks through a de Bruijn graph\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -p, --prefix=NAME                load the FM-index with prefix NAME\n"
"      -s, --start-kmer=K               First kmer size used to attempt to resolve each gap (default: 91)\n"
"      -e, --end-kmer=K                 Last kmer size used to attempt to resolve each gap (default: 51)\n"
"      -x, --kmer-threshold=T           only use kmers seen at least T times\n"
"      -t, --threads=NUM                use NUM computation threads\n"
"      -d, --sample-rate=N              use occurrence array sample rate of N in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

//static const char* PROGRAM_IDENT = PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static int startKmer = 91;
    static int endKmer = 51;
    static int stride = 10;
    static int kmerThreshold = 3;
    static int sampleRate = 128;
    static int cacheLength = 10;

    static std::string scaffoldFile;
    static std::string prefix;
    static std::string outFile = "scaffolds.gapfill.fa";
}

static const char* shortopts = "o:s:e:t:x:p:s:d:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "threads",       required_argument, NULL, 't' },
    { "prefix",        required_argument, NULL, 'p' },
    { "outfile",       required_argument, NULL, 'o' },
    { "start-kmer",    required_argument, NULL, 's' },
    { "end-kmer",      required_argument, NULL, 'e' },
    { "kmer-threshold",required_argument, NULL, 'x' },
    { "sample-rate",   required_argument, NULL, 'd' },
    { "help",          no_argument,       NULL, OPT_HELP },
    { "version",       no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int gapfillMain(int argc, char** argv)
{
    parseGapFillOptions(argc, argv);

    // In the BWTs and create interval caches
    assert(!opt::prefix.empty());
    BWT* pBWT = new BWT(opt::prefix + BWT_EXT, opt::sampleRate);
    pBWT->printInfo();

    BWTIntervalCache* pBWTCache = new BWTIntervalCache(opt::cacheLength, pBWT);

    GapFillParameters parameters;
    parameters.pBWT = pBWT;
    parameters.pRevBWT = NULL;
    parameters.pBWTCache = pBWTCache;
    parameters.pRevBWTCache = NULL;
    parameters.startKmer = opt::startKmer;
    parameters.endKmer = opt::endKmer;
    parameters.stride = opt::stride;
    parameters.kmerThreshold = opt::kmerThreshold;
    parameters.verbose = opt::verbose;

    GapFillProcess processor(parameters);

    std::ostream* pWriter = createWriter(opt::outFile);

    SeqReader reader(opt::scaffoldFile, SRF_NO_VALIDATION | SRF_KEEP_CASE);
    SeqRecord record;
    while(reader.get(record))
    {
        GapFillResult result = processor.processScaffold(record.seq.toString());
        record.seq = result.scaffold;
        record.write(*pWriter);
    }

    // Cleanup
    delete pWriter;
    delete pBWT;
    delete pBWTCache;

    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

// 
// Handle command line arguments
//
void parseGapFillOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 's': arg >> opt::startKmer; break;
            case 'e': arg >> opt::endKmer; break;
            case 'p': arg >> opt::prefix; break;
            case 'o': arg >> opt::outFile; break;
            case 't': arg >> opt::numThreads; break;
            case 'd': arg >> opt::sampleRate; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << GAPFILL_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << GAPFILL_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    // Validate parameters
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

    if(opt::prefix.empty())
    {
        std::cerr << SUBPROGRAM ": error a --prefix for the FM-index must be supplied\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << GAPFILL_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    opt::scaffoldFile = argv[optind++];
}
