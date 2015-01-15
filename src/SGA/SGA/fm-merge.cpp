//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// fm-merge - Merge reads/sequences using the FM-index
// without explicitly constructing the full graph
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "fm-merge.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "ASQG.h"
#include "gzstream.h"
#include "SequenceProcessFramework.h"
#include "OverlapProcess.h"
#include "ReadInfoTable.h"
#include "FMMergeProcess.h"

//
// Getopt
//
#define SUBPROGRAM "fm-merge"
static const char *FMMERGE_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *FMMERGE_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Merge unambiguously sequences from the READSFILE using the FM-index.\n"
"This program requires rmdup to be run before it.\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -p, --prefix=PREFIX              use PREFIX for the names of the index files (default: prefix of the input file)\n"
"      -t, --threads=NUM                use NUM worker threads (default: no threading)\n"
"      -m, --min-overlap=LEN            minimum overlap required between two reads to merge (default: 45)\n"
"      -o, --outfile=FILE               write the merged sequences to FILE (default: basename.merged.fa)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static std::string readsFile;
    static std::string outFile;
    static std::string prefix;
    static unsigned int minOverlap = DEFAULT_MIN_OVERLAP;
}

static const char* shortopts = "p:m:d:e:t:l:s:o:vix";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "prefix",      required_argument, NULL, 'p' },
    { "verbose",     no_argument,       NULL, 'v' },
    { "threads",     required_argument, NULL, 't' },
    { "min-overlap", required_argument, NULL, 'm' },
    { "outfile",     required_argument, NULL, 'o' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int FMMergeMain(int argc, char** argv)
{
    parseFMMergeOptions(argc, argv);

    BWT* pBWT = new BWT(opt::prefix + BWT_EXT);
    BWT* pRBWT = new BWT(opt::prefix + RBWT_EXT);
    OverlapAlgorithm* pOverlapper = new OverlapAlgorithm(pBWT, pRBWT,0.0f, 0,0,true); 
    pOverlapper->setExactModeOverlap(true);
    pOverlapper->setExactModeIrreducible(true);
    Timer* pTimer = new Timer(PROGRAM_IDENT);
    pBWT->printInfo();

    // Construct a bitvector indicating what reads have been used
    // All the processes read from this vector and only the post processor
    // writes to it.
    BitVector markedReads(pBWT->getNumStrings());

    std::ostream* pWriter = createWriter(opt::outFile);
    FMMergePostProcess postProcessor(pWriter, &markedReads);

    if(opt::numThreads <= 1)
    {
        printf("[%s] starting serial-mode read merging\n", PROGRAM_IDENT);
        FMMergeProcess processor(pOverlapper, opt::minOverlap, &markedReads);
        SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
                                                         FMMergeResult, 
                                                         FMMergeProcess, 
                                                         FMMergePostProcess>(opt::readsFile, &processor, &postProcessor);
    }
    else
    {
        printf("[%s] starting parallel-mode read merging computation with %d threads\n", PROGRAM_IDENT, opt::numThreads);
        
        std::vector<FMMergeProcess*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            FMMergeProcess* pProcessor = new FMMergeProcess(pOverlapper, opt::minOverlap, &markedReads);
            processorVector.push_back(pProcessor);
        }

        SequenceProcessFramework::processSequencesParallel<SequenceWorkItem,
                                                         FMMergeResult, 
                                                         FMMergeProcess, 
                                                         FMMergePostProcess>(opt::readsFile, processorVector, &postProcessor);
        
        for(size_t i = 0; i < processorVector.size(); ++i)
        {
            delete processorVector[i];
            processorVector[i] = NULL;
        }
    }

    // Check that every bit was set in the bit vector
    size_t numSet = 0;
    size_t numTotal = pBWT->getNumStrings();
    for(size_t i = 0; i < numTotal; ++i)
    {
        if(markedReads.test(i))
            ++numSet;
    }

    // Get the number of strings in the BWT, this is used to pre-allocated the read table
    delete pOverlapper;
    delete pBWT; 
    delete pRBWT;
    delete pWriter;

    // Cleanup
    delete pTimer;
    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

// 
// Handle command line arguments
//
void parseFMMergeOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'm': arg >> opt::minOverlap; break;
            case 'p': arg >> opt::prefix; break;
            case 'o': arg >> opt::outFile; break;
            case 't': arg >> opt::numThreads; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << FMMERGE_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << FMMERGE_VERSION_MESSAGE;
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
        std::cout << "\n" << FMMERGE_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::readsFile = argv[optind++];

    if(opt::prefix.empty())
    {
        opt::prefix = stripFilename(opt::readsFile);
    }

    if(opt::outFile.empty())
    {
        opt::outFile = opt::prefix + ".merged.fa";
    }
 
}

