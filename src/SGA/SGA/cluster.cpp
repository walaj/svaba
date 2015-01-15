//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// cluster - find connected components in an 
// overlap graph
//
#include <iostream>
#include <fstream>
#include "cluster.h"
#include "Util.h"
#include "SGUtil.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"
#include "SGSearch.h"
#include "Timer.h"
#include "OverlapCommon.h"
#include "BitVector.h"
#include "ClusterProcess.h"
#include "ClusterReader.h"
#include <unistd.h>

// Local functions
void readLimitKmers(std::set<std::string>* limit_kmers);

// Defines to clarify awful template function calls
#define PROCESS_CLUSTER_SERIAL SequenceProcessFramework::processSequencesSerial<SequenceWorkItem, ClusterResult, \
                                                                                ClusterProcess, ClusterPostProcess>

#define PROCESS_CLUSTER_PARALLEL SequenceProcessFramework::processSequencesParallel<SequenceWorkItem, ClusterResult, \
                                                                                   ClusterProcess, ClusterPostProcess>


#define PROCESS_EXTEND_SERIAL SequenceProcessFramework::processWorkSerial<ClusterVector, ClusterResult, \
                                                                         ClusterReader, ClusterProcess, \
                                                                         ClusterPostProcess>

#define PROCESS_EXTEND_PARALLEL SequenceProcessFramework::processWorkParallelPthread<ClusterVector, ClusterResult, \
                                                                                     ClusterReader, ClusterProcess, \
                                                                                     ClusterPostProcess>
//
// Getopt
//
#define SUBPROGRAM "cluster"

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

static const char *CLUSTER_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *CLUSTER_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] READS\n"
"Extract clusters of reads belonging to the same connected components.\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -o, --out=FILE                   write the clusters to FILE (default: clusters.txt)\n"
"      -c, --min-cluster-size=N         only write clusters with at least N reads (default: 2)\n"
"      -x, --max-cluster-size=N         abort the search if the cluster size exceeds N\n"
"      -m, --min-overlap=N              require an overlap of at least N bases between reads (default: 45)\n"
"      -e, --error-rate                 the maximum error rate allowed to consider two sequences aligned (default: exact matches only)\n"
"      -t, --threads=NUM                use NUM worker threads to compute the overlaps (default: no threading)\n"
"      -i, --iterations=NUM             limit cluster extension to NUM iterations\n"
"          --extend=FILE                extend previously existing clusters in FILE\n"
"          --limit=FILE                 do not extend through the sequences in FILE\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static size_t minSize = 2;
    static size_t maxSize = (size_t)-1;
    static std::string outFile;
    static std::string readsFile;
    static std::string prefix;
    static std::string extendFile;
    static std::string limitFile;
    static int seedLength = 0;
    static int seedStride = 0;//seedLength;
    static int numThreads = 1;
    static size_t limitKmer = 21;
    static int maxIterations = -1;
    static double errorRate = 0.0f;
    static unsigned int minOverlap = DEFAULT_MIN_OVERLAP;
}

static const char* shortopts = "o:m:c:t:e:p:x:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_EXTEND, OPT_LIMIT };

static const struct option longopts[] = {
    { "verbose",          no_argument,       NULL, 'v' },
    { "out",              required_argument, NULL, 'o' },
    { "min-cluster-size", required_argument, NULL, 'c' },
    { "max-cluster-size", required_argument, NULL, 'x' },
    { "min-overlap",      required_argument, NULL, 'm' },
    { "error-rate",       required_argument, NULL, 'e' },
    { "threads",          required_argument, NULL, 't' },
    { "iterations",       required_argument, NULL, 'i' },
    { "prefix",           required_argument, NULL, 'p' },
    { "extend",           required_argument, NULL, OPT_EXTEND },
    { "limit",            required_argument, NULL, OPT_LIMIT },
    { "help",             no_argument,       NULL, OPT_HELP },
    { "version",          no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int clusterMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("sga cluster");
    parseClusterOptions(argc, argv);
    cluster();
    delete pTimer;

    return 0;
}

void cluster()
{
    BWT* pBWT = new BWT(opt::prefix + BWT_EXT);
    BWT* pRBWT = new BWT(opt::prefix + RBWT_EXT);
    OverlapAlgorithm* pOverlapper = new OverlapAlgorithm(pBWT, pRBWT,opt::errorRate, opt::seedLength, opt::seedStride, true);

    pOverlapper->setExactModeOverlap(opt::errorRate < 0.001f);
    pOverlapper->setExactModeIrreducible(opt::errorRate < 0.001f);

    BitVector markedReads(pBWT->getNumStrings());

    std::string preclustersFile = opt::outFile + ".preclusters";
    std::ostream* pPreWriter = createWriter(preclustersFile);
    ClusterPostProcess postProcessor(pPreWriter, opt::minSize, &markedReads);
    
    // Set the cluster parameters
    ClusterParameters parameters;
    parameters.pOverlapper = pOverlapper;
    parameters.minOverlap = opt::minOverlap;
    parameters.maxClusterSize = opt::maxSize;
    parameters.maxIterations = opt::maxIterations;
    parameters.pMarkedReads = &markedReads;

    // Read the limit kmer sequences, if provided
    std::set<std::string>* pLimitKmers = NULL;

    if(!opt::limitFile.empty())
    {
        // Read in the limit sequences
        pLimitKmers = new std::set<std::string>;
        readLimitKmers(pLimitKmers);
        parameters.pLimitKmers = pLimitKmers;
        parameters.limitK = opt::limitKmer;
    }
    else
    {
        parameters.pLimitKmers = NULL;
        parameters.limitK = 0;
    }

    // Make pre-clusters from the reads
    if(opt::numThreads <= 1)
    {
        printf("[%s] starting serial-mode read clustering\n", PROGRAM_IDENT);
        ClusterProcess processor(parameters);
        
        // If the extend file is empty, build new clusters
        if(opt::extendFile.empty())
        {
            PROCESS_CLUSTER_SERIAL(opt::readsFile, &processor, &postProcessor);
        }
        else
        {
            // Process a set of preexisting clusters
            ClusterReader clusterReader(opt::extendFile);
            PROCESS_EXTEND_SERIAL(clusterReader, &processor, &postProcessor);
        }
    }
    else
    {
        printf("[%s] starting parallel-mode read clustering computation with %d threads\n", PROGRAM_IDENT, opt::numThreads);
        
        std::vector<ClusterProcess*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            ClusterProcess* pProcessor = new ClusterProcess(parameters);
            processorVector.push_back(pProcessor);
        }
        
        if(opt::extendFile.empty())
        {
            PROCESS_CLUSTER_PARALLEL(opt::readsFile, processorVector, &postProcessor);
        }
        else
        {
            ClusterReader clusterReader(opt::extendFile);
            PROCESS_EXTEND_PARALLEL(clusterReader, processorVector, &postProcessor);
        }
        
        for(size_t i = 0; i < processorVector.size(); ++i)
        {
            delete processorVector[i];
            processorVector[i] = NULL;
        }
    }
    delete pPreWriter;
    delete pBWT;
    delete pRBWT;
    delete pOverlapper;

    // Deallocate limit kmers
    if(pLimitKmers != NULL)
        delete pLimitKmers;

    // Open the preclusters file and convert them to read names
    SuffixArray* pFwdSAI = new SuffixArray(opt::prefix + SAI_EXT);
    ReadInfoTable* pRIT = new ReadInfoTable(opt::readsFile, pFwdSAI->getNumStrings());

    size_t seedIdx = 0;
    std::istream* pPreReader = createReader(preclustersFile);
    std::ostream* pClusterWriter = createWriter(opt::outFile);
    std::string line;
    while(getline(*pPreReader,line))
    {
        std::stringstream parser(line);
        std::string clusterName;
        std::string readSequence;
        size_t clusterSize;
        int64_t lowIdx;
        int64_t highIdx;
        parser >> clusterName >> clusterSize >> readSequence >> lowIdx >> highIdx;

        if(lowIdx > highIdx)
        {
            // This is an extra read that is not present in the FM-index
            // Output a record with a fake read ID
            *pClusterWriter << clusterName << "\t" << clusterSize << "\tseed-" << seedIdx++ << "\t" << readSequence << "\n";
        }
        else
        {
            for(int64_t i = lowIdx; i <= highIdx; ++i)
            {
                const ReadInfo& targetInfo = pRIT->getReadInfo(pFwdSAI->get(i).getID());
                std::string readName = targetInfo.id;
                *pClusterWriter << clusterName << "\t" << clusterSize << "\t" << readName << "\t" << readSequence << "\n";
            }
        }
    }
    unlink(preclustersFile.c_str());

    delete pFwdSAI;
    delete pRIT;
    delete pPreReader;
    delete pClusterWriter;
}

//
void readLimitKmers(std::set<std::string>* limit_kmers)
{
    assert(!opt::limitFile.empty());
    assert(opt::limitKmer > 0);
    ClusterReader limit_reader(opt::limitFile);
    ClusterRecord record;
    while(limit_reader.readCluster(record))
    {
        // Insert kmers into the output set

        // Skip short sequences
        if(record.sequence.size() < opt::limitKmer)
            continue;

        size_t nk = record.sequence.size() - opt::limitKmer + 1;
        for(size_t i = 0; i < nk; ++i) 
        {
            std::string kmer = record.sequence.substr(i, opt::limitKmer);
            limit_kmers->insert(kmer);
        }
    }
}

// 
// Handle command line arguments
//
void parseClusterOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'o': arg >> opt::outFile; break;
            case 'c': arg >> opt::minSize; break;
            case 'x': arg >> opt::maxSize; break;
            case 'e': arg >> opt::errorRate; break;
            case 'm': arg >> opt::minOverlap; break;
            case 't': arg >> opt::numThreads; break;
            case 'i': arg >> opt::maxIterations; break;
            case 'p': arg >> opt::prefix; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_EXTEND: arg >> opt::extendFile; break;
            case OPT_LIMIT: arg >> opt::limitFile; break;
            case OPT_HELP:
                std::cout << CLUSTER_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << CLUSTER_VERSION_MESSAGE;
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

    if (die) 
    {
        std::cout << "\n" << CLUSTER_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filename
    opt::readsFile = argv[optind++];
    opt::prefix = stripFilename(opt::readsFile);

    if(opt::outFile.empty())
        opt::outFile = opt::prefix + ".clusters";
}
