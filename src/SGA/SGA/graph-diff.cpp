//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// graph-diff - Find strings that are only present
// in one of two input graphs
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "SuffixArray.h"
#include "SampledSuffixArray.h"
#include "BWT.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "SequenceProcessFramework.h"
#include "SGACommon.h"
#include "GraphCompare.h"
#include "VCFTester.h"
#include "DindelRealignWindow.h"
#include "QualityTable.h"
#include "BloomFilter.h"
#include "Verbosity.h"
#include "graph-diff.h"

// Functions
void runVCFTester(GraphCompareParameters& parameters);
void runGraphDiff(GraphCompareParameters& parameters);
void runDebug(GraphCompareParameters& parameters);
void runInteractive(GraphCompareParameters& parameters);
void preloadBloomFilter(const ReadTable* pReadTable, size_t k, BloomFilter* pBloomFilter);

// Defines to clarify awful template function calls
#define PROCESS_GDIFF_SERIAL SequenceProcessFramework::processSequencesSerial<SequenceWorkItem, GraphCompareResult, \
                                                                              GraphCompare, GraphCompareAggregateResults>

#define PROCESS_GDIFF_PARALLEL SequenceProcessFramework::processSequencesParallelOpenMP<SequenceWorkItem, GraphCompareResult, \
                                                                                        GraphCompare, GraphCompareAggregateResults>

   
//
// Getopt
//
#define SUBPROGRAM "graph-diff"
static const char *GRAPH_DIFF_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2011 Wellcome Trust Sanger Institute\n";

static const char *GRAPH_DIFF_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] --base BASE.fa --variant VARIANT.fa --ref REFERENCE.fa\n"
"Find and report strings only present in the graph of VARIANT when compared to BASE\n"
"\n"
"General options:\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -p, --prefix=NAME                prefix the output files with NAME\n"
"      -t, --threads=NUM                use NUM computation threads\n"
"          --genome-size=N              (optional) set the size of the genome to be N bases\n"
"                                       this is used to determine the number of bits to use in the bloom filter\n"
"                                       if unset, it will be calculated from the reference genome FASTA file\n"
"          --precache-reference=STR     precache the named chromosome of the reference genome\n"
"                                       If STR is \"all\" the entire reference will be cached\n"
//"          --test=VCF                   test the variants in the provided VCF file\n"
"\n"
"Index options:\n"
"      -r, --variant=FILE               call variants present in the read set in FILE\n"
"      -b, --base=FILE                  use the read set in FILE as the base line for comparison\n"
"                                       if this option is not given, reference-based calls will be made\n"
"          --reference=FILE             use the reference sequence in FILE\n"
"\n"
"Algorithm options:\n"
"      -k, --kmer=K                     use K-mers to discover variants\n"
"      -x, --min-discovery-count=T      require a variant k-mer to be seen at least T times\n"
"      -a, --algorithm=STR              select the assembly algorithm to use from: debruijn, string\n"
"      -m, --min-overlap=N              require at least N bp overlap when assembling using a string graph\n" 
"          --min-dbg-count=T            only use k-mers seen T times when assembling using a de Bruijn graph\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    // These parameters control run time and memory usage
    static unsigned int verbose = 0;
    static int numThreads = 1;
    static int cacheLength = 10;
    static int sampleRate = 128;
    static int bloomGenomeSize = -1;
    static std::string precacheReference;

    // Calling parameters
    static int kmer = 55;
    static int minDiscoveryCount = 2;
    static int minOverlap = 61;
    static int minDBGCount = 2;

    // Calling modes
    static GraphCompareAlgorithm algorithm = GCA_DEBRUIJN_GRAPH;
    static bool lowCoverage = false;
    static bool referenceMode = false;
    static bool useQualityScores = false;
    static bool interactiveMode = false;

    // I/O
    static std::string outPrefix = "sgavariants";
    static std::string indexPrefix;
    static std::string debugFile;
    static std::string referenceFile;
    static std::string baseFile;
    static std::string variantFile;
    static std::string inputVCFFile;
}

static const char* shortopts = "b:r:o:k:d:t:x:y:p:m:a:v";

enum { OPT_HELP = 1, 
       OPT_VERSION, 
       OPT_REFERENCE, 
       OPT_TESTVCF, 
       OPT_DEBUG, 
       OPT_MIN_DBG_COUNT, 
       OPT_INDEX, 
       OPT_LOWCOVERAGE, 
       OPT_QUALSCORES,
       OPT_BLOOM_GENOME,
       OPT_PRECACHE_REFERENCE,
       OPT_INTERACTIVE };

static const struct option longopts[] = {
    { "verbose",              no_argument,       NULL, 'v' },
    { "threads",              required_argument, NULL, 't' },
    { "base",                 required_argument, NULL, 'b' },
    { "variant",              required_argument, NULL, 'r' },
    { "kmer",                 required_argument, NULL, 'k' },
    { "min-discovery-count",  required_argument, NULL, 'x' },
    { "sample-rate",          required_argument, NULL, 'd' },
    { "prefix",               required_argument, NULL, 'p' },
    { "min-overlap",          required_argument, NULL, 'm' },
    { "algorithm",            required_argument, NULL, 'a' },
    { "low-coverage",         no_argument,       NULL, OPT_LOWCOVERAGE },
    { "use-quality-scores",   no_argument,       NULL, OPT_QUALSCORES},
    { "interactive",          no_argument,       NULL, OPT_INTERACTIVE},
    { "index",                required_argument, NULL, OPT_INDEX },
    { "min-dbg-count",        required_argument, NULL, OPT_MIN_DBG_COUNT },
    { "debug",                required_argument, NULL, OPT_DEBUG },
    { "reference",            required_argument, NULL, OPT_REFERENCE },
    { "genome-size",          required_argument, NULL, OPT_BLOOM_GENOME },
    { "precache-reference",   required_argument, NULL, OPT_PRECACHE_REFERENCE },
    { "test"      ,           required_argument, NULL, OPT_TESTVCF },
    { "help",                 no_argument,       NULL, OPT_HELP },
    { "version",              no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int graphDiffMain(int argc, char** argv)
{
    parseGraphDiffOptions(argc, argv);

    // Set the verbosity level for the entire package
    Verbosity::Instance().setPrintLevel(opt::verbose);

    if(opt::lowCoverage)
        std::cout << "Initializing population calling\n";
    else if(opt::referenceMode)
        std::cout << "Initializing single genome calling\n";
    else
        std::cout << "Initializing contrast calling (" << opt::variantFile << " vs " << opt::baseFile << ")\n";

    // Create indices for the variant reads
    std::string variantPrefix = stripGzippedExtension(opt::variantFile);

    // Use debug index prefix if specified
    if(!opt::indexPrefix.empty())
        variantPrefix = opt::indexPrefix;

    //
    // Initialize indices
    //

    // Variant reads
    std::cout << "Loading variant read index... " << std::flush;
    BWTIndexSet variantIndex;
    variantIndex.pBWT = new BWT(variantPrefix + BWT_EXT, opt::sampleRate);
    variantIndex.pSSA = new SampledSuffixArray(variantPrefix + SAI_EXT, SSA_FT_SAI);
    variantIndex.pCache = new BWTIntervalCache(opt::cacheLength, variantIndex.pBWT);
    if(opt::lowCoverage)
        variantIndex.pPopIdx = new PopulationIndex(variantPrefix + POPIDX_EXT);

    // Generate a quality table for the variant sequences. If --use-quality
    // is not set, this will just return default quality scores for all reads
    QualityTable* variantQuals = new QualityTable;
    if(opt::useQualityScores)
        variantQuals->loadQualities(opt::variantFile);
    variantIndex.pQualityTable = variantQuals;
    std::cout << "done" << std::endl; 

    // Reference genome
    std::cout << "Loading reference index... " << std::flush;
    BWTIndexSet referenceIndex;
    std::string refPrefix = stripExtension(opt::referenceFile);

    referenceIndex.pBWT = new BWT(refPrefix + BWT_EXT, opt::sampleRate);
    referenceIndex.pSSA = new SampledSuffixArray(refPrefix + SSA_EXT);
    
    // Read in the reference sequences
    ReadTable refTable(opt::referenceFile, SRF_NO_VALIDATION);
    refTable.indexReadsByID();
    std::cout << "done" << std::endl;

    // Base reads, only if not in reference mode
    BWTIndexSet baseIndex;

    if(!opt::referenceMode)
    {
        std::cout << "Loading base read index index... " << std::flush;
        std::string basePrefix = stripGzippedExtension(opt::baseFile);
        baseIndex.pBWT = new BWT(basePrefix + BWT_EXT, opt::sampleRate);
        baseIndex.pSSA = new SampledSuffixArray(basePrefix + SAI_EXT, SSA_FT_SAI);
        baseIndex.pCache = new BWTIntervalCache(opt::cacheLength, baseIndex.pBWT);
        baseIndex.pQualityTable = new QualityTable();

        QualityTable* baseQuals = new QualityTable;
        if(opt::useQualityScores)
            baseQuals->loadQualities(opt::baseFile);
        baseIndex.pQualityTable = baseQuals;
        std::cout << "done" << std::endl;
    }
    else
    {
        baseIndex.pBWT = referenceIndex.pBWT;
        baseIndex.pSSA = referenceIndex.pSSA;
        baseIndex.pCache = NULL;
    }

    // Set the parameters shared between all threads
    GraphCompareParameters sharedParameters;

    sharedParameters.baseIndex = baseIndex;
    sharedParameters.variantIndex = variantIndex;
    sharedParameters.referenceIndex = referenceIndex;
    sharedParameters.pRefTable = &refTable;
    sharedParameters.bReferenceMode = opt::referenceMode;
    sharedParameters.algorithm = opt::algorithm;
    sharedParameters.kmer = opt::kmer;
    sharedParameters.minDiscoveryCount = opt::minDiscoveryCount;
    sharedParameters.minDBGCount = opt::minDBGCount;
    sharedParameters.minOverlap = opt::minOverlap;
    sharedParameters.maxHaplotypes = 5;
    sharedParameters.maxReads = 10000;
    sharedParameters.maxExtractionIntervalSize = 500;
    sharedParameters.maxDiscoveryCount = 500;

    // Set population variant calling parameters
    if(opt::lowCoverage)
    {
        sharedParameters.maxDiscoveryCount = 10000;
        sharedParameters.maxExtractionIntervalSize = 5000;
        sharedParameters.maxReads = 100000;
        sharedParameters.dindelRealignParameters.multiSample = 1;
    }

    if(!opt::debugFile.empty())
    {
        runDebug(sharedParameters);
    }
    else if(opt::interactiveMode)
    {
        runInteractive(sharedParameters);
    }
    else
    {
        // If a VCF file was provided, just test the variants in the file without any discovery
        if(!opt::inputVCFFile.empty())
            runVCFTester(sharedParameters);
        else
            runGraphDiff(sharedParameters);
    }

    // Cleanup
    if(!opt::referenceMode)
    {
        delete baseIndex.pBWT;
        delete baseIndex.pSSA;
        delete baseIndex.pCache;
        delete baseIndex.pQualityTable;
    }

    // Cleanup indices
    delete variantIndex.pBWT;
    delete variantIndex.pSSA;
    delete variantIndex.pCache;
    delete variantIndex.pQualityTable;
    if(opt::lowCoverage)
        delete variantIndex.pPopIdx;

    delete referenceIndex.pBWT;
    delete referenceIndex.pSSA;

    return 0;
}

void runGraphDiff(GraphCompareParameters& parameters)
{
    // Initialize a bloom filter
    size_t expected_bits = 0;
    if(opt::bloomGenomeSize == -1)
    {
        // Calculate the expected bits using the number of bases in the reference
        for(size_t i = 0; i < parameters.pRefTable->getCount(); ++i)
            expected_bits += parameters.pRefTable->getReadLength(i);
    }
    else 
    {
        expected_bits = opt::bloomGenomeSize;
    }

    size_t occupancy_factor = 20;
    size_t bloom_size = occupancy_factor * expected_bits;

    BloomFilter* pBloomFilter = new BloomFilter(bloom_size, 5);
    parameters.pBloomFilter = pBloomFilter;
    preloadBloomFilter(parameters.pRefTable, parameters.kmer, pBloomFilter);

    Timer gdbenchmark("benchmark");

    // This call can throw via dindel
    GraphCompareAggregateResults* pSharedResults;
    try 
    {
        StringVector samples;

        // If in multi-sample mode, write sample names in the VCF header
        if(parameters.variantIndex.pPopIdx != NULL)
            samples = parameters.variantIndex.pPopIdx->getSamples();
        pSharedResults = new GraphCompareAggregateResults(opt::outPrefix, samples, *parameters.pRefTable, stripDirectories(opt::referenceFile));
    }
    catch(std::string e)
    {
        std::cout << "Exception: " << e << "\n";
        exit(EXIT_FAILURE);
    }


    if(opt::numThreads <= 1)
    {
        GraphCompare graphCompare(parameters); 
        PROCESS_GDIFF_SERIAL(opt::variantFile, &graphCompare, pSharedResults);
        graphCompare.updateSharedStats(pSharedResults);
    }
    else
    {
        std::vector<GraphCompare*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            GraphCompare* pProcessor = new GraphCompare(parameters);
            processorVector.push_back(pProcessor);
        }
        
        PROCESS_GDIFF_PARALLEL(opt::variantFile, processorVector, pSharedResults);
        
        for(size_t i = 0; i < processorVector.size(); ++i)
        {
            // Update the shared stats
            processorVector[i]->updateSharedStats(pSharedResults);

            delete processorVector[i];
            processorVector[i] = NULL;
        }
    }

    pSharedResults->printStats();

    delete pBloomFilter;
    delete pSharedResults;
}

//
void preloadBloomFilter(const ReadTable* pReadTable, size_t k, BloomFilter* pBloomFilter)
{
    std::cout << "Initializing bloom filter... " << std::flush;
    for(size_t i = 0; i < pReadTable->getCount(); ++i)
    {
        const SeqItem& si = pReadTable->getRead(i);
        if(si.id == opt::precacheReference || opt::precacheReference == "all")
        {
            const DNAString& seq = si.seq;
            for(size_t j = 0; j < seq.length() - k + 1; ++j)
            {       
                std::string kmer = seq.substr(j, k);
                std::string rc_kmer = reverseComplement(kmer);
                std::string& key_kmer = kmer < rc_kmer ? kmer : rc_kmer;
                pBloomFilter->add(key_kmer.c_str(), k);
            }
        }
    }
    std::cout << "done" << std::endl;
}

//
// DEBUG functions
//
void runVCFTester(GraphCompareParameters& parameters)
{
    std::cout << "Testing variants in " << opt::inputVCFFile << "\n";
    std::string line;

    VCFTester tester(parameters);
    try
    {
        VCFFile input(opt::inputVCFFile, "r");
        VCFFile::VCFEntry record;

        while(1)
        {
            record = input.getNextEntry();
            if(record.isEmpty())
                break;
            else
                tester.process(record);
        }
    }
    catch(std::string e)
    {
        std::cerr << "Exception: " << e << "\n";
        exit(EXIT_FAILURE);
    }
}


// Run in debug mode
void runDebug(GraphCompareParameters& parameters)
{
    GraphCompare graphCompare(parameters); 
    graphCompare.testKmersFromFile(opt::debugFile);
}

// Run interactively, parse a sequence from cin, process it then wait
void runInteractive(GraphCompareParameters& parameters)
{
    // Initialize a bloom filter
    size_t expected_bits = 0;
    if(opt::bloomGenomeSize == -1)
    {
        // Calculate the expected bits using the number of bases in the reference
        for(size_t i = 0; i < parameters.pRefTable->getCount(); ++i)
            expected_bits += parameters.pRefTable->getReadLength(i);
    }
    else 
    {
        expected_bits = opt::bloomGenomeSize;
    }

    size_t occupancy_factor = 20;
    size_t bloom_size = occupancy_factor * expected_bits;

    BloomFilter* pBloomFilter = new BloomFilter(bloom_size, 5);
    parameters.pBloomFilter = pBloomFilter;
    preloadBloomFilter(parameters.pRefTable, parameters.kmer, pBloomFilter);
    
    GraphCompare graphCompare(parameters); 
    while(1) 
    {
        std::string sequence;
        std::cout << "Input a sequence: ";
        std::cin >> sequence;

        if(sequence == "done")
            break;
        SequenceWorkItem item;
        item.read.seq = sequence;
        item.read.id = "input";

        std::cout << "processing: " << sequence << "\n";

        graphCompare.process(item);
    }
}

// 
// Handle command line arguments
//
void parseGraphDiffOptions(int argc, char** argv)
{
    std::string algorithmString;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'k': arg >> opt::kmer; break;
            case 'x': arg >> opt::minDiscoveryCount; break;
            case 'b': arg >> opt::baseFile; break;
            case 'r': arg >> opt::variantFile; break;
            case 't': arg >> opt::numThreads; break;
            case 'd': arg >> opt::sampleRate; break;
            case 'p': arg >> opt::outPrefix; break;
            case 'm': arg >> opt::minOverlap; break;
            case 'a': arg >> algorithmString; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_REFERENCE: arg >> opt::referenceFile; break;
            case OPT_LOWCOVERAGE: opt::lowCoverage = true; break;
            case OPT_MIN_DBG_COUNT: arg >> opt::minDBGCount; break;
            case OPT_BLOOM_GENOME: arg >> opt::bloomGenomeSize; break;
            case OPT_PRECACHE_REFERENCE: arg >> opt::precacheReference; break;
            case OPT_DEBUG: arg >> opt::debugFile; break;
            case OPT_TESTVCF: arg >> opt::inputVCFFile; break;
            case OPT_INDEX: arg >> opt::indexPrefix; break;
            case OPT_QUALSCORES:  opt::useQualityScores = true; break;
            case OPT_INTERACTIVE:  opt::interactiveMode = true; break;
            case OPT_HELP:
                std::cout << GRAPH_DIFF_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << GRAPH_DIFF_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    // Validate parameters
    if (argc - optind > 1) 
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if(opt::numThreads <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
        die = true;
    }

    if(opt::variantFile.empty())
    {
        std::cerr << SUBPROGRAM ": error a --base and --variant file must be provided\n";
        die = true;
    }

    if(opt::baseFile.empty())
        opt::referenceMode = true;

    if(algorithmString == "paired") {
        opt::algorithm = GCA_PAIRED_DEBRUIJN_GRAPH;
    } else if(algorithmString == "debruijn") {
        opt::algorithm = GCA_DEBRUIJN_GRAPH;
    } else if(algorithmString == "string") {
        opt::algorithm = GCA_STRING_GRAPH;
    } else {
        std::cerr << "Error: unrecognized algorithm string: " << algorithmString << "\n";
        die = true;
    }

    if(opt::referenceFile.empty())
    {
        std::cerr << SUBPROGRAM ": error, a --reference file must be provided\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << GRAPH_DIFF_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

}
