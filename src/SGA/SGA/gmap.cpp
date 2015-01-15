//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// gmap - Map sequences to the vertices of a graph
//
#include <iostream>
#include <fstream>
#include <algorithm>
#include "Util.h"
#include "gmap.h"
#include "overlap.h"
#include "Timer.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "SequenceProcessFramework.h"
#include "RmdupProcess.h"
#include "BWTDiskConstruction.h"
#include <unistd.h>

struct GmapData
{
    std::string id;
    bool isRC;
};

typedef std::vector<GmapData> GmapVector;

//
size_t computeGmapHitsSerial(const std::string& prefix, const std::string& readsFile, 
                              const OverlapAlgorithm* pOverlapper, StringVector& filenameVec);

size_t computeGmapHitsParallel(int numThreads, const std::string& prefix, const std::string& readsFile, 
                               const OverlapAlgorithm* pOverlapper, StringVector& filenameVec);
//
// Getopt
//
#define SUBPROGRAM "gmap"
static const char *GMAP_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *GMAP_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... TARGETS QUERYS\n"
"Map the reads in QUERYS to the reads in TARGETS. TARGETS must be indexed\n"
"Used to place reads that were rmdup'd onto a graph.\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -o, --out=FILE                   write the output to FILE (default: READFILE.gmap)\n"
"      -p, --prefix=PREFIX              use PREFIX instead of the prefix of the reads filename for the input/output files\n"
"      -e, --error-rate                 the maximum error rate allowed to consider two sequences identical (default: exact matches required)\n"
"      -t, --threads=N                  use N threads (default: 1)\n"
"      -d, --sample-rate=N              sample the symbol counts every N symbols in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static std::string prefix;
    static std::string outFile;
    static std::string readsFile;
    static std::string targetsFile;
    static unsigned int numThreads;
    static double errorRate;
    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
}

static const char* shortopts = "p:o:e:t:d:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_VALIDATE };

static const struct option longopts[] = {
    { "verbose",        no_argument,       NULL, 'v' },
    { "prefix",         required_argument, NULL, 'p' },
    { "out",            required_argument, NULL, 'o' },
    { "error-rate",     required_argument, NULL, 'e' },
    { "threads",        required_argument, NULL, 't' },
    { "sample-rate",    required_argument, NULL, 'd' },
    { "help",           no_argument,       NULL, OPT_HELP },
    { "version",        no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int gmapMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("sga gmap");
    parseGmapOptions(argc, argv);
    gmap();
    delete pTimer;

    if(opt::numThreads > 1)
        pthread_exit(NULL);
    return 0;
}

void gmap()
{
    StringVector hitsFilenames;
    BWT* pBWT = new BWT(opt::prefix + BWT_EXT, opt::sampleRate);
    BWT* pRBWT = new BWT(opt::prefix + RBWT_EXT, opt::sampleRate);
    OverlapAlgorithm* pOverlapper = new OverlapAlgorithm(pBWT, pRBWT, 
                                                         opt::errorRate, 0, 
                                                         0, false);
    Timer* pTimer = new Timer(PROGRAM_IDENT);
    if(opt::numThreads <= 1)
    {
        printf("[%s] starting serial-mode overlap computation\n", PROGRAM_IDENT);
        computeGmapHitsSerial(opt::prefix, opt::readsFile, pOverlapper, hitsFilenames);
    }
    else
    {
        printf("[%s] starting parallel-mode overlap computation with %d threads\n", PROGRAM_IDENT, opt::numThreads);
        computeGmapHitsParallel(opt::numThreads, opt::prefix, opt::readsFile, pOverlapper, hitsFilenames);
    }

    delete pOverlapper;
    delete pBWT; 
    delete pRBWT;
    delete pTimer;

    parseGmapHits(hitsFilenames);
}

// Compute the hits for each read in the input file without threading
// Return the number of reads processed
size_t computeGmapHitsSerial(const std::string& prefix, const std::string& readsFile, 
                              const OverlapAlgorithm* pOverlapper, StringVector& filenameVec)
{
    std::string filename = prefix + GMAPHITS_EXT + GZIP_EXT;
    filenameVec.push_back(filename);

    RmdupProcess processor(filename, pOverlapper);
    RmdupPostProcess postProcessor;

    size_t numProcessed = 
           SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
                                                            OverlapResult, 
                                                            RmdupProcess, 
                                                            RmdupPostProcess>(readsFile, &processor, &postProcessor);
    return numProcessed;
}

// Compute the hits for each read in the SeqReader file with threading
// The way this works is we create a vector of numThreads OverlapProcess pointers and 
// pass this to the SequenceProcessFragmework which wraps the processes
// in threads and distributes the reads to each thread.
// The number of reads processsed is returned
size_t computeGmapHitsParallel(int numThreads, const std::string& prefix, const std::string& readsFile, 
                                const OverlapAlgorithm* pOverlapper, StringVector& filenameVec)
{
    std::string filename = prefix + RMDUPHITS_EXT + GZIP_EXT;

    std::vector<RmdupProcess*> processorVector;
    for(int i = 0; i < numThreads; ++i)
    {
        std::stringstream ss;
        ss << prefix << "-thread" << i << GMAPHITS_EXT << GZIP_EXT;
        std::string outfile = ss.str();
        filenameVec.push_back(outfile);
        RmdupProcess* pProcessor = new RmdupProcess(outfile, pOverlapper);
        processorVector.push_back(pProcessor);
    }

    // The post processing is performed serially so only one post processor is created
    RmdupPostProcess postProcessor;
    
    size_t numProcessed = 
           SequenceProcessFramework::processSequencesParallel<SequenceWorkItem,
                                                              OverlapResult, 
                                                              RmdupProcess, 
                                                              RmdupPostProcess>(readsFile, processorVector, &postProcessor);
    for(int i = 0; i < numThreads; ++i)
        delete processorVector[i];
    return numProcessed;
}

void parseGmapHits(const StringVector& hitsFilenames)
{
    // Load the suffix array index and the reverse suffix array index
    // Note these are not the full suffix arrays
    SuffixArray* pFwdSAI = new SuffixArray(opt::prefix + SAI_EXT);
    SuffixArray* pRevSAI = new SuffixArray(opt::prefix + RSAI_EXT);

    // Load the read table and output the initial vertex set, consisting of all the reads
    ReadInfoTable* pRIT = new ReadInfoTable(opt::targetsFile, pFwdSAI->getNumStrings());

    std::ostream* pWriter = createWriter(opt::outFile);
    int numRead = 0;
    int numMapped = 0;
    size_t num_files = hitsFilenames.size();

    for(size_t i = 0; i < num_files; ++i)
    {
        std::cout << "Opening " << hitsFilenames[i] << "\n";
        std::istream* pReader = createReader(hitsFilenames[i]);
        std::string line;

        while(getline(*pReader, line))
        {
           ++numRead;
            // Parse the hit
            std::string id;
            std::string sequence;

            std::stringstream parser(line);
            parser >> id;
            parser >> sequence;

            // Parse the overlap blocks into the list of IDs that this read matches
            std::string hitsStr;
            getline(parser, hitsStr);

            GmapVector matchedReads;
            size_t readIdx;
            size_t numBlocks;
            int isSubstring;
            std::istringstream convertor(hitsStr);
            convertor >> readIdx >> isSubstring >> numBlocks;

            // Convert the blocks to IDs
            for(size_t i = 0; i < numBlocks; ++i)
            {
                // Read the block
                OverlapBlock record;
                convertor >> record;

                // Iterate through the range and write the overlaps
                for(int64_t j = record.ranges.interval[0].lower; j <= record.ranges.interval[0].upper; ++j)
                {
                    const SuffixArray* pCurrSAI = (record.flags.isTargetRev()) ? pRevSAI : pFwdSAI;
                    int64_t saIdx = j;

                    // The index of the second read is given as the position in the SuffixArray index
                    const ReadInfo& targetInfo = pRIT->getReadInfo(pCurrSAI->get(saIdx).getID());

                    // Avoid self-matches to the opposite strand for palindromes
                    GmapData data = {targetInfo.id, record.flags.isReverseComplement()};
                    if(data.id == id && data.isRC == true)
                        continue;

                    matchedReads.push_back(data);
                }
            }

            // Build the output record
            GmapRecord record;
            record.readID = id;
            record.readSeq = sequence;
            record.mappedID = "-";
            record.isRC = false;

            if(matchedReads.size() == 1)
            {
                record.mappedID = matchedReads.front().id;
                if(record.mappedID == id)
                    record.isRC = false;
                else
                    record.isRC = matchedReads.front().isRC;

                numMapped += 1;
            }
            else if(matchedReads.size() > 1)
            {
                record.mappedID = "MM";
            }

            *pWriter << record << "\n";
        }
        delete pReader;

        // Delete the hits file
        unlink(hitsFilenames[i].c_str());
    }
    
    std::cout << "Read: " << numRead << "\n";
    std::cout << "Mapped: " << numMapped << "\n";

    // Delete allocated data
    delete pFwdSAI;
    delete pRevSAI;
    delete pRIT;
    delete pWriter;
}

// 
// Handle command line arguments
//
void parseGmapOptions(int argc, char** argv)
{
    // Set defaults
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'p': arg >> opt::prefix; break;
            case 'o': arg >> opt::outFile; break;
            case 'e': arg >> opt::errorRate; break;
            case 'd': arg >> opt::sampleRate; break;
            case 't': arg >> opt::numThreads; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << GMAP_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << GMAP_VERSION_MESSAGE;
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

    if (die) 
    {
        std::cerr << "Try `" << SUBPROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::targetsFile = argv[optind++];
    opt::readsFile = argv[optind++];

    if(opt::prefix.empty())
    {
        opt::prefix = stripFilename(opt::targetsFile);
    }

    if(opt::outFile.empty())
    {
        opt::outFile = stripFilename(opt::readsFile) + ".gmap.gz";
    }
}

