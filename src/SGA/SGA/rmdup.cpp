//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// rmdup - remove duplicated reads from the data set
//
#include <iostream>
#include <fstream>
#include "Util.h"
#include "rmdup.h"
#include "overlap.h"
#include "Timer.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "SequenceProcessFramework.h"
#include "RmdupProcess.h"
#include "BWTDiskConstruction.h"
#include <unistd.h>

// functions
size_t computeRmdupHitsSerial(const std::string& prefix, const std::string& readsFile, 
                              const OverlapAlgorithm* pOverlapper, StringVector& filenameVec);
 
size_t computeRmdupHitsParallel(int numThreads, const std::string& prefix, const std::string& readsFile, 
                                const OverlapAlgorithm* pOverlapper, StringVector& filenameVec);

//
// Getopt
//
#define SUBPROGRAM "rmdup"
static const char *RMDUP_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *RMDUP_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READFILE\n"
"Remove duplicate reads from the data set.\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -o, --out=FILE                   write the output to FILE (default: READFILE.rmdup.fa)\n"
"      -p, --prefix=PREFIX              use PREFIX instead of the prefix of the reads filename for the input/output files\n"
"      -e, --error-rate                 the maximum error rate allowed to consider two sequences identical (default: exact matches required)\n"
"      -t, --threads=N                  use N threads (default: 1)\n"
"      -d, --sample-rate=N              sample the symbol counts every N symbols in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 256)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static std::string prefix;
    static std::string outFile;
    static std::string readsFile;
    static unsigned int numThreads;
    static double errorRate;
    static bool bReindex = true;
    static int sampleRate = 256;
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
int rmdupMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("sga rmdup");
    parseRmdupOptions(argc, argv);
    rmdup();
    delete pTimer;

    if(opt::numThreads > 1)
        pthread_exit(NULL);
    return 0;
}

void rmdup()
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
        computeRmdupHitsSerial(opt::prefix, opt::readsFile, pOverlapper, hitsFilenames);
    }
    else
    {
        printf("[%s] starting parallel-mode overlap computation with %d threads\n", PROGRAM_IDENT, opt::numThreads);
        computeRmdupHitsParallel(opt::numThreads, opt::prefix, opt::readsFile, pOverlapper, hitsFilenames);
    }

    delete pOverlapper;
    delete pBWT; 
    delete pRBWT;
    delete pTimer;
    
    std::string out_prefix = stripFilename(opt::outFile);
    std::string dupsFile = parseDupHits(hitsFilenames, out_prefix);

    // Rebuild the indices without the duplicated sequences
    if(opt::bReindex)
    {
        std::cout << "Rebuilding indices without duplicated reads\n";
        removeReadsFromIndices(opt::prefix, dupsFile, out_prefix, BWT_EXT, SAI_EXT, false, opt::numThreads);
        removeReadsFromIndices(opt::prefix, dupsFile, out_prefix, RBWT_EXT, RSAI_EXT, true, opt::numThreads);
    }
}

// Compute the hits for each read in the input file without threading
// Return the number of reads processed
size_t computeRmdupHitsSerial(const std::string& prefix, const std::string& readsFile, 
                              const OverlapAlgorithm* pOverlapper, StringVector& filenameVec)
{
    std::string filename = prefix + RMDUPHITS_EXT + GZIP_EXT;
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
size_t computeRmdupHitsParallel(int numThreads, const std::string& prefix, const std::string& readsFile, 
                                const OverlapAlgorithm* pOverlapper, StringVector& filenameVec)
{
    std::string filename = prefix + RMDUPHITS_EXT + GZIP_EXT;

    std::vector<RmdupProcess*> processorVector;
    for(int i = 0; i < numThreads; ++i)
    {
        std::stringstream ss;
        ss << prefix << "-thread" << i << RMDUPHITS_EXT << GZIP_EXT;
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

std::string parseDupHits(const StringVector& hitsFilenames, const std::string& out_prefix)
{
    // Load the suffix array index and the reverse suffix array index
    // Note these are not the full suffix arrays
    SuffixArray* pFwdSAI = new SuffixArray(opt::prefix + SAI_EXT);
    SuffixArray* pRevSAI = new SuffixArray(opt::prefix + RSAI_EXT);

    // Load the read table to look up the lengths of the reads and their ids.
    // When rmduping a set of reads, the ReadInfoTable can actually be larger than the
    // BWT if the names of the reads are very long. Previously, when two reads
    // are duplicated, the read with the lexographically lower read name was chosen
    // to be kept. To save memory here, we break ties using the index in the ReadInfoTable
    // instead. This allows us to avoid loading the read names.
    ReadInfoTable* pRIT = new ReadInfoTable(opt::readsFile, pFwdSAI->getNumStrings(), RIO_NUMERICID);

    std::string outFile = out_prefix + ".fa";
    std::string dupFile = out_prefix + ".dups.fa";
    std::ostream* pWriter = createWriter(outFile);
    std::ostream* pDupWriter = createWriter(dupFile);

    size_t substringRemoved = 0;
    size_t identicalRemoved = 0;
    size_t kept = 0;
    size_t buffer_size = SequenceProcessFramework::BUFFER_SIZE;

    // The reads must be output in their original ordering.
    // The hits are in the blocks of buffer_size items. We read
    // buffer_size items from the first hits file, then buffer_size
    // from the second and so on until all the hits have been processed.
    size_t num_files = hitsFilenames.size();
    std::vector<std::istream*> reader_vec(num_files, 0);

    for(size_t i = 0; i < num_files; ++i)
    {
        std::cout << "Opening " << hitsFilenames[i] << "\n";
        reader_vec[i] = createReader(hitsFilenames[i]);
    }

    bool done = false;
    size_t currReaderIdx = 0;
    size_t numRead = 0;
    size_t numReadersDone = 0;
    std::string line;

    while(!done)
    {
        // Parse a line from the current file
        bool valid = getline(*reader_vec[currReaderIdx], line);
        ++numRead;
        // Deal with switching the active reader and the end of files
        if(!valid || numRead == buffer_size)
        {
            // Switch the reader
            currReaderIdx = (currReaderIdx + 1) % num_files;
            numRead = 0;

            // Break once all the readers are invalid
            if(!valid)
            {
                ++numReadersDone;
                if(numReadersDone == num_files)
                {
                    done = true;
                    break;
                }
            }
        }

        // Parse the data
        if(valid)
        {
            std::string id;
            std::string sequence;
            std::string hitsStr;
            size_t readIdx;
            size_t numCopies;
            bool isSubstring;

            std::stringstream parser(line);
            parser >> id;
            parser >> sequence;
            getline(parser, hitsStr);

            OverlapVector ov;
            OverlapCommon::parseHitsString(hitsStr, pRIT, pRIT, pFwdSAI, pRevSAI, true, readIdx, numCopies, ov, isSubstring);
            
            bool isContained = false;
            if(isSubstring)
            {
                ++substringRemoved;
                isContained = true;
            }
            else
            {
                for(OverlapVector::iterator iter = ov.begin(); iter != ov.end(); ++iter)
                {
                    if(iter->isContainment() && iter->getContainedIdx() == 0)
                    {
                        // This read is contained by some other read
                        ++identicalRemoved;
                        isContained = true;
                        break;
                    }
                }
            }

            SeqItem item = {id, sequence};
            std::stringstream meta;
            meta << id << " NumDuplicates=" << numCopies;

            if(isContained)
            {
                // The read's index in the sequence data base
                // is needed when removing it from the FM-index.
                // In the output fasta, we set the reads ID to be the index
                // and record its old id in the fasta header.
                std::stringstream newID;
                newID << item.id << ",seqrank=" << readIdx;
                item.id = newID.str();

                // Write some metadata with the fasta record
                item.write(*pDupWriter, meta.str());
            }
            else
            {
                ++kept;
                // Write the read
                item.write(*pWriter, meta.str());
            }
        }
    }

    for(size_t i = 0; i < num_files; ++i)
    {
        delete reader_vec[i];
        unlink(hitsFilenames[i].c_str());
    }

    
    printf("[%s] Removed %zu substring reads\n", PROGRAM_IDENT, substringRemoved);
    printf("[%s] Removed %zu identical reads\n", PROGRAM_IDENT, identicalRemoved);
    printf("[%s] Kept %zu reads\n", PROGRAM_IDENT, kept);

    // Delete allocated data
    delete pFwdSAI;
    delete pRevSAI;
    delete pRIT;
    delete pWriter;
    delete pDupWriter;

    return dupFile;
}

// 
// Handle command line arguments
//
void parseRmdupOptions(int argc, char** argv)
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
                std::cout << RMDUP_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << RMDUP_VERSION_MESSAGE;
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
        std::cerr << "Try `" << SUBPROGRAM << " --help' for more information.\n";
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
        opt::outFile = stripFilename(opt::readsFile) + ".rmdup.fa";
    }
}

