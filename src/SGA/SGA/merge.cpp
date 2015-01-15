//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// merge - merge read files and their indices
//
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include "SGACommon.h"
#include "Util.h"
#include "merge.h"
#include "SACAInducedCopying.h"
#include "BWTDiskConstruction.h"
#include "BWT.h"
#include "PopulationIndex.h"
#include <unistd.h>

//
void removeFiles(const std::string& inFile);

//
// Getopt
//
#define SUBPROGRAM "merge"

static const char *MERGE_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *MERGE_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READS1 READS2\n"
"Merge the sequence files READS1, READS2 into a single file/index\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"  -t, --threads=NUM                    use NUM threads to merge the indices (default: 1)\n"
"  -p, --prefix=PREFIX                  write final index to files starting with PREFIX (the default is to concatenate the input filenames)\n"
"  -r, --remove                         remove the original BWT, SAI and reads files after the merge\n"
"  -g, --gap-array=N                    use N bits of storage for each element of the gap array. Acceptable values are 4,8,16 or 32. Lower\n"
"                                       values can substantially reduce the amount of memory required at the cost of less predictable memory usage.\n"
"                                       When this value is set to 32, the memory requirement is essentially deterministic and requires ~5N bytes where\n"
"                                       N is the size of the FM-index of READS2.\n"
"                                       The default value is 4.\n"
"      --no-sequence                    Suppress merging of the sequence files. Use this option when merging the index(es) separate e.g. in parallel\n"
"      --no-forward                     Suppress merging of the forward index. Use this option when merging the index(es) separate e.g. in parallel\n"
"      --no-reverse                     Suppress merging of the reverse index. Use this option when merging the index(es) separate e.g. in parallel\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string prefix;
    static int numThreads = 1;
    static bool bRemove;
    static int gapArrayStorage = 4;
	static bool bMergeSequence = true;
	static bool bMergeForward = true;
	static bool bMergeReverse = true;
}

static const char* shortopts = "p:m:t:g:vr";

enum { OPT_HELP = 1, OPT_VERSION, OPT_NO_SEQUENCE, OPT_NO_FWD, OPT_NO_REV };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "prefix",      required_argument, NULL, 'p' },
    { "remove",      no_argument,       NULL, 'r' },
    { "threads",     required_argument, NULL, 't' },
    { "gap-array",   required_argument, NULL, 'g' },
    { "no-sequence", no_argument,       NULL, OPT_NO_SEQUENCE },
    { "no-forward", no_argument,       NULL, OPT_NO_FWD },
    { "no-reverse", no_argument,       NULL, OPT_NO_REV },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

int mergeMain(int argc, char** argv)
{
    parseMergeOptions(argc, argv);
    StringVector inFiles;
    while(optind < argc)
    {
        inFiles.push_back(argv[optind++]);
    }

    assert(inFiles.size() == 2);
    if(inFiles[0] == inFiles[1])
        return 0; // avoid self-merge

    if(opt::prefix.empty())
    {
        std::string basename1 = stripFilename(inFiles[0]);
        std::string basename2 = stripFilename(inFiles[1]);
        opt::prefix = basename1 + "." + basename2;
    }

    // Merge the indices
	if(opt::bMergeForward)
	{
		mergeIndependentIndices(inFiles[0], inFiles[1], opt::prefix, BWT_EXT, SAI_EXT, false, opt::numThreads, opt::gapArrayStorage);
	}
    
    std::string prefix1 = stripGzippedExtension(inFiles[0]);
    std::string prefix2 = stripGzippedExtension(inFiles[1]);

    // Skip merging the reverse indices if the reverse bwt file does not exist. 
    std::string rbwt_filename_1 = prefix1 + RBWT_EXT;
    std::string rbwt_filename_2 = prefix2 + RBWT_EXT;

    struct stat file_s_1;
    struct stat file_s_2;
    int ret1 = stat(rbwt_filename_1.c_str(), &file_s_1);
    int ret2 = stat(rbwt_filename_2.c_str(), &file_s_2);

    if((ret1 == 0 || ret2 == 0) && opt::bMergeReverse)
	{
		mergeIndependentIndices(inFiles[0], inFiles[1], opt::prefix, RBWT_EXT, RSAI_EXT, true, opt::numThreads, opt::gapArrayStorage);
	}
		
    // Merge the read files
	if(opt::bMergeSequence)
	{
		mergeReadFiles(inFiles[0], inFiles[1], opt::prefix);
	}

    // Merge any population index files
    std::string popidx_filename_1 = prefix1 + POPIDX_EXT;
    std::string popidx_filename_2 = prefix2 + POPIDX_EXT;
    ret1 = stat(popidx_filename_1.c_str(), &file_s_1);
    ret2 = stat(popidx_filename_2.c_str(), &file_s_2);
    if(ret1 == 0 && ret2 == 0)
        PopulationIndex::mergeIndexFiles(popidx_filename_1, popidx_filename_2, opt::prefix + POPIDX_EXT);

    if(opt::bRemove)
    {
        // Delete the original reads, bwt and sai files
        removeFiles(inFiles[0]);
        removeFiles(inFiles[1]);
    }
    return 0;
}

//
void removeFiles(const std::string& inFile)
{
    std::string prefix = stripFilename(inFile);
    std::string bwt_file = prefix + BWT_EXT;
    std::string rbwt_file = prefix + RBWT_EXT;

    std::string sai_file = prefix + SAI_EXT;
    std::string rsai_file = prefix + RSAI_EXT;

    unlink(bwt_file.c_str());
    unlink(rbwt_file.c_str());
    unlink(sai_file.c_str());
    unlink(rsai_file.c_str());
    unlink(inFile.c_str());
}

// 
// Handle command line arguments
//
void parseMergeOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'p': arg >> opt::prefix; break;
            case 'r': opt::bRemove = true; break;
            case '?': die = true; break;
            case 't': arg >> opt::numThreads; break;
            case 'g': arg >> opt::gapArrayStorage; break;
            case 'v': opt::verbose++; break;
			case OPT_NO_SEQUENCE: opt::bMergeSequence = false; break;
			case OPT_NO_FWD: opt::bMergeForward = false; break;
			case OPT_NO_REV: opt::bMergeReverse = false; break;
            case OPT_HELP:
                std::cout << MERGE_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << MERGE_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if(opt::gapArrayStorage != 4 && opt::gapArrayStorage != 8 &&
       opt::gapArrayStorage != 16 && opt::gapArrayStorage != 32)
    {
        std::cerr << SUBPROGRAM ": invalid argument, --gap-array,-g must be one of 4,8,16,32 (found: " << opt::gapArrayStorage << ")\n";
        die = true;
    }

    if(opt::numThreads <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
        die = true;
    }    

    if (argc - optind < 1) 
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } 

    if (argc - optind > 2)
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << MERGE_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}
