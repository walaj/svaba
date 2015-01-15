//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// extract - extract all sequences of a given length from a bwt
//
#include <iostream>
#include <fstream>
#include <iterator>
#include "extract.h"
#include "Util.h"
#include "BWT.h"
#include "SGACommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "BWTTraverse.h"

//
// Getopt
//
#define SUBPROGRAM "extract"
static const char *EXTRACT_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *EXTRACT_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Extract all sequences of a given length from a BWT\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -l, --length=LEN                 extract sequences of length LEN\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static unsigned int length;
    static std::string prefix;
    static std::string readsFile;
}

static const char* shortopts = "p:l:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "length",      required_argument, NULL, 'l' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int extractMain(int argc, char** argv)
{
    parseExtractOptions(argc, argv);

    // Load the BWT for this read set
    BWT* pBWT = new BWT(opt::prefix + BWT_EXT);
    BWT* pRevBWT = new BWT(opt::prefix + RBWT_EXT);

    BWTTraverse::extractSG(pBWT, pRevBWT, opt::length);

    delete pBWT;
    delete pRevBWT;
    return 0;
}

// 
// Handle command line arguments
//
void parseExtractOptions(int argc, char** argv)
{
    // Set defaults
    opt::length = DEFAULT_EXTRACT_LEN;

    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'l': arg >> opt::length; break;
            case 'p': arg >> opt::prefix; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << EXTRACT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << EXTRACT_VERSION_MESSAGE;
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
        std::cout << "\n" << EXTRACT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::readsFile = argv[optind++];

    if(opt::prefix.empty())
    {
        opt::prefix = stripFilename(opt::readsFile);
    }
}


