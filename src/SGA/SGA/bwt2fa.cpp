//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// bwt2fa - Transform a bwt back into a set of sequences
//
#include <iostream>
#include <fstream>
#include "SGACommon.h"
#include "Util.h"
#include "bwt2fa.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "Timer.h"
#include "BWTAlgorithms.h"

//
// Getopt
//
#define SUBPROGRAM "bwt2fa"

static const char *BWT2FA_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2011 Wellcome Trust Sanger Institute\n";

static const char *BWT2FA_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... BWTFILE\n"
"Transform the bwt BWTFILE back into a set of sequences\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -o,--outfile=FILE                write the sequences to FILE\n"
"      -p,--prefix=STR                  prefix the names of the reads with STR\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string bwtFile;
    static std::string outFile;
    static std::string readPrefix;
    static int sampleRate = 256;
}

static const char* shortopts = "p:o:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_NO_REVERSE };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "prefix",      required_argument, NULL, 'p' },
    { "outfile",     required_argument, NULL, 'o' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

int bwt2faMain(int argc, char** argv)
{
    Timer t("sga bwt2fa");
    parseBWT2FAOptions(argc, argv);

    BWT* pBWT = new BWT(opt::bwtFile, opt::sampleRate);
    pBWT->printInfo();

    std::ostream* pWriter = createWriter(opt::outFile);

    SeqItem outItem;
    outItem.id = "";
    size_t n = pBWT->getNumStrings();
    for(size_t i = 0; i < n; ++i)
    {
        std::stringstream nameSS;
        nameSS << opt::readPrefix << "-" << i;
        outItem.id = nameSS.str();
        outItem.seq = BWTAlgorithms::extractString(pBWT, i);
        outItem.write(*pWriter);
    }

    delete pBWT;
    delete pWriter;

    return 0;
}

// 
// Handle command line arguments
//
void parseBWT2FAOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'p': arg >> opt::readPrefix; break;
            case '?': die = true; break;
            case 'o': arg >> opt::outFile; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << BWT2FA_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << BWT2FA_VERSION_MESSAGE;
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
        std::cout << "\n" << BWT2FA_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::bwtFile = argv[optind++];
    if(opt::readPrefix.empty())
        opt::readPrefix = stripFilename(opt::bwtFile);

    if(opt::outFile.empty())
        opt::outFile = stripFilename(opt::bwtFile) + ".fa";
}
