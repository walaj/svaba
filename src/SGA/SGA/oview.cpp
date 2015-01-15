//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// oview - view overlap alignments
//
#include <iostream>
#include <fstream>
#include <iterator>
#include <math.h>
#include "Util.h"
#include "oview.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGUtil.h"
#include "MultiOverlap.h"

void detectMisalignments(const ReadTable* pRT, const OverlapMap* pOM);
void detect(const SeqItem& read, const ReadTable* pRT, const OverlapMap* pOM);


//
// Getopt
//
#define SUBPROGRAM "oview"
static const char *OVIEW_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *OVIEW_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... ASQGFILE\n"
"Draw overlaps in ASQGFILE\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -i, --id=ID                      only show overlaps for read with ID\n"
"      -m, --max-overhang=D             only show D overhanging bases of the alignments (default: 6)\n"
"      -d, --default-padding=D          pad the overlap lines with D characters (default: 20)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static int max_overhang = 6;
    static int padding = 20;
    static std::string prefix;
    static std::string asqgFile;
    static std::string readFilter;
}

static const char* shortopts = "p:m:i:d:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",         no_argument,       NULL, 'v' },
    { "id",              required_argument, NULL, 'i' },
    { "prefix",          required_argument, NULL, 'p' },
    { "max-overhang",    required_argument, NULL, 'm' },
    { "default-padding", required_argument, NULL, 'd' },
    { "help",            no_argument,       NULL, OPT_HELP },
    { "version",         no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int oviewMain(int argc, char** argv)
{
    parseOviewOptions(argc, argv);

    // parse reads and index them
    ReadTable* pRT = new ReadTable();
    OverlapMap* pOM = new OverlapMap;

    parseASQG(opt::asqgFile, pRT, pOM);
    pRT->indexReadsByID();

    // draw mode
    if(!opt::readFilter.empty())
    {
        drawAlignment(opt::readFilter, pRT, pOM);
    }
    else
    {
        // Output each overlap
        for(size_t i = 0; i < pRT->getCount(); ++i)
            drawAlignment(pRT->getRead(i).id, pRT, pOM);
    }

    delete pRT;
    delete pOM;
    return 0;
}

//
void drawAlignment(std::string rootID, const ReadTable* pRT, const OverlapMap* pOM)
{
    std::string rootSeq = pRT->getRead(rootID).seq.toString();
    MultiOverlap multi_overlap(rootID, rootSeq);

    // Get all the overlaps for this read
    OverlapMap::const_iterator finder = pOM->find(rootID);
    if(finder == pOM->end())
        return;

    const OverlapVector& overlaps = finder->second;
    for(size_t j = 0; j < overlaps.size(); ++j)
    {
        Overlap curr = overlaps[j];
        // Swap root read into first position if necessary
        if(curr.id[0] != rootID)
            curr.swap();
        assert(curr.id[0] == rootID);
        std::string otherSeq = pRT->getRead(curr.id[1]).seq.toString();
        multi_overlap.add(otherSeq, curr);
    }
    multi_overlap.print(opt::padding, opt::max_overhang);
}

void parseASQG(std::string filename, ReadTable* pRT, OverlapMap* pOM)
{
    std::istream* pReader = createReader(filename);
    int stage = 0;
    int line = 0;
    std::string recordLine;
    while(getline(*pReader, recordLine))
    {
        ASQG::RecordType rt = ASQG::getRecordType(recordLine);
        switch(rt)
        {
            case ASQG::RT_HEADER:
            {
                if(stage != 0)
                {
                    std::cerr << "Error: Unexpected header record found at line " << line << "\n";
                    exit(EXIT_FAILURE);
                }

                ASQG::HeaderRecord headerRecord(recordLine);
                (void)headerRecord; // do nothing with the header at this point
                break;
            }
            case ASQG::RT_VERTEX:
            {
                // progress the stage if we are done the header
                if(stage == 0)
                    stage = 1;

                if(stage != 1)
                {
                    std::cerr << "Error: Unexpected vertex record found at line " << line << "\n";
                    exit(EXIT_FAILURE);
                }

                ASQG::VertexRecord vertexRecord(recordLine);
                SeqItem si = { vertexRecord.getID(), vertexRecord.getSeq() };
                pRT->addRead(si);
                break;
            }
            case ASQG::RT_EDGE:
            {
                if(stage == 1)
                    stage = 2;
                
                if(stage != 2)
                {
                    std::cerr << "Error: Unexpected edge record found at line " << line << "\n";
                    exit(EXIT_FAILURE);
                }

                ASQG::EdgeRecord edgeRecord(recordLine);
                const Overlap& ovr = edgeRecord.getOverlap();
                if(opt::readFilter.empty() || ovr.id[0] == opt::readFilter || ovr.id[1] == opt::readFilter)
                {
                    (*pOM)[ovr.id[0]].push_back(ovr);
                    (*pOM)[ovr.id[1]].push_back(ovr);
                }
            }
        }
        ++line;
    }
    delete pReader;
}


// 
// Handle command line arguments
//
void parseOviewOptions(int argc, char** argv)
{
    bool die = false;

    // Defaults
    opt::max_overhang = 6;

    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'p': arg >> opt::prefix; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case 'i': arg >> opt::readFilter; break;
            case 'm': arg >> opt::max_overhang; break;
            case 'd': arg >> opt::padding; break;
            case OPT_HELP:
                std::cout << OVIEW_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << OVIEW_VERSION_MESSAGE;
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
        std::cout << "\n" << OVIEW_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::asqgFile = argv[optind++];

    if(opt::prefix.empty())
        opt::prefix = stripFilename(opt::asqgFile);
}
