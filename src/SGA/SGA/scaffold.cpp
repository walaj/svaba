//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// scaffold - Scaffold a set of contigs using
// distances estimates and/or reference information.
//
#include <getopt.h>
#include "Util.h"
#include "config.h"
#include "ScaffoldGraph.h"
#include "ScaffoldVisitors.h"
#include "ScaffoldSearch.h"
#include "ScaffoldAlgorithms.h"

#define SUBPROGRAM "scaffold"
static const char *SCAFFOLD_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *SCAFFOLD_USAGE_MESSAGE =
"Usage: " SUBPROGRAM " [OPTION] ... [--pe PE-DE] [--mate-pair MATEPAIR-DE] CONTIGSFILE\n"
"Construct scaffolds from CONTIGSFILE using distance estimates. \n"
"The distance estimates are read from the --pe and --matepair parameters\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"          --pe=FILE                    load links derived from paired-end (short insert) libraries from FILE\n"
"          --mate-pair=FILE             load links derived from mate-pair (long insert) libraries from FILE\n"
"      -m, --min-length=N               only use contigs at least N bp in length to build scaffolds (default: no minimun).\n"
"      -g, --asqg-file=FILE             optionally load the sequence graph from FILE\n"
"      -a, --astatistic-file=FILE       load Myers' A-statistic values from FILE. This is used to\n"
"                                       determine unique and repetitive contigs with the -u/--unique-astat\n"
"                                       and -r/--repeat-astat parameters (required)\n"
"      -u, --unique-astat=FLOAT         Contigs with an a-statitic value about FLOAT will be considered unique (default: 20.0)\n"
"      -c, --min-copy-number=FLOAT      remove vertices with estimated copy number less than FLOAT (default: 0.5f)\n"
"      -s, --max-sv-size=N              collapse heterozygous structural variation if the event size is less than N (default: 0)\n"
"      -o, --outfile=FILE               write the scaffolds to FILE (default: CONTIGSFILE.scaf\n"
"          --remove-conflicting         if two contigs have multiple distance estimates between them and they do not agree, break the scaffold\n"
"                                       at this point\n"
"          --strict                     perform strict consistency checks on the scaffold links. If a vertex X has multiple edges, a path will\n"
"                                       be searched for that contains every vertex linked to X. If no such path can be found, the edge of X are removed.\n"
"                                       This builds very conservative scaffolds that should be highly accurate.\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string contigsFile;
    static StringVector peDistanceEstFiles;
    static StringVector mateDistanceEstFiles;
    static std::string astatFile;
    static std::string outFile;
    static std::string asqgFile;
    static bool strict = false;
    static bool removeConflicting = false;
    static double uniqueAstatThreshold = 20.0f;
    static double minEstCopyNumber = 0.3f;
    static int maxSVSize = 0;
    static int minContigLength = 0;
}

static const char* shortopts = "vm:a:u:r:o:g:s:c:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_PE, OPT_MATEPAIR, OPT_CUTCONFLICT, OPT_STRICT };

static const struct option longopts[] = {
    { "verbose",            no_argument,       NULL, 'v' },
    { "min-length",         required_argument, NULL, 'm' },
    { "asgq-file",          required_argument, NULL, 'g' }, 
    { "astatistic-file",    required_argument, NULL, 'a' },
    { "unique-astat",       required_argument, NULL, 'u' },
    { "repeat-astat",       required_argument, NULL, 'r' },
    { "outfile",            required_argument, NULL, 'o' },
    { "max-sv-size",        required_argument, NULL, 's' },
    { "min-copy-number",    required_argument, NULL, 'c' },
    { "remove-conflicting", no_argument,       NULL, OPT_CUTCONFLICT },
    { "strict",             no_argument,       NULL, OPT_STRICT },
    { "pe",                 required_argument, NULL, OPT_PE },
    { "mate-pair",          required_argument, NULL, OPT_MATEPAIR },
    { "help",               no_argument,       NULL, OPT_HELP },
    { "version",            no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
void parseScaffoldOptions(int argc, char** argv);

//
int scaffoldMain(int argc, char** argv)
{
    parseScaffoldOptions(argc, argv);

    int maxOverlap = 100;

    ScaffoldGraph graph;
    
    graph.loadVertices(opt::contigsFile, opt::minContigLength);

    for(size_t i = 0; i < opt::peDistanceEstFiles.size(); ++i)
        graph.loadDistanceEstimateEdges(opt::peDistanceEstFiles[i], false, opt::verbose);
    
    for(size_t i = 0; i < opt::mateDistanceEstFiles.size(); ++i)
        graph.loadDistanceEstimateEdges(opt::mateDistanceEstFiles[i], true, opt::verbose);

    // Load the a-stat data and mark vertices as unique and repeat
    if(!opt::astatFile.empty())
    {
        graph.loadAStatistic(opt::astatFile);
        ScaffoldAStatisticVisitor astatVisitor(opt::uniqueAstatThreshold, 
                                               opt::minEstCopyNumber);
        graph.visit(astatVisitor);
    }
    else
    {
        std::cout << "=====\nWARNING -- no a-stat file provided, assuming all vertices are not repeats\n";
        std::cout << "It is highly suggested that an a-stat file is provided " <<
                     "to give copy number estimates of the input sequences\n=====\n\n";
    }

    std::cout << "[sga-scaffold] Removing non-unique vertices from scaffold graph\n";
    graph.deleteVertices(SVC_REPEAT);
    if(opt::removeConflicting)
    {
        ScaffoldConflictingVisitor conflictVisitor;
        graph.visit(conflictVisitor);
        graph.deleteVertices(SVC_REPEAT);
    }

    if(opt::strict)
    {
        std::cout << "Performing strict resolutions\n";
        ScaffoldTransitiveReductionVisitor trVisit;
        graph.visit(trVisit);
    
        // Check for cycles in the graph
        ScaffoldAlgorithms::destroyStrictCycles(&graph, "scaffold.cycles.out");

        ScaffoldMultiEdgeRemoveVisitor meVisit;
        graph.visit(meVisit);
    }
    else
    {
        // Remove polymorphic nodes from the graph
        ScaffoldPolymorphismVisitor polyVisitor(maxOverlap);
        while(graph.visit(polyVisitor)) {}

        // If requested, collapse structural variation in the graph
        if(opt::maxSVSize > 0)
        {
            ScaffoldSVVisitor svVisit(opt::maxSVSize);
            graph.visit(svVisit);
        }

        // Break any links in the graph that are inconsistent
        ScaffoldLinkValidator linkValidator(100, 0.05f, opt::verbose);
        graph.visit(linkValidator);
        graph.deleteVertices(SVC_REPEAT);
        
        // Check for cycles in the graph using the old cycle finding algorithm
        ScaffoldAlgorithms::removeInternalCycles(&graph);
    }


    // Linearize the scaffolds
    ScaffoldAlgorithms::makeScaffolds(&graph);

    // TODO Place floating contigs and repeats into the gaps.

    // Break any remaining multi-edge contigs scaffolds
    ScaffoldMultiEdgeRemoveVisitor cutVisitor;
    graph.visit(cutVisitor);

    // Write out the scaffolds
    ScaffoldWriterVisitor writer(opt::outFile);
    graph.visit(writer);

    return 0;
}

//
void parseScaffoldOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case 'm': arg >> opt::minContigLength; break;
            case 'a': arg >> opt::astatFile; break;
            case 'g': arg >> opt::asqgFile; break;
            case 'u': arg >> opt::uniqueAstatThreshold; break;
            case 'o': arg >> opt::outFile; break;
            case 's': arg >> opt::maxSVSize; break;
            case 'c': arg >> opt::minEstCopyNumber; break;
            case OPT_CUTCONFLICT: opt::removeConflicting = true; break;
            case OPT_STRICT: opt::strict = true; break;
            case OPT_PE: 
                opt::peDistanceEstFiles.push_back(arg.str()); 
                break;
            case OPT_MATEPAIR: 
                opt::mateDistanceEstFiles.push_back(arg.str()); 
                break;
            case OPT_HELP:
                std::cout << SCAFFOLD_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << SCAFFOLD_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 1) 
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
        std::cout << "\n" << SCAFFOLD_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // 
    opt::contigsFile = argv[optind++];

    if(opt::contigsFile.empty())
    {
        std::cerr << SUBPROGRAM ": a contigs file must be provided\n";
        exit(1);
    }

    if(opt::peDistanceEstFiles.empty() && opt::mateDistanceEstFiles.empty())
    {
        std::cerr << SUBPROGRAM ": at least one distance estimation file must be provided using the --pe and/or --mate-pair options\n";
        exit(1);
    }

    /*
    if(opt::astatFile.empty())
    {
        std::cerr << SUBPROGRAM ": an a-statistic file must be provided\n";
        exit(1);
    }
    */

    if(opt::outFile.empty())
    {
        opt::outFile = opt::contigsFile + ".scaf";
    }

    
}

