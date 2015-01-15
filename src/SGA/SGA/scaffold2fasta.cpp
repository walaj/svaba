//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// scaffold2fasta - Convert a scaffold file into a fasta
// file representing the scaffold with gaps
//
#include <getopt.h>
#include "Util.h"
#include "config.h"
#include "ScaffoldGraph.h"
#include "ScaffoldVisitors.h"
#include "ScaffoldRecord.h"
#include "SGUtil.h"
#include "OverlapTools.h"
#include "ScaffoldSequenceCollection.h"

//
void writeUnplaced(std::ostream* pWriter, StringGraph* pGraph, int minLength);


#define SUBPROGRAM "scaffold2fasta"
static const char *SCAFFOLD2FASTA_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *SCAFFOLD2FASTA_USAGE_MESSAGE =
"Usage: " SUBPROGRAM " [OPTION] ... [-f CONTIGSFILE | -a ASQGFILE] SCAFFOLDFILE\n"
"Write out a fasta file for the scaffolds indicated in SCAFFOLDFILE\n"
"One of -f CONTIGSFILE or -a ASQGFILE must be provided. If an asqg file is provided,\n"
"the program can attempt to determine the sequence linking the scaffold components by\n"
"walking the graph/\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -f, --contig-file=FILE           read the contig sequences from FILE\n"
"      -a, --asqg-file=FILE             read the contig string graph from FILE. This supercedes --contig-file\n"
"                                       this is usually the output from the sga-assemble step\n"
"          --no-singletons              do not output scaffolds that consist of a single contig\n"
"      -o, --outfile=FILE               write the scaffolds to FILE (default: scaffolds.fa)\n"
"      -m, --min-length=N               only output scaffolds longer than N bases\n"
"          --write-unplaced             output unplaced contigs that are larger than minLength\n"
"          --write-names                write the name of contigs contained in the scaffold in the FASTA header\n"
"          --min-gap-length=N           separate contigs by at least N bases. All predicted gaps less\n"
"                                       than N will be extended to N (default: 25)\n"
"          --use-overlap                attempt to merge contigs using predicted overlaps.\n"
"                                       This can help close gaps in the scaffolds but comes\n"
"                                       with a small risk of collapsing tandem repeats.\n"
"      -g, --graph-resolve=MODE         if an ASQG file is present, attempt to resolve the links\n"
"                                       between contigs using walks through the graph. The MODE parameter\n"
"                                       is a string describing the algorithm to use.\n"
"                                       The MODE parameter must be one of best-any|best-unique|unique|none.\n"
"                                       best-any: The walk with length closest to the estimated\n"
"                                       distance between the contigs will be chosen to resolve the gap.\n"
"                                       If multiple best walks are found, the tie is broken arbitrarily.\n"
"                                       best-unique: as above but if there is a tie no walk will be chosen.\n"
"                                       unique: only resolve the gap if there is a single walk between the contigs\n"
"                                       none: do not resolve gaps using the graph\n"
"                                       The most conservative most is unique, then best-unique with best-any being the most\n"
"                                       aggressive. The default is unique\n"
"      -d, --distanceFactor=T           Accept a walk as correctly resolving a gap if the walk length is within T standard \n"
"                                       deviations from the estimated distance (default: 3.0f)\n" 
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string contigFile;
    static std::string asqgFile;
    static std::string outFile = "scaffolds.fa";
    static std::string scafFile;
    static int minGapLength = 25;
    static int minOverlap = 20;
    static int maxOverlap = 100;
    static int resolveMask = RESOLVE_GRAPH_UNIQUE;
    static double maxErrorRate = 0.05f;
    static bool bNoSingletons = false;
    static bool bWriteUnplaced = false;
    static bool bWriteNames = false;
    static int minScaffoldLength = 200;
    static float distanceFactor = 3.0f;
}

static const char* shortopts = "vm:o:f:a:g:d:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_NOSINGLETON, OPT_USEOVERLAP, OPT_MINGAPLENGTH, OPT_WRITEUNPLACED, OPT_WRITENAMES };

static const struct option longopts[] = {
    { "verbose",        no_argument,       NULL, 'v' },
    { "min-length",     required_argument, NULL, 'm' },
    { "outfile",        required_argument, NULL, 'o' },
    { "contig-file",    required_argument, NULL, 'f' },
    { "asqg-file",      required_argument, NULL, 'a' },
    { "graph-resolve",  required_argument, NULL, 'g' },
    { "distanceFactor", required_argument, NULL, 'd' },
    { "min-gap-length", required_argument, NULL, OPT_MINGAPLENGTH },
    { "write-unplaced", no_argument,       NULL, OPT_WRITEUNPLACED },
    { "write-names",    no_argument,       NULL, OPT_WRITENAMES },
    { "no-singleton",   no_argument,       NULL, OPT_NOSINGLETON },
    { "use-overlap",    no_argument,       NULL, OPT_USEOVERLAP },
    { "help",           no_argument,       NULL, OPT_HELP },
    { "version",        no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
void parseScaffold2fastaOptions(int argc, char** argv);

//
int scaffold2fastaMain(int argc, char** argv)
{
    parseScaffold2fastaOptions(argc, argv);

    std::cout << "Graph best: " << (opt::resolveMask & RESOLVE_GRAPH_BEST) << "\n";
    std::cout << "Graph unique: " << (opt::resolveMask & RESOLVE_GRAPH_UNIQUE) << "\n";
    std::cout << "Find overlaps: " << (opt::resolveMask & RESOLVE_OVERLAP) << "\n";

    assert((opt::asqgFile.empty() || opt::contigFile.empty()) && !(opt::asqgFile.empty() && opt::contigFile.empty()));
    ResolveParams resolveParams;

    if(!opt::asqgFile.empty())
    {
        resolveParams.pGraph = SGUtil::loadASQG(opt::asqgFile, 0, true);
        resolveParams.pSequenceCollection = new GraphSequenceCollection(resolveParams.pGraph);
    }
    else
    {
        resolveParams.pGraph = NULL;
        resolveParams.pSequenceCollection = new MapSequenceCollection(opt::contigFile);
    }
    std::istream* pReader = new std::ifstream(opt::scafFile.c_str());
    std::ostream* pWriter = createWriter(opt::outFile);

    // Statistics tracking object
    ResolveStats stats;

    // Set up the parameters for the gap resolution function
    resolveParams.minOverlap = opt::minOverlap;
    resolveParams.maxOverlap = opt::maxOverlap;
    resolveParams.maxErrorRate = opt::maxErrorRate;
    resolveParams.resolveMask = opt::resolveMask;
    resolveParams.minGapLength = opt::minGapLength;
    resolveParams.distanceFactor = opt::distanceFactor;
    resolveParams.pStats = &stats;

    std::string line;
    size_t idx = 0;
    while(getline(*pReader, line))
    {
        ScaffoldRecord record;
        record.parse(line);
        if(record.getNumComponents() > 1 || !opt::bNoSingletons)
        {
            StringVector ids;
            std::string sequence = record.generateString(resolveParams, ids);
            
            // Write out the sequence of contigs to a stringstream
            std::stringstream contig_ss;
            contig_ss << "Contigs=";
            std::copy(ids.begin(), ids.end(), std::ostream_iterator<std::string>(contig_ss, ";"));

            //
            std::stringstream id_ss;
            id_ss << "scaffold-" << idx;

            if(opt::bWriteNames)
                id_ss << "\t" << contig_ss.str();
            
            writeFastaRecord(pWriter, id_ss.str(), sequence);
            ++idx;
        }
    }

    if(opt::bWriteUnplaced)
        resolveParams.pSequenceCollection->writeUnplaced(pWriter, opt::minScaffoldLength);

    delete resolveParams.pSequenceCollection;
    if(resolveParams.pGraph != NULL)
        delete resolveParams.pGraph;

    delete pReader;
    delete pWriter;
    return 0;
}

//
void parseScaffold2fastaOptions(int argc, char** argv)
{
    bool die = false;
    std::string modeStr;

    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case 'm': arg >> opt::minScaffoldLength; break;
            case 'f': arg >> opt::contigFile; break;
            case 'a': arg >> opt::asqgFile; break;
            case 'o': arg >> opt::outFile; break;
            case 'g': arg >> modeStr; break;
            case 'd': arg >> opt::distanceFactor; break;
            case OPT_WRITEUNPLACED: opt::bWriteUnplaced = true; break;
            case OPT_WRITENAMES: opt::bWriteNames = true; break;
            case OPT_MINGAPLENGTH: arg >> opt::minGapLength; break;
            case OPT_NOSINGLETON: opt::bNoSingletons = true; break;
            case OPT_USEOVERLAP: opt::resolveMask |= RESOLVE_OVERLAP; break;
            case OPT_HELP:
                std::cout << SCAFFOLD2FASTA_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << SCAFFOLD2FASTA_VERSION_MESSAGE;
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
        std::cout << "\n" << SCAFFOLD2FASTA_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    if(!modeStr.empty())
    {
        // Clear the graph bits of the mode string
        opt::resolveMask &= ~(RESOLVE_GRAPH_UNIQUE | RESOLVE_GRAPH_BEST);

        if(modeStr == "best-any")
        {
            opt::resolveMask |= RESOLVE_GRAPH_BEST;
        }
        else if(modeStr == "best-unique")
        {
            opt::resolveMask |= RESOLVE_GRAPH_BEST;
            opt::resolveMask |= RESOLVE_GRAPH_UNIQUE;
        }
        else if(modeStr == "unique")
        {
            opt::resolveMask |= RESOLVE_GRAPH_UNIQUE;
        }
        else if(modeStr == "none")
        {
            opt::resolveMask &= ~(RESOLVE_GRAPH_UNIQUE | RESOLVE_GRAPH_BEST);
        }
        else
        {
            std::cerr << "Unknown graph resolve mode string: " << modeStr << "\n"; 
            exit(1);
        }
    }

    // If a graph is not provided, turn off graph resolution
    if(opt::asqgFile.empty())
    {
        if(opt::resolveMask & RESOLVE_GRAPH_UNIQUE || 
           opt::resolveMask & RESOLVE_GRAPH_BEST)
        {
            std::cout << SUBPROGRAM << " WARNING: graph resolution was enabled but no graph provided\n";
            std::cout << "graph resolution will not be performed\n\n";
            opt::resolveMask &= ~(RESOLVE_GRAPH_UNIQUE | RESOLVE_GRAPH_BEST);
        }
             
    }

    // 
    opt::scafFile = argv[optind++];

    if(opt::contigFile.empty() && opt::asqgFile.empty())
    {
        std::cerr << SUBPROGRAM ": a contigs file or asqg file must be provided\n";
        exit(1);
    }
}

