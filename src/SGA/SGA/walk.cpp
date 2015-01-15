//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// walk - extract a walks between vertices from the graph
//
#include <iostream>
#include <fstream>
#include "walk.h"
#include "Util.h"
#include "SGUtil.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"
#include "Timer.h"
#include "SGSearch.h"
#include "Bigraph.h"

// functions
void pairWalk(StringGraph* pGraph, SGWalkVector& outWalks);
void explicitWalk(StringGraph* pGraph, SGWalkVector& outWalks);
void componentWalk(StringGraph* pGraph, SGWalkVector& outWalks);

//
// Getopt
//
#define SUBPROGRAM "walk"
static const char *WALK_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *WALK_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... ASQGFILE\n"
"Extract walks from ID1 to ID2 from the graph.\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -d, --distance=NUM               the maximum length of the walks to find.\n"
"      -s,--start ID                    start the walk at vertex with ID\n"
"      -e,--end ID                      end the walk at vertex with ID\n"
"      -w,walk-str STR                  extract the defined walk given by STR\n"
"                                       STR should be a comma delimited list of IDs\n"
"         --component-walks             find all possible walks through the largest connected component\n"
"                                       of the graph.\n"
"         --longest-n=NUM               when using --component-walks, only output the longest NUM walks\n"
"      -o,--out-file=FILE               write the walks to FILE in FASTA format (default: walks.fa)\n"
"      -p, --prefix=PREFIX              write final walks sequence ids starting with PREFIX (default: walk-N)\n"
"         --sam=FILE                    write the walk descriptions to FILE in SAM format\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string asqgFile;
    static std::string walkStr;
    static std::string id1;
    static std::string id2;
    static std::string outFile;
    static std::string prefix;
    static std::string samFile;
    static int maxDistance = 500;
    static bool componentWalks = false;
    static int numOutputWalks = -1;
}

static const char* shortopts = "o:d:s:e:w:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_COMPONENT, OPT_SAM, OPT_LONGEST_N };

static const struct option longopts[] = {
    { "verbose",           no_argument,       NULL, 'v' },
    { "out-file",          required_argument, NULL, 'o' },
    { "prefix",            required_argument, NULL, 'p' },
    { "sam",               required_argument, NULL, OPT_SAM },
    { "distance",          required_argument, NULL, 'd' },
    { "start",             required_argument, NULL, 's' },
    { "end",               required_argument, NULL, 'e' },
    { "walk-str",          required_argument, NULL, 'w' },
    { "component-walks",   no_argument,       NULL, OPT_COMPONENT },
    { "longest-n",         required_argument, NULL, OPT_LONGEST_N },
    { "help",              no_argument,       NULL, OPT_HELP },
    { "version",           no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int walkMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("sga walk");
    parseWalkOptions(argc, argv);
    walk();
    delete pTimer;

    return 0;
}

void walk()
{
    StringGraph* pGraph = SGUtil::loadASQG(opt::asqgFile, 0, true);
    pGraph->printMemSize();
    std::ostream* pWriter = createWriter(opt::outFile);
    std::ostream* pSAMWriter = NULL;
    if(!opt::samFile.empty())
        pSAMWriter = createWriter(opt::samFile);

    // Build walks
    SGWalkVector walkVector;
    if(!opt::id1.empty() && !opt::id2.empty())
        pairWalk(pGraph, walkVector);
    else if(!opt::walkStr.empty())
        explicitWalk(pGraph, walkVector);
    else if(opt::componentWalks)
        componentWalk(pGraph, walkVector);
    else
        assert(false);
    
    // Output walks
    // The first pass over the vector writes the walk sequences to a FASTA
    // file and builds a vector with their lengths. The second pass
    // writes a SAM file of the reads making up the walk
    StringVector walkNames;
    for(size_t i = 0; i < walkVector.size(); ++i)
    {
        SGWalk& walk = walkVector[i];

        std::stringstream idSS;
        if(opt::prefix.empty())
            idSS << "walk-" << i;
        else
            idSS << opt::prefix << "-walk-" << i;

        std::string walkID = idSS.str();
        walkNames.push_back(walkID);

        std::string str = walk.getString(SGWT_START_TO_END);
        if(opt::verbose > 0)
        {
            std::cout << walkID << "\n";
            std::cout << "Str: " << str << "\n";
            std::cout << "Len: " << str.length() << "\n";
            walk.print();
        }
        writeFastaRecord(pWriter, walkID, str);
        
        // Print a header line to the SAM file
        if(pSAMWriter != NULL)
            *pSAMWriter << "@SQ\tSN:" << walkID << "\tLN:" << str.length() << "\n";
    }
    
    // Write the placement of the sequences in the walk as SAM
    if(pSAMWriter != NULL)
    {
        for(size_t i = 0; i < walkVector.size(); ++i)
        {
            SGWalk& walk = walkVector[i];
            SGWalkVertexPlacementVector placementVector;
            walk.getString(SGWT_START_TO_END, &placementVector);

            for(size_t j = 0; j < placementVector.size(); ++j)
            {
                SGWalkVertexPlacement& placement = placementVector[j];
                *pSAMWriter << placement.pVertex->getID() << "\t"; // read name
                *pSAMWriter << (placement.isRC ? "16" : "0") << "\t"; // flag
                *pSAMWriter << walkNames[i] << "\t"; // ref name
                *pSAMWriter << placement.position + 1 << "\t"; // ref postion, 1-based
                *pSAMWriter << 99 << "\t"; // map quality
                *pSAMWriter << placement.pVertex->getSeqLen() << "M\t"; // CIGAR
                *pSAMWriter << "*\t"; // mate reference
                *pSAMWriter << "0\t"; // mate position
                *pSAMWriter << "0\t"; // template length
                
                std::string seq = placement.pVertex->getSeq().toString();
                if(placement.isRC)
                    seq = reverseComplement(placement.pVertex->getSeq().toString());

                *pSAMWriter << seq << "\t"; // sequence
                *pSAMWriter << "*"; // quality
                *pSAMWriter << "\n";
            }
        }
        delete pSAMWriter;
    }

    delete pGraph;
    delete pWriter;
}

// Find walks between a pair of vertices
void pairWalk(StringGraph* pGraph, SGWalkVector& outWalks)
{
    // Search the the walk between id1 and id2
    Vertex* pX = pGraph->getVertex(opt::id1);
    Vertex* pY = pGraph->getVertex(opt::id2);
    SGSearch::findWalks(pX, pY, ED_SENSE, opt::maxDistance, 1000, false, outWalks);
    SGSearch::findWalks(pX, pY, ED_ANTISENSE, opt::maxDistance, 1000, false, outWalks);
}

// Find walks using an explicitly provided string
void explicitWalk(StringGraph* pGraph, SGWalkVector& outWalks)
{
    std::cout << "Building walk from string: " << opt::walkStr << "\n";
    StringVector vertIDs = split(opt::walkStr, ',');
    assert(vertIDs.size() > 0);
    Vertex* pCurr = pGraph->getVertex(vertIDs.front());
    SGWalk out(pCurr);
    
    for(size_t i = 1; i < vertIDs.size(); ++i)
    {
        Vertex* pNext = pGraph->getVertex(vertIDs[i]);
        EdgePtrVec epv = pCurr->findEdgesTo(pNext->getID());
        if(epv.size() == 0)
        {
            std::cerr << "[sga walk] Error: edge between " << pCurr->getID() << " and " << pNext->getID() << " not found\n";
            exit(EXIT_FAILURE);
        }
        out.addEdge(epv.front());
        pCurr = pNext;
    }

    outWalks.push_back(out);
}

// Find all walks through the largest component of the graph
void componentWalk(StringGraph* pGraph, SGWalkVector& outWalks)
{
    // Perform some simplifications on the graph
    
    std::cout << "Component walK: Removing contained vertices and lingering transitive edges\n";
    SGContainRemoveVisitor containVisit;
    pGraph->visit(containVisit);

    SGTransitiveReductionVisitor trVisit;
    pGraph->visit(trVisit);
    pGraph->stats();

    typedef std::vector<VertexPtrVec> ComponentVector;
    VertexPtrVec allVertices = pGraph->getAllVertices();
    ComponentVector components;
    SGSearchTree::connectedComponents(allVertices, components);

    // Select the largest component
    int selectedIdx = -1;
    size_t largestSize = 0;

    for(size_t i = 0; i < components.size(); ++i)
    {
        std::cout << "Component " << i << ": " << components[i].size() << " vertices\n";
        if(components[i].size() > largestSize)
        {
            selectedIdx = i;
            largestSize = components[i].size();
        }
    }

    assert(selectedIdx != -1);
    VertexPtrVec selectedComponent = components[selectedIdx];

    std::cout << "component-walk: selected component of size " << selectedComponent.size() << "\n";

    // Build a vector of the terminal vertices
    VertexPtrVec terminals;
    for(size_t i = 0; i < selectedComponent.size(); ++i)
    {
        Vertex* pVertex = selectedComponent[i];
        size_t asCount = pVertex->getEdges(ED_ANTISENSE).size();
        size_t sCount = pVertex->getEdges(ED_SENSE).size();

        if(asCount == 0 || sCount == 0)
            terminals.push_back(pVertex);
    }

    std::cout << "selected component has " << terminals.size() << " terminal vertices\n";

    // Find walks between all-pairs of terminal vertices
    SGWalkVector tempWalks;
    for(size_t i = 0; i < terminals.size(); ++i)
    {
        for(size_t j = i + 1; j < terminals.size(); j++)
        {
            Vertex* pX = terminals[i];
            Vertex* pY = terminals[j];
            SGSearch::findWalks(pX, pY, ED_SENSE, opt::maxDistance, 1000000, false, tempWalks);
            SGSearch::findWalks(pX, pY, ED_ANTISENSE, opt::maxDistance, 1000000, false, tempWalks);   
        }
    }

    // Remove duplicate walks
    std::map<std::string, SGWalk> walkMap;
    for(size_t i = 0; i < tempWalks.size(); ++i)
    {
        std::string walkString = tempWalks[i].getString(SGWT_START_TO_END);
        walkMap.insert(std::make_pair(walkString, tempWalks[i]));
    }

    // Copy unique walks to the output
    for(std::map<std::string, SGWalk>::iterator mapIter = walkMap.begin(); mapIter != walkMap.end(); ++mapIter)
        outWalks.push_back(mapIter->second);

    // Sort the walks by string length
    std::sort(outWalks.begin(), outWalks.end(), SGWalk::compareByTotalLength);

    if(opt::numOutputWalks > 0 && opt::numOutputWalks < (int)outWalks.size())
    {
        assert(outWalks.begin() + opt::numOutputWalks < outWalks.end());
        outWalks.erase(outWalks.begin() + opt::numOutputWalks, outWalks.end());
    }
}

// 
// Handle command line arguments
//
void parseWalkOptions(int argc, char** argv)
{
    // Set defaults
    opt::outFile = "walks.fa";

    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'o': arg >> opt::outFile; break;
            case 'p': arg >> opt::prefix; break;
            case 'd': arg >> opt::maxDistance; break;
            case 's': arg >> opt::id1; break;
            case 'e': arg >> opt::id2; break;
            case 'w': arg >> opt::walkStr; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_SAM: arg >> opt::samFile; break;
            case OPT_COMPONENT: opt::componentWalks = true; break;
            case OPT_LONGEST_N: arg >> opt::numOutputWalks; break;
            case OPT_HELP:
                std::cout << WALK_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << WALK_VERSION_MESSAGE;
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

    if(!opt::componentWalks && opt::walkStr.empty() && (opt::id1.empty() || opt::id2.empty()))
    {
        std::cerr << SUBPROGRAM ": one of a start/end vertex pair (-s/-e), a walk string (-w) or --component-walks must be provided\n";
        exit(EXIT_FAILURE);
    }

    if(opt::numOutputWalks == 0)
    {
        std::cerr << SUBPROGRAM ": error --longest-n set to zero. Nothing to do\n";
        exit(EXIT_FAILURE);
    }

    if(opt::id1.empty() != opt::id2.empty())
    {
        std::cerr << SUBPROGRAM ": both a start and end vertex must be provided\n";
        exit(EXIT_FAILURE);
    }

    // Parse the input filename
    opt::asqgFile = argv[optind++];
}
