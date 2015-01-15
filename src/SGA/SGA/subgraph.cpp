//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// subgraph - extract a subgraph from an assembly graph
//
#include <iostream>
#include <fstream>
#include "subgraph.h"
#include "Util.h"
#include "SGUtil.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"
#include "Timer.h"

// functions
void addNeighborsToSubgraph(Vertex* pCurrVertex, StringGraph* pSubgraph, int span);
void copyVertexToSubgraph(StringGraph* pSubgraph, const Vertex* pVertex);


//
// Getopt
//
#define SUBPROGRAM "subgraph"
static const char *SUBGRAPH_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *SUBGRAPH_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... ID ASQGFILE\n"
"Extract the subgraph around the sequence with ID from an asqg file.\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -o, --out=FILE                   write the subgraph to FILE (default: subgraph.asqg.gz)\n"
"      -s, --size=N                     the size of the subgraph to extract, all vertices that are at most N hops\n"
"                                       away from the root will be included (default: 5)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string asqgFile;
    static std::string outFile;
    static std::string rootID;
    static unsigned int span = 5;
}

static const char* shortopts = "o:s:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",        no_argument,       NULL, 'v' },
    { "out",            required_argument, NULL, 'o' },
    { "size",           required_argument, NULL, 's' },
    { "help",           no_argument,       NULL, OPT_HELP },
    { "version",        no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int subgraphMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("sga subgraph");
    parseSubgraphOptions(argc, argv);
    subgraph();
    delete pTimer;

    return 0;
}

void subgraph()
{
    StringGraph* pGraph = SGUtil::loadASQG(opt::asqgFile, 0, true);
    pGraph->printMemSize();

    /*
    // Remove containments from the graph
    SGContainRemoveVisitor containVisit;
    std::cout << "Removing contained vertices\n";
    while(pGraph->hasContainment())
    {
        pGraph->visit(containVisit);
    }
    */
    StringGraph* pSubgraph = new StringGraph;

    // Set the graph parameters to match the main graph
    pSubgraph->setContainmentFlag(pGraph->hasContainment());
    pSubgraph->setTransitiveFlag(pGraph->hasTransitive());
    pSubgraph->setMinOverlap(pGraph->getMinOverlap());
    pSubgraph->setErrorRate(pGraph->getErrorRate());

    // Get the root vertex
    Vertex* pRootVertex = pGraph->getVertex(opt::rootID);
    if(pRootVertex == NULL)
    {
        std::cout << "Vertex " << opt::rootID << " not found in the graph.\n";
    }
    else
    {
        copyVertexToSubgraph(pSubgraph, pRootVertex);
        pRootVertex->setColor(GC_BLACK);

        // Recursively add neighbors
        addNeighborsToSubgraph(pRootVertex, pSubgraph, opt::span);

        // Write the subgraph
        pSubgraph->writeASQG(opt::outFile);
    }

    delete pSubgraph;
    delete pGraph;
}

//
void addNeighborsToSubgraph(Vertex* pCurrVertex, StringGraph* pSubgraph, int span)
{
    if(span <= 0)
        return;

    // These are the edges in the main graph
    EdgePtrVec edges = pCurrVertex->getEdges();
    for(size_t i = 0; i < edges.size(); ++i)
    {
        if(edges[i]->getColor() != GC_BLACK)
        {
            Vertex* pY = edges[i]->getEnd();
            copyVertexToSubgraph(pSubgraph, pY);
            Overlap ovr = edges[i]->getOverlap();
            SGAlgorithms::createEdgesFromOverlap(pSubgraph, ovr, true);
            edges[i]->setColor(GC_BLACK);
            edges[i]->getTwin()->setColor(GC_BLACK);

            // Recurse
            addNeighborsToSubgraph(pY, pSubgraph, span - 1);
        }
    }
}

// We do not directly add the vertex pointer to the graph, as each graph
// manages its set of vertices and we do not want to double-free the vertices
void copyVertexToSubgraph(StringGraph* pSubgraph, const Vertex* pVertex)
{
    // Make sure the vertex hasn't been added yet
    if(pSubgraph->getVertex(pVertex->getID()) == NULL)
    {
        Vertex* pCopy = new(pSubgraph->getVertexAllocator()) Vertex(pVertex->getID(), pVertex->getSeq().toString());
        pSubgraph->addVertex(pCopy);
    }
}

// 
// Handle command line arguments
//
void parseSubgraphOptions(int argc, char** argv)
{
    // Set defaults
    opt::outFile = "subgraph.asqg.gz";

    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'o': arg >> opt::outFile; break;
            case 's': arg >> opt::span; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << SUBGRAPH_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << SUBGRAPH_VERSION_MESSAGE;
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
        std::cout << "\n" << SUBGRAPH_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the vertex id
    opt::rootID = argv[optind++];

    if(opt::rootID.empty())
    {
        std::cerr << SUBPROGRAM ": missing root ID\n";
        exit(EXIT_FAILURE);
    }

    // Parse the input filename
    opt::asqgFile = argv[optind++];
}
