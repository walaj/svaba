//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// connect - Determine the complete sequence of a 
// paired end fragment by finding a walk that
// connects the ends.
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "connect.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "ASQG.h"
#include "gzstream.h"
#include "SequenceProcessFramework.h"
#include "ConnectProcess.h"
#include "SGUtil.h"
#include "gmap.h"
#include "SGSearch.h"
#include "api/BamReader.h"

//#define DEBUG_CONNECT 1

// Structs

// Functions
void markWalkVertices(SGWalk& walk, GraphColor color);
void writeWalk(const std::string& name, int walkIdx, const std::string& fragment, std::ostream* pWriter);

//
// Getopt
//
#define SUBPROGRAM "connect"
static const char *CONNECT_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2010 Wellcome Trust Sanger Institute\n";

static const char *CONNECT_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... ASQGFILE BAMFILE\n"
"Resolve the complete sequence of a paired end fragment by finding a walk through the graph connecting the ends.\n"
"The graph is read from the ASQGFILE. The reads alignments (to the contigs that make up the graph graph) are\n"
"read from the BAMFILE.\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
//"      -t, --threads=NUM                use NUM threads to compute the overlaps (default: 1)\n"
"      -l, --min-distance=LEN           minimum expected distance between the PE reads (start to end). Default: 150.\n"
"      -m, --max-distance=LEN           maximum expected distance between the PE reads (start to end). This option specifies\n"
"                                       how long the search should proceed for. Default: 250\n"
"      -o, --outfile=FILE               write the connected reads to FILE\n"
"          --connected-only             only write the connected read pairs, not the unresolved vertices in the graph\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;

    static int minDistance = 150;
    static int maxDistance = 500;
    
    static bool bWriteCovered = true;
    static bool bWriteUnresolved = false;

    // In hetSVMode, we will output up to two paths between the pairs
    static bool hetSVMode = false;

    static std::string outFile;
    static std::string unconnectedFile = "unconnected.fa";
    static std::string asqgFile;
    static std::string bamFile;
}

static const char* shortopts = "p:m:e:t:l:s:o:d:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_METRICS, OPT_HETSV, OPT_CONNECTED_ONLY, };

static const struct option longopts[] = {
    { "verbose",         no_argument,       NULL, 'v' },
    { "threads",         required_argument, NULL, 't' },
    { "max-distance",    required_argument, NULL, 'd' },
    { "outfile",         required_argument, NULL, 'o' },
    { "min-distance",    required_argument, NULL, 'l' },
    { "connected-only",  no_argument,       NULL, OPT_CONNECTED_ONLY },
    { "het-sv",          no_argument,       NULL, OPT_HETSV },
    { "help",            no_argument,       NULL, OPT_HELP },
    { "version",         no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

// Visit every vertex in the graph and write out
// the vertices that were covered by a path but not
// resolved
struct CoveredVertexVisitor
{

    CoveredVertexVisitor(std::ostream* pWriter) : m_pWriter(pWriter), m_notUsedWriter("notused.fa"), resolved(0), internal(0), not_used(0) {}
    ~CoveredVertexVisitor()
    {
        std::cout << "Resolved: " << resolved << "\n";
        std::cout << "Internal: " << internal << "\n";
        std::cout << "Not used: " << not_used << "\n";
    }

    // not used
    void previsit(StringGraph*) {}
    void postvisit(StringGraph*) {}

    //
    bool visit(StringGraph*, Vertex* pX)
    {
        if(pX->getColor() == GC_BLACK)
        {
            resolved += 1;
        }

        if(pX->getColor() == GC_WHITE)
        {
            not_used += 1;
            SeqRecord record;
            record.id = pX->getID();
            record.seq = pX->getSeq().toString();
            record.write(m_notUsedWriter);
        }

        if(pX->getColor() == GC_RED)
        {
            SeqRecord record;
            record.id = pX->getID();
            record.seq = pX->getSeq().toString();
            record.write(*m_pWriter);
            internal += 1;
        }

        return false;
    }

    int getNumWrote() const { return internal; }

    //
    std::ostream* m_pWriter;
    std::ofstream m_notUsedWriter;
    int resolved;
    int internal;
    int not_used;

};

//
// Main
//
int connectMain(int argc, char** argv)
{
    parseConnectOptions(argc, argv);

    // Read the graph and compute walks
    StringGraph* pGraph = SGUtil::loadASQG(opt::asqgFile, 0, false);

    Timer* pTimer = new Timer(PROGRAM_IDENT);    

    // Open the bam file for reading
    BamTools::BamReader* pBamReader = new BamTools::BamReader;
    pBamReader->Open(opt::bamFile);

    std::ostream* pWriter = createWriter(opt::outFile);
    std::cout << "NUM REFS: " << pBamReader->GetReferenceCount() << "\n";
    
    int numPairsAttempted = 0;
    int numPairsResolved = 0;
    int numUnresolvedWrote = 0;
    int numCoveredWrote = 0;
    int numFailedNoPath = 0;
    int numFailedMultiPaths = 0;
    int numFailedUnaligned = 0;
    int numPathsRejectLow = 0;
    int numPathsRejectHigh = 0;
    int numPathsRejectOrientation = 0;    

    // In heterozygous SV mode, write up to 2 paths
    size_t maxPaths = (opt::hetSVMode ? 2 : 1);

    BamTools::BamAlignment record1;
    BamTools::BamAlignment record2;
    bool done = false;

    const BamTools::RefVector& referenceVector = pBamReader->GetReferenceData();
    while(!done)
    {
        // Read a pair from the BAM
        // Read record 1. Skip secondary alignments of the previous pair
        do
        {
            if(!pBamReader->GetNextAlignment(record1))
            {
                done = true;
                break;
            }
        } while(!record1.IsPrimaryAlignment());

        // Stop if no alignment could be parsed from the stream
        if(done)
            break;

        // Read record 2. Skip any 
        do
        {
            if(!pBamReader->GetNextAlignment(record2))
            {
                done = true;
                break;
            }
        } while(!record2.IsPrimaryAlignment());

        // If this read failed, there is a mismatch between the pairing
        if(done)
        {
            std::cout << "Could not read pair for read: " << record1.Name << "\n";
        }
    
        if(!record1.IsMapped() || !record2.IsMapped())
        {
            numFailedUnaligned += 1;
            continue;
        }

        // Ensure the pairing is correct
        assert(record1.Name == record2.Name);
        
        std::string vertexID1 = referenceVector[record1.RefID].RefName;
        std::string vertexID2 = referenceVector[record2.RefID].RefName;

        // Get the vertices for this pair using the mapped IDs
        Vertex* pX = pGraph->getVertex(vertexID1);
        Vertex* pY = pGraph->getVertex(vertexID2);


        // Ensure that the vertices are found
        assert(pX != NULL && pY != NULL);

#ifdef DEBUG_CONNECT
        std::cout << "Finding path from " << vertexID1 << " to " << vertexID2 << "\n";
#endif

        EdgeDir walkDirectionXOut = ED_SENSE;
        EdgeDir walkDirectionYIn = ED_SENSE;

        // Flip walk directions if the alignment is to the reverse strand
        if(record1.IsReverseStrand())
            walkDirectionXOut = !walkDirectionXOut;
        
        if(record2.IsReverseStrand())
            walkDirectionYIn = !walkDirectionYIn;

        int fromX = walkDirectionXOut == ED_SENSE ? record1.Position : record1.GetEndPosition();
        int toY = walkDirectionYIn == ED_SENSE ? record2.Position : record2.GetEndPosition();

        // Calculate the amount of contig X that already covers the fragment
        // Using this number, we calculate how far we should search
        int coveredX = walkDirectionXOut == ED_SENSE ? pX->getSeqLen() - fromX : fromX;
        int maxWalkDistance = opt::maxDistance - coveredX;

        SGWalkVector walks;
        SGSearch::findWalks(pX, pY, walkDirectionXOut, maxWalkDistance, 10000, true, walks);

        // Mark used vertices in the graph
        // If the entire path was resolved, mark black
        // otherwise mark as red
        GraphColor used_color = (walks.size() <= maxPaths) ? GC_BLACK : GC_RED;

        for(size_t i = 0; i < walks.size(); i +=1 )
            markWalkVertices(walks[i], used_color);

        if(!walks.empty() && walks.size() <= maxPaths)
        {
            for(size_t i = 0; i < walks.size(); ++i)
            {
                // Validate that the path is as expected
                // This has 2 conditions:
                // 1) The inferred fragment is orientated correctly
                // 2) The fragment size is within the expected range
                bool correctOrientation = true;
                WARN_ONCE("check orientation of result");

                std::string fragment = walks[i].getFragmentString(pX, 
                                                                  pY, 
                                                                  fromX,
                                                                  toY,
                                                                  walkDirectionXOut,
                                                                  walkDirectionYIn);

                // Calculate the seqcoord on the path string representing the paired end fragment
                int fragSize = fragment.length();
                bool correctSize = !fragment.empty();

                if(fragSize < opt::minDistance)
                {
                    correctSize = false;
                    numPathsRejectLow += 1;
                }

                if(fragSize > opt::maxDistance)
                {
                    correctSize = false;
                    numPathsRejectHigh += 1;
                }


                if(correctOrientation && correctSize)
                {                    
                    writeWalk(getPairBasename(record1.Name), i, fragment, pWriter);
                    
                    // Mark all the vertices in this walk as resolved
                    markWalkVertices(walks[i], GC_BLACK);
                    numPairsResolved += 1;
                }
            }
        }
        else
        {
            if(walks.empty())
            {
                /*
                std::cout << "No walk found for pair: \n";
                std::cout << record1.Name << " " << vertexID1 << " -> " << vertexID2 << "\n";
                std::cout << "Map qual: " << record1.MapQuality << " " << record2.MapQuality << "\n";
                int xa1 = 0;
                int xa2 = 0;
                record1.GetTag("X1", xa1);
                record2.GetTag("X1", xa2);
                std::cout << "Num sub1: " << xa1 << "\n";
                std::cout << "Num sub2: " << xa2 << "\n";
                std::cout << "Total alts: " << (xa1 + xa2) << "\n";
                */
                numFailedNoPath += 1;
            }
            else if(walks.size() > maxPaths)
            {
                /*
                for(size_t i = 0; i < walks.size(); ++i)
                {
                    std::string fragment = walks[i].getFragmentString(pX, 
                                                                      pY, 
                                                                      fromX,
                                                                      toY,
                                                                      walkDirectionXOut,
                                                                      walkDirectionYIn);

                    //std::cout << i << " " << walks[i].pathSignature() << " size: " << fragment.size() << "\n";
                    //std::cout << ">" << i << "\n" << fragment << "\n";
                }
                */

                numFailedMultiPaths += 1;
            }
            if(opt::bWriteUnresolved)
            {
                // Write the unconnected reads
                SeqRecord unresolved1;
                unresolved1.id = record1.Name;
                unresolved1.seq = record1.QueryBases;
                
                SeqRecord unresolved2;
                unresolved2.id = record2.Name;
                unresolved2.seq = record2.QueryBases;

                unresolved1.write(*pWriter);
                unresolved2.write(*pWriter);

                numUnresolvedWrote += 2;
            }
        }
        numPairsAttempted += 1;
        
        if(numPairsAttempted % 50000 == 0)
            printf("[sga connect] Processed %d pairs\n", numPairsAttempted);
    }

    //
    if(opt::bWriteCovered)
    {
        CoveredVertexVisitor cvv(pWriter);
        pGraph->visit(cvv);
        numCoveredWrote += cvv.getNumWrote();
    }

    double proc_time_secs = pTimer->getElapsedWallTime();
    printf("connect: Resolved %d out of %d pairs (%lf) in %lfs (%lf pairs/s)\n",
            numPairsResolved, numPairsAttempted, 
            (double)numPairsResolved / numPairsAttempted,
            proc_time_secs,
            numPairsAttempted / proc_time_secs);

    printf("Num failed due to no valid path: %d\n", numFailedNoPath);
    printf("Num failed due to multiple valid paths: %d\n", numFailedMultiPaths);
    printf("Num failed due to part of the pair not aligning to the graph: %d\n", numFailedUnaligned);
    printf("Num paths rejected because they are shorter than the minimum distance: %d\n", numPathsRejectLow);
    printf("Num paths rejected because they are longer than the maximum distance: %d\n", numPathsRejectHigh);
    printf("Num paths rejected because they are do not have the correct orientation: %d\n", numPathsRejectOrientation);

    printf("Wrote %d unconnected pairs\n", numUnresolvedWrote);
    printf("Wrote %d vertices that were covered by a path but not full resolved\n", numCoveredWrote);

    delete pTimer;
    delete pGraph;
    delete pWriter;
    delete pBamReader;
    return 0;
}

// Write the given walk out to the file
void writeWalk(const std::string& name, int walkIdx, const std::string& fragment, std::ostream* pWriter)
{
    std::stringstream idSS;
    idSS << name;

    if(opt::hetSVMode)
        idSS << "-walk:" << walkIdx;

    SeqRecord resolved;
    resolved.id = idSS.str();
    resolved.seq = fragment;
    resolved.write(*pWriter);
}

// Read all the alignments for 
bool readPairAlignments()
{
    return false;
}

// Mark all the vertices in the walk as color 
// unless they are already black (we never change from black->red)
void markWalkVertices(SGWalk& walk, GraphColor color)
{
    VertexPtrVec verts = walk.getVertices();
    for(size_t i = 0; i < verts.size(); ++i)
    {
        if(verts[i]->getColor() != GC_BLACK)
            verts[i]->setColor(color);
    }
}

// 
// Handle command line arguments
//
void parseConnectOptions(int argc, char** argv)
{
    std::string algo_str;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'o': arg >> opt::outFile; break;
            case 'd': arg >> opt::maxDistance; break;
            case 't': arg >> opt::numThreads; break;
            case 'l': arg >> opt::minDistance; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_CONNECTED_ONLY: 
                opt::bWriteCovered = false; 
                opt::bWriteUnresolved = false; 
                break;
            case OPT_HETSV: opt::hetSVMode = true; break;
            case OPT_HELP:
                std::cout << CONNECT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << CONNECT_VERSION_MESSAGE;
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

    if(opt::numThreads <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << CONNECT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::asqgFile = argv[optind++];
    opt::bamFile = argv[optind++];

    if(opt::outFile.empty())
    {
        std::string prefix = stripFilename(opt::bamFile);
        opt::outFile = prefix + ".connect.fa";
        opt::unconnectedFile = prefix + ".single.fa";
    }
}
