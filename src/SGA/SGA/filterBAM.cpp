//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// filterBAM - remove erroneous paired end connections
// from a bam file.
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "filterBAM.h"
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
#include "api/BamWriter.h"

// Structs

// Functions
bool filterByGraph(StringGraph* pGraph, 
                   const BamTools::RefVector& referenceVector, 
                   BamTools::BamAlignment& record1, 
                   BamTools::BamAlignment& record2);

double getErrorRate(BamTools::BamAlignment& record);

bool readAlignmentPair(BamTools::BamReader* pReader, 
                       BamTools::BamAlignment& record1,
                       BamTools::BamAlignment& record2);

int64_t getMaxKmerDepth(const std::string& w, const BWT* pBWT, const BWT* pRBWT);

//
// Getopt
//
#define SUBPROGRAM "filterBAM"
static const char *FILTERBAM_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2011 Wellcome Trust Sanger Institute\n";

static const char *FILTERBAM_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... ASQGFILE BAMFILE\n"
"Discard mate-pair alignments from a BAM file that are potentially erroneous\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -a, --asqg=FILE                  load an asqg file and filter pairs that are shorter than --max-distance\n"
"      -d, --max-distance=LEN           search the graph for a path completing the mate-pair fragment. If the path is less than LEN\n"
"                                       then the pair will be discarded.\n"
"      -e, --error-rate=F               filter out pairs where one read has an error rate higher than F (default: no filter)\n"
"      -q, --min-quality=F              filter out pairs where one read has mapping quality less than F (default: 10)\n"
"      -o, --out-bam=FILE               write the filtered reads to FILE\n"
"      -p, --prefix=STR                 load the FM-index with prefix STR\n"
"      -x, --max-kmer-depth=N           filter out pairs that contain a kmer that has been seen in the FM-index more than N times\n"
"      -c, --mate-contamination         filter out pairs aligning with FR orientation, which may be contiminates in a mate pair library\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;

    static int maxDistance = 700;
    static int minQuality = 10;
    static double maxError = 1.0f;

    static std::string outFile;
    static std::string asqgFile;
    static std::string bamFile;
    static std::string fmIndexPrefix;

    bool filterFRContamination = true;
    int minDistanceToEnd = 500;

    static int kmerSize = 31;
    static int maxKmerDepth = -1;
    static int sampleRate = 256;
}

static const char* shortopts = "d:t:o:q:e:a:p:x:t:c:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",            no_argument,       NULL, 'v' },
    { "threads",            required_argument, NULL, 't' },
    { "asqg-file",          required_argument, NULL, 'a' },
    { "max-distance",       required_argument, NULL, 'd' },
    { "error-rate",         required_argument, NULL, 'e' },
    { "min-quality",        required_argument, NULL, 'q' },
    { "outfile",            required_argument, NULL, 'o' },
    { "fmIndexPrefix",      required_argument, NULL, 'p' },
    { "end-distance",       required_argument, NULL, 't' },
    { "mate-contamination", required_argument, NULL, 'c' },
    { "help",               no_argument,       NULL, OPT_HELP },
    { "version",            no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int filterBAMMain(int argc, char** argv)
{
    parseFilterBAMOptions(argc, argv);

    // Read the graph if distance-filtering mode is enabled
    StringGraph* pGraph = NULL;
    if(!opt::asqgFile.empty())
        pGraph = SGUtil::loadASQG(opt::asqgFile, 0, false);

    // Read the BWTs if depth-filtering mode is enabled
    BWT* pBWT = NULL;
    BWT* pRBWT = NULL;
    if(!opt::fmIndexPrefix.empty())
    {
        pBWT = new BWT(opt::fmIndexPrefix + BWT_EXT, opt::sampleRate);
        pRBWT = new BWT(opt::fmIndexPrefix + RBWT_EXT, opt::sampleRate);
    }

    Timer* pTimer = new Timer(PROGRAM_IDENT);    

    // 
    int numPairsTotal = 0;
    int numPairsFilteredByDistance = 0;
    int numPairsFilteredByER = 0;
    int numPairsFilteredByQuality = 0;
    int numPairsFilteredByDepth = 0;
    int numPairsUnmapped = 0;
    int numPairsWrote = 0;
    int numPairsFilteredFRContamination = 0;
    int numPairsTooCloseToEnd = 0;

    // Open the bam files for reading/writing
    BamTools::BamReader* pBamReader = new BamTools::BamReader;
    pBamReader->Open(opt::bamFile);

    BamTools::BamWriter* pBamWriter = new BamTools::BamWriter;
    pBamWriter->Open(opt::outFile, pBamReader->GetHeaderText(), pBamReader->GetReferenceData());
    const BamTools::RefVector& referenceVector = pBamReader->GetReferenceData();


    BamTools::BamAlignment record1;
    BamTools::BamAlignment record2;
    bool done = false;

    while(!done)
    {
        if(numPairsTotal++ % 200000 == 0)
            printf("[sga filterBAM] Processed %d pairs\n", numPairsTotal);

        done = !readAlignmentPair(pBamReader, record1, record2);
        if(done)
            break;

        if(!record1.IsMapped() || !record2.IsMapped())
        {
            numPairsUnmapped += 1;
            continue;
        }

        // Ensure the pairing is correct
        if(record1.Name != record2.Name)
            std::cout << "NAME FAIL: " << record1.Name << " " << record2.Name << "\n";
        assert(record1.Name == record2.Name);
        bool bPassedFilters = true;

        // Check if the error rate is below the max
        double er1 = getErrorRate(record1);
        double er2 = getErrorRate(record2);

        if(er1 > opt::maxError || er2 > opt::maxError)
        {
            bPassedFilters = false;
            numPairsFilteredByER += 1;
        }

        if(record1.MapQuality < opt::minQuality || record2.MapQuality < opt::minQuality)
        {
            bPassedFilters = false;
            numPairsFilteredByQuality += 1;
        }

        // Perform depth check for pairs aligning to different contigs
        if(bPassedFilters && (pBWT != NULL && pRBWT != NULL && opt::maxKmerDepth > 0) && (record1.RefID != record2.RefID))
        {
            int maxDepth1 = getMaxKmerDepth(record1.QueryBases, pBWT, pRBWT);
            int maxDepth2 = getMaxKmerDepth(record1.QueryBases, pBWT, pRBWT);
            if(maxDepth1 > opt::maxKmerDepth || maxDepth2 > opt::maxKmerDepth)
            {
                bPassedFilters = false;
                numPairsFilteredByDepth += 1;
            }
        }

        // Filter forward-reverse contimating pairs in a mate pair library
        if(opt::filterFRContamination)
        {
            if(record1.RefID == record2.RefID)
            {
                // Check the orientation of the pairs
                // We discard the pair if they are like this:
                //  ------1---->
                //                <------2------
                BamTools::BamAlignment* pUpstream;
                BamTools::BamAlignment* pDownstream;
                if(record1.Position < record2.Position)
                {
                    pUpstream = &record1;
                    pDownstream = &record2;
                }
                else
                {
                    pUpstream = &record2;
                    pDownstream = &record1;
                }
                
                // Upstream half of the pair (more 5') should be forward, downstream should be reverse
                if(!pUpstream->IsReverseStrand() && pDownstream->IsReverseStrand())
                {
                    numPairsFilteredFRContamination += 1;
                    bPassedFilters = false;
                }
            }

            if(bPassedFilters && record1.RefID != record2.RefID)
            {
                int distanceToLeftEnd1 = record1.Position;
                int distanceToRightEnd1 = referenceVector[record1.RefID].RefLength - record1.GetEndPosition();
                int distance1 = std::min(distanceToLeftEnd1, distanceToRightEnd1);
                
                int distanceToLeftEnd2 = record2.Position;
                int distanceToRightEnd2 = referenceVector[record2.RefID].RefLength - record2.GetEndPosition();
                int distance2 = std::min(distanceToLeftEnd2, distanceToRightEnd2);
                if(distance1 < opt::minDistanceToEnd || distance2 < opt::minDistanceToEnd)
                {
                    bPassedFilters = false;
                    numPairsTooCloseToEnd += 1;
                }
            }
        }

        // Perform short-insert pair check
        if(pGraph != NULL)
        {
            bPassedFilters = bPassedFilters && filterByGraph(pGraph, referenceVector, record1, record2);
            numPairsFilteredByDistance += 1;
        }

        if(bPassedFilters)
        {
            pBamWriter->SaveAlignment(record1);
            pBamWriter->SaveAlignment(record2);
            numPairsWrote += 1;
        }
    }

    std::cout << "Total pairs: " << numPairsTotal << "\n";
    std::cout << "Total pairs output: " << numPairsWrote << "\n";
    std::cout << "Total filtered because one pair is unmapped: " << numPairsUnmapped << "\n";
    std::cout << "Total filtered by distance: " << numPairsFilteredByDistance << "\n";
    std::cout << "Total filtered by error rate: " << numPairsFilteredByER << "\n";
    std::cout << "Total filtered by quality: " << numPairsFilteredByQuality << "\n";
    std::cout << "Total filtered by depth: " << numPairsFilteredByDepth << "\n";
    std::cout << "Total filtered by FR orientation: " << numPairsFilteredFRContamination << "\n";
    std::cout << "Total filtered by alignment too close to contig end: " << numPairsTooCloseToEnd << "\n";
    
    if(pGraph != NULL)
        delete pGraph;

    if(pBWT != NULL)
        delete pBWT;

    if(pRBWT != NULL)
        delete pRBWT;

    pBamWriter->Close();
    pBamReader->Close();

    delete pTimer;
    delete pBamReader;
    delete pBamWriter;
    return 0;
}

// Returns true if the paired reads are a short-insert pair
bool filterByGraph(StringGraph* pGraph, 
                   const BamTools::RefVector& referenceVector, 
                   BamTools::BamAlignment& record1, 
                   BamTools::BamAlignment& record2)
{
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

    bool bShortInsertPair = false;
    if(pX == pY)
    {
        if(abs(record1.InsertSize) < opt::maxDistance)
            bShortInsertPair = true;
    }
    else
    {

        SGWalkVector walks;
        SGSearch::findWalks(pX, pY, walkDirectionXOut, maxWalkDistance, 10000, true, walks);

        if(!walks.empty())
        {
            for(size_t i = 0; i < walks.size(); ++i)
            {
                std::string fragment = walks[i].getFragmentString(pX, 
                                                                  pY, 
                                                                  fromX,
                                                                  toY,
                                                                  walkDirectionXOut,
                                                                  walkDirectionYIn);
                if((int)fragment.size() < opt::maxDistance)
                {
                    bShortInsertPair = true;
                    //std::cout << "Found completing fragment (" << pX->getID() << " -> " << pY->getID() << ": " << fragment.size() << "\n";
                    break;
                }
            }
        }
    }
    
    return bShortInsertPair;
}

// Calculate the error rate between the read and the reference
double getErrorRate(BamTools::BamAlignment& record)
{
    int nm = 0;
    bool hasNM = record.GetTag("NM", nm);
    if(hasNM)
        return (double)nm / record.Length;
    else
        return 0.0f;
}

int64_t getMaxKmerDepth(const std::string& w, 
                        const BWT* pBWT, 
                        const BWT* /*pRBWT*/)
{
    int64_t max = 0;
    int nk = w.size() - opt::kmerSize + 1;
    if(w.find_first_not_of("ACGT") != std::string::npos)
        return 0;

    for(int i = 0; i < nk; ++i)
    {
        std::string kmer = w.substr(i, opt::kmerSize);
        int count = BWTAlgorithms::countSequenceOccurrences(kmer, pBWT);
        if(count > max)
            max = count;
    }
    return max;
}

// Read an alignment pair from the BamReader.
// Returns false if the read fails
bool readAlignmentPair(BamTools::BamReader* pReader, 
                       BamTools::BamAlignment& record1,
                       BamTools::BamAlignment& record2)
{
    // Read a pair from the BAM
    // Read record 1. Skip secondary alignments of the previous pair
    do
    {
        if(!pReader->GetNextAlignment(record1))
            return false;
    } while(!record1.IsPrimaryAlignment());

    // Read record 2.
    do
    {
        if(!pReader->GetNextAlignment(record2))
            return false;
    } while(!record2.IsPrimaryAlignment());
    return true;
}

// 
// Handle command line arguments
//
void parseFilterBAMOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'o': arg >> opt::outFile; break;
            case 'd': arg >> opt::maxDistance; break;
            case 't': arg >> opt::numThreads; break;
            case 'e': arg >> opt::maxError; break;
            case 'q': arg >> opt::minQuality; break;
            case 'a': arg >> opt::asqgFile; break;
            case 'p': arg >> opt::fmIndexPrefix; break;
            case 'x': arg >> opt::maxKmerDepth; break;
            case 'c': opt::filterFRContamination = true; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << FILTERBAM_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << FILTERBAM_VERSION_MESSAGE;
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

    if(opt::numThreads <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << FILTERBAM_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::bamFile = argv[optind++];
    if(opt::outFile.empty())
        opt::outFile = "filtered." + opt::bamFile;
}
