//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// assemble - convert read overlaps into contigs
//
#include <iostream>
#include <fstream>
#include "Util.h"
#include "assemble.h"
#include "SGUtil.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"
#include "Timer.h"
#include "EncodedString.h"
#include <unistd.h>
#include <string>

//
// Getopt
//
#define SUBPROGRAM "assemble"
static const char *ASSEMBLE_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *ASSEMBLE_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... ASQGFILE\n"
"Create contigs from the assembly graph ASQGFILE.\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -o, --out-prefix=NAME            use NAME as the prefix of the output files (output files will be NAME-contigs.fa, etc)\n"
"      -m, --min-overlap=LEN            only use overlaps of at least LEN. This can be used to filter\n"
"                                       the overlap set so that the overlap step only needs to be run once.\n"
"          --transitive-reduction       remove transitive edges from the graph. Off by default.\n"
"          --max-edges=N                limit each vertex to a maximum of N edges. For highly repetitive regions\n"
"                                       this helps save memory by culling excessive edges around unresolvable repeats (default: 128)\n"
"\nBubble/Variation removal parameters:\n"
"      -b, --bubble=N                   perform N bubble removal steps (default: 3)\n"
"      -d, --max-divergence=F           only remove variation if the divergence between sequences is less than F (default: 0.05)\n"
"      -g, --max-gap-divergence=F       only remove variation if the divergence between sequences when only counting indels is less than F (default: 0.01)\n"
"                                       Setting this to 0.0 will suppress removing indel variation\n"
"          --max-indel=D                do not remove variation that is an indel of length greater than D (default: 20)\n"
"\n"
"\nTrimming parameters:\n"
"      -x, --cut-terminal=N             cut off terminal branches in N rounds (default: 10)\n"
"      -l, --min-branch-length=LEN      remove terminal branches only if they are less than LEN bases in length (default: 150)\n"

"\nSmall repeat resolution parameters:\n"
"      -r,--resolve-small=LEN           resolve small repeats using spanning overlaps when the difference between the shortest\n"
"                                       and longest overlap is greater than LEN (default: not performed)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string asqgFile;
    static std::string outContigsFile;
    static std::string outVariantsFile;
    static std::string outGraphFile;

    static unsigned int minOverlap;
    static bool bEdgeStats = false;
    static bool bSmoothGraph = false;
    static int resolveSmallRepeatLen = -1;
    static size_t maxEdges = 128;

    // Trim parameters
    static int numTrimRounds = 10;
    static size_t trimLengthThreshold = 300;
    
    // Bubble parameters
    static int numBubbleRounds = 3;
    static double maxBubbleDivergence = 0.05f;
    static double maxBubbleGapDivergence = 0.01f;
    static int maxIndelLength = 20;

    // 
    static bool bValidate;
    static bool bExact = true;
    static bool bPerformTR = false;
}

static const char* shortopts = "p:o:m:d:g:b:a:r:x:l:sv";

enum { OPT_HELP = 1, OPT_VERSION, OPT_VALIDATE, OPT_EDGESTATS, OPT_EXACT, OPT_MAXINDEL, OPT_TR, OPT_MAXEDGES };

static const struct option longopts[] = {
    { "verbose",               no_argument,       NULL, 'v' },
    { "out-prefix",            required_argument, NULL, 'o' },
    { "min-overlap",           required_argument, NULL, 'm' },
    { "bubble",                required_argument, NULL, 'b' },
    { "cut-terminal",          required_argument, NULL, 'x' },
    { "min-branch-length",     required_argument, NULL, 'l' },
    { "resolve-small",         required_argument, NULL, 'r' },
    { "max-divergence",        required_argument, NULL, 'd' },
    { "max-gap-divergence",    required_argument, NULL, 'g' },
    { "max-indel",             required_argument, NULL, OPT_MAXINDEL },
    { "max-edges",             required_argument, NULL, OPT_MAXEDGES },
    { "smooth",                no_argument,       NULL, 's' },
    { "transitive-reduction",  no_argument,       NULL, OPT_TR },
    { "edge-stats",            no_argument,       NULL, OPT_EDGESTATS },
    { "exact",                 no_argument,       NULL, OPT_EXACT },
    { "help",                  no_argument,       NULL, OPT_HELP },
    { "version",               no_argument,       NULL, OPT_VERSION },
    { "validate",              no_argument,       NULL, OPT_VALIDATE},
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int assembleMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("sga assemble");
    parseAssembleOptions(argc, argv);
    //assemble();
    delete pTimer;

    return 0;
}

void assemble(std::stringstream& asqg_stream, int minOverlap, int maxEdges, bool bExact, 
	      int trimLengthThreshold, bool bPerformTR, bool bValidate, int numTrimRounds, 
              int resolveSmallRepeatLen, int numBubbleRounds, double maxBubbleGapDivergence, 
              double maxBubbleDivergence, int maxIndelLength, int cutoff, std::string prefix, 
              ContigVector &contigs)
{

  /*std::cerr << "ASSEMBLE  MinOverlap: " << minOverlap << " MaxEdges: " << maxEdges << " bExact: " << bExact << 
      " TrimLengthThreshold: " << trimLengthThreshold << " bPerformTR: " << bPerformTR << " bValidate: " << bValidate << 
      " numTrimRounds: " << numTrimRounds << " resolveSmallRepeatLen: " << resolveSmallRepeatLen << " numBubbleRounds: " << 
      numBubbleRounds << " maxBubbleGapDivergence: " << maxBubbleGapDivergence << " maxBubbleDivergence: " << maxBubbleDivergence <<
    " opt::maxIndelLength: " << maxIndelLength << " cutoff: " << cutoff << " name: " << prefix << std::endl;
  */
  cout << "HERHERHERHEHREHRHERHEHRH" << endl;
  exit(1);
    //Timer t("sga assemble");
    //StringGraph* pGraph = SGUtil::loadASQG(opt::asqgFile, opt::minOverlap, true, opt::maxEdges);
    StringGraph* pGraph = SGUtil::loadASQG(asqg_stream, minOverlap, true, maxEdges);

    if(bExact)
        pGraph->setExactMode(true);
    //pGraph->printMemSize();

    // Visitor functors
    SGTransitiveReductionVisitor trVisit;
    SGGraphStatsVisitor statsVisit;
    SGTrimVisitor trimVisit(trimLengthThreshold);
    SGContainRemoveVisitor containVisit;
    SGValidateStructureVisitor validationVisit;

    // Pre-assembly graph stats
    //std::cout << "[Stats] Input graph:\n";
    //pGraph->visit(statsVisit);    

    // Remove containments from the graph
    //std::cout << "Removing contained vertices from graph\n";
    //std::cout << "Containments: " << pGraph->hasContainment() << std::endl;
    //std::cout << prefix << std::endl;

    //debug
    //pGraph->writeASQG("after_build.asqg");

    while(pGraph->hasContainment())
        pGraph->visit(containVisit);

    //debug
    //pGraph->writeASQG("after_contain.asqg");
    /*std::cout << "After containments" << endl;
      pGraph->visit(statsVisit);    
      pGraph->writeASQG("/home/unix/jwala/tmp.graph.aftercontainments.asqg");
      pGraph->writeASQG("/home/unix/jwala/tmp.graph.aftercontainments.asqg");*/
    //std::cout << asqg_stream.str();
    //    std::cerr << prefix << std::endl;
    /*VertexPtrMap vt = pGraph->getVertexMap();
    for (VertexPtrMap::const_iterator it = vt.begin(); it != vt.end(); it++)
      for (int jt = 0; jt < it->second->getEdges().size(); jt++) 
    	std::cerr << it->second->getEdges()[jt]->getStartID() << " " << std::endl;
    */

    //std::cerr << it->second->getID() << " Num edges: " << it->second->countEdges() << std::endl;
    // Pre-assembly graph stats
    //std::cout << "[Stats] After removing contained vertices:\n";
    //  pGraph->visit(statsVisit);    

    // Remove any extraneous transitive edges that may remain in the graph
    if(bPerformTR)
    {
      //std::cout << "Removing transitive edges\n";
        pGraph->visit(trVisit);
    }

    // Compact together unbranched chains of vertices
    pGraph->simplify();

    //debug
    //pGraph->writeASQG("after_simplify.asqg");
    
    if(bValidate)
    {
      //std::cout << "Validating graph structure\n";
        pGraph->visit(validationVisit);
    }

    // Remove dead-end branches from the graph
    //numTrimRounds = 1;
    if(numTrimRounds > 0)
    {
      //std::cout << "Trimming bad vertices\n"; 
        int numTrims = numTrimRounds;
        while(numTrims-- > 0)
           pGraph->visit(trimVisit);
	std::cout << "\n[Stats] Graph after trimming:\n";
        pGraph->visit(statsVisit);
    }

    // Resolve small repeats
    if(resolveSmallRepeatLen > 0 && false)
    {
        SGSmallRepeatResolveVisitor smallRepeatVisit(resolveSmallRepeatLen);
	// std::cout << "Resolving small repeats\n";

        //int totalSmallRepeatRounds = 0;
        //while(pGraph->visit(smallRepeatVisit))
        //    std::cout << "Finished small repeat resolve round " << totalSmallRepeatRounds++ << "\n";
        
        //std::cout << "\n[Stats] After small repeat resolution:\n";
        //pGraph->visit(statsVisit);
    }

    /*std::cerr << prefix << " NumBubble: " << numBubbleRounds << "\n\n\n\n";
    vt = pGraph->getVertexMap();
    std::cerr << "MapSize: " << vt.size() << std::endl;
    for (VertexPtrMap::const_iterator it = vt.begin(); it != vt.end(); it++)
    	std::cerr << it->second->getSeq() << " " << std::endl;
    */

    // Peform another round of simplification
    //pGraph->simplify();
    //pGraph->writeASQG("/home/unix/jwala/tmp.graph.aftersimp2.asqg");
    //cout << "NumBUBBLGE: " << numBubbleRounds << endl;
    //SGVisitorContig avtmp;
    //pGraph->visit(avtmp);
    //for (ContigVector::const_iterator it = avtmp.m_ct.begin(); it != avtmp.m_ct.end(); it++) {
     // cout << "PRE BUBBLE LEN: " << it->getLength() << endl;
    //if (it->getLength() >= cutoff) {
	//std::string new_name = it->getID() + "_L" + it->getLength();
	//it->setID(new_name);
	//contigs.push_back(*it);
    //}
    //}


    //debug
    SGVisitorContig av_TEST;
    pGraph->visit(av_TEST);
    cout << "checking before bubble " << endl;
    for (auto& i : av_TEST.m_ct) {
      cout << "BEFRE BUBBLEin assembly " << i.getID() << " " << i.getSeq() << " len " << i.getSeq().length() << endl;
    }


    if(numBubbleRounds > 0)
    {
        //std::cout << "\nPerforming variation smoothing\n";
        SGSmoothingVisitor smoothingVisit(opt::outVariantsFile, maxBubbleGapDivergence, maxBubbleDivergence, maxIndelLength);
        int numSmooth = numBubbleRounds;
        while(numSmooth-- > 0)
            pGraph->visit(smoothingVisit);
        pGraph->simplify();
    }

    //debug
    //pGraph->writeASQG("after_bubble.asqg");

    /* std::cerr << prefix << " NumBubble2: " << numBubbleRounds << "\n\n\n\n";
    vt = pGraph->getVertexMap();
    std::cerr << "MapSize: " << vt.size() << std::endl;
    for (VertexPtrMap::const_iterator it = vt.begin(); it != vt.end(); it++)
    	std::cerr << it->second->getSeq() << " " << std::endl;
    */
    //pGraph->writeASQG("/home/unix/jwala/tmp.graph.final.asqg");
    //pGraph->visit(statsVisit);
    //std::string cmd = "sga-asqg2dot.pl /home/unix/jwala/tmp.graph.final.asqg > /home/unix/jwala/tmp.graph.final.dot; dot -Tpdf /home/unix/jwala/tmp.graph.aftersimple.dot -o /home/unix/jwala/tmp.graph.final.pdf";
    //system(cmd.c_str());
    //pGraph->renameVertices("contig-");
    pGraph->renameVertices(prefix);

    //cout << "FINAL" << endl;
    //pGraph->visit(statsVisit);

    //std::vector<std::string> idvec;
    //pGraph->getToKeep(idvec);
    //std::cout << "Len: " << idvec.size() << " Olen: " << olen << std::endl;
    //for (int i = 0; i < idvec.size(); i++)
    //  std::cout << idvec[i] << std::endl;

    // Rename the vertices to have contig IDs instead of read IDs
    //pGraph->renameVertices("contig-");

    //JEREMIAH commented below
    // Write the results
    //SGFastaVisitor av(opt::outContigsFile); 
    //pGraph->visit(av);

    //JEREMIAH
    //SGFastaVisitorSeqRecord av(cutoff); 
    //SGVisitorContig av(cutoff);
    SGVisitorContig av;
    pGraph->visit(av);

    /*   std::cerr << prefix << "ASFASDFSDFSDFSDF\n\n\n\n";
    vt = pGraph->getVertexMap();
    std::cerr << "MapSize: " << vt.size() << std::endl;
    for (VertexPtrMap::const_iterator it = vt.begin(); it != vt.end(); it++)
    	std::cerr << it->second->getID() << " " << std::endl;
    */
 
    //std::cerr << "NUM CONTIGS: " << av.m_ct.size() << std::endl;
    //for (int i = 0; i < av.m_ct.size(); i++)
    // std::cerr << "LEN: " << av.m_ct[i].getLength() << std::endl;

    //debug
    /*
    SGVisitorContig av_TEST;
    pGraph->visit(av_TEST);
    for (auto& i : avTEST.m_ct) {
      cout << "in assembly " << i.getID() << " " << i.getSeq() << " len " << i.getSeq().length() << endl;
      }*/
    
    ContigVector tmp = av.m_ct;
    for (ContigVector::const_iterator it = tmp.begin(); it != tmp.end(); it++) {
      if (it->getLength() >= cutoff) {
	//std::string new_name = it->getID() + "_L" + it->getLength();
	//it->setID(new_name);
	contigs.push_back(*it);
      }
    }

      //for (unsigned i = 0; i < av.m_ct.size(); i++)
      /*	if (av.m_ct[i].getLength() > cutoff) {
	  //append the length string
	  std::string new_name = av.m_ct[i].getID() + "_L" + av.m_ct[i].getLength();
	  av.m_ct[i].setID()
	    contigs.push_back(av.m_ct[i]);
	  //}
	  }*/

    /*    for (int i = 0; i < av.m_ct.size(); i++)
      if (av.m_ct[i].getLength() > cutoff) {
	std::cerr << av.m_ct[i].getID() << std::endl;
	std::cerr << av.m_ct[i].printAlignments();
	}*/
    
    //VertexPtrVec vecr = pGraph->getAllKeepVertices();
    //std::cerr << "ToKeepLen: " << vecr.size() << std::endl;

    /*
    // match contigs to reads
    ReadTable* pRT_C = new ReadTable(srv);
    SuffixArray* pSAf = new SuffixArray(pRT_, 1, false); //1 is num threads. false is isReverse
    RLBWT *pBWT= new RLBWT(pSAf, pRT_C);
    pRT_C->reverseAll();
    SuffixArray * pSAr = new SuffixArray(pRT_C, 1, true); //1 is num threads. false is isReverse
    RLBWT *pRBWT = new RLBWT(pSAr, pRT_C);
    pRT_C->reverseAll();

    pSAf->writeIndex();
    pSAr->writeIndex();
    double errorRate = 0.05;
    int seedLength = 60;
    bool bIrreducibleOnly = true; // default    
    OverlapAlgorithm* pOverlapper = new OverlapAlgorithm(pBWT, pRBWT, 
                                                         errorRate, seedLength, 
                                                         seedStride, bIrreducibleOnly); */



    //pGraph->writeASQG("/home/unix/jwala/testgraph.asqg");
    

    delete pGraph;
   
}

// 
// Handle command line arguments
//
void parseAssembleOptions(int argc, char** argv)
{
    // Set defaults
    opt::minOverlap = 0;
    std::string prefix = "default";
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'o': arg >> prefix; break;
            case 'm': arg >> opt::minOverlap; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case 'l': arg >> opt::trimLengthThreshold; break;
            case 'b': arg >> opt::numBubbleRounds; break;
            case 'd': arg >> opt::maxBubbleDivergence; break;
            case 'g': arg >> opt::maxBubbleGapDivergence; break;
            case 's': opt::bSmoothGraph = true; break;
            case 'x': arg >> opt::numTrimRounds; break;
            case 'r': arg >> opt::resolveSmallRepeatLen; break;
            case OPT_MAXEDGES: arg >> opt::maxEdges; break;
            case OPT_TR: opt::bPerformTR = true; break;
            case OPT_MAXINDEL: arg >> opt::maxIndelLength; break;
            case OPT_EXACT: opt::bExact = true; break;
            case OPT_EDGESTATS: opt::bEdgeStats = true; break;
            case OPT_VALIDATE: opt::bValidate = true; break;
            case OPT_HELP:
                std::cout << ASSEMBLE_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << ASSEMBLE_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
                
        }
    }

    // Build the output names
    opt::outContigsFile = prefix + "-contigs.fa";
    opt::outVariantsFile = prefix + "-variants.fa";
    opt::outGraphFile = prefix + "-graph.asqg.gz";

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
        std::cout << "\n" << ASSEMBLE_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filename
    opt::asqgFile = argv[optind++];
}
