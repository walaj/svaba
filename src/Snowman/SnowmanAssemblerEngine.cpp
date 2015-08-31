#include "SnowmanAssemblerEngine.h"

#include "SnowmanAssemble.h"
#include "SnowmanOverlapAlgorithm.h"
#include "SnowmanASQG.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "CorrectionThresholds.h"
#include <map>

#define MAX_OVERLAPS_PER_ASSEMBLY 20000
//#define DEBUG_ENGINE 1
#define NUM_ROUNDS 3

static std::string POLYA = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
static std::string POLYT = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
static std::string POLYC = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
static std::string POLYG = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
static std::string POLYAT = "ATATATATATATATATATATATATATATATATATATATAT";
static std::string POLYTC = "TCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTC";
static std::string POLYAG = "AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG";
static std::string POLYCG = "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG";
static std::string POLYTG = "TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG";
static std::string POLYCA = "CACACACACACACACACACACACACACACACACACACACA";

void SnowmanAssemblerEngine::fillReadTable(const std::vector<std::string>& r) {

  int count = 0;
  for (auto& i : r) {
    
    SeqItem si;
    
    assert(i.length());
    si.id = "read_" + std::to_string(++count); 
    si.seq = i;
    
    m_pRT.addRead(si);
    
  }
  
}

void SnowmanAssemblerEngine::fillReadTable(SnowTools::BamReadVector& r)
{
  
  m_reads = r;

  // make the reads tables
  for (auto& i : r) {

    SeqItem si;
    string sr, seq = "";

    // get the sequence
    int dum = 0;
    sr = i.GetZTag("SR");
    if (!sr.length())
      sr = i.Qname();

    seq = i.GetZTag("KC");
    if (!seq.length()) {
      seq = i.QualityTrimmedSequence(4, dum);
    } 
    assert(sr.length());
    assert(seq.length());

    if (hasRepeat(seq))
      continue;

    si.id = sr;
    si.seq = seq;
    m_pRT.addRead(si);

  }
  
}

bool SnowmanAssemblerEngine::hasRepeat(const std::string& seq) {
  if (seq.length() < 40)
    return false;
  if ((seq.find(POLYT) == std::string::npos) && 
      (seq.find(POLYA) == std::string::npos) && 
      (seq.find(POLYC) == std::string::npos) && 
      (seq.find(POLYG) == std::string::npos) && 
      (seq.find(POLYCG) == std::string::npos) && 
      (seq.find(POLYAT) == std::string::npos) && 
      (seq.find(POLYTC) == std::string::npos) && 
      (seq.find(POLYAG) == std::string::npos) && 
      (seq.find(POLYCA) == std::string::npos) && 
      (seq.find(POLYTG) == std::string::npos) && 
      (seq.find("N") == std::string::npos))
    return false;
  
  return true;

}

bool SnowmanAssemblerEngine::performAssembly() 
{
  if (m_pRT.getCount() < 3)
    return false;

#ifdef DEBUG_ENGINE
  std::cout << "Doing assembly on: " << m_id << " with " << m_pRT.getCount() << " reads" << std::endl; 
#endif
  

  ContigVector contigs0;

#ifdef DEBUG_ENGINE
    std::cout << "...round 1" << std::endl;
#endif    
 
  // do the first round (on raw reads)
  doAssembly(&m_pRT, contigs0, 0);
  
  //m_contigs = contigs0; return true; //debug

  for (size_t yy = 1; yy != (NUM_ROUNDS); yy++) {

    if (contigs0.size() < 2) {
      for (auto& c: contigs0)
	if (c.getSeq().length() >= m_readlen + 20)
	  m_contigs.push_back(c);
      continue;
    }

#ifdef DEBUG_ENGINE
    std::cout << "...round " << (yy+1) << std::endl;
#endif    
    
    // do the second round (on assembled contigs)
    ReadTable pRTc0(contigs0);
    m_contigs.clear();
    doAssembly(&pRTc0, m_contigs, yy);      
    contigs0 = m_contigs;
    
  }

  return true;
}


// call the assembler
void SnowmanAssemblerEngine::doAssembly(ReadTable *pRT, ContigVector &contigs, int pass) {
  
  if (pRT->getCount() == 0)
    return;

  //#ifdef CLOCK_COUNTER
  //clock_t ct1_before = clock();
  //#endif
    
  // forward
  SuffixArray* pSAf = new SuffixArray(pRT, 1, false); //1 is num threads. false is silent/no
  RLBWT *pBWT= new RLBWT(pSAf, pRT);

  // reverse
  pRT->reverseAll();
  SuffixArray * pSAr = new SuffixArray(pRT, 1, false);
  RLBWT *pRBWT = new RLBWT(pSAr, pRT);
  pRT->reverseAll();

  //if (pass == 0) {
    
  //} // end read correction
  
  // rmdup attempt
  SnowmanOverlapAlgorithm* pRmDupOverlapper = new SnowmanOverlapAlgorithm(pBWT, pRBWT, 
                                                       0, 0, 
                                                       0, false);
  pRT->setZero();
  ReadTable * pRT_nd = new ReadTable();
  SeqItem sir;
  while (pRT->getRead(sir)) {
    OverlapBlockList OBout;
    SeqRecord read;
    read.id = sir.id;
    read.seq = sir.seq;
    OverlapBlockList obl;
    OverlapResult rr = pRmDupOverlapper->alignReadDuplicate(read, &OBout);

    if (!rr.isSubstring)
      pRT_nd->addRead(sir);
    //else
    //  std::cerr << std::endl << "Read " << sir.id << " is substring with seq " << sir.seq.toString() << std::endl << std::endl;;
  }

  // forward
  SuffixArray* pSAf_nd = new SuffixArray(pRT_nd, 1, false); //1 is num threads. false is silent/no
  RLBWT *pBWT_nd = new RLBWT(pSAf_nd, pRT_nd);

  // reverse
  pRT_nd->reverseAll();
  SuffixArray * pSAr_nd = new SuffixArray(pRT_nd, 1, false);
  RLBWT *pRBWT_nd = new RLBWT(pSAr_nd, pRT_nd);
  pRT_nd->reverseAll();

  //delete pRT;
  //delete pRmDupOverlapper;
  //pRT = pRT_nodup;
  
 
  /*
  pRT->setZero();
  SeqItem si2;
  ReadTable * pRT_nodup = new ReadTable();
  while(pRT->getRead(si2)) {
    if (!findOverlapBlocksExactSnow(si2.seq.toString(), pBWT, pRBWT))
      pRT_nodup->addRead(si2);
  }
  if (pRT_nodup->getCount() == 0)
    return;

  // forward
  SuffixArray* pSAf_rd = new SuffixArray(pRT_nodup, 1, false); //1 is num threads. false is silent/no
  RLBWT *pBWT_rd= new RLBWT(pSAf_rd, pRT_nodup);

  // reverse
  pRT->reverseAll();
  SuffixArray * pSAr_rd = new SuffixArray(pRT_nodup, 1, false);
  RLBWT *pRBWT_rd = new RLBWT(pSAr_rd, pRT_nodup);
  pRT_nodup->reverseAll();
  */

  pSAf_nd->writeIndex();
  pSAr_nd->writeIndex();

  //#ifdef CLOCK_COUNTER  
  //(*clock_counter) << pass << "," << "suffix" << "," << (clock() - ct1_before)<< "," << pRT->getCount() << std::endl;
  //#endif

  double errorRate = m_error_rate;
  int min_overlap = m_min_overlap;

  int cutoff = 0;
  if (pass == 0)
    //cutoff = m_readlen + 10;
    cutoff = 0; // debug
  if (pass > 0) {
    min_overlap = 35;
    errorRate = 0.05;
    cutoff = m_readlen + 20;
  }

  int seedLength = min_overlap;
  int seedStride = seedLength;
  seedLength = 20; seedStride = 20; // debug
  bool bIrreducibleOnly = true; // default

  //std::cerr << "min overlap " << min_overlap << " seed length " << seedLength << std::endl;
  
  SnowmanOverlapAlgorithm* pOverlapper = new SnowmanOverlapAlgorithm(pBWT_nd, pRBWT_nd, 
                                                       errorRate, seedLength, 
                                                       seedStride, bIrreducibleOnly);

  bool exact = errorRate < 0.001f;
  //pOverlapper->setExactModeOverlap(opt::assemb::error_rate < 0.001f/*false*/);
  //pOverlapper->setExactModeIrreducible(opt::assemb::error_rate < 0.001f/*false*/);
  pOverlapper->setExactModeOverlap(exact);
  pOverlapper->setExactModeIrreducible(exact);

  stringstream hits_stream;
  stringstream asqg_stream;

  SnowmanASQG::HeaderRecord headerRecord;
  headerRecord.setOverlapTag(min_overlap);
  headerRecord.setErrorRateTag(errorRate);
  headerRecord.setInputFileTag("");
  headerRecord.setContainmentTag(true); // containments are always present
  headerRecord.setTransitiveTag(!bIrreducibleOnly);
  headerRecord.write(asqg_stream);    

  pRT_nd->setZero();

  size_t workid = 0;
  SeqItem si;
  //#ifdef CLOCK_COUNTER
  //clock_t ct_before = clock();
  //#endif

  size_t ocount = 0;
  while (pRT_nd->getRead(si) && (++ocount < MAX_OVERLAPS_PER_ASSEMBLY)) {

    SeqRecord read;
    read.id = si.id;
    read.seq = si.seq;
    OverlapBlockList obl;
    
    OverlapResult rr = pOverlapper->overlapRead(read, min_overlap, &obl);
    
    pOverlapper->writeOverlapBlocks(hits_stream, workid, rr.isSubstring, &obl);

    SnowmanASQG::VertexRecord record(read.id, read.seq.toString());
    record.setSubstringTag(rr.isSubstring);
    record.write(asqg_stream);

    workid++;

  }
  //#ifdef CLOCK_COUNTER  
  //(*clock_counter) << pass << "," <<"overlaps" << "," << (clock() - ct_before)  << "," << pRT->getCount() << std::endl;
  //#endif

  string line;
  bool bIsSelfCompare = true;
  ReadInfoTable* pQueryRIT = new ReadInfoTable(pRT_nd);

  while(getline(hits_stream, line)) {
    size_t readIdx;
    size_t totalEntries;
    bool isSubstring; 
    OverlapVector ov;
    OverlapCommon::parseHitsString(line, pQueryRIT, pQueryRIT, pSAf_nd, pSAr_nd, bIsSelfCompare, readIdx, totalEntries, ov, isSubstring);
    for(OverlapVector::iterator iter = ov.begin(); iter != ov.end(); ++iter)
    {
       SnowmanASQG::EdgeRecord edgeRecord(*iter);
       edgeRecord.write(asqg_stream);
    }

  }
  
  // optionally output the graph structure
  if (m_write_asqg) {

    // write ASQG to file for visualization
    std::stringstream asqgfile;
    asqgfile << m_id << "pass_" << pass << ".asqg";
    std::ofstream ofile(asqgfile.str(), std::ios::out);
    ofile << asqg_stream.str();
    ofile.close();

    // write the hits stream file
    std::stringstream hitsfile;
    hitsfile << m_id << "pass_" << pass << ".hits";
    std::ofstream ofile_hits(hitsfile.str(), std::ios::out);
    ofile_hits << hits_stream.str();
    ofile_hits.close();
    
  }
  
  // Get the number of strings in the BWT, this is used to pre-allocated the read table
  delete pOverlapper;
  delete pBWT; 
  delete pRBWT;
  delete pSAf;
  delete pSAr;

  delete pRmDupOverlapper;
  delete pBWT_nd; 
  delete pRBWT_nd;
  delete pSAf_nd;
  delete pSAr_nd;
  delete pRT_nd;
  
  //#ifdef CLOCK_COUNTER
  //clock_t ctA_before = clock();
  //#endif
  
  // PERFORM THE ASSMEBLY
  assemble(asqg_stream, min_overlap, maxEdges, bExact, 
	   trimLengthThreshold, bPerformTR, bValidate, numTrimRounds, 
	   resolveSmallRepeatLen, numBubbleRounds, gap_divergence, 
	   divergence, maxIndelLength, cutoff, m_id + "_", contigs);

  // remove exact dups
  std::set<std::string> ContigDeDup;
  ContigVector cvec;
  for (auto& i : contigs) {
    if (!ContigDeDup.count(i.getSeq())) {
      ContigDeDup.insert(i.getSeq());
      cvec.push_back(i);
    } else {
      std::cerr << "Filtered out a contig" << std::endl;
    }
  }

  contigs = cvec;
  
  //#ifdef CLOCK_COUNTER  
  //(*clock_counter) << pass << "," << "assembly" << "," << (clock() - ctA_before) << "," << pRT->getCount() << std::endl;
  //#endif
  
  delete pQueryRIT;
  asqg_stream.str("");
  hits_stream.str("");
  
  // print out some results
#ifdef DEBUG_ENGINE
  if (contigs.size() >= 1) {
    std::cout << "Contig Count: " << contigs.size() << " at " << m_id << std::endl;
    //if (opt::verbose > 3)
    for (auto& i : contigs) 
    	std::cout << "   " << i.getID() << " " << i.getSeq().length() << " " << i.getSeq() << std::endl;
  }
#endif
  
  return;
}
