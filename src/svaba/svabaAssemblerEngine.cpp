#include "svabaAssemblerEngine.h"

#include <map>
#include <algorithm>

#include "SGACommon.h"

#include "svabaASQG.h"
#include "svabaAssemble.h"
#include "svabaOverlapAlgorithm.h"

#include "OverlapCommon.h"
#include "CorrectionThresholds.h"

#define MAX_OVERLAPS_PER_ASSEMBLY 20000
//#define DEBUG_ENGINE 1

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

void svabaAssemblerEngine::fillReadTable(const std::vector<std::string>& r) {

  int count = 0;
  for (auto& i : r) {

    if (i.length() < m_min_overlap)
      continue;
    
    SeqItem si;
    
    assert(i.length() && i.length() >= m_min_overlap);
    si.id = "read_" + std::to_string(++count); 
    si.seq = i;
    
    m_pRT.addRead(si);
    
  }
  
}

void svabaAssemblerEngine::fillReadTable(SeqLib::BamRecordVector& r)
{
  
  m_reads = r;

  // make the reads tables
  for (auto& i : r) {

    SeqItem si;
    std::string sr, seq; 

    // get the sequence
    sr = i.GetZTag("SR");
    if (!sr.length())
      sr = i.Qname();

    seq = i.GetZTag("KC");
    if (!seq.length()) {
      seq = i.QualitySequence(); //i.QualityTrimmedSequence(4, dum);
    } 
    assert(sr.length());
    assert(seq.length());

    if (hasRepeat(seq) || seq.length() < m_min_overlap)
      continue;
      
    si.id = sr;

    // put onto the foward strand if not
    if (!i.MappedFlag() && !i.MateReverseFlag())
      SeqLib::rcomplement(seq);

    si.seq = seq;
    m_pRT.addRead(si);

  }
  
}

bool svabaAssemblerEngine::hasRepeat(const std::string& seq) {

  if (seq.find("N") != std::string::npos)
    return true;
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

bool svabaAssemblerEngine::performAssembly(int num_assembly_rounds) 
{
  if (m_pRT.getCount() < 2)
    return false;

#ifdef DEBUG_ENGINE
  std::cout << "Doing assembly on: " << m_id << " with " << m_pRT.getCount() << " reads" << std::endl; 
#endif

  SeqLib::UnalignedSequenceVector contigs0;
  //ContigVector contigs0;

#ifdef DEBUG_ENGINE
  std::cout << "...round 1" << std::endl;
#endif    
  
  // do the first round (on raw reads)
  doAssembly(&m_pRT, m_contigs, 0);
  
  for (int yy = 1; yy != num_assembly_rounds; yy++) {
    
    if (m_contigs.size() < 2) 
      continue; // break because too few contigs to assemle
    
#ifdef DEBUG_ENGINE
    std::cout << "...round " << (yy+1) << std::endl;
#endif    
    
    // do the second round (on assembled contigs)
    SeqLib::UnalignedSequenceVector tmpc;
    for (auto& j : m_contigs)
      if (j.Seq.length() > m_readlen)
	tmpc.push_back(j);
    
    ReadTable pRTc0(tmpc);
    m_contigs.clear();
    doAssembly(&pRTc0, m_contigs, yy);      
    
  }
  return true;
}


// call the assembler
void svabaAssemblerEngine::doAssembly(ReadTable *pRT, SeqLib::UnalignedSequenceVector &contigs, int pass) {
  
  if (pRT->getCount() == 0)
    return;

  // clear the hits stream
  std::stringstream hits_stream, asqg_stream;

  // set the paramters for this run
  double errorRate = m_error_rate;
  int min_overlap = m_min_overlap;
  int cutoff = 0;

  if (pass > 1) {
    //errorRate = 0.015f; // up to one bp mismatch
    errorRate = 0;
    m_error_rate = errorRate;
    min_overlap = m_min_overlap * 2;
  }

  bool exact = errorRate < 0.001f;

  // remove duplicates if running in exact mode
  ReadTable * pRT_nd = exact ? removeDuplicates(pRT) : pRT;    

  // forward
  SuffixArray* pSAf_nd = new SuffixArray(pRT_nd, 1, false); //1 is num threads. false is silent/no
  RLBWT *pBWT_nd = new RLBWT(pSAf_nd, pRT_nd);

  // reverse
  pRT_nd->reverseAll();
  SuffixArray * pSAr_nd = new SuffixArray(pRT_nd, 1, false);
  RLBWT *pRBWT_nd = new RLBWT(pSAr_nd, pRT_nd);
  pRT_nd->reverseAll();

  pSAf_nd->writeIndex();
  pSAr_nd->writeIndex();

  bool bIrreducibleOnly = true; // default
  int seedLength = 0;
  int seedStride = 0;
  if (!exact)
    calculateSeedParameters(m_readlen, min_overlap, seedLength, seedStride);

  svabaOverlapAlgorithm* pOverlapper = new svabaOverlapAlgorithm(pBWT_nd, pRBWT_nd, 
								     errorRate, seedLength,
								     seedStride, bIrreducibleOnly);
  
  pOverlapper->setExactModeOverlap(exact);
  pOverlapper->setExactModeIrreducible(exact);

  svabaASQG::HeaderRecord headerRecord;
  headerRecord.setOverlapTag(min_overlap);
  headerRecord.setErrorRateTag(errorRate);
  headerRecord.setInputFileTag("");
  headerRecord.setContainmentTag(true); // containments are always present
  headerRecord.setTransitiveTag(!bIrreducibleOnly);
  headerRecord.write(asqg_stream);    

  pRT_nd->setZero();

  size_t workid = 0;
  SeqItem si;

  size_t ocount = 0;
  while (pRT_nd->getRead(si) && (++ocount < MAX_OVERLAPS_PER_ASSEMBLY)) {
    
    SeqRecord read;
    read.id = si.id;
    read.seq = si.seq;
    OverlapBlockList obl;
    
    OverlapResult rr = pOverlapper->overlapRead(read, min_overlap, &obl);
    pOverlapper->writeOverlapBlocks(hits_stream, workid, rr.isSubstring, &obl);

    svabaASQG::VertexRecord record(read.id, read.seq.toString());
    record.setSubstringTag(rr.isSubstring);
    record.write(asqg_stream);

    ++workid;

  }

  std::string line;
  bool bIsSelfCompare = true;
  ReadInfoTable* pQueryRIT = new ReadInfoTable(pRT_nd);

  while(std::getline(hits_stream, line)) {
    size_t readIdx;
    size_t totalEntries;
    bool isSubstring; 
    OverlapVector ov;
    OverlapCommon::parseHitsString(line, pQueryRIT, pQueryRIT, pSAf_nd, pSAr_nd, bIsSelfCompare, readIdx, totalEntries, ov, isSubstring);

    for(OverlapVector::iterator iter = ov.begin(); iter != ov.end(); ++iter) {
       svabaASQG::EdgeRecord edgeRecord(*iter);
       edgeRecord.write(asqg_stream);
    }

  }
  
  // Get the number of strings in the BWT, this is used to pre-allocated the read table
  delete pOverlapper;
  delete pBWT_nd; 
  delete pRBWT_nd;
  delete pSAf_nd;
  delete pSAr_nd;
  if (exact) // only in exact mode did we actually allocate for pRT_nd, otherwise just pRT which we want to keep
    delete pRT_nd; 

  //#ifdef CLOCK_COUNTER
  //clock_t ctA_before = clock();
  //#endif
  
  // PERFORM THE ASSMEBLY
  trimLengthThreshold = 100; 
  numTrimRounds = 1; 
  StringGraph * oGraph = assemble(asqg_stream, min_overlap, maxEdges, bExact, 
	   trimLengthThreshold, bPerformTR, bValidate, numTrimRounds, 
	   resolveSmallRepeatLen, numBubbleRounds, gap_divergence, 
				  divergence, maxIndelLength, cutoff, m_id + "_", contigs, (pass > 0), m_write_asqg);
  
  // optionally output the graph structure
  if (m_write_asqg)
    write_asqg(oGraph, asqg_stream, hits_stream, pass);

  // this was allocated in assemble
  delete oGraph;
  delete pQueryRIT;

  // remove exact dups
  remove_exact_dups(contigs);
  
  // print out some results
#ifdef DEBUG_ENGINE
  std::cerr << " PASS " << pass << std::endl;
  print_results(contigs);
#endif
  
  return;
}

// not totally sure this works...
ReadTable* svabaAssemblerEngine::removeDuplicates(ReadTable* pRT) {

  // forward
  SuffixArray* pSAf = new SuffixArray(pRT, 1, false); //1 is num threads. false is silent/no
  RLBWT *pBWT= new RLBWT(pSAf, pRT);

  // reverse
  pRT->reverseAll();
  SuffixArray * pSAr = new SuffixArray(pRT, 1, false);
  RLBWT *pRBWT = new RLBWT(pSAr, pRT);
  pRT->reverseAll();

  svabaOverlapAlgorithm* pRmDupOverlapper = new svabaOverlapAlgorithm(pBWT, pRBWT, 
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
  }

  delete pRmDupOverlapper;
  delete pBWT; 
  delete pRBWT;
  delete pSAf;
  delete pSAr;

  return pRT_nd;
}

void svabaAssemblerEngine::write_asqg(const StringGraph* oGraph, std::stringstream& asqg_stream, std::stringstream& hits_stream, int pass) const {

  // write ASQG to file for visualization
  std::ofstream ofile(m_id + "pass_" + std::to_string(pass) + ".asqg", std::ios::out);
  ofile << asqg_stream.str();
  ofile.close();
  
  // write the hits stream file
  std::ofstream ofile_hits(m_id + "pass_" + std::to_string(pass) + ".hits", std::ios::out);
  ofile_hits << hits_stream.str();
  ofile_hits.close();

  // write the connected components file
  std::ofstream ofile_cc(m_id + "pass_" + std::to_string(pass) + ".cc", std::ios::out);
  for (auto& i : oGraph->m_connected_components)
    ofile_cc << i.first << "\t" << i.second << std::endl;
  ofile_cc.close();
  
  // write the labels
  std::ofstream ofile_lb(m_id + "pass_" + std::to_string(pass) + ".lb", std::ios::out);
  for (auto& i : oGraph->m_reads_on_contigs)
    for (auto& r : i.second)
      ofile_lb << i.first << "\t" << r << std::endl;
  ofile_lb.close();

}

void svabaAssemblerEngine::remove_exact_dups(SeqLib::UnalignedSequenceVector& cc) const {

  std::set<std::string> ContigDeDup;
  SeqLib::UnalignedSequenceVector cvec;
  for (auto& i : cc) {
    if (!ContigDeDup.count(i.Seq)) {
      ContigDeDup.insert(i.Seq);
      cvec.push_back(i);
    } else {
      //std::cerr << "Filtered out a contig for having exact duplicate with another contig" << std::endl;
    }
  }
  cc = cvec;
}

void svabaAssemblerEngine::print_results(const SeqLib::UnalignedSequenceVector& cc) const {

  if (cc.size() >= 1) {
    std::cout << "Contig Count: " << cc.size() << " at " << m_id << std::endl;
    for (auto& i : cc) 
    	std::cout << "   " << i.Name << " " << i.Seq.length() << " " << i.Seq << std::endl;
  }

}

// lifted from OverlapAlgorithm
void svabaAssemblerEngine::calculateSeedParameters(int read_len, int minOverlap, int& seed_length, int& seed_stride) const {

    seed_length = 0;
    
    // The maximum possible number of differences occurs for a fully-aligned read
    int max_diff_high = static_cast<int>(m_error_rate * read_len);

    // Calculate the seed length to use
    // If the error rate is so low that no differences are possible just seed
    // over the entire minOverlap region
    if(max_diff_high > 0)
    {
        // Calculate the maximum number of differences between two sequences that overlap
        // by minOverlap
        int max_diff_low = static_cast<int>(m_error_rate * minOverlap);

         if(max_diff_low == 0)
            max_diff_low = 1;
         
         int seed_region_length = static_cast<int>(ceil(max_diff_low / m_error_rate));
         int num_seeds_low = max_diff_low + 1;
         seed_length = static_cast<int>(seed_region_length / num_seeds_low);
         if(seed_length > static_cast<int>(minOverlap))
            seed_length = minOverlap;
    }
    else
    {
        seed_length = minOverlap;
    }
    seed_stride = seed_length;    

}
