#include "SnowmanAssemblerEngine.h"

#include "SnowmanAssemble.h"
#include "SnowmanOverlapAlgorithm.h"
#include "SnowmanASQG.h"
#include "SGACommon.h"
#include "OverlapCommon.h"

#define MAX_OVERLAPS_PER_ASSEMBLY 200000

void SnowmanAssemblerEngine::fillReadTable(ReadVec& r)
{
  
  m_reads = r;

  // make the reads tables
  for (auto& i : r) {

    SeqItem si;
    string sr, seq = "";

    r_get_SR(i, sr);
    r_get_trimmed_seq(i, seq); 
    assert(sr.length());
    assert(seq.length());
    
    si.id = sr;
    si.seq = seq;
    m_pRT.addRead(si);

  }
  
}

bool SnowmanAssemblerEngine::performAssembly() 
{
  if (m_pRT.getCount() < 3)
    return false;

#ifdef DEBUG_ENGINE
  std::cout << "Doing assembly on: " << m_id << " with " << m_pRT.getCount() << " reads" << std::endl; 
#endif
  

  ContigVector contigs0;
  
  // do the first round (on raw reads)
  doAssembly(&m_pRT, contigs0, 0);
  
  for (size_t yy = 1; yy != 2; yy++) {
    
    // do the second round (on assembled contigs)
    ReadTable pRTc0(contigs0);
    m_contigs.clear();
    doAssembly(&pRTc0, m_contigs, yy);      
    contigs0 = m_contigs;
    
  }

#ifdef DEBUG_ENGIN
  std::cout << "Did assembly on: " << m_id << " with " << pRT.getCount() << " reads and " << m_contigs.size() << " contigs " << std::endl;
#endif

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
  
  // rmdup attempt
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

  pSAf->writeIndex();
  pSAr->writeIndex();

  //#ifdef CLOCK_COUNTER  
  //(*clock_counter) << pass << "," << "suffix" << "," << (clock() - ct1_before)<< "," << pRT->getCount() << std::endl;
  //#endif

  double errorRate = m_error_rate;
  int min_overlap = m_min_overlap;

  int cutoff = 0;
  if (pass > 0) {
    min_overlap = 50;
    errorRate = 0.05;
    cutoff = m_readlen + 30;
  }

  int seedLength = min_overlap;
  int seedStride = seedLength;
  bool bIrreducibleOnly = true; // default

  SnowmanOverlapAlgorithm* pOverlapper = new SnowmanOverlapAlgorithm(pBWT, pRBWT, 
                                                       errorRate, seedLength, 
                                                       seedStride, bIrreducibleOnly);


  bool exact = false;
  //exact = errorRate < 0.001f;
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

  pRT->setZero();

  size_t workid = 0;
  SeqItem si;
  //#ifdef CLOCK_COUNTER
  //clock_t ct_before = clock();
  //#endif

  size_t ocount = 0;
  while (pRT->getRead(si) && (++ocount < MAX_OVERLAPS_PER_ASSEMBLY)) {

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
  ReadInfoTable* pQueryRIT = new ReadInfoTable(pRT);

  while(getline(hits_stream, line)) {
    size_t readIdx;
    size_t totalEntries;
    bool isSubstring; 
    OverlapVector ov;
    OverlapCommon::parseHitsString(line, pQueryRIT, pQueryRIT, pSAf, pSAr, bIsSelfCompare, readIdx, totalEntries, ov, isSubstring);
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
  
  //#ifdef CLOCK_COUNTER
  //clock_t ctA_before = clock();
  //#endif
  
  // PERFORM THE ASSMEBLY
  assemble(asqg_stream, min_overlap, maxEdges, bExact, 
	   trimLengthThreshold, bPerformTR, bValidate, numTrimRounds, 
	   resolveSmallRepeatLen, numBubbleRounds, gap_divergence, 
	   divergence, maxIndelLength, cutoff, m_id + "_", contigs);
  
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
    if (opt::verbose > 3)
      for (auto& i : contigs) 
	std::cout << "   " << i.getID() << " " << i.getSeq().length() << " " << i.getSeq() << std::endl;
  }
#endif
  
  return;
}
