#include "SnowmanBamWalker.h"
#include <sstream>

//#define DEBUG_SNOWMAN_BAMWALKER 1
#define MIN_MAPQ_FOR_MATE_LOOKUP 1

static const std::string FWD_ADAPTER_A = "AGATCGGAAGAGC";
static const std::string FWD_ADAPTER_B = "AGATCGGAAAGCA";
static const std::string REV_ADAPTER = "GCTCTTCCGATCT";

void SnowmanBamWalker::addCigar(BamRead &r)
{
  // this is a 100% match
  if (r.CigarSize() == 0)
    return;

  std::stringstream ss;
  int pos = r.Position(); // position ON REFERENCE
  
  for (auto& i : r.GetCigar()) {
      // if it's a D or I, add it to the list
      if (i.Type == 'D' || i.Type == 'I') {	
	//ss << r_id(r) << "_" << pos << "_" << /*r_cig_len(r,i) <<*/ r_cig_type(r, i);
	ss << r.ChrID() << "_" << pos << "_" << i.Length << i.Type;
	cigmap[ss.str()]++;
	ss.str("");
      }
      
    // move along the REFERENCE
      if (!(i.Type == 'I') && !(i.Type == 'S') && !(i.Type == 'H'))
      pos += i.Length;
  }
  
  return;
  
}

bool SnowmanBamWalker::isDuplicate(BamRead &r)
{

  // deduplicate by query-bases / position
  std::string sname = std::to_string(r.Position()) + "_" + std::to_string(r.MatePosition()); // + r.Sequence();    
  // deduplicate by Name
  //std::string uname = r.Qname() + "_" + std::to_string(r.FirstFlag());
  
  // its not already add, insert
  bool sname_pass = false;
  //if (name_map.count(uname) == 0) { // && seq_map.count(sname) == 0) {  
  //  uname_pass = true;
  //  name_map.insert(std::pair<std::string, int>(uname, true));
  //} 
  if (seq_set.count(sname) == 0) {
    sname_pass = true;
    seq_set.insert(sname);
    //seq_set.insert(std::pair<std::string, int>(sname, true));
  }
  
  return !sname_pass;
  //return (!uname_pass || !sname_pass);

}

void SnowmanBamWalker::readBam()
{
  
  BamRead r;

  //BamReadVector all_reads;

  bool rule_pass;
  int reads_to_start = reads.size();

  //std::cerr << "**Starting read for " << (prefix == "n" ? "NORMAL" : "TUMOR") << (get_mate_regions ? "**" : " MATE REGIONS**");
  //if (get_mate_regions && m_region.size())
  //   std::cerr << " on region " << m_region[0] << std::endl;
  //else if (get_mate_regions && m_region.size() == 0)
  //  std::cerr << " on WHOLE GENOME" << std::endl;
  //else
  //  std::cerr << " on " << m_region.size() << " regions " << std::endl;

  SnowTools::BamReadVector mate_reads;

  size_t countr = 0;
  while (GetNextRead(r, rule_pass))
    {

      bool qcpass = !r.DuplicateFlag() && !r.QCFailFlag() && !r.SecondaryFlag();

      // add to all reads pile for kmer correction
      if (qcpass) {
	cov.addRead(r);
	//all_reads.push_back(r);
      }
      
      // check if it passed blacklist
      bool blacklisted = false;
      if (blacklist.size() && blacklist.findOverlapping(r.asGenomicRegion()))
	blacklisted = true;
      rule_pass = rule_pass && !blacklisted;

      // check if has adapter
      if (adapter_trim) 
	rule_pass = rule_pass && !hasAdapter(r);

      // add to weird coverage
      if (rule_pass)
	weird_cov.addRead(r);

      // add to the cigar map for all non-duplicate reads
      if (qcpass)
	addCigar(r);

      bool is_dup = isDuplicate(r);
      //bool is_dup = false;

      if (!rule_pass)
	r.AddIntTag("VR", -1); 
      else if (rule_pass && !is_dup)

	{
	  // optional tag processing
	  //r.RemoveAllTags(); // cut down on memory
	  //r_remove_tag(r, "R2");
	  //r_remove_tag(r, "Q2");
	  //r_remove_tag(r, "OQ");
	  //r_add_Z_tag(r, "RL", rule_pass);

	  ++countr;
	  if (countr % 10000 == 0 && m_region.size() == 0)
	    std::cerr << "...read in " << SnowTools::AddCommas<size_t>(countr) << " weird reads for whole genome read-in. At pos " << r.Brief() << std::endl;
	  
	  if (disc_only) {
	    //	    r.SetSequence("A");
	    //r.RemoveAllTags();
	  }

	  // add the ID tag
	  std::string srn =  prefix+std::to_string(r.AlignmentFlag()) + "_" + r.Qname();
	  r.AddZTag("SR", srn);
	  r.AddIntTag("VR", 1);
	  
	  // get the weird coverage
	  weird_cov.addRead(r);

	  reads.push_back(r); // adding later because of kmer correction

	  // check that we're not above the read limit
	  if (reads.size() > m_limit && m_limit != 0)
	    return;
	  
	}
    } // end the read loop

  // clear out the added reads if hit limit on mate lookup
  if (m_keep_limit > 0 && m_num_reads_kept*1.05 > m_keep_limit) {
    BamReadVector tmpr;
    for (int i = 0; i < reads_to_start; ++i)
      tmpr.push_back(reads[i]);
    reads = tmpr;
  }

  if (reads.size() < 3)
    return;
      
  // get rid of repats
  //removeRepeats();

  // clean out the buffer
  subSampleToWeirdCoverage(max_cov);

  // calculate the mate region
  if (get_mate_regions && m_region.size() /* don't get mate regions if reading whole bam */) {
    calculateMateRegions();
  }
  
}

void SnowmanBamWalker::KmerCorrect() {

  // do the kmer filtering
  KmerFilter kmer;
  /*int kcor = */kmer.correctReads(/*all_reads*/reads);
  //std::cerr << "     kmer corrected " << SnowTools::AddCommas<int>(kcor) << " reads of " << SnowTools::AddCommas<size_t>(reads.size()) << std::endl; 

}


void SnowmanBamWalker::subSampleToWeirdCoverage(double max_coverage) {
  
  BamReadVector new_reads;

  for (auto& r : reads)
    {
      double this_cov1 = weird_cov.getCoverageAtPosition(r.ChrID(), r.Position());
      double this_cov2 = weird_cov.getCoverageAtPosition(r.ChrID(), r.PositionEnd());
      double this_cov = std::max(this_cov1, this_cov2);
      double sample_rate = 1; // dummy, always set if max_coverage > 0
      if (this_cov > 0) 
	sample_rate = 1 - (this_cov - max_coverage) / this_cov; // if cov->inf, sample_rate -> 0. if cov -> max_cov, sample_rate -> 1
      
      // this read should be randomly sampled, cov is too high
      if (this_cov > max_coverage) 
	{
	  uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(r.Qname().c_str()) ^ m_seed);
	  if ((double)(k&0xffffff) / 0x1000000 <= sample_rate) // passed the random filter
	    new_reads.push_back(r);
	}
      else // didn't have a coverage problems
	{
	  new_reads.push_back(r);
	}
      
    }

  reads = new_reads;
}

void SnowmanBamWalker::calculateMateRegions() {

  if (m_region.size() == 0) {
    std::cerr << "Attempting to calculate mate region with no BamWalker region defined." << std::endl;
    return;
  }
  SnowTools::GenomicRegion main_region = m_region.at(0);

  MateRegionVector tmp_mate_regions;

  for (auto& r : reads)
    {
      
      if (r.MateChrID() > 22) // no Y or M
	continue;
      
      MateRegion mate(r.MateChrID(), r.MatePosition(), r.MatePosition());
      mate.pad(500);
      mate.partner = main_region;
      
      // if mate not in main interval, add a padded version
      if (!main_region.getOverlap(mate) && r.MapQuality() >= MIN_MAPQ_FOR_MATE_LOOKUP) 
	tmp_mate_regions.add(mate);
       
    }

  // merge it down
  tmp_mate_regions.mergeOverlappingIntervals();

  // create an interval tree for fast lookup
  tmp_mate_regions.createTreeMap();

#ifdef DEBUG_SNOWMAN_BAMWALKER
  std::cerr << "Mate regions are" << std::endl;
  for (auto& i : tmp_mate_regions) 
    std::cerr << i << std::endl;
#endif
  
  // get the counts
  for (auto& r : reads)
    {
      SnowTools::GenomicRegion mate(r.MateChrID(), r.MatePosition(), r.MatePosition());

      // skip if it is unmapped
      if (r.MatePosition() < 0)
	continue; 

      //debug
      //for (auto& kk : tmp_mate_regions)
      //	if (kk.getOverlap(mate))
      //  std::cerr << "tmp_region " << kk << " overlaps with " << mate << std::endl;

      // if mate not in main interval, check which mate regions it's in
      if (!main_region.getOverlap(mate) && r.MapQuality() > 0) 
	{
	  for (auto& k : tmp_mate_regions)
	    if (k.getOverlap(mate)) {
	      k.count++;
	      continue;
	    }
	  //std::vector<int32_t> query_id, subject_id;
	  //GRC tmp(mate);
	  //tmp.createTreeMap();
	  //tmp_mate_regions.findOverlaps(tmp, query_id, subject_id);

	  // we should re-find it. If not, probably because it's right on the edge
	  // and we dont count mate regions in the original window. It's OK
	  //if (query_id.size() != 1) {
	    //std::cerr << "SnowmanBamWalker::calculateMateRegions error. Read mate does not match a mate region for mate read " << mate << std::endl;
	  // continue;
	  //}
	      
	  // update the count
	  //tmp_mate_regions[query_id[0]].count++;
	}


    }

  // keep only ones with 2+ reads
  for (auto& i : tmp_mate_regions)
    {
      if (i.count >= 2)
	mate_regions.add(i);
    }
  
  #ifdef DEBUG_SNOWMAN_BAMWALKER
  std::cerr << "Final mate regions are:" << std::endl;
  for (auto& i : mate_regions)
    std::cerr << i << " read count " << i.count << std::endl;
  #endif

}

void SnowmanBamWalker::removeRepeats()
{
#ifdef DEBUG_SNOWMAN_BAMWALKER
  std::cerr << "...removing repeats on " << reads.size() << " reads" << std::endl;
#endif
  std::string POLYA = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
  std::string POLYT = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
  std::string POLYC = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
  std::string POLYG = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
  std::string POLYAT = "ATATATATATATATATATATATATATATATATATATATAT";
  std::string POLYTC = "TCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTC";
  std::string POLYAG = "AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG";
  std::string POLYCG = "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG";
  std::string POLYTG = "TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG";
  std::string POLYCA = "CACACACACACACACACACACACACACACACACACACACA";
  
  BamReadVector new_reads;

  for (auto& r : reads)
    {
      int dum = 0;
      std::string seq = r.QualityTrimmedSequence(4, dum);

      if (seq.length() >= 40)
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
	  new_reads.push_back(r);
    }

#ifdef DEBUG_SNOWMAN_BAMWALKER
  std::cerr << "...removed " << (reads.size() - new_reads.size()) << " repeats" << std::endl;
#endif


  reads = new_reads;
}

void SnowmanBamWalker::filterMicrobial(SnowTools::BWAWrapper * b) {

  BamReadVector new_reads;

  for (auto& r : reads) {
    BamReadVector micro_alignments;
    b->alignSingleSequence(r.Sequence(), r.Qname(), micro_alignments, true);
    if (micro_alignments.size() == 0)
      new_reads.push_back(r);
    else if (micro_alignments[0].NumMatchBases() <= r.Length() * 0.8)
      new_reads.push_back(r);
    
    // store the microbe alignments
    //for (auto& a : micro_alignments) {
    //  if (a.NumMatchBases() > opt::readlen * 0.8) 
    //	bb_microbe.push_back(a);//b_microbe_writer.WriteAlignment(a);
  }
  //std::cerr << "...filtered out " << (reads.size() - new_reads.size()) << " microbial reads "  << std::endl;
  reads = new_reads;
}

bool SnowmanBamWalker::hasAdapter(const BamRead& r) const {

  // keep it if it has indel
  if (r.MaxDeletionBases() || r.MaxInsertionBases() || !r.InsertSize() || r.NumClip() < 5)
    return false;
  
  // toss it then if isize explans clip
  int exp_ins_size = r.Length() - r.NumClip(); // expected isize if has adapter
  if ((exp_ins_size - 4) < std::abs(r.InsertSize()) && (exp_ins_size+4) > std::abs(r.InsertSize()))
    return true;

  return false;
  /*
  if (std::abs(r.InsertSize()) < 300 && 
	   r.PairMappedFlag() && (r.ChrID() == r.MateChrID())) {
    std::string seqr = r.Sequence();
    if (seqr.find("AGATCGGAAGAGC") != std::string::npos || seqr.find("AGATCGGAAAGCA") != std::string::npos || 
	seqr.find("GCTCTTCCGATCT") != std::string::npos) {
      has_adapter = true;
    }
  */  

}
