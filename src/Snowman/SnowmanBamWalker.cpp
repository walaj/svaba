#include "SnowmanBamWalker.h"
#include <sstream>

void SnowmanBamWalker::addCigar(Read &r)
{
  // this is a 100% match
  if (r_cig_size(r) == 0)
    return;

  std::stringstream ss;
  int pos = r_pos(r); // position ON REFERENCE
  
  for (size_t i = 0; i < r_cig_size(r); i++) {

    // if it's a D or I, add it to the list
    if (r_cig_type(r,i) == 'D' || r_cig_type(r,i) == 'I') {	
      //ss << r_id(r) << "_" << pos << "_" << /*r_cig_len(r,i) <<*/ r_cig_type(r, i);
      ss << r_id(r) << "_" << pos << "_" << r_cig_len(r,i) << r_cig_type(r, i);
      cigmap[ss.str()]++;
      ss.str("");
    }

    // move along the REFERENCE
    if (!(r_cig_type(r, i) == 'I') && !(r_cig_type(r, i) == 'S') && !(r_cig_type(r,i) == 'H'))
      pos += r_cig_len(r, i);
  }
  
  return;
  
}

bool SnowmanBamWalker::checkIfDuplicate(Read &r)
{

  // deduplicate by query-bases / position
  std::string sname = std::to_string(r_id(r)) + "_" + std::to_string(r_pos(r)) + "_" + std::to_string(r_mid(r)) + "_" + std::to_string(r_mpos(r));    
  // deduplicate by Name
  std::string uname = r_qname(r) + "_" + std::to_string(r_is_first(r));
  
  // its not already add, insert
  bool uname_pass = false, sname_pass = false;
  if (name_map.count(uname) == 0) { // && seq_map.count(sname) == 0) {  
    uname_pass = true;
    name_map.insert(pair<string, int>(uname, true));
  } 
  if (seq_map.count(sname) == 0) {
    sname_pass = true;
    seq_map.insert(pair<string, int>(sname, true));
  }
  
  return (!uname_pass || !sname_pass);

}

void SnowmanBamWalker::readBam()
{
  
  Read r;

  // start the timer
  //clock_t startr = clock();

  std::string rule_pass;
  while (GetNextRead(r, rule_pass))
    {
      bool valid = rule_pass != "";
      bool qcpass = !r_is_dup(r) && r_is_qc_fail(r);

      // add to the coverage calcualtion
      if (get_coverage && qcpass)
	cov.addRead(r);
      
      // add to the cigar map for all non-duplicate reads
      if (qcpass)
	addCigar(r);

      //bool is_dup = checkIfDuplicate(r);
      bool is_dup = false;
      
      if (valid && !is_dup)
	{
	  // optional tag processing
	  //r_remove_tag(r, "R2");
	  //r_remove_tag(r, "Q2");
	  //r_remove_tag(r, "OQ");
	  //r_add_Z_tag(r, "RL", rule_pass);
	  
	  // add the ID tag
	  std::string srn =  prefix+to_string(r_flag(r)) + "_" + r_qname(r);
	  r_add_Z_tag(r, "SR", srn);
	  
	  // get the weird coverage
	  weird_cov.addRead(r);

	  // add it to the final buffer
	  reads.push_back(r);

	  // check that we're not above the read limit
	  if (reads.size() > m_limit && m_limit != 0)
	    return;
	  
	}
    } // end the read loop

  // get rid of repats
  removeRepeats();

  // clean out the buffer
  subSampleToWeirdCoverage(max_cov);

  // calculate the mate regions
  if (get_mate_regions)
    calculateMateRegions();
}

void SnowmanBamWalker::subSampleToWeirdCoverage(double max_coverage)
{
  
  ReadVec new_reads;

  for (auto& r : reads)
    {
      double this_cov = weird_cov.getCoverageAtPosition(r_pos(r));
      double sample_rate = 1; // dummy, always set if max_coverage > 0
      if (this_cov > 0) 
	sample_rate = 1 - (this_cov - max_coverage) / this_cov; // if cov->inf, sample_rate -> 0. if cov -> max_cov, sample_rate -> 1

      // this read should be randomly sampled, cov is too high
      if (this_cov > max_coverage) 
	{
	  uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(bam_get_qname(r.get())) ^ m_seed);
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

void SnowmanBamWalker::calculateMateRegions() 
{

  SnowTools::GenomicRegion main_region = m_region.at(0);

  MateRegionVector tmp_mate_regions;
  //SnowTools::GenomicRegionCollection<MateRegion> tmp_mate_regions;

  for (auto& r : reads)
    {
      
      MateRegion mate(r_mid(r), r_mpos(r), r_mpos(r));
      mate.pad(500);
      
      // if mate not in main interval, add a padded version
      if (!main_region.getOverlap(mate) && r_mapq(r) > 0) 
	tmp_mate_regions.add(mate);
       
    }

  // merge it down
  tmp_mate_regions.mergeOverlappingIntervals();

  // create an interval tree for fast lookup
  tmp_mate_regions.createTreeMap();

  // get the counts
  for (auto& r : reads)
    {
      SnowTools::GenomicRegion mate(r_mid(r), r_mpos(r), r_mpos(r));

      // if mate not in main interval, check which mate regions it's in
      if (!main_region.getOverlap(mate) && r_mapq(r) > 0) 
	{
	  std::vector<size_t> query_id, subject_id;
	  GRC tmp(mate);
	  tmp_mate_regions.findOverlaps(tmp, query_id, subject_id);
	  
	  // we should re-find it
	  assert(query_id.size() == 1);

	  // update the count
	  tmp_mate_regions[query_id[0]].count++;
	}


    }

  // keep only ones with 3+ reads
  for (auto& i : tmp_mate_regions)
    {
      if (i.count > 2)
	mate_regions.add(i);
    }
}

void SnowmanBamWalker::removeRepeats()
{

  std::string POLYA = "AAAAAAAAAAAAAAAAAAAAAAAAA";
  std::string POLYT = "TTTTTTTTTTTTTTTTTTTTTTTTT";
  std::string POLYC = "CCCCCCCCCCCCCCCCCCCCCCCCC";
  std::string POLYG = "GGGGGGGGGGGGGGGGGGGGGGGGG";
  std::string POLYAT = "ATATATATATATATATATATATATATATATATATATATAT";
  std::string POLYCG = "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG";
  
  std::vector<Read> new_reads;

  for (auto& r : reads)
    {
      std::string seq;
      r_get_trimmed_seq(r, seq);

      if (seq.length() >= 40)
	if ((seq.find(POLYT) == std::string::npos) && 
	    (seq.find(POLYA) == std::string::npos) && 
	    (seq.find(POLYC) == std::string::npos) && 
	    (seq.find(POLYG) == std::string::npos) && 
	    (seq.find(POLYCG) == std::string::npos) && 
	    (seq.find(POLYAT) == std::string::npos) && 
	    (seq.find("N") == std::string::npos))
	  new_reads.push_back(r);
    }


  reads = new_reads;
}
