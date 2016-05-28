#include "SnowmanBamWalker.h"

//#define DEBUG_SNOWMAN_BAMWALKER 1
#define MIN_MAPQ_FOR_MATE_LOOKUP 0

//#define QNAME "D0ENMACXX111207:3:2105:3186:55517"

// dont read in more reads than this at once
#define FAIL_SAFE 8000

static const std::string ILLUMINA_PE_PRIMER_2p0 = "CAAGCAGAAGACGGCAT";
static const std::string FWD_ADAPTER_A = "AGATCGGAAGAGC";
static const std::string FWD_ADAPTER_B = "AGATCGGAAAGCA";
static const std::string REV_ADAPTER = "GCTCTTCCGATCT";

bool SnowmanBamWalker::addCigar(BamRead &r, const SnowTools::DBSnpFilter* d) {

  // this is a 100% match
  if (r.CigarSize() == 1)
    return true;
  std::stringstream cigar_ss;
  cigar_ss.str(std::string());
  int pos = r.Position(); // position ON REFERENCE

  bool has_dbsnp_hit = false;
  
  for (auto& i : r.GetCigar()) {
      // if it's a D or I, add it to the list
      if (i.Type() == 'D' || i.Type() == 'I') {	
	//cigar_ss << r_id(r) << "_" << pos << "_" << /*r_cig_len(r,i) <<*/ r_cig_type(r, i);

	cigar_ss << r.ChrID() << "_" << pos << "_" << i.Length() << i.Type();
	++cigmap[cigar_ss.str()];

	//cigar_ss << r.ChrID() << "_" << pos; // << "_" << i.Length() << i.Type();
	if (prefix == "n")
	  for (int kk = -2; kk <= 2; ++kk) { // add this and neighbors
	    uint32_t cigs = r.ChrID() * 1e9 + pos;
	    ++cig_pos[cigs];
	  }
	
	/*if (d) {
	  cigar_ss.str(std::string());
	  cigar_ss << r.ChrID() << "_" << pos;
	  has_dbsnp_hit = d->queryHash(cigar_ss.str());
	  }*/
	//std::cerr << d << " has_db " << has_dbsnp_hit << " hash " << cigar_ss.str() << std::endl;
	cigar_ss.str(std::string());
      }
      
    // move along the REFERENCE
      if (!(i.Type() == 'I') && !(i.Type() == 'S') && !(i.Type() == 'H'))
      pos += i.Length();
  }
  
  return has_dbsnp_hit;
  
}

bool SnowmanBamWalker::isDuplicate(BamRead &r)
{

  // deduplicate by query-bases / position
  //int pos_key = r.Position() + r.MatePosition()*3 + 1332 * r.NumMatchBases() + 12345 * r.AlignmentFlag() + 221*r.InsertSize();
  //std::string key = std::to_string(pos_key) + r.Sequence();
  std::string key = r.QualitySequence() + std::to_string(r.Position()) + "_" + std::to_string(r.MatePosition());

#ifdef QNAME
  if (r.Qname() == QNAME)
    std::cerr << "dedupe key " << key << std::endl;
#endif

  // its not already add, insert
  if (!seq_set.count(key)) {
    seq_set.insert(key);
    return false;
  }
  
  return true;
}

SnowTools::GRC SnowmanBamWalker::readBam(std::ofstream * log, const SnowTools::DBSnpFilter* dbs)
{
  
  BamRead r;

  bool rule_pass;
  int reads_to_start = reads.size();

  //std::cerr << "**Starting read for " << (prefix == "n" ? "NORMAL" : "TUMOR") << (get_mate_regions ? "**" : " MATE REGIONS**");
  //if (get_mate_regions && m_region.size())
  //   std::cerr << " on region " << m_region[0] << std::endl;
  //else if (get_mate_regions && m_region.size() == 0)
  //  std::cerr << " on WHOLE GENOME" << std::endl;
  //else
  //  std::cerr << " on " << m_region.size() << " regions " << std::endl;

  SnowTools::GRC bad_regions;

  SnowTools::BamReadVector mate_reads;

  size_t countr = 0;
  int curr_reads = reads.size();
  while (GetNextRead(r, rule_pass)) {

#ifdef QNAME
      if (r.Qname() == QNAME) 
	std::cerr << " read seen " << r << std::endl;
#endif

      if ((reads.size() - curr_reads) > m_limit && log && m_limit > 0) {
	(*log) << "...breaking at " << r.Brief() << 
	  (get_mate_regions ? " in main window " : " in mate window " ) 
	       << m_region[m_region_idx]
	       << " with " << SnowTools::AddCommas(reads.size() - curr_reads) 
	       << " weird reads. Limit: " << SnowTools::AddCommas(m_limit) << std::endl;
	
	curr_reads = reads.size();
	bad_regions.add(m_region[m_region_idx]);
	
	// try next region, return if no others to try
	++m_region_idx; // increment to next region
	if (m_region_idx >= m_region.size()) /// no more regions left
	  break;
	else { // move to next region
	  countr = 0;
	  __set_region(m_region[m_region_idx]);
	  continue;
	}
	break;
      }

      bool qcpass = !r.DuplicateFlag() && !r.QCFailFlag(); 

      // add to all reads pile for kmer correction
      if (!r.DuplicateFlag() && !r.QCFailFlag()) {
	cov.addRead(r); 
      }
      
      // check if it passed blacklist
      bool blacklisted = false;
      if (blacklist.size() && blacklist.findOverlapping(r.asGenomicRegion()))
	blacklisted = true;
      rule_pass = rule_pass && !blacklisted;

      // check if has adapter
      //if (adapter_trim && rule_pass)  //debug
	rule_pass = rule_pass && !hasAdapter(r);

      // add to weird coverage
      if (rule_pass)
	weird_cov.addRead(r);
      
      if (prefix == "n" && !r.DuplicateFlag() && !r.QCFailFlag())  {
	if (r.CigarSize() >= 6 || r.MapQuality() <= 5 || (r.QualitySequence().length() < 0.6 * readlen))
	  bad_cov.addRead(r, true);
	if (r.NumClip()) {
	  clip_cov.addRead(r,true);
	}
      }

      //bool dbsnp_fail = false;
      // add to the cigar map for all non-duplicate reads
      if (qcpass && get_mate_regions && rule_pass) // only add cigar for non-mate regions
	addCigar(r, dbs);
      //dbsnp_fail = dbsnp_fail && (r.MaxInsertionBases() || r.MaxDeletionBases());
      //dbsnp_fail = false; // turn this off;

      bool is_dup = false;
      if (rule_pass)
	is_dup = isDuplicate(r);

#ifdef QNAME
      if (r.Qname() == QNAME) 
	std::cerr << " read has qcpass " << qcpass << " blacklisted " << blacklisted << " rule_pass " << rule_pass << " isDup " << is_dup << " " << r << std::endl;
      if (r.Qname() == QNAME)
	std::cerr << " qual " << r.Qualities() << " trimmed seq " << r.QualitySequence() << std::endl;
#endif

      //if (!rule_pass)
      //r.AddIntTag("VR", -1); 
      if (rule_pass && !is_dup/* && (!dbsnp_fail || !(r.MaxInsertionBases() || r.MaxDeletionBases()))*/)
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

	  // add the ID tag
	  std::string srn =  prefix + "_" + std::to_string(r.AlignmentFlag()) + "_" + r.Qname();
	  assert(srn.length());
	  r.AddZTag("SR", srn);
	  //r.AddIntTag("VR", 1);
	  
	  // get the weird coverage
	  weird_cov.addRead(r);
	  
#ifdef QNAME
	  if (r.Qname() == QNAME) 
	std::cerr << " read added " << r << std::endl;
#endif
      reads.push_back(r); // adding later because of kmer correction

	}
  } // end the read loop

#ifdef QNAME
  for (auto& j : reads)
      if (j.Qname() == QNAME) 
	std::cerr << " read KEPT right after loop " << j << std::endl;
#endif


  // clear out the added reads if hit limit on mate lookup
  if (m_keep_limit > 0 && m_num_reads_kept*1.05 > m_keep_limit) {
    BamReadVector tmpr;
    for (int i = 0; i < reads_to_start; ++i)
      tmpr.push_back(reads[i]);
    reads = tmpr;
  }

#ifdef QNAME
  for (auto& j : reads)
      if (j.Qname() == QNAME) 
	std::cerr << " read KEPT-PREFILTER " << j << std::endl;
#endif

  if (reads.size() < 3)
    return bad_regions;
      
  // get rid of repats
  //removeRepeats();

  // clean out the buffer
  subSampleToWeirdCoverage(max_cov);

#ifdef QNAME
  for (auto& j : reads)
      if (j.Qname() == QNAME) 
	std::cerr << " read KEPT-POST SUBSAMPLE " << j << std::endl;
#endif

  // realign discordants
  //if (main_bwa)
  //  realignDiscordants(reads);

  // calculate the mate region
  if (get_mate_regions && m_region.size() /* don't get mate regions if reading whole bam */) {
    calculateMateRegions();
  }

#ifdef QNAME
  for (auto& j : reads)
      if (j.Qname() == QNAME) 
	std::cerr << " read KEPT " << j << std::endl;
#endif

 
  return bad_regions;
 
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
	  #ifdef QNAME
	  if (r.Qname() == QNAME) {
	    std::cerr << "subsampling because this_cov is " << this_cov << " and max cov is " << max_coverage << " at position " << r.Position() << " and end position " << r.PositionEnd() << std::endl;
	    std::cerr << " this cov 1 " << this_cov1 << " this_cov2 " << this_cov2 << std::endl;
	  }
	  #endif
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
	  //std::cerr << r << " getOverlap " << tmp_mate_regions[0].getOverlap(mate) << std::endl;
	  for (auto& k : tmp_mate_regions)
	    if (k.getOverlap(mate)) {
	      ++k.count;
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


// defunct
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
      //bool was_trimmed = false;
      //std::string seq = r.QualityTrimmedSequence(4, dum, was_trimmed);
      std::string seq = ""; // dummy

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
    bool hardclip = false;
    b->alignSingleSequence(r.Sequence(), r.Qname(), micro_alignments, hardclip, true, 1000);
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
  if (r.MaxDeletionBases() || r.MaxInsertionBases() || !r.InsertSize() || !r.NumClip())
    return false;
  
  // toss it then if isize explans clip
  int exp_ins_size = r.Length() - r.NumClip(); // expected isize if has adapter
  if ((exp_ins_size - 4) < std::abs(r.InsertSize()) && (exp_ins_size+4) > std::abs(r.InsertSize()))
    return true;

  std::string seq = r.Sequence();
  if (seq.find(ILLUMINA_PE_PRIMER_2p0) != std::string::npos)
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

void SnowmanBamWalker::realignDiscordants(SnowTools::BamReadVector& reads) {

  const int secondary_cap = 20;
  SnowTools::BamReadVector brv;

  for (auto& r : reads) {

    if (std::abs(r.InsertSize()) >= 1000) {
      SnowTools::BamReadVector als;
      main_bwa->alignSingleSequence(r.Sequence(), r.Qname(), als, false, 0.95, secondary_cap);
      if (als.size() < 10) {
	brv.push_back(r);
	r.AddIntTag("DD", als.size());
      }

    } else {
      brv.push_back(r);
    }
  }
  
  //std::cerr << "...removed " << (reads.size()-brv.size()) << " of " << reads.size() << std::endl;
  reads = brv;
  
}
