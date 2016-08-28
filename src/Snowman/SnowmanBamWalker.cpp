#include "SnowmanBamWalker.h"

//#define DEBUG_SNOWMAN_BAMWALKER 1
#define MIN_MAPQ_FOR_MATE_LOOKUP 0

#define ALL_READS_FAIL_SAFE 50000
//#define QNAME "H01PEALXX140819:2:2104:25705:6987"
//#define QFLAG 145

#define MIN_ISIZE_FOR_DISCORDANT_REALIGNMENT 1000

// how much to search aroudn mate region for read alignent
#define DISC_REALIGN_MATE_PAD 100

// trim this many bases from front and back of read when determining coverage
// this should be synced with the split-read buffer in BreakPoint2 for more accurate 
// representation of covearge of INFORMATIVE reads (eg ones that could be split)
#define INFORMATIVE_COVERAGE_BUFFER 4

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
  if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1))
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
      if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1)) 
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
      if (!r.DuplicateFlag() && !r.QCFailFlag() && !r.SecondaryFlag()) {
	cov.addRead(r, INFORMATIVE_COVERAGE_BUFFER, false); 
      }
      
      // check if it passed blacklist
      bool blacklisted = false;
      if (blacklist.size() && blacklist.findOverlapping(r.asGenomicRegion()))
	blacklisted = true;

      // check if in simple-seq
      if (!blacklisted && simple_seq->size()) {

	// check simple sequence overlaps
	SnowTools::GRC ovl = simple_seq->findOverlaps(r.asGenomicRegion(), true);
	
	int msize = 0;
	for (auto& j: ovl) {
	  int nsize = j.width() - r.MaxDeletionBases() - 1;
	  if (nsize > msize && nsize > 0)
	    msize = nsize;
	}

	if (msize > 30)
	  blacklisted = true;
      }

      rule_pass = rule_pass && !blacklisted;
      
      // check if has adapter
      rule_pass = rule_pass && !hasAdapter(r);

      // add to weird coverage
      if (rule_pass)
	weird_cov.addRead(r, 0, false);
      
      if (prefix == "n" && !r.DuplicateFlag() && !r.QCFailFlag())  {
	if (r.CigarSize() >= 6 || r.MapQuality() <= 5 || (r.QualitySequence().length() < 0.6 * readlen))
	  bad_cov.addRead(r, 0, true);
	if (r.NumClip()) {
	  clip_cov.addRead(r, 0, true);
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
      if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1)) 
	std::cerr << " read has qcpass " << qcpass << " blacklisted " << blacklisted << " rule_pass " << rule_pass << " isDup " << is_dup << " " << r << std::endl;
      if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1))
	std::cerr << " qual " << r.Qualities() << " trimmed seq " << r.QualitySequence() << std::endl;
#endif

      // add all the reads for kmer correction
      if (do_kmer_filtering && all_seqs.size() < ALL_READS_FAIL_SAFE && qcpass && !r.NumHardClip() && !blacklisted) {
	uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(r.Qname().c_str()) ^ m_seed);
	if ((double)(k&0xffffff) / 0x1000000 <= kmer_subsample) {
	  std::string qq = r.QualitySequence();
	  char * tmp_loc = (char*) malloc(qq.length() + 1); // + 1 for /0 terminator
	  all_seqs.push_back(strcpy(tmp_loc, qq.c_str()));
	}
      }
      
      if (rule_pass && !is_dup/* && (!dbsnp_fail || !(r.MaxInsertionBases() || r.MaxDeletionBases()))*/)
	{
	  ++countr;
	  if (countr % 10000 == 0 && m_region.size() == 0)
	    std::cerr << "...read in " << SnowTools::AddCommas<size_t>(countr) << " weird reads for whole genome read-in. At pos " << r.Brief() << std::endl;

	  // for memory conservation
	  r.RemoveTag("BQ");
	  r.RemoveTag("OQ");
	  r.RemoveTag("XT");
	  r.RemoveTag("XA");
	  r.RemoveTag("SA");

	  // add the ID tag
	  std::string srn =  prefix + "_" + std::to_string(r.AlignmentFlag()) + "_" + r.Qname();
	  assert(srn.length());
	  r.AddZTag("SR", srn);
	  
	  // get the weird coverage
	  weird_cov.addRead(r, 0, false);
	  
#ifdef QNAME
	  if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1)) 
	    std::cerr << " read added " << r << std::endl;
#endif
	  

	  reads.push_back(r); // adding later because of kmer correction
	  
	}
  } // end the read loop
  
#ifdef QNAME
  for (auto& j : reads)
      if (j.Qname() == QNAME && (j.AlignmentFlag() == QFLAG || QFLAG == -1)) 
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
      if (j.Qname() == QNAME && (j.AlignmentFlag() == QFLAG || QFLAG == -1)) 
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
      if (j.Qname() == QNAME && (j.AlignmentFlag() == QFLAG || QFLAG == -1)) 
	std::cerr << " read KEPT-POST SUBSAMPLE " << j << std::endl;
#endif

  // realign discordants
  if (main_bwa) {
    realignDiscordants(reads);

    // filter out the bad qnames
    SnowTools::BamReadVector tmp_vec;
    for (auto& j : reads) {
      if (j.NumClip() < 5 && !j.MaxInsertionBases() && !j.MaxDeletionBases() && bad_qnames.count(j.Qname())) {
#ifdef QNAME
	if (j.Qname() == QNAME && (j.AlignmentFlag() == QFLAG || QFLAG == -1))
	  std::cerr << " LOST TO QNAME FILTER AFTER REALIGN DISCORDANTS " << std::endl; 
#endif
	;
      } else {
	tmp_vec.push_back(j);
      }
    }
    reads = tmp_vec;
    

  }


  // calculate the mate region
  if (get_mate_regions && m_region.size() /* don't get mate regions if reading whole bam */) {
    calculateMateRegions();
  }

#ifdef QNAME
  for (auto& j : reads)
      if (j.Qname() == QNAME && (j.AlignmentFlag() == QFLAG || QFLAG == -1)) 
	std::cerr << " read KEPT FINAL " << j << std::endl;
#endif

 
  return bad_regions;
 
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
	  if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1)) {
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
      
      if (r.MateChrID() > 22 || r.GetIntTag("DD") < 0) // no Y or M
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

      // skip if it is unmapped or is secondary alignment
      if (r.MatePosition() < 0 || r.SecondaryFlag())
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

  // keep it if it has indel or unmapped read
  if (r.MaxDeletionBases() || r.MaxInsertionBases() || !r.InsertSize()) // || !r.NumClip())
    return false;
  
  // toss if isize is basically read length (completely overlaping)
  if (std::abs(r.FullInsertSize() - r.Length()) < 5)
    return true;

  // now only consider reads where clip can be explained by isize
  if (!r.NumClip())
    return false;

  // toss it then if isize explans clip
  int exp_ins_size = r.Length() - r.NumClip(); // expected isize if has adapter
  if ((exp_ins_size - 6) < std::abs(r.InsertSize()) && (exp_ins_size + 6) > std::abs(r.InsertSize()))
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

    // if already re-aligned, move on, or if part of unmapped
    if (r.GetIntTag("DD") != 0 || !r.MappedFlag() || !r.MateMappedFlag())
      continue;

    int fi = r.FullInsertSize();

    // if read is big disc or interchromosomal, double check its mapq
    if (fi >= MIN_ISIZE_FOR_DISCORDANT_REALIGNMENT || r.Interchromosomal()) { //r.PairOrientation() != FRORIENTATION) {

      SnowTools::BamReadVector als;
      main_bwa->alignSingleSequence(r.Sequence(), r.Qname(), als, false, 0.60, secondary_cap);

      // no alignments, so label as bad
      if (als.size() == 0) {
	r.AddIntTag("DD", -3);
	continue;
      }

      // discordant alignments
      SnowTools::GenomicRegion gr  = r.asGenomicRegion();
      SnowTools::GenomicRegion grm = r.asGenomicRegionMate();
      grm.pad(DISC_REALIGN_MATE_PAD);

      if (als.size() >= 1) {
	bool has_orig = false;
	bool maps_near_mate = false;
	for (auto& i : als) {

#ifdef QNAME
	  if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1))
	    std::cerr << " REALIGN HIT " << i << " OF TOTAL " << als.size() << std::endl;
#endif
	  
	  // if the realignment has a better mapq at differnt location
	  // then take that, and say that it is too uncertain to be a 
	  // reliable alignment
	  if (i.MapQuality() > r.MapQuality() && gr.getOverlap(i.asGenomicRegion())) {
	    i.SetChrIDMate(r.MateChrID());
	    i.SetPositionMate(r.MatePosition());
	    i.SetPairMappedFlag();
	    if (r.MateReverseFlag())
	      i.SetMateReverseFlag();
	    i.AddZTag("SR", r.GetZTag("SR"));
	    r = i; // reassign the read
	    r.AddIntTag("DD", -5); // read is re-assigned, so too worrisome for discordant
	    continue;
	  }
	
	  // if there is another alignment that overlaps with the original
	  // but has a lower MAPQ, take the new mapq
	  if (gr.getOverlap(i.asGenomicRegion())) {
	    has_orig = true;
	    r.SetMapQuality(std::min(i.MapQuality(), r.MapQuality()));
	  }


#ifdef QNAME
	  if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1)) {
	    std::cerr << " CHECKING AGAINST MATE? " << (r.MateMappedFlag() && (fi >= MIN_ISIZE_FOR_DISCORDANT_REALIGNMENT || r.Interchromosomal())) << std::endl;
	    std::cerr << " HAS OVERLAP WITH MATE REGION? " << " GR " << i.asGenomicRegion() << " MATE " << grm << " OVERLAP " << (grm.getOverlap(i.asGenomicRegion())) << " overlap count  " << r.OverlappingCoverage(i)  << std::endl;
	  }
#endif
	  
	  // if mate is mapped, and read maps to near mate region wo clips, it's not disc
	  if (grm.getOverlap(i.asGenomicRegion())) {
	    // if not clipped or overlaps same covered region as original, 
	    // then it has a secondary mapping
	    if (r.NumClip() < 20 || r.OverlappingCoverage(i) >= 20) 
	      maps_near_mate = true;
#ifdef QNAME	      
	    if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1))
	      std::cerr << " FOUND MATCH NEAR MATE. MATCH ALIGNMENT IS " << i << std::endl;
#endif
	  }
	}
	
	// add tags
	if (!has_orig) 
	  r.AddIntTag("DD", -1);
	else if (maps_near_mate) 
	  r.AddIntTag("DD", -2);
	else
	  r.AddIntTag("DD", als.size());
	
      }
    }
    
  }

  // remove discordant reads without re-creatable mappings, 
  // or with better mappings at mate region
  // ACTUALLY lets keep them all, but don't call them discordantx
  for (auto& r : reads) {
    //if (r.GetIntTag("DD") >= 0)
    brv.push_back(r);
    //else
    //bad_qnames.insert(r.Qname());
#ifdef QNAME    
    if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1)) 
      std::cerr << " AFTER DISC REALIGN " << r << std::endl;
#endif
  }
  reads = brv;

}
