#include "svabaBamWalker.h"
#include "svabaRead.h"

//#define QNAME "H01PEALXX140819:3:2218:11657:19504"
//#define QFLAG -1

#ifdef QNAME
#define DEBUG(msg, read)				\
  { if (read.Qname() == QNAME && (read.AlignmentFlag() == QFLAG || QFLAG == -1)) { std::cerr << (msg) << " read " << (read) << std::endl; } }
#else
#define DEBUG(msg, read)
#endif

//#define DEBUG_SVABA_BAMWALKER 1
#define MIN_MAPQ_FOR_MATE_LOOKUP 0
//#define TRAIN_READS_FAIL_SAFE 50000
#define MIN_ISIZE_FOR_DISCORDANT_REALIGNMENT 1000
#define DISC_REALIGN_MATE_PAD 100
#define MAX_SECONDARY_HIT_DISC 10
#define MATE_REGION_PAD 250

// trim this many bases from front and back of read when determining coverage
// this should be synced with the split-read buffer in BreakPoint2 for more accurate 
// representation of covearge of INFORMATIVE reads (eg ones that could be split)
#define INFORMATIVE_COVERAGE_BUFFER 4

static const std::string ILLUMINA_PE_PRIMER_2p0 = "CAAGCAGAAGACGGCAT";
static const std::string FWD_ADAPTER_A = "AGATCGGAAGAGC";
static const std::string FWD_ADAPTER_B = "AGATCGGAAAGCA";
static const std::string REV_ADAPTER = "GCTCTTCCGATCT";

void svabaBamWalker::addCigar(SeqLib::BamRecord &r) {

  // this is a 100% match
  if (r.CigarSize() == 1)
    return;
  std::stringstream cigar_ss;
  cigar_ss.str(std::string());
  int pos = r.Position(); // position ON REFERENCE
  
  for (auto& i : r.GetCigar()) {

       // if it's a D or I, add it to the list
      if (i.Type() == 'D' || i.Type() == 'I') {	

	cigar_ss << r.ChrID() << "_" << pos << "_" << i.Length() << i.Type();
	++cigmap[cigar_ss.str()];

	cigar_ss.str(std::string());
      }
      
      // move along the REFERENCE
      if (!(i.Type() == 'I') && !(i.Type() == 'S') && !(i.Type() == 'H'))
	pos += i.Length();
  }
  
}

bool svabaBamWalker::isDuplicate(const SeqLib::BamRecord &r) {

  // deduplicate by query-bases / position
  std::string key = r.QualitySequence() + std::to_string(r.Position()) + "_" + std::to_string(r.MatePosition());

  // its not already add, insert
  if (!seq_set.count(key)) {
    seq_set.insert(key);
    return false;
  }
  
  return true;
}

SeqLib::GRC svabaBamWalker::readBam(std::ofstream * log) {

  // these are setup to only use one bam, so just shortcut it
  SeqLib::_Bam * tb = &m_bams.begin()->second;

  SeqLib::BamRecord r;

  // keep track of regions where we hit the weird-read limit
  SeqLib::GRC bad_regions;

  // if we have a LOT of reads (aka whole genome run), keep track for printing
  size_t countr = 0;

  // loop the reads
  while (GetNextRecord(r)) {

    // check if it passed blacklist
    if (blacklist.size() && blacklist.CountOverlaps(r.AsGenomicRegion())) {
      continue;
    }

    // dont even mess with them
    if (r.CountNBases())
      continue;

    // set some things to check later
    bool is_dup = false;
    bool rule_pass = false;
    bool qcpass = !r.DuplicateFlag() && !r.QCFailFlag();
    bool pass_all = true;

    svabaRead s(r, prefix);

    // quality score trim read. Stores in GV tag
    QualityTrimRead(s);
    
    // if its less than 20, dont even mess with it
    //if (s.QualitySequence().length() < 20)
    if (s.SeqLength() < 20)
      //if (s.QualitySequence().length() < 20)
      continue;

    //debug
    s.AddZTag("GV", s.Seq());

    rule_pass = m_mr->isValid(r);
      
    DEBUG("SBW read seen", r);

    // if hit the limit of reads, log it and try next region
    if (countr > m_limit && m_limit > 0) {

      if (log)
	(*log) << "\tbreaking at " << r.Brief() << " in window " 
	       << (m_region.size() ? m_region[tb->m_region_idx].ToString() : " whole BAM")
	       << " with " << SeqLib::AddCommas(countr) 
	       << " weird reads. Limit: " << SeqLib::AddCommas(m_limit) << std::endl;

      if (m_region.size()) 
	bad_regions.add(m_region[tb->m_region_idx]);

      // clear these reads out
      if ((int)reads.size() - countr > 0)
	reads.erase(reads.begin(), reads.begin() + countr);
      
      // force it to try the next region, or return if none left
      ++tb->m_region_idx; // increment to next region
      if (tb->m_region_idx >= m_region.size()) /// no more regions left
	break;
      else { // move to next region
	countr = 0;
	tb->SetRegion(m_region[tb->m_region_idx]);
	continue;
      }
      break;
    }
    
    // add to all reads pile for kmer correction
    if (qcpass && get_coverage) 
      cov.addRead(r, INFORMATIVE_COVERAGE_BUFFER, false); 
    
    // check if in simple-seq
    if (simple_seq->size()) {
      
      // check simple sequence overlaps
      SeqLib::GRC ovl = simple_seq->FindOverlaps(r.AsGenomicRegion(), true);
      
      int msize = 0;
      for (auto& j: ovl) {
	int nsize = j.Width() - r.MaxDeletionBases() - 1;
	if (nsize > msize && nsize > 0)
	  msize = nsize;
      }
      
      if (msize > 30)
        qcpass = false;
    }
    
    pass_all = pass_all && qcpass && rule_pass;
    
    // check if has adapter
    pass_all = pass_all && !hasAdapter(r);

    // check if duplicated
    is_dup = isDuplicate(r);
    pass_all = !is_dup && pass_all;
    
    // add to weird coverage
    if (pass_all && get_coverage)
      weird_cov.addRead(r, 0, false);
    
    // add to the cigar map for all non-duplicate reads
    if (pass_all && get_mate_regions) // only add cigar for non-mate regions
      addCigar(r);

    DEBUG("SBW read has qcpass?: " + std::to_string(qcpass), r);
    DEBUG("SBW duplicated? " + std::to_string(is_dup), r);
    DEBUG("SBW rule pass? " + std::to_string(rule_pass), r); 
    DEBUG("SBW pass all? " + std::to_string(pass_all), r);
   
    // add all the reads for kmer correction
    if (qcpass && do_kmer_filtering && all_seqs.size() < (m_limit * 5) && qcpass && !r.NumHardClip()) {
      std::string qq = r.QualitySequence();
      bool train = pass_all && qq.length() > 40;
      // if not 
      if (!pass_all) {
	uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(r.Qname().c_str()) ^ m_seed);
	if ((double)(k&0xffffff) / 0x1000000 <= kmer_subsample) 
	  train = true;
      }
      
      // in bfc addsequence, memory is copied. for all_seqs,i copy explicitly
      if (train) {
	if (bfc)
	  bfc->AddSequence(qq.c_str(), r.Qualities().c_str(), r.Qname().c_str()); // for BFC correciton
	else {
	  all_seqs.push_back(strdup(qq.c_str()));
	}
      }

    }
     
    if (!pass_all)
      continue;
    
    ++countr;
    if (countr % 10000 == 0 && m_region.size() == 0 && log)
      (*log) << "...read in " << SeqLib::AddCommas<size_t>(countr) << " weird reads for whole genome read-in. At pos " << r.Brief() << std::endl;
    
    // for memory conservation
    s.RemoveTag("BQ");
    s.RemoveTag("OQ");
    //r.RemoveTag("XT");
    //r.RemoveTag("XA");
    //r.RemoveTag("SA");
    
    // add the ID tag
    //std::string srn =  prefix + "_" + std::to_string(r.AlignmentFlag());// + "_" + r.Qname();
    //r.AddZTag("SR", srn);
    
    DEBUG("SBW read added ", r);

    //
    //std::string g = r.QualitySequence();
    //svabaRead s(r, prefix);
    //s.AddSeq(r.QualitySequence()); // copies the sequence
    //r.ClearSeqQualAndTags();       // also clears same r within svabaRead
    //s.ClearSeqQualAndTags();
    s.SetSequence(std::string()); // clear out the sequence & qual in the htslib
    s.RemoveTag("GV");

    //reads.push_back(r); // adding later because of kmer correction
    reads.push_back(s); // adding later because of kmer correctiona
    
  } // end the read loop

#ifdef QNAME
  for (auto& j : reads) { DEBUG("SBW read kept pre-filter", j); }
#endif
  
  if (reads.size() < 3)
    return bad_regions;
  
  // clean out the buffer
  if (get_coverage)
    subSampleToWeirdCoverage(max_cov);

#ifdef QNAME
  for (auto& j : reads) {   DEBUG("SBW read kept post subsample", j);  }
#endif

  // realign discordants
  if (main_bwa) 
    realignDiscordants(reads);
  
  // calculate the mate region
  if (get_mate_regions && m_region.size() /* don't get mate regions if reading whole bam */) 
    calculateMateRegions();

#ifdef QNAME
  for (auto& j : reads) { DEBUG("SBW read kept FINAL FINAL", j); }
#endif
  
  return bad_regions;
  
}

void svabaBamWalker::subSampleToWeirdCoverage(double max_coverage) {
  
  svabaReadVector new_reads;

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

void svabaBamWalker::calculateMateRegions() {

  assert(m_region.size());

  SeqLib::GenomicRegion main_region = m_region.at(0);

  // hold candidate regions. Later trim based on count
  MateRegionVector tmp_mate_regions;

  // loop the reads and add mate reads to MateRegionVector
  for (auto& r : reads) {

    //int dd = r.GetIntTag("DD");
    int dd = r.GetDD(); 
    
    // throw away reads that have too many different discordant mappings, or 
    // are otherwise bad (dd < 0)
    if (!r.PairedFlag() || r.MateChrID() > 22 || r.ChrID() > 22 || !r.MappedFlag() 
	|| !r.MateMappedFlag() || dd < 0 || dd > MAX_SECONDARY_HIT_DISC) // no Y or M, no bad discordants
      continue;

    MateRegion mate(r.MateChrID(), r.MatePosition(), r.MatePosition());
    mate.Pad(MATE_REGION_PAD);
    mate.partner = main_region;
    
    // if mate not in main interval, add a padded version
    if (!main_region.GetOverlap(mate) && r.MapQuality() >= MIN_MAPQ_FOR_MATE_LOOKUP) {
      tmp_mate_regions.add(mate);
    }
    
  }

  // merge it down to get the mate regions
  tmp_mate_regions.MergeOverlappingIntervals();

  // get the counts by overlapping mate-reads with the newly defined mate regions
  for (auto& r : reads) {
    
    SeqLib::GenomicRegion mate(r.MateChrID(), r.MatePosition(), r.MatePosition());

    // if mate not in main interval, check which mate regions it's in
    if (!main_region.GetOverlap(mate) && r.MapQuality() > 0) {

      for (auto& k : tmp_mate_regions)
	if (k.GetOverlap(mate)) {
	  ++k.count;
	  continue;
	}
    }
  }

#ifdef DEBUG_SVABA_BAMWALKER
  std::cerr << "SBW: Mate regions are" << std::endl;
  for (auto& i : tmp_mate_regions) 
    std::cerr << "    " << i << " count " << i.count << std::endl;
#endif

  // keep only ones with 2+ reads
  for (auto& i : tmp_mate_regions) {
    if (i.count >= 2)
      mate_regions.add(i);
  }
  
#ifdef DEBUG_SVABA_BAMWALKER
  std::cerr << "SBW: Final mate regions are:" << std::endl;
  for (auto& i : mate_regions) 
    std::cerr << "    " << i << " read count " << i.count << std::endl;
#endif
  
}

bool svabaBamWalker::hasAdapter(const SeqLib::BamRecord& r) const {

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

}

//void svabaBamWalker::realignDiscordants(SeqLib::BamRecordVector& reads) {
void svabaBamWalker::realignDiscordants(svabaReadVector& reads) {

  size_t realigned_count = 0;
  for (auto& r : reads) {
    
    if (!dr.ShouldRealign(r))
      continue;

    // modifies the read in place
    if(dr.RealignRead(r, main_bwa)) {
      ++realigned_count;
      bad_discordant.insert(r.Qname());
    }
  }

  // loop through and set DD < 0 for all read-pairs that have a bad discordant
  for (auto& r : reads)
    if (r.GetDD() >= 0 && bad_discordant.count(r.Qname()))
      r.SetDD(DiscordantRealigner::MATE_BAD_DISC);

}

void svabaBamWalker::QualityTrimRead(svabaRead& r) const {

  int32_t startpoint = 0, endpoint = 0;
  r.QualityTrimmedSequence(3, startpoint, endpoint);
  int32_t new_len = endpoint - startpoint;
  if (endpoint != -1 && new_len < r.Length() && new_len > 0 && new_len - startpoint >= 0 && startpoint + new_len <= r.Length()) { 
    try { 
      //r.AddZTag("GV", r.Sequence().substr(startpoint, new_len));
      r.SetSeq(r.Sequence().substr(startpoint, new_len));
      //assert(r.GetZTag("GV").length());
    } catch (...) {
      std::cerr << "Subsequence failure with sequence of length "  
		<< r.Sequence().length() << " and startpoint "
		<< startpoint << " endpoint " << endpoint 
		<< " newlen " << new_len << std::endl;
    }
    // read is fine
  } else {
    //r.AddZTag("GV", r.Sequence());
    r.SetSeq(r.Sequence()); // copies the sequence
  }
  

}
