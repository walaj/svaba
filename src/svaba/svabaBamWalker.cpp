#include "svabaBamWalker.h"
#include "svabaRead.h"
#include "svaba_params.h"
#include "SvabaSharedConfig.h"
#include "svabaOptions.h"
#include "svabaLogger.h"

//#define QNAME "H01PEALXX140819:3:2218:11657:19504"
//#define QFLAG -1

using SeqLib::AddCommas;

#ifdef QNAME
#define DEBUG(msg, read)						\
  { if (read.Qname() == QNAME && (read.AlignmentFlag() == QFLAG || QFLAG == -1)) { std::cerr << (msg) << " read " << (read) << std::endl; } }
#else
#define DEBUG(msg, read)
#endif

//static const std::string ILLUMINA_PE_PRIMER_2p0 = "CAAGCAGAAGACGGCAT";
//static const std::string FWD_ADAPTER_A = "AGATCGGAAGAGC";
//static const std::string FWD_ADAPTER_B = "AGATCGGAAAGCA";
//static const std::string REV_ADAPTER = "GCTCTTCCGATCT";

svabaBamWalker::svabaBamWalker(SvabaSharedConfig& sc_) : sc(sc_) {}

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

/*bool svabaBamWalker::isDuplicate(const SeqLib::BamRecord &r) {
  
  // deduplicate by query-bases / position
  std::string key = r.Sequence() + std::to_string(r.Position()) + "_" + std::to_string(r.MatePosition()); 
  
  // its not already add, insert
  if (!seq_set.count(key)) {
    seq_set.insert(key);
    return false;
  }
  
  return true;
  }*/

SeqLib::GRC svabaBamWalker::readBam() {

  // keep track of regions where we hit the weird-read limit
  SeqLib::GRC bad_regions;
  
  // buffer to store reads from a region
  // if we don't hit the limit, then add them to be assembled
  svabaReadVector read_buffer;
  
  // if we have a LOT of reads (aka whole genome run), keep track for printing
  size_t countr = 0;
  
  // keep track of which region we are in 
  int current_region = this->region_idx_;
  
  // loop the reads
  while ( auto bamrec = this->Next()) {

    // makes a deep copy. So when bamrec goes out at end of this iteration, s is still safe
    svabaRead s{ *bamrec, prefix_ };
    
    // when we move regions, clear the buffer
    if (this->region_idx_ != current_region) {
      
      current_region = region_idx_;
      // cache the reads from the last region into BamWalker cache
      reads.insert(reads.end(), read_buffer.begin(), read_buffer.end());
      read_buffer.clear();
    }
    
    // check if this read passes the rules for potential SV reads
    if (s.DuplicateFlag() ||
	s.QCFailFlag() ||
	s.NumHardClip() ||
	s.CountNBases() ||
	sc.blacklist.CountOverlaps(s.AsGenomicRegion()))
      continue;
    
    // quality score trim read
    s.QualityTrimRead(); // copies sequence into svabaRead.seq_corrected
    
    // if its less than 40, dont even mess with it
    if (s.CorrectedSeqLength() < 40)
      continue;
    
    // this is the main rule check 
    bool rule_pass = sc.mr.isValid(s);
    
    // add to all reads pile for kmer correction
    if (get_coverage) { 
      cov.addRead(s, 0); //, false); 
    }
    
    // add to weird coverage
    if (rule_pass && get_coverage) 
      weird_cov.addRead(s, 0); 
    
    // add to the cigar map for all non-duplicate reads
    if (rule_pass && get_mate_regions) // only add cigar for non-mate regions
      addCigar(s); 
    
    // add all the reads for kmer correction
    if (!sc.opts.ecCorrectType.empty() && sc.opts.ecCorrectType != "0") {
      
      // always include weird reads into correction
      s.train = rule_pass;
      
      // if its not a weird read, include for training purposes
      // but can subsample for speed
      if (!rule_pass) { 
	uint32_t hashed  = __ac_Wang_hash(__ac_X31_hash_string(s.Qname().c_str()));
	if ((double)(hashed&0xffffff) / 0x1000000 <= kmer_subsample) 
	  s.train = true;
      }
      
      // in bfc addsequence, memory is copied. for all_seqs (SGA correction), copy explicitly
      if (s.train) {
	assert(bfc.AddSequence(s.CorrectedSeq(),
			       ""/*r.Qualities().c_str()*/,
			       s.UniqueName())); // for BFC correciton
      }
    }
    
    // if not a weird read, move on
    if (!rule_pass)
      continue;
    
    // message logging
    ++countr;
    if (countr % 10000 == 0 && regions_.size() == 0 && log) {
      sc.logger.log(sc.opts.verbose > 1, false/*TODO*/,
		    "...read in ", AddCommas<size_t>(countr),
		    " weird reads for whole genome read-in. At pos ",
		    s.Brief());
    }
    
    // adding first to a temp buffer, in case hit limit and then
    // just don't add any reads
    read_buffer.push_back(s);
    
    // if hit the limit of reads, log it and try next region
    if (read_buffer.size() > m_limit && m_limit > 0) { 

      std::string regstr = (regions_.size() ? regions_[region_idx_].ToString(this->Header()) : " whole BAM");
      sc.logger.log(sc.opts.verbose > 1,
		    false /*TODO*/,
		    "\tstopping read lookup at ",
		    s.Brief(), " in window ", regstr, 
		    "with ", AddCommas(read_buffer.size()),
		    " weird reads. Limit: ",
		    AddCommas(m_limit));
      
      // add this region to the bad regions list
      if (regions_.size())
	bad_regions.add(regions_[region_idx_]);
      
      // clear these reads out, not a good region
      read_buffer.clear();
      
      // force it to try the next region, or return if none left
      ++region_idx_; // increment to next region
      if (region_idx_ >= regions_.size()) {/// no more regions left and last one was bad
	return bad_regions; // totally done
      } 
      break;
    }
  } // end the read loop
  
    // "reads" is total cache for this walker.
    // "read_buffer" was temporary cache in readBam for last region
  reads.insert(reads.end(), read_buffer.begin(), read_buffer.end());
  
  // subsamples reads if too high of coverage
  if (get_coverage)
    subSampleToWeirdCoverage(max_cov);
  
  // realign discordants
  realignDiscordants(reads);
  
  // calculate the mate region
  if (get_mate_regions &&
      regions_.size() /* don't get mate regions if reading whole bam */) 
    calculateMateRegions();
  
  read_buffer.clear(); // should do this, but be extra sure
  return bad_regions;
  
}

void svabaBamWalker::AddBackReadsToCorrect() {
  for (const auto& r : reads) {
    bfc.AddSequence(r.CorrectedSeq(), "", r.UniqueName());
  }
}

void svabaBamWalker::ClearTraining() {
  bfc.ClearReads();
}

void svabaBamWalker::Train() {
  bfc.Train();
}

void svabaBamWalker::ErrorCorrect() {
  bfc.ErrorCorrect();
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

  assert(regions_.size());

  SeqLib::GenomicRegion main_region = regions_.at(0);

  // hold candidate regions. Later trim based on count
  MateRegionVector tmp_mate_regions;

  // loop the reads and add mate reads to MateRegionVector
  for (auto& r : reads) {

    int dd = r.GetDD(); 
    
    // throw away reads that have too many different discordant mappings
    // or are otherwise bad (dd < 0)
    if (!r.PairedFlag() ||
	r.MateChrID() > 22 ||
	r.ChrID() > 22 ||
	!r.MappedFlag() ||
	!r.MateMappedFlag() ||
	dd < 0 ||
	dd > MAX_SECONDARY_HIT_DISC) // no Y or M, no bad discordants
      continue;

    // get the mate region position, will check next if different
    // from current main region
    MateRegion mate(r.MateChrID(), r.MatePosition(), r.MatePosition());
    mate.Pad(MATE_REGION_PAD);
    mate.partner = main_region;
    
    // if mate not in main interval, add a padded version
    if (!main_region.GetOverlap(mate) &&
	r.MapQuality() >= MIN_MAPQ_FOR_MATE_LOOKUP) {
      tmp_mate_regions.add(mate);
    }
    
  }

  // merge it down to get the mate regions
  tmp_mate_regions.MergeOverlappingIntervals();

  // get the counts by overlapping mate-reads
  // with the newly defined mate regions
  for (auto& r : reads) {
    
    SeqLib::GenomicRegion mate(r.MateChrID(), r.MatePosition(), r.MatePosition());

    // if mate not in main interval, check which mate regions it's in
    if (!main_region.GetOverlap(mate) && r.MapQuality() > 0) {

      // loop all of the mate regions -slow?
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

// bool svabaBamWalker::hasAdapter(const SeqLib::BamRecord& r) const {

//   // keep it if it has indel or unmapped read
//   if (r.MaxDeletionBases() || r.MaxInsertionBases() || !r.InsertSize()) // || !r.NumClip())
//     return false;
  
//   // toss if isize is basically read length (completely overlaping)
//   if (std::abs(r.FullInsertSize() - r.Length()) < 5)
//     return true;

//   // now only consider reads where clip can be explained by isize
//   if (!r.NumClip())
//     return false;

//   // toss it then if isize explans clip
//   int exp_ins_size = r.Length() - r.NumClip(); // expected isize if has adapter
//   if ((exp_ins_size - 6) < std::abs(r.InsertSize()) && (exp_ins_size + 6) > std::abs(r.InsertSize()))
//     return true;

//   std::string seq = r.Sequence();
//   if (seq.find(ILLUMINA_PE_PRIMER_2p0) != std::string::npos)
//     return true;

//   return false;

// }

 void svabaBamWalker::realignDiscordants(svabaReadVector& reads) {

   SeqLib::BWAAligner bwa_aligner(sc.bwa_idx);
   
   size_t realigned_count = 0;
   std::unordered_set<std::string> bad_discordant;   
   for (auto& r : reads) {
     
     if (!discordantRealigner.ShouldRealign(r))
       continue;
     
     // modifies the read in place
     if(discordantRealigner.RealignRead(r, bwa_aligner)) {
       ++realigned_count;
       bad_discordant.insert(r.Qname());
     }
   }
   
   // loop through and set DD < 0 for all read-pairs that have a bad discordant
   for (auto& r : reads)
     if (r.GetDD() >= 0 && bad_discordant.count(r.Qname()))
       r.SetDD(DiscordantRealigner::MATE_BAD_DISC);

 }

 void svabaBamWalker::TagDiscordantReads() {
   
   // only warn about missing information once per read grouip
   //static std::unordered_set<std::string> warned_rgs;

   // find the bamStats
   const auto it = sc.bamStats.find(this->prefix_);
   
   // cache RG cutoffs so don't have to recalculate
   std::unordered_map<std::string, int> isize_cutoff_per_rg;
   
   for (auto& r : reads) {
     
     // if suspicious as discordant, remove from clustering
     if (r.GetDD() < 0)
       continue;

     // has to be mapped and have mostly not-hardclip
     if (!r.PairMappedFlag() || r.NumMatchBases() > r.NumHardClip())
       continue;
     
     // get read group
     std::string RG;
     if (!r.GetZTag("RG", RG))
       RG = "NA";
     
     auto cc = isize_cutoff_per_rg.find(RG);
     if (cc == isize_cutoff_per_rg.end()) {
       
       // check that we learned this bam
       if (it == sc.bamStats.end()) {
	 sc.logger.log(true, true, "SBW: clusterDiscordantReads - Can't find",
		       "LearnBamParams for ", this->prefix_);
	 isize_cutoff_per_rg[RG] = DEFAULT_ISIZE_THRESHOLD;
       } 

       else {
	 const auto rgg = it->second.bam_read_groups.find(RG);
	 // check that we learned this read group
	 if (rgg == it->second.bam_read_groups.end()) {
	   sc.logger.log(true, true, "SBW: clusterDiscordantReads - Can't find",
		       "leared read group: ", RG, " for bam ",
			 prefix_);
	   isize_cutoff_per_rg[RG] = DEFAULT_ISIZE_THRESHOLD;
	 } else {
	   isize_cutoff_per_rg[RG] =
	     rgg->second.isize_mean + rgg->second.sd_isize * sc.opts.sdDiscCutoff;
	 }

	 // re-check for calculated cutoff now that we placed it above
	 cc = isize_cutoff_per_rg.find(RG);
	 assert(cc != isize_cutoff_per_rg.end());
	 
       }
     } // end isize cutoff calculation
     
     // accept as discordant if not FR, has large enough isize, is inter-chromosomal, 
     // and has both mates mapping. Also dont cluster on weird chr
     bool weird_orientation = r.PairOrientation() != FRORIENTATION;
     if ( weird_orientation ||
	  abs(r.FullInsertSize()) >= cc->second || // cc-> second is isize cutoff
	  r.Interchromosomal()
	  )
       r.dd = 1; // set as discordant
     
   } // end read loop
 }
 
