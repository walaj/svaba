#include "SvabaBamWalker.h"
#include "SvabaRead.h"
#include "SvabaSharedConfig.h"
#include "SvabaOptions.h"
#include "SvabaLogger.h"
#include "SvabaThreadUnit.h"
#include "SvabaDebug.h"
#include <string_view>

#define QNAME "LH00306:129:227V5CLT4:6:1268:41599:27601"
#define QFLAG 147

using SeqLib::AddCommas;
using SeqLib::UnalignedSequenceVector;
using SeqLib::UnalignedSequence;

inline void debug_read(std::string_view msg,
                       const svabaReadPtr read,
                       std::string_view qname_filter = {},
                       int qflag_filter = -1) {

    if (!read) return;
    if (!qname_filter.empty() && read->Qname() != qname_filter) return;
    if (qflag_filter != -1 && read->AlignmentFlag() != qflag_filter) return;

    std::cerr << "DEBUG: " << msg << " read " << read << '\n';
}

static const std::string ILLUMINA_PE_PRIMER_2p0 = "CAAGCAGAAGACGGCAT";

// Detect adapter read-through: the fragment is shorter than the read
// length, so the tail of the read sequences into the adapter.  Returns
// true if the read looks like adapter contamination and should be
// skipped.  Reads with indels or unmapped mates are exempted — those
// might genuinely support an SV even if the geometry looks adapter-ish.
static bool hasAdapter(const svabaRead& r) {

  // keep reads with indels or unmapped mate — could be real SV signal
  if (r.MaxDeletionBases() || r.MaxInsertionBases() || !r.InsertSize())
    return false;

  // toss if insert size ≈ read length (completely overlapping pair)
  if (std::abs(r.FullInsertSize() - r.Length()) < 5)
    return true;

  // from here, only consider reads that have soft-clips
  if (!r.NumClip())
    return false;

  // toss if insert size explains the clip length (adapter read-through)
  int exp_ins_size = r.Length() - r.NumClip();
  if (std::abs(exp_ins_size - std::abs(r.InsertSize())) <= 6)
    return true;

  // toss if literal Illumina adapter sequence is present
  if (r.Sequence().find(ILLUMINA_PE_PRIMER_2p0) != std::string::npos)
    return true;

  return false;
}
//static const std::string REV_ADAPTER = "GCTCTTCCGATCT";

svabaBamWalker::svabaBamWalker(SvabaSharedConfig& sc_) : sc(sc_) {}

void svabaBamWalker::addCigar(const svabaReadPtr& r) {
  
  if (r->CigarSize() == 1)
    return;

  // NB: the -1 is NOT a 0- vs 1-based thing.
  // It's that VCF (and BP) store location of an indel as the base PRIOR to the event.
  // but accumulating positions here has "pos" starting at the first base of the next
  // position - which would be the first base of the indel, so we minus 1 to keep
  // convention c/w BPs file. It's all 0-based as usual
  int pos = r->Position() - 1;  // position ON REFERENCE

  for (const auto& i : r->GetCigar()) {
    if (i.Type() == 'D' || i.Type() == 'I') {
      std::string key = std::to_string(r->ChrID()) + "_" +
	std::to_string(pos) + "_" +
	std::to_string(i.Length()) + i.Type();
      ++cigmap[key];
      // //debug
      // if (pos > 79347219 & pos < 79347403)
      // 	std::cerr << r->Qname() << " pos " << r->Position() << " pos " << pos <<
      // 	  " " << i.Length() << i.Type() << std::endl;
    }
    
    // advance along the reference for alignment types that consume reference
    if (i.Type() != 'I' && i.Type() != 'S' && i.Type() != 'H')
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

SeqLib::GRC svabaBamWalker::readBam(svabaThreadUnit& unit) {

  // keep track of regions where we hit the weird-read limit
  SeqLib::GRC bad_regions;
  
  // buffer to store reads from a region
  // if we don't hit the limit, then add them to be assembled
  svabaReadPtrVector read_buffer;
  UnalignedSequenceVector train_buffer;
  train_reads.clear();  // member — cleared per readBam call
  
  // if we have a LOT of reads (aka whole genome run), keep track for printing
  size_t countr = 0;
  size_t seen = 0;
  
  // keep track of which region we are in 
  int current_region = this->region_idx_;

  // should just skip regions that have excessive coverage
  // cov = readlength + reads / region_width. So 150 * X / 25000 = 10000 to set max coverage to 10000
  // so 25'000 * 10'000 / 150 = 1'666'667 for max cov 10'000
  const size_t max_coverage_before_ignore_region = 10'000;
  size_t seen_cutoff = regions_.at(current_region).Width() * max_coverage_before_ignore_region / std::max(1, sc.readlen);
  sc.logger.log(sc.opts.verbose > 4, false, "...seen cutoff ", SeqLib::AddCommas(seen_cutoff),
		" region width ", regions_.at(current_region), " read len ", sc.readlen);
  
  // loop the reads
  while ( auto bamrec = this->Next()) {
    
    // makes a deep copy. So when bamrec goes out at end of this iteration, s is still safe
    svabaReadPtr s = std::make_shared<svabaRead>(*bamrec, prefix_);

    ++seen;
    if (seen % 100'000 == 0)
      sc.logger.log(sc.opts.verbose > 3, false, "...bamwlker reading read ",
		    s->Brief(), " - ",
		    SeqLib::AddCommas(seen));
    
    // when we move regions, clear the buffer
    if (this->region_idx_ != current_region) {

      seen = 0;
      sc.logger.log(sc.opts.verbose > 3, false, "...svabaBamWalker switched region from ",
		    regions_.at(current_region), " to ", regions_.at(region_idx_));
      
      current_region = region_idx_;
      
      // recalculate seen limit for this region
      if (current_region < regions_.size()) {
	seen_cutoff = regions_.at(current_region).Width() * max_coverage_before_ignore_region / std::max(1, sc.readlen);
	sc.logger.log(sc.opts.verbose > 4, false, "...seen cutoff ", SeqLib::AddCommas(seen_cutoff));	
      } else { // shouldn't really get here
	seen_cutoff = 0;
      }
      
      // cache the reads from the last region into BamWalker cache
      reads.insert(reads.end(), read_buffer.begin(), read_buffer.end());
      read_buffer.clear();
    }

    // if hit the limit of reads, log it and try next region
    if ( (read_buffer.size() > m_limit && m_limit > 0) ||
	 seen > seen_cutoff) { 
      
      // don't train on these
      train_buffer.clear();

      // clear the buffer for this region
      seen = 0;
      
      std::string regstr =
	(regions_.size() ? regions_[region_idx_].ToString(this->Header()) : " whole BAM");
      sc.logger.log(sc.opts.verbose > 1,
		    sc.opts.verbose_log,
		    "\tstopping read lookup at ",
		    s->Brief(), " in window ", regstr, 
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
    }

    SVABA_READ_TRACE(s->Qname(),
      "SEEN flag=" << s->AlignmentFlag()
      << " mapq=" << s->MapQuality()
      << " pos=" << s->ChrName(sc.header) << ":" << s->Position()
      << " isize=" << s->InsertSize()
      << " cigar=" << s->CigarString()
      << " len=" << s->Length());

    // check if this read passes the rules for potential SV reads
    if (s->DuplicateFlag() ||
	s->QCFailFlag() ||
	s->NumHardClip() ||
	sc.blacklist.CountOverlaps(s->AsGenomicRegionMate()) ||
	sc.blacklist.CountOverlaps(s->AsGenomicRegion()))
      {
	SVABA_READ_TRACE(s->Qname(),
	  "SKIP dup=" << s->DuplicateFlag()
	  << " qcfail=" << s->QCFailFlag()
	  << " hardclip=" << s->NumHardClip()
	  << " blacklist_mate=" << sc.blacklist.CountOverlaps(s->AsGenomicRegionMate())
	  << " blacklist_self=" << sc.blacklist.CountOverlaps(s->AsGenomicRegion()));
	continue;
      }

    // quality score trim read
    s->QualityTrimRead(); // copies sequence into svabaRead.seq_corrected

    SVABA_READ_TRACE(s->Qname(),
      "AFTER_TRIM corrected_len=" << s->CorrectedSeqLength());

    // this is the main rule check
    bool rule_pass = sc.mr.isValid(*s);

    SVABA_READ_TRACE(s->Qname(),
      "RULE_CHECK rule_pass=" << rule_pass);

    // special case to also check against high mismatch
    // --no-nm disables this path entirely for maximum efficiency
    if (!rule_pass && !sc.opts.noNmSalvage) {
      int32_t nm_tag = 0;
      s->GetIntTag("NM", nm_tag);
      float nm_frac = static_cast<float>(nm_tag) / static_cast<float>(s->Length());
      bool is_adapter = hasAdapter(*s);
      SVABA_READ_TRACE(s->Qname(),
        "NM_SALVAGE nm=" << nm_tag
        << " nm_frac=" << nm_frac
        << " adapter=" << is_adapter
        << " pass=" << (nm_frac > 0.02 && !is_adapter));
      if (nm_frac > 0.02 && !is_adapter) {
	s->to_assemble = false;
	read_buffer.push_back(s);
	SVABA_READ_TRACE(s->Qname(), "ADDED via NM salvage (to_assemble=false)");
      }
    }

    // if its less than 40, dont even mess with it
    if (s->CorrectedSeqLength() < 40) {
      SVABA_READ_TRACE_IF(rule_pass, s->Qname(),
        "REJECT corrected_len=" << s->CorrectedSeqLength() << " < 40");
      rule_pass = false;
    }

    // special rule to filter "RF" "discordants" that are just overlaps
    if (s->PairOrientation() == SeqLib::Orientation::RF &&
	std::abs(s->InsertSize()) < std::max(200, sc.readlen * 2)) {
      SVABA_READ_TRACE_IF(rule_pass, s->Qname(),
        "REJECT RF overlap isize=" << s->InsertSize()
        << " threshold=" << std::max(200, sc.readlen * 2));
      rule_pass = false;
    }

    // for mate regions, be more strict
    if (!get_mate_regions && s->MapQuality() == 0) {
      SVABA_READ_TRACE_IF(rule_pass, s->Qname(),
        "REJECT mapq=0 in non-mate region");
      rule_pass = false;
    }

    // add to all reads pile for kmer correction
    if (get_coverage) {
      cov.addRead(*s, 0);
    }

    // add to weird coverage
    if (rule_pass && get_coverage)
      weird_cov.addRead(*s, 0);

    // add to the cigar map for all non-duplicate reads
    if (get_mate_regions) // only add cigar for non-mate regions
      addCigar(s);

    // add all the reads for kmer correction
    // TODO this will add reads to training even if this entire buffer is rejected, so need to make a training buffer as well before adding them
    if (!sc.opts.ecCorrectType.empty() && sc.opts.ecCorrectType != "0") {

      // if its not a weird read, include for training purposes
      // but can subsample for speed
      if (!rule_pass) {
	uint32_t hashed  = __ac_Wang_hash(__ac_X31_hash_string(s->Qname().c_str()));
	if ((double)(hashed&0xffffff) / 0x1000000 <= kmer_subsample)
	  train_buffer.push_back(UnalignedSequence(s->UniqueName(),
							   s->CorrectedSeq()));
      }
    }

    // if not a weird read, move on
    if (!rule_pass) {
      SVABA_READ_TRACE(s->Qname(), "SKIP not weird (rule_pass=false after all gates)");
      continue;
    }

    // skip adapter read-through (fragment shorter than read length)
    // unless the read has an indel that could indicate real SV signal
    if (hasAdapter(*s)) {
      SVABA_READ_TRACE(s->Qname(),
        "SKIP adapter isize=" << s->InsertSize()
        << " len=" << s->Length()
        << " nclip=" << s->NumClip());
      continue;
    }

    // label as discordant or not
    // modifies r in place
    TagDiscordant(s);

    // message logging
    ++countr;
    if (countr % 25000 == 0) {
      sc.logger.log(sc.opts.verbose > 1, false,
		    "...read in ", AddCommas<size_t>(countr),
		    " weird reads for whole genome read-in. At pos ",
		    s->Brief());
    }

    // adding first to a temp buffer, in case hit limit and then
    // just don't add any reads
    read_buffer.push_back(s);
    SVABA_READ_TRACE(s->Qname(),
      "ADDED as weird read (to_assemble=true) buffer_size=" << read_buffer.size());
  } // end the read loop
  
    // "reads" is total cache for this walker.

    // "read_buffer" was temporary cache in readBam for last region
  /*reads.insert(reads.end(),
    std::make_move_iterator(read_buffer.begin()),
    std::make_move_iterator(read_buffer.end()));
    read_buffer.clear();  */ // shouldn't need to, but avoid "lingering state"
  reads.insert(reads.end(), read_buffer.begin(), read_buffer.end());
  read_buffer.clear();
  // stash training reads as member for the pooled BFC in SvabaRegionProcessor
  train_reads.insert(train_reads.end(), train_buffer.begin(), train_buffer.end());
  train_buffer.clear();
  
  // subsamples reads if too high of coverage
  //  if (get_coverage)
  //  subSampleToWeirdCoverage(max_cov);

  // calculate the mate region
  if (get_mate_regions && 
      regions_.size() /* don't get mate regions if reading whole bam */) 
    calculateMateRegions();

  return bad_regions;
  
}

void svabaBamWalker::AddBackReadsToCorrect() {
  for (const auto& r : reads) {
    if (!r->to_assemble) continue;   // same filter as Train() side
    bfc.AddSequence(r->CorrectedSeq(), "", r->UniqueName());
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
  auto keep = [&](const std::shared_ptr<svabaRead>& r) {
    double cov1 = weird_cov.getCoverageAtPosition(r->ChrID(), r->Position());
    double cov2 = weird_cov.getCoverageAtPosition(r->ChrID(), r->PositionEnd());
    double this_cov = std::max(cov1, cov2);
    
    if (this_cov <= max_coverage)
      return true;
    
    double sample_rate = 1 - (this_cov - max_coverage) / this_cov;
    uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(r->Qname().c_str()) ^ m_seed);
    return ((double)(k & 0xffffff) / 0x1000000) <= sample_rate;
  };
  
  reads.erase(std::remove_if(reads.begin(), reads.end(),
                             [&](const std::shared_ptr<svabaRead>& r) { return !keep(r); }),
              reads.end());
}

void svabaBamWalker::calculateMateRegions() {

  assert(regions_.size());

  SeqLib::GenomicRegion main_region = regions_.at(0);

  // hold candidate regions. Later trim based on count
  MateRegionVector tmp_mate_regions;

  // loop the reads and add mate reads to MateRegionVector
  for (const auto& r : reads) {
    
    int dd = r->GetDD(); 
    
    // throw away reads that have too many different discordant mappings
    // or are otherwise bad (dd < 0)
    if (!r->PairedFlag() ||
	r->MateChrID() > 22 ||
	r->MapQuality() == 0 || 
	r->ChrID() > 22 ||
	!r->MappedFlag() ||
	!r->MateMappedFlag() ||
	dd < 0 ||
	dd > MAX_SECONDARY_HIT_DISC) // no Y or M, no bad discordants
      continue;

    // get the mate region position, will check next if different
    // from current main region
    MateRegion mate(r->MateChrID(), r->MatePosition(), r->MatePosition());
    mate.Pad(MATE_REGION_PAD);
    mate.partner = main_region;
    
    // if mate not in main interval, add a padded version
    if (!main_region.GetOverlap(mate) &&
	r->MapQuality() >= MIN_MAPQ_FOR_MATE_LOOKUP &&
	!local_blacklist.CountOverlaps(mate) && 
	!sc.blacklist.CountOverlaps(mate)) {
      tmp_mate_regions.add(mate);
    }
    
  }

  // merge it down to get the mate regions
  tmp_mate_regions.MergeOverlappingIntervals();



  // get the counts by overlapping mate-reads
  // with the newly defined mate regions
  for (const auto& r : reads) {
    
    SeqLib::GenomicRegion mate(r->MateChrID(), r->MatePosition(), r->MatePosition());

    // if mate not in main interval, check which mate regions it's in
    if (!main_region.GetOverlap(mate) && r->MapQuality() > 0) {

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
  
  // hard cutoff on mate regions TODO make more elegant
  if (mate_regions.size() > 6) {
    sc.logger.log(sc.opts.verbose > 2, true, "...cutting off mate lookups for region ",
		  regions_.at(0), " with too many mate regions of: ", mate_regions.size());
    mate_regions.clear();
  }
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

void svabaBamWalker::RealignDiscordants(svabaThreadUnit& unit) {

   size_t realigned_count = 0;
   size_t realign_attempted = 0;
   std::unordered_set<std::string> bad_discordant;   
   for (auto& r : reads) {
     
     if (!discordantRealigner.ShouldRealign(r))
       continue;

     ++realign_attempted;
     if (realign_attempted % 10000 == 0)
       sc.logger.log(sc.opts.verbose > 3, false, "Realigning read ",
		     r->Brief(), " - ", SeqLib::AddCommas(realign_attempted), " - ",
		     realigned_count);
     
     // modifies the read in place
     if(discordantRealigner.RealignRead(r, unit.bwa_aligner)) {
       ++realigned_count;
       bad_discordant.insert(r->Qname());
     }
   }
   
   // loop through and set DD < 0 for all read-pairs that have a bad discordant
   for (auto& r : reads)
     if (r->GetDD() >= 0 && bad_discordant.count(r->Qname()))
       r->SetDD(DiscordantRealigner::MATE_BAD_DISC);
   
}

void svabaBamWalker::TagDiscordant(svabaReadPtr& r) {
  
  // find the bamStats
  const auto it = sc.bamStats.find(this->prefix_);
  
  // get read group
  std::string RG;
  if (!r->GetZTag("RG", RG))
    RG = "NA";

  // this is a local cache
  auto cc = isize_cutoff_per_rg.find(RG);

  // if not found, go find it 
  if (cc == isize_cutoff_per_rg.end()) {
    
    // check that we learned this bam
    // emit a warning if didn't find this RG
    // by putting it here after caching, emits warnings
    // once per RG. Then get a default size
    if (it == sc.bamStats.end()) {
      if (!sc.warned.count(this->prefix_)) {
	//debug
	sc.logger.log(true, true, "SBW: TagDiscordant - Can't find ",
		      "LearnBamParams for ", this->prefix_);
	sc.warned.insert(this->prefix_);
      }
      isize_cutoff_per_rg[RG] = DEFAULT_ISIZE_THRESHOLD;
      
      // re-check for calculated cutoff now that we placed it above
      cc = isize_cutoff_per_rg.find(RG);
      assert(cc != isize_cutoff_per_rg.end());
    } 
    
    // the bamStats exists for this 
    else {
      const auto rgg = it->second.bam_read_groups.find(RG);
      
      // check that we learned this read group
      if (rgg == it->second.bam_read_groups.end()) {
	if (!sc.warned.count(RG)) {
	  //debug
	  //sc.logger.log(true, true, "SBW: TagDiscordant - Can't find ",
	  //		"read group: ", RG, " for bam ",
	  //		prefix_, " - setting default isize cutoff ", DEFAULT_ISIZE_THRESHOLD);
	  sc.warned.insert(RG);
	} 
	isize_cutoff_per_rg[RG] = DEFAULT_ISIZE_THRESHOLD;

	// RG is found in this bamStats. Cache the cutoff
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
  bool weird_orientation = r->PairOrientation() != SeqLib::Orientation::FR;
  if ( weird_orientation ||
       abs(r->FullInsertSize()) >= cc->second || // cc-> second is isize cutoff
       r->Interchromosomal()
       )
    r->dd = 1; // set as discordant

  // not discordant if read or mate in blacklist
  if (sc.blacklist.CountOverlaps(r->AsGenomicRegion()) ||
      sc.blacklist.CountOverlaps(r->AsGenomicRegionMate()) ||
      local_blacklist.CountOverlaps(r->AsGenomicRegion()) ||
      local_blacklist.CountOverlaps(r->AsGenomicRegionMate()))
    r->dd = 0; 
  
  //debug
  // if (r->Qname() == "LH00306:129:227V5CLT4:6:2114:6074:25724") {
  //   std::cerr << " TAG " << r->dd << " -- " << *r << " weird orientation " <<
  //     weird_orientation << " cutoff " << cc->second << " INS " <<
  //     abs(r->FullInsertSize()) << std::endl;
  // }
      
  return;
}
 
