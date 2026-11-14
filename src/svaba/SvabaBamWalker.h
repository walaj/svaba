#pragma once

#include <vector>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
//#include "malloc.h" //debug

#include "SeqLib/BamReader.h"
#include "SeqLib/ReadFilter.h"
#include "STCoverage.h"
#include "SeqLib/BWAAligner.h"
#include "DiscordantRealigner.h"

#include "SeqLib/BFC.h"

class svabaThreadUnit;
class SvabaSharedConfig;

// storage container for mate lookup-regions
class MateRegion: public SeqLib::GenomicRegion
{
 public:
  MateRegion() {}
  MateRegion (int32_t c, uint32_t p1, uint32_t p2, char s = '*') : SeqLib::GenomicRegion(c, p1, p2, s) {}
  size_t count = 0;// read count
  SeqLib::GenomicRegion partner;

};

typedef SeqLib::GenomicRegionCollection<MateRegion> MateRegionVector;

class svabaBamWalker: public SeqLib::BamReader {
  
 public:
  
  svabaBamWalker(SvabaSharedConfig& sc_);

  // read in the reads
  SeqLib::GRC readBam(svabaThreadUnit& unit); 

  // clear it out
  void clear() { 
    cov.clear();
    weird_cov.clear();
    
    bfc.ClearReads();
    cigmap.clear();
    
    mate_regions.clear();
    reads.clear();
    excluded_reads.clear();

    get_coverage = true;
    get_mate_regions = true;
    
    local_blacklist.clear();
    train_reads.clear();

    //malloc_trim(0);//debug
  }

  void AddBackReadsToCorrect();
  
  //
  void Train();

  void ErrorCorrect();

  void ClearTraining();
  
  // set the id for this bam e.g. t001
  void SetPrefix(std::string_view pref) { prefix_ = pref; }
  
  void RealignDiscordants(svabaThreadUnit& unit);
  
  ///bool hasAdapter(const SeqLib::BamRecord& r) const;
  
  void addCigar(const svabaReadPtr& r);
  
  //bool isDuplicate(const SeqLib::BamRecord &r);
  
  void subSampleToWeirdCoverage(double max_coverage);
  
  void calculateMateRegions();

  void TagDiscordant(svabaReadPtr& r);

  /// Return the per-RG isize discordant cutoff, lazily cached.
  int getIsizeCutoff(const std::string& RG);

  // should we store the mate regions?
  bool get_mate_regions = true;

  // place to store reads when we get them
  svabaReadPtrVector reads; //c

  // SvABA2 somatic-safety net: reads EXCLUDED from assembly/r2c (adapter
  // read-through, blacklist-self) but which may still carry a breakpoint
  // junction kmer. Stored as (qname, raw seq) and scanned post-assembly
  // against each BP's jxn_kmer so normal junction evidence can never be
  // silently lost (a missed normal read => false somatic). Cleared per
  // region in clear(). Bounded by MAX_EXCLUDED_POOL to cap memory in
  // dense regions.
  static constexpr size_t MAX_EXCLUDED_POOL = 100000;
  std::vector<std::pair<std::string, std::string>> excluded_reads;
  void stashExcluded(const svabaRead& r);

  // cov is the all-read coverage tracker
  // weird-cov just tracks coverage of accepted (clip, disc, etc reads)
  // buffered cov is coverage minus first 8 and last 8 bp. Why?
  //    because when looking for variant-supporting reads, we require
  //    alignment of a read to the variant to overlap it by 8 bp
  //    to reduce false-positive alt reads. So we use the buffered
  //    coverage to compare against this buffered alt cov.
  STCoverage cov, weird_cov; //c

  // hash of cigars for indels
  SeqLib::CigarMap cigmap; //c

  // mate regions to lookup
  MateRegionVector mate_regions; //c

  // object for realigning discordant reads
  DiscordantRealigner discordantRealigner; //c
  
  // maximum coverage of accepted reads, before subsampling
  size_t max_cov = 100;

  // should we keep reads for learning correction 
  double kmer_subsample = 0.5;
  
  // should we subsample the learning reads?
  bool do_kmer_filtering = true;

  // should we get the read coverage
  bool get_coverage = true;

  // set a hard limit on how many reads to accept
  size_t m_limit = 0;

  // an extra learned region-specific blacklist
  SeqLib::GRC local_blacklist;
  
  SeqLib::BFC bfc;

  // subsampled non-weird reads kept for BFC training. Populated during
  // readBam, consumed by the pooled BFC in SvabaRegionProcessor.
  SeqLib::UnalignedSequenceVector train_reads;

 private:

  // for setting the SR tag
  std::string prefix_; // eg. tumor, normal

  // might want these in case we are looking for duplicates
  //std::unordered_set<std::string> seq_set; //c

  // seed for the kmer-learning subsampling
  uint32_t m_seed = 1337;

  // for logging to console, options etc
  SvabaSharedConfig& sc;

  // cache RG cutoffs so don't have to recalculate
  std::unordered_map<std::string, int> isize_cutoff_per_rg;
  
  
};
