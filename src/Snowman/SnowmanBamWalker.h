#ifndef SNOWMAN_SNOWMAN_BAM_WALKER_H__
#define SNOWMAN_SNOWMAN_BAM_WALKER_H__

#include <vector>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include "SeqLib/BamReader.h"
#include "SeqLib/ReadFilter.h"
#include "STCoverage.h"
#include "SeqLib/BWAWrapper.h"
#include "DiscordantRealigner.h"

#include "SeqLib/BFC.h"

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

class SnowmanBamWalker: public SeqLib::BamReader {
  
 public:
  
  SnowmanBamWalker() {}

  // for discordant read realignments
  SeqLib::BWAWrapper * main_bwa = nullptr;

  // for setting the SR tag
  std::string prefix; // eg. tumor, normal

  // regions to blacklist
  SeqLib::GRC blacklist;

  // read in the reads
  SeqLib::GRC readBam(std::ofstream* log = nullptr);

  void realignDiscordants(SeqLib::BamRecordVector& reads);
  
  bool hasAdapter(const SeqLib::BamRecord& r) const;
  
  void addCigar(SeqLib::BamRecord &r);
  
  bool isDuplicate(const SeqLib::BamRecord &r);
  
  void subSampleToWeirdCoverage(double max_coverage);
  
  void calculateMateRegions();

  // should we store the mate regions?
  bool get_mate_regions = true;

  // place to store reads when we get them
  SeqLib::BamRecordVector reads;

  // store raw sequences for kmer correction learning
  std::vector<char*> all_seqs;

  // cov is the all-read coverage tracker
  // weird-cov just tracks coverage of accepted (clip, disc, etc reads)
  STCoverage cov, weird_cov;

  // hash of cigars for indels
  SeqLib::CigarMap cigmap;

  // mate regions to lookup
  MateRegionVector mate_regions;

  // filter out reads at simple repeats?
  SeqLib::GRC * simple_seq;
  
  // object for realigning discordant reads
  DiscordantRealigner dr; 
  
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

  // set a read filter
  SeqLib::Filter::ReadFilterCollection * m_mr;

  // 
  SeqLib::BFC * bfc = nullptr;

 private:

  // might want these in case we are looking for duplicates
  std::unordered_set<std::string> seq_set;

  // keep track of which reads were flagged for being bad discordant
  std::unordered_set<std::string> bad_discordant;

  // seed for the kmer-learning subsampling
  uint32_t m_seed = 1337;

  // quality trim the readd
  void QualityTrimRead(SeqLib::BamRecord& r) const;
  
};

#endif
