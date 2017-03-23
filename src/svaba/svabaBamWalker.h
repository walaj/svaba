#ifndef SVABA_BAM_WALKER_H__
#define SVABA_BAM_WALKER_H__

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

class svabaBamWalker: public SeqLib::BamReader {
  
 public:
  
  svabaBamWalker() {}

  // for discordant read realignments
  SeqLib::BWAWrapper * main_bwa = nullptr;

  // for setting the SR tag
  std::string prefix; // eg. tumor, normal

  // regions to blacklist
  SeqLib::GRC blacklist;

  // read in the reads
  SeqLib::GRC readBam(std::ofstream* log = nullptr);

  // clear it out
  void clear() { 
    cov.clear();
    cigmap.clear();
    weird_cov.clear();
    mate_regions.clear();
    reads.clear();
    get_coverage = true;
    get_mate_regions = true;
    all_seqs.clear();
    seq_set.clear();
    bad_discordant.clear();
  }

  void realignDiscordants(svabaReadVector& reads);
  
  bool hasAdapter(const SeqLib::BamRecord& r) const;
  
  void addCigar(SeqLib::BamRecord &r);
  
  bool isDuplicate(const SeqLib::BamRecord &r);
  
  void subSampleToWeirdCoverage(double max_coverage);
  
  void calculateMateRegions();

  // should we store the mate regions?
  bool get_mate_regions = true;

  // place to store reads when we get them
  svabaReadVector reads; //c

  // store raw sequences for kmer correction learning
  std::vector<char*> all_seqs; //c

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

  // filter out reads at simple repeats?
  SeqLib::GRC * simple_seq;
  
  // object for realigning discordant reads
  DiscordantRealigner dr; //c
  
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
  SeqPointer<SeqLib::BFC> bfc;
  //  SeqLib::BFC * bfc = nullptr;

 private:

  // might want these in case we are looking for duplicates
  std::unordered_set<std::string> seq_set; //c

  // keep track of which reads were flagged for being bad discordant
  std::unordered_set<std::string> bad_discordant; //c

  // seed for the kmer-learning subsampling
  uint32_t m_seed = 1337;

  // quality trim the readd
  void QualityTrimRead(svabaRead& r) const;
  
};

#endif
