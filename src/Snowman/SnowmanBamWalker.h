#ifndef SNOWMAN_SNOWMAN_BAM_WALKER_H__
#define SNOWMAN_SNOWMAN_BAM_WALKER_H__

#include <vector>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include "htslib/khash.h"

#include "SeqLib/BamReader.h"
#include "SeqLib/ReadFilter.h"
//#include "SnowTools/Fractions.h"
#include "STCoverage.h"
#include "SeqLib/BWAWrapper.h"


#include "DiscordantRealigner.h"

#include "KmerFilter.h"

#define GET_COVERAGE 1

//typedef std::unordered_map<std::string, size_t> CigarMap;

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

  SeqLib::BWAWrapper * main_bwa = nullptr;
  int readlen;  
  std::string prefix; // eg. tumor, normal
  SeqLib::GRC blacklist;
  
  SeqLib::GRC readBam(std::ofstream* log = nullptr);

  void realignDiscordants(SeqLib::BamRecordVector& reads);
  
  void filterMicrobial(SeqLib::BWAWrapper * b);
  
  void KmerCorrect();
  
  bool hasAdapter(const SeqLib::BamRecord& r) const;
  
  bool addCigar(SeqLib::BamRecord &r);
  
  bool isDuplicate(const SeqLib::BamRecord &r);
  
  void subSampleToWeirdCoverage(double max_coverage);
  
  void calculateMateRegions();
  
  void removeRepeats();
  
  bool get_coverage = true;
  bool get_mate_regions = true;
 
  // bad qnames, as found from discordant realignment
  std::unordered_set<std::string> bad_qnames;
  
  //ReadVec reads;
  SeqLib::BamRecordVector reads;

  SeqLib::BamRecordVector all_reads;

  std::vector<char*> all_seqs;

  //std::unordered_map<uint32_t, size_t> cig_pos;

  // cov is the all-read coverage tracker
  // weird-cov just tracks coverage of accepted (clip, disc, etc reads)
  STCoverage cov, weird_cov;
  
  SeqLib::CigarMap cigmap;

  MateRegionVector mate_regions;
  
  //ReadCount rc;
  
  SeqLib::GRC * simple_seq;
  
  // object for realigning discordant reads
  DiscordantRealigner dr; 
  
  // maximum coverage of accepted reads, before subsampling
  size_t max_cov = 100;
  
  // should we trim adapter sequences
  bool adapter_trim = true;

  // should we keep reads for learning correction 
  double kmer_subsample = 0.5;
  // should we subsample the learning reads?
  bool do_kmer_filtering = true;

  //size_t m_keep_limit = 0;

  // set a hard limit on how many reads to accept
  size_t m_limit = 0;

  //size_t m_num_reads_kept = 0;

  // set a read filter
  SeqLib::Filter::ReadFilterCollection * m_mr;

 private:
  
  // private string stream. Initialize once here for speed

  
  // might want these in case we are looking for duplicates
  std::unordered_set<std::string> seq_set;
  
  uint32_t m_seed = 1337;
  
};

#endif
