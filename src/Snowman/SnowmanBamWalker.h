#ifndef SNOWMAN_SNOWMAN_BAM_WALKER_H__
#define SNOWMAN_SNOWMAN_BAM_WALKER_H__

#include <vector>
#include <unordered_map>

#include "KmerFilter.h"
#include "SnowTools/BamWalker.h"
#include "SnowTools/STCoverage.h"
#include "SnowTools/BWAWrapper.h"

#include "htslib/khash.h"

#define GET_COVERAGE 1

typedef std::unordered_map<std::string, size_t> CigarMap;
using SnowTools::GRC;
using SnowTools::BamReadVector;
using SnowTools::STCoverage;
using SnowTools::BamRead;

class MateRegion: public SnowTools::GenomicRegion
{
 public:
  MateRegion() {}
  MateRegion (int32_t c, uint32_t p1, uint32_t p2, bool s=true) : SnowTools::GenomicRegion(c, p1, p2, s) {}
  size_t count = 0;// read count

};

typedef SnowTools::GenomicRegionCollection<MateRegion> MateRegionVector;

class SnowmanBamWalker: public SnowTools::BamWalker
{
  
 public:

  SnowmanBamWalker() {}

  SnowmanBamWalker(const std::string& in) : SnowTools::BamWalker(in) {}

    void readBam();

    void filterMicrobial(SnowTools::BWAWrapper * b);

  void KmerCorrect();

  void addCigar(BamRead &r);

  bool checkIfDuplicate(BamRead &r);

  void subSampleToWeirdCoverage(double max_coverage);

  void calculateMateRegions();

  void removeRepeats();

  bool get_coverage = true;
  bool get_mate_regions = true;

  SnowTools::GRC blacklist;

  //ReadVec reads;
  BamReadVector reads;

  STCoverage cov, weird_cov;

  CigarMap cigmap;
  
  //SnowTools::GenomicRegion coverage_region;

  size_t max_weird_cov = 100;

  MateRegionVector mate_regions;

  SnowTools::ReadCount rc;

  std::string prefix = ""; // eg. tumor, normal

  size_t max_cov = 100;

  bool do_kmer_filtering = true;

  bool disc_only = false;
 private:

  // might want these in case we are looking for duplicates
  std::unordered_map<std::string, bool> name_map;
  std::unordered_map<std::string, bool> seq_map;

  size_t m_limit = 0;

  uint32_t m_seed = 1337;



};

#endif
