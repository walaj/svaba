#ifndef SNOWMAN_SNOWMAN_BAM_WALKER_H__
#define SNOWMAN_SNOWMAN_BAM_WALKER_H__

#include <vector>
#include <unordered_map>

#include "SnowTools/BamWalker.h"
#include "SnowTools/GenomicRegionCollection.h"
#include "SnowTools/GenomicRegion.h"
#include "SnowTools/HTSTools.h"

#include "Coverage.cpp"
#include "htslib/khash.h"

typedef std::unordered_map<std::string, size_t> CigarMap;
using SnowTools::GRC;

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

  void addCigar(Read &r);

  bool checkIfDuplicate(Read &r);

  void subSampleToWeirdCoverage(double max_coverage);

  void calculateMateRegions();

  void removeRepeats();

  bool get_coverage = true;
  bool get_mate_regions = true;

  ReadVec reads;

  Coverage cov, weird_cov;
  
  CigarMap cigmap;
  
  size_t max_weird_cov = 100;

  MateRegionVector mate_regions;

 private:

  // might want these in case we are looking for duplicates
  std::unordered_map<std::string, bool> name_map;
  std::unordered_map<std::string, bool> seq_map;

  std::string prefix = ""; // eg. tumor, normal

  size_t m_limit = 0;

  uint32_t m_seed = 1337;
  
  size_t max_cov = 100;
    
};

#endif
