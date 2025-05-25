#ifndef SVABA_LEARN_BAM_PARAMS_H__
#define SVABA_LEARN_BAM_PARAMS_H__

#include <string>
#include <ostream>
#include <unordered_map>
#include <vector>

#include "SeqLib/BamRecord.h"

struct BamParams {
  
  BamParams() {}

  BamParams(int r, double fc, double fd, double fb, double m) : readlen(r), frac_clip(fc), frac_disc(fd), frac_bad(fb), max_mapq(m) {}

  BamParams(const std::string& rg) : read_group(rg) {}

  void collectStats();

  friend std::ostream& operator<<(std::ostream& out, const BamParams& p);
  
  int visited = 0;
  int num_clip = 0;
  int num_disc = 0;
  int num_bad = 0;

  std::vector<int> isize_vec;

  int lp = 0; // 0.025%
  int hp = 0; // 97.5%

  int readlen = 0;
  double frac_clip = 0;
  double frac_disc = 0;
  double frac_bad = 0;
  int max_mapq = 0;
  double mean_cov = 0;

  double mean_isize = 0;
  double median_isize = 0;
  double sd_isize = 0; 

  std::string read_group;
  
};

typedef std::unordered_map<std::string, BamParams> BamParamsMap;

class LearnBamParams {

 public:
  LearnBamParams(const std::string& b) : bam(b) { };
  
  void learnParams(BamParamsMap& p);

  // universal limit to number reads to learn from per RG
  size_t per_rg_limit = 1000000;
  
 private:

  std::string bam;
  
  void process_read(const SeqLib::BamRecord& r, BamParamsMap& p,
		    size_t& satisfied); 

  std::unordered_map<std::string, size_t> rg_counts;
  size_t num_reads_seen = 0;
};



#endif
