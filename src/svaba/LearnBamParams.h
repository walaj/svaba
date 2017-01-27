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
  
  void learnParams(BamParams& p, int max_count);

  void learnParams(BamParamsMap& p, int max_count);

 private:
  std::string bam;

  void process_read(const SeqLib::BamRecord& r, size_t count, 
		    BamParamsMap& p, double& pos1, double& pos2, double& chr, int& wid) const;

};



#endif
