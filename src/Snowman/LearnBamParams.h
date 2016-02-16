#ifndef SNOWMAN_LEARN_BAM_PARAMS_H__
#define SNOWMAN_LEARN_BAM_PARAMS_H__

#include <string>
#include <ostream>

struct BamParams {
  
  BamParams() : readlen(0), frac_clip(0), frac_disc(0), frac_bad(0), max_mapq(0) {}
  BamParams(int r, double fc, double fd, double fb, double m) : readlen(r), frac_clip(fc), frac_disc(fd), frac_bad(fb), max_mapq(m) {}

  friend std::ostream& operator<<(std::ostream& out, const BamParams& p);
  
  int readlen;
  double frac_clip;
  double frac_disc;
  double frac_bad;
  int max_mapq;
  double mean_cov;

};

class LearnBamParams {

 public:
 LearnBamParams(const std::string& b) : bam(b) { };
  
  void learnParams(BamParams& p, int max_count);

 private:
  std::string bam;

};

#endif
