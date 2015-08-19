#ifndef SNOWMAN_READSIM_H__
#define SNOWMAN_READSIM_H__

#include <string>
#include <vector>
#include "SnowTools/GenomicRegion.h"

struct Indel {
  
  Indel() : len(0), type('N'), pos(0) {}

  Indel(int l, char t, uint32_t p) : len(l), type(t), pos(p) {}

  int len;
  char type;
  uint32_t pos;
  SnowTools::GenomicRegion gr;
  int frag_id;

  friend std::ostream& operator<<(std::ostream& out, const Indel& i);


};

class ReadSim {


 public:

  ReadSim() : m_error(0), m_imean(350), m_isd(50) {}

  /** Give the sequene string SNP errors at rate of er */
  void makeSNVErrors(std::string& s, double er);

  Indel placeDel(std::string& s, uint32_t rpos, uint32_t ds, const std::string& refseq);

  void sampleReadsToCoverage(std::vector<std::string>& reads,  int cov, double error_rate, 
			     double ins_error_rate, double del_error_rate, int readlen);

  void samplePairedEndReadsToCoverage(std::vector<std::string>& reads1, std::vector<std::string>& reads2, 
				      int cov, double error_rate, double ins_error_rate, double del_error_rate,
				      int readlen, double mean_isize, double sd_isize);


  int getRandomIndelSize() const;

  void addAllele(const std::string& s, double af);

  Indel makeInsErrors(std::string& s, bool keep_size = true);

  Indel makeDelErrors(std::string& s, int sstart, const std::string& refseq);

  Indel makeDelErrors(std::string& s);

  void makeClipErrors(std::string& s, double er, int min_clip_len, int max_clip_len);

 private:

  // get the cumulative sum
  double __get_cumsum(std::vector<double>& cs) const; 

  int __random_allele(std::vector<double>& cs) const;

  // sequences to simulate from
  std::vector<std::string> m_seq;

  // allelic fractions for each
  std::vector<double> m_frac;

  // read SNP error rate
  double m_error;

  // insert size distribution
  double m_imean;
  double m_isd;

  

};

#endif

