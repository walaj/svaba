#ifndef SNOWMAN_READSIM_H__
#define SNOWMAN_READSIM_H__

#include <string>
#include <vector>
#include <cassert>
#include "SnowTools/GenomicRegion.h"

struct Indel {
  
  Indel() : len(0), type('N') {}

  Indel(size_t l, char t, const std::string& rseq, const std::string& aseq, const std::string& lseq) : len(l), type(t) {
    ref_seq = rseq;
    alt_seq = aseq;
    lead_base = lseq;
    assert(lseq.length() == 1);
    assert(t == 'D' || alt_seq.length() == l); 
    assert(t == 'I' || ref_seq.length() == l); 
  }

  size_t len;
  char type;
  std::string ref_seq, alt_seq, lead_base;
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

  void baseQualityRelevantErrors(std::string& s, const std::string& bq);

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

