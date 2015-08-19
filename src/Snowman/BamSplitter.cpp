#include "BamSplitter.h"

#include "htslib/khash.h"

void BamSplitter::splitBam() {

  assert(m_writers.size() == m_frac.size());

  std::vector<double> csum;
  double cs = 0;
  for (auto& i : m_frac) {
    cs += i;
    csum.push_back(cs);
  }
    

  SnowTools::BamRead r;
  
  bool rule_pass;

  std::cerr << "**Starting to split BAM "; 
  if (m_region.size() == 1)
    std::cerr << " on region " << m_region[0] << std::endl;
  else if (m_region.size() == 0)
    std::cerr << " on WHOLE GENOME" << std::endl;
  else
    std::cerr << " on " << m_region.size() << " regions " << std::endl;

  size_t countr = 0;
  std::vector<size_t> all_counts(m_writers.size(), 0);
  while (GetNextRead(r, rule_pass))
    {

      // print a message
      ++countr;
      if (countr % 1000000 == 0) 
	std::cerr << "...at position " << r.Brief() << std::endl; 

      // put in one BAM or another
      uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(r.Qname().c_str()) ^ m_seed);
      double hash_val = (double)(k&0xffffff) / 0x1000000;

      for (size_t i = 0; i < m_writers.size(); ++i) {
	// decide whether to keep
	if (hash_val <= csum[i]/*sample_rate*/) {
	  r.RemoveAllTags();
	  m_writers[i].WriteAlignment(r);
	  ++all_counts[i];
	  break;
	}
      }
      
      
    }
  
}


void BamSplitter::setWriters(const std::vector<std::string>& writers, const std::vector<double>& fracs) {

  assert(fracs.size() == writers.size());  

  for (auto& i : writers) {

    std::cerr << "...setting up BAM" << i << std::endl;
    SnowTools::BamWalker w;
    bam_hdr_t * h = bam_hdr_dup(br.get()/*header()*/);
    w.SetWriteHeader(h);
    w.OpenWriteBam(i);
    m_writers.push_back(w);

  }

  m_frac = fracs;


}
