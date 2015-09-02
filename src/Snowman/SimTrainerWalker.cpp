#include "SimTrainerWalker.h"

void SimTrainerWalker::train() {

  SnowTools::BamRead r;
  bool rule;
  size_t count = 0;
  while (GetNextRead(r, rule)) {
   
    ++count;
    
    // check occasionally that we still have read groups
    if (count % 10000 == 0)
      assert(r.GetZTag("RG").length());

    if (count % 1000000 == 0)
      std::cerr << "...training on read " << SnowTools::AddCommas(count) << " at read at " << r.Brief(br.get()) << std::endl;

    m_bam_stats.addRead(r);
    
  }
   

}

std::string SimTrainerWalker::printBamStats() const {

  std::stringstream ss;
  ss << m_bam_stats;
  return ss.str();
}
