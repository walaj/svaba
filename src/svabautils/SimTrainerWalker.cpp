#include "SimTrainerWalker.h"

void SimTrainerWalker::train() {

  SeqLib::BamRecord r;
  size_t count = 0;
  while (GetNextRecord(r)) {
   
    ++count;
    
    // check occasionally that we still have read groups
    if (count % 10000 == 0)
      assert(r.GetZTag("RG").length());

    if (count % 1000000 == 0)
      std::cerr << "...training on read " << SeqLib::AddCommas(count) << " at read at " << r.Brief() << std::endl;

    m_bam_stats.addRead(r);
    
  }
   

}

std::string SimTrainerWalker::printBamStats() const {

  std::stringstream ss;
  ss << m_bam_stats;
  return ss.str();
}
