#ifndef SNOWMAN_SIMTRAINER_WALKER_H__
#define SNOWMAN_SIMTRAINER_WALKER_H__

#include "SeqLib/BamReader.h"
#include "BamStats.h"

class SimTrainerWalker : public SeqLib::BamReader {

 public:
  
  SimTrainerWalker() {}

  void train();

  std::string printBamStats() const;
 
 private:

  BamStats m_bam_stats;
    
};


#endif
