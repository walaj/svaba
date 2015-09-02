#ifndef SNOWMAN_SIMTRAINER_WALKER_H__
#define SNOWMAN_SIMTRAINER_WALKER__

#include "SnowTools/BamWalker.h"
#include "SnowTools/BamStats.h"

class SimTrainerWalker : public SnowTools::BamWalker {

 public:
  
  SimTrainerWalker(const std::string& bam) : SnowTools::BamWalker(bam) {}

  void train();

  std::string printBamStats() const;
 
 private:

  SnowTools::BamStats m_bam_stats;
    
};


#endif
