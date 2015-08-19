#ifndef SNOWMAN_BAM_SPLITTER_H__
#define SNOWMAN_BAM_SPLITTER_H__

#include "SnowTools/BamWalker.h"

class BamSplitter: public SnowTools::BamWalker 
{

 public:

  BamSplitter() {}

  BamSplitter(const std::string& in, uint32_t seed) : SnowTools::BamWalker(in), m_seed(seed) {}
    
    void setWriters(const std::vector<std::string>& writers, const std::vector<double>& fracs);

  void splitBam();

 private:

    uint32_t m_seed = 0;

    std::vector<double> m_frac;

    std::vector<SnowTools::BamWalker> m_writers;

};

#endif
