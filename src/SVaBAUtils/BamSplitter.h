#ifndef SNOWMAN_BAM_SPLITTER_H__
#define SNOWMAN_BAM_SPLITTER_H__

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "Fractions.h"

class BamSplitter: public SeqLib::BamReader
{

 public:

  BamSplitter() {}

  BamSplitter(uint32_t seed) : m_seed(seed) {}
    
    void setWriters(const std::vector<std::string>& writers, const std::vector<double>& fracs);

  void splitBam();

  void fractionateBam(const std::string& outbam, Fractions& f);

 private:

  std::string __startMessage() const;

    uint32_t m_seed = 0;

    std::vector<double> m_frac;

    std::vector<SeqLib::BamWriter> m_writers;

};

#endif
