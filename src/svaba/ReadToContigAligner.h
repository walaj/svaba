#ifndef SNOWMAN_READ_TO_CONTIG_ALIGNER_H__
#define SNOWMAN_READ_TO_CONTIG_ALIGNER_H__

#include "SnowTools/AlignedContig.h"
#include "SnowTools/BamRead.h"
#include "SnowTools/BWAWrapper.h"

class ReadToContigAligner {

 public: 

  ReadToContigAligner() {}

  ReadToContigAligner(const SnowTools::USeqVector& usv, BamReadVector& bav_this, std::vector<SnowTools::AlignedContig>& this_alc) {
  
 private:

  BamReadVector m_reads;

  SnowTools::USeqVector m_usv;

  SnowTools::BWAWrapper& m_bw;

  std::vector<SnowTools::AlignedContig> m_alc;

};

#endif
