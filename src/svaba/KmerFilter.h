#pragma once

#include <map>
#include <string>

#include "SeqLib/BamRecord.h"

#include "CorrectionThresholds.h"
#include "OverlapCommon.h"
#include "svabaRead.h"

typedef std::map<std::string, int> KmerCountMap;

class KmerFilter {

 public:
  
  KmerFilter() : pBWT(nullptr), pSAf(nullptr) {}

    ~KmerFilter() { delete pBWT; delete pSAf; }

    int correctReads(SeqLib::BamRecordVector& vec);

    int correctReads(svabaReadVector& vec);

    void makeIndex(const std::vector<char*>& v);

    //void makeIndex(SeqLib::BamRecordVector& vec);
  
 private: 

  RLBWT* pBWT;
  SuffixArray* pSAf;

  int m_kmer_len = 31;

  bool attemptKmerCorrection(size_t i, size_t k_idx, size_t minCount, std::string& readSequence, BWTIndexSet& inds);



};
