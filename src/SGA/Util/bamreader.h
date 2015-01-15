#ifndef BAMREADER_H
#define BAMREADER_h

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "Util.h"

class SVBamReader {

 public:

  Bamreader(const std::string filename, const std::string prefix, const int isize = 2000) 
    : m_bam(filename), m_prefix(prefix), m_isize(isize) {}
  ~Bamreader() { }

  void bamToSeqRecordVector(SeqRecordVector &srv);

  bool setBamRegion(int refid1, int refid2, int pos1, int pos2);

  bool findBamIndex();

 private:
  std::string m_bam;
  std::string m_prefix;
  int m_isize;
  std::string m_index;
  BamTools::BamReader m_reader;
  BamTools::BamRegion m_region;
}


#endif
