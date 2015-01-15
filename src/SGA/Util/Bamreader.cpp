#include "bamreader.h"

bool Bamreader::bamToSeqRecordVector(SeqRecordVector &srv) {
 
  if (!m_reader.Open(m_bam) || !reader.SetRegion(m_region)) 
    return false;

  BamTools::BamAlignment a;

  while (m_reader.GetNextAlignment())

}

bool Bamreader::findBamIndex() {

  m_bai = m_bam.substr(0, m_bam.size()-3.append("bai"));
  return m_reader.OpenIndex(m_bai);

}

bool Bamreader::setBamRegion(int refid1, int refid2, int pos1, int pos2);
