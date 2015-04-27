#ifndef SNOW_SUFFIX_H__
#define SNOW_SUFFIX_H__

#include "reads.h"
#include <vector>

struct SuffArrayElement {

  uint16_t id;
  uint16_t pos;

};

typedef std::vector<SuffArrayElement> SuffArrayElementVector;

class SuffArray {

  //
  SuffArray() {}

  SuffArray(ReadVec &rv);

  void sacaInducedCopying();

 private:
  
  size_t m_num_reads;
  SuffArrayElementVector m_reads;
};

#endif
