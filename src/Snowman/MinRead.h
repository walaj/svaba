#ifndef MINREAD_H
#define MINREAD_H

#define READLEN 101

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include "EncodedString.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"

using namespace std;

// minimal representation of a sequencing read
struct MinRead {

  MinRead(string mseq = "", int mpos = 0, int mrefID = 0, int mstart = 0, int mstop = 0) : 
    seq(mseq), pos(mpos), refID(mrefID), start(mstart), stop(mstop) { }

  DNAString seq; // 2-bit encoded
  uint16_t pos;
  uint16_t refID;
  uint16_t start; // after q-trim start
  uint16_t stop;  // after q-trim stop
  uint16_t bin;
  int32_t is;
  uint32_t flag;
  vector<CigarOp> cig;

};

typedef shared_ptr<MinRead> MinReadSP;
typedef vector<MinReadSP> MinReadSPVector;
typedef vector<MinRead> MinReadVector;


#endif
