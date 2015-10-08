#ifndef SNOWMAN_READ_H__
#define SNOWMAN_READ_H__

#include "SnowTools/BamRead.h"

class SnowmanRead {
  
 public:

  SnowmanRead() {}

  std::string qseq;
  std::string sr;
  SnowTools::BamRead r;

};


#endif
