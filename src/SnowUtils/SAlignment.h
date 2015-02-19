#ifndef SNOW_ALIGNMENT_H
#define SNOW_ALIGNMENT_H

struct SAlignment {

  char seq[101];
  char qual[101];
  uint32_t flag;
  int32_t pos;
  int32_t len;
  int32_t mapq;
  

};

#endif
