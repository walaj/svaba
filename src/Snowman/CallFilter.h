#ifndef SNOWMAN_CALL_FILTER
#define SNOWMAN_CALL_FILTER

#include "SnowTools/BreakPoint.h"

using SnowTools::BreakPoint;

class IndelCallFilter {

 public:
  
  IndelCallFilter(int mint, int minn, int minr) : m_min_tsplit(mint), m_min_nsplit(minn), m_min_som_ratio(minr) {}
    
    bool filterBreakpoint(BreakPoint& bp);
    
 private:

  int m_min_tsplit;
  int m_min_nsplit;
  int m_min_som_ratio;

};

#endif
