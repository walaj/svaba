#ifndef SVABA_WORKITEM_H__
#define SVABA_WORKITEM_H__

#include "SnowTools/GenomicRegionCollection.h"
#include "run_snowman.h"

/** @brief p-thread work item that calls svaba on a small region

    Detailed description follows here.
    @author X. XYZ, DESY
    @date March 2008
*/

class svabaWorkItem {

 private:
  SnowTools::GenomicRegion m_gr;
  int m_number;  

 public:
  svabaWorkItem(const SnowTools::GenomicRegion& gr, int number)  
    : m_gr(gr), m_number(number) {}
    ~svabaWorkItem() {}
    
    int getNumber() { return m_number; }
    
    bool run() { return grabReads(m_gr); }
    
};


#endif
