#pragma once

#include <memory>
#include <unordered_map>
#include <vector>
#include <cstdint>
#include <vector>
#include <cassert> 
#include <iostream>
#include <memory>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"

#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"

namespace SeqLib {
  class BamRecord;
}

typedef std::unordered_map<int32_t,uint32_t> CovMap;

/** Hold base-pair or binned coverage across an interval or genome
 *
 * Currently stores coverage as an unordered_map.
 */
class STCoverage {
  
private:
  
  std::unordered_map<int32_t, CovMap> m_map; //[chr [pos,cov] ] 

 public:

  STCoverage() = default;
  
  /** Clear the coverage map */
  void clear();
      
  /** Add a read to this coverage track 
   * @param reserve_size Upper bound estimate for size of map. Not a hard
   *   cutoff but improves performance if total number of positions is less than this, 
   *   as it will not rehash. */
  void addRead(const SeqLib::BamRecord &r, int buff); //, bool full_length);

  /** Print the entire data */
  friend std::ostream& operator<<(std::ostream &out, const STCoverage &c);

  /** Return the coverage count at a position */
  int getCoverageAtPosition(int chr, int pos) const;
  
};

