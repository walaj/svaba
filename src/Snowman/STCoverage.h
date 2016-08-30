#ifndef SNOWMAN_SEQLIB_COVERAGE_H__
#define SNOWMAN_SEQLIB_COVERAGE_H__

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

#include "SeqLib/BamRecord.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"

typedef std::shared_ptr<std::vector<uint16_t>> uint16_sp;
typedef std::unordered_map<int,int> CovMap;
//typedef std::unordered_map<int,CovMap> CovMapMap;

  /** Hold base-pair or binned coverage across an interval or genome
   *
   * Currently stores coverage as an unordered_map.
   */
class STCoverage {
  
 private:

  SeqLib::GRC m_grc;
  SeqLib::GenomicRegion m_gr;

  //CovMapMap m_map;
  std::vector<CovMap> m_map;

  uint16_sp v;


 public:

  /** Clear the coverage map */
  void clear();


  /** */
  void settleCoverage();
      
  /** Add a read to this coverage track 
   * @param reserve_size Upper bound estimate for size of map. Not a hard
   *   cutoff but improves performance if total number of positions is less than this, 
   *   as it will not rehash. */
  void addRead(const SeqLib::BamRecord &r, int buff, bool full_length);

  /** Make a new coverage object at interval gr */
  STCoverage(const SeqLib::GenomicRegion& gr);

  uint16_t maxCov() const;

  /** Make an empty coverage */
  STCoverage() {}

  /*! Add to coverage objects together to get total coverge 
   * 
   * @param cov Coverage object to add to current object
   */
  //void combineCoverage(Coverage &cov);

  /** Return a short summary string of this coverage object */
  void ToBedgraph(std::ofstream * o, const bam_hdr_t * h) const;
  
  /** Print the entire data */
  friend std::ostream& operator<<(std::ostream &out, const STCoverage &c);

  /** Return the coverage count at a position */
  int getCoverageAtPosition(int chr, int pos) const;
  
};

#endif
