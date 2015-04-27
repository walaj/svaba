#ifndef SNOWMAN_COVERAGE_H__
#define SNOWMAN_COVERAGE_H__

#include <memory>
#include <unordered_map>
#include <vector>
#include "bedFile.h"
#include <cstdint>
#include <vector>
#include <cassert> 
#include <iostream>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"

typedef std::shared_ptr<bam1_t> Read;

struct cfree_delete {
  void operator()(void* x) { free (x); }
};

typedef std::shared_ptr<std::vector<uint16_t>> uint16_sp;
//typedef std::unordered_map<int32_t, uint16_sp> CovMap;

/*! Class to hold base-pair or binned coverage across the genome
 */
class Coverage {
  
 private:

  size_t id = 0;  
  size_t start = 0;
  size_t end = 0;

  uint16_sp v;
  
 public:

  Coverage() : id(0), start(0), end(0) {}
  Coverage(int id, int start, int end);
  void addRead(Read &r);
  
  /*! Add to coverage objects together to get total coverge 
   * 
   * @param cov Coverage object to add to current object
   */
  //void combineCoverage(Coverage &cov);

  friend std::ostream& operator<<(std::ostream &out, const Coverage &c);

  uint16_t getCoverageAtPosition(size_t pos) const;
  
};

#endif
