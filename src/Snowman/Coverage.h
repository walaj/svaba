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
  //CovMap m_cov_map;;

  const uint32_t m_chr_len[25] = {249250621, 243199373, 198022430, 191154276, //1-4
				 180915260, 171115067, //5-6
				 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, //7-12
				 115169878, 107349540, 102531392, 90354753,  81195210,  78077248, //13-18
				 59128983,  63025520,  48129895,  51304566,  155270560, 59373566, //19-24
                                 16571}; //25
  const uint32_t m_chr_len_10[25] = {24925062, 24319937, 19802243, 19115427, //1-4
				 18091526, 17111506, //5-6
				 15913866, 14636402, 14121343, 13553474, 13500651, 13385189, //7-12
				 11516987, 10734954, 10253139, 9035475,  8119521,  7807724, //13-18
				 5912898,  6302552,  4812989,  5130456,  15527056, 5937356, //19-24
                                 1657}; //25

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
  void combineCoverage(Coverage &cov);

  friend std::ostream& operator<<(std::ostream &out, const Coverage &c);

  uint16_t getCoverageAtPosition(int pos) const;
  
};

#endif
