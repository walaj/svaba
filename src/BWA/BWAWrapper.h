#ifndef BWA_WRAPPER_H
#define BWA_WRAPPER_H

#define MAX_READS 100

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <memory>

extern "C" {
  #include "bwa.h"
  #include <stdlib.h>
}

typedef std::pair<std::string, std::string> BWARead; // name, seq
typedef std::vector<BWARead> BWAReadVec;

struct SamRecord {
  
  SamRecord() {}
  ~SamRecord() {}

  SamRecord(std::string sam);

  friend std::ostream& operator<<(std::ostream &out, const SamRecord &sam);
  
  static void parseXPString(std::string xp, std::vector<int> &mapq_vec, std::vector<int> &nm_vec);

  std::string name;
  std::string seq;
  int chr;
  int pos;
  int mapq  = -1;
  //char strand;
  int flag;
  std::string cigar;
  int as;
  int xs;
  int nm = -1;
  std::string xp = "";

  std::vector<int> mapq_vec;
  std::vector<int> nm_vec;

  std::string record;
  
  bool hasXP() const { return xp != ""; } 
  
  int minMapq() const {return *std::min_element(mapq_vec.begin(), mapq_vec.end()); };

  int maxMapq() const { return *std::max_element(mapq_vec.begin(), mapq_vec.end()); };

  int minNM() const {  return *std::min_element(nm_vec.begin(), nm_vec.end()); };

  int maxNM() const { return *std::max_element(nm_vec.begin(), nm_vec.end()); };

};

typedef std::vector<SamRecord> SamRecordVec;

class BWAWrapper {

 public:
  BWAWrapper() { 
    m_bseqs = (bseq1_t*)malloc(sizeof(bseq1_t)); 
  }

  ~BWAWrapper() { 
    free (m_bseqs); 
  }

  bseq1_t t;

  void addSequences(const BWAReadVec &seqs, std::unique_ptr<bwaidx_t>* idx, SamRecordVec &sam);

  //void makebseq1(const std::string &name, const std::string &seq, bseq1_t *s);

  void memProcess(bwaidx_t *idx);
  
 private:

  bseq1_t * m_bseqs;
  int m_blen = 0;
  
};



/*
typedef struct {
  int l_seq, id;
  char *name, *comment, *seq, *qual, *sam;
  } bseq1_t;
*/

/*
class Bseq : public bseq1_t {
  
 public:
  Bseq() : l_seq(), id(0), name(0), comment(0), seq(0), qual(0), sam(0) { }
    //void clear() { buf_clear(this); }
    //bool print() { return buf_print(this); }
    //bool append(const char* p, unsigned c)
    //{ return buf_append(this, p, c); }
};
*/

#endif
