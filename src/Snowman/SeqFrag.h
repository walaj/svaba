#ifndef SNOWMAN_SEQFRAG_H__
#define SNOWMAN_SEQFRAG_H__

#include <string>

#include "SnowTools/GenomicRegion.h"
#include "SnowTools/BWAWrapper.h"
#include "ReadSim.h"

using SnowTools::GenomicRegion;

class SeqFrag {

 public:

  SeqFrag(const GenomicRegion& gr) : m_gr(gr) {}

  void getSeqFromRef(faidx_t * findex);

  int getLeftSide() const;

  int getRightSide() const;

  friend std::ostream& operator<<(std::ostream& out, const SeqFrag& s);

  std::string m_seq;

  std::string m_parent;

  char getStrand() const { return m_gr.strand; }

  void addDels(size_t n, faidx_t * findex);

  void addIns();

  std::vector<Indel> m_indels;  

  int frag_id;
 private:
  
  GenomicRegion m_gr;



};


#endif
