#ifndef SNOWMAN_SEQFRAG_H__
#define SNOWMAN_SEQFRAG_H__

#include <string>

#include "SnowTools/GenomicRegion.h"
#include "SnowTools/BWAWrapper.h"
#include "ReadSim.h"

using SnowTools::GenomicRegion;

class SeqFrag {

 public:

 SeqFrag(const GenomicRegion& gr, faidx_t * findex) : m_gr(gr), m_index(findex) {}

  void getSeqFromRef(faidx_t * findex);

  int getLeftSide() const;

  int getRightSide() const;

  friend std::ostream& operator<<(std::ostream& out, const SeqFrag& s);

  std::string m_seq;

  void addScrambledEnds(size_t left_len, size_t right_len);

  char getStrand() const { return m_gr.strand; }

  void addIndels(size_t n);

  void addIns();

  void spikeMicrobe();

  std::vector<Indel> m_indels;  

  int frag_id;

  std::string left_scramble = "";
  std::string right_scramble = "";
 
  GenomicRegion m_gr;   
  
  int32_t phage_site = -1;
  std::string phage_string = "";

 private:

  faidx_t * m_index;

};


#endif
