#include "svabaRead.h"

svabaRead::svabaRead(const SeqLib::BamRecord r, const std::string& prefix) {
    b = r.shared_pointer();
    seq = nullptr;
    assert(prefix.length() >= 4);
    memcpy(p, prefix.data(), 4);
  }

std::string svabaRead::Seq() const {

  assert(seq);
  if (!seq)
    return Sequence();
  else
    return std::string(seq.get());

}

void svabaRead::SetSeq(const std::string& nseq) {
  seq = SeqPointer<char>(strdup(nseq.c_str()));
}

std::string svabaRead::SR() const {
  return(std::string(p, 4) + "_" + std::to_string(AlignmentFlag()) + "_" + Qname());
}

/*void svabaRead::Reassign(const svabaRead& s) {

  //std::string sr = SRTAG(r);
  sr = r.SR();
  
  // make a deep copy
  // now if s is deleted, it doesn't affect r 
  // (t is new memory location)
  bam1_t* t = bam_dup1(s.raw());

  // should carry with it everything (tags etc)
  assign(t);

  SetChrIDMate(s.MateChrID());
  SetPositionMate(s.MatePosition());
  SetPairMappedFlag();

  if (s.MateReverseFlag())
    SetMateReverseFlag();

  r = SeqPointer<char>(s.seq);

  }*/
