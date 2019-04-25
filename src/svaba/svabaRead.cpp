#include "svabaRead.h"

svabaRead::svabaRead() : seq(nullptr) { }


void svabaRead::AddR2C(const std::string& contig_name, const r2c& r) {
    if (!m_r2c) 
      m_r2c = SeqPointer<R2CMap>(new R2CMap());

    auto it = m_r2c->find(contig_name);
    if(it != m_r2c->end()) 
      it->second = r;
    else
      m_r2c->insert({contig_name, r});
  }

r2c& svabaRead::GetR2C(const std::string& contig_name) const {
    assert(m_r2c);
    R2CMap::const_iterator ff = m_r2c->find(contig_name);
    assert(ff != m_r2c->end());
    return m_r2c->find(contig_name)->second;
  }

int svabaRead::SeqLength() const { 
  if (!seq)
    return Sequence().length();
  else
    return strlen(seq.get());
}

std::ostream& operator<<(std::ostream& out, const r2c& a) {
  out << "[" << a.start_on_contig << "," << a.end_on_contig << "] -- " 
      << "[" << a.start_on_read << "]" 
      << a.cig;
  return out;
}

svabaRead::svabaRead(const SeqLib::BamRecord r, const std::string& prefix) {

    b = r.shared_pointer();
    seq = nullptr;
    assert(prefix.length() >= 4);
    memcpy(p, prefix.data(), 4);
  }

std::string svabaRead::Prefix() const { 
  assert(p[0]=='t' || p[0] == 'n');
  return std::string(p, 4); 
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
