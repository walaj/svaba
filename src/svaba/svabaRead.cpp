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

void svabaRead::QualityTrimRead() {

  int32_t startpoint = 0, endpoint = 0;
  QualityTrimmedSequence(3, startpoint, endpoint);
  int32_t new_len = endpoint - startpoint;
  if (endpoint != -1 && new_len < r.Length() && new_len > 0 && new_len - startpoint >= 0 && startpoint + new_len <= r.Length()) { 
    try { 
      SetSeq(Sequence().substr(startpoint, new_len));
    } catch (...) {
      std::cerr << "Subsequence failure with sequence of length "  
		<< Sequence().length() << " and startpoint "
		<< startpoint << " endpoint " << endpoint 
		<< " newlen " << new_len << std::endl;
    }

  } else {
    SetSeq(r.Sequence()); // copies the sequence
  }

  // remove the HTSlib version of qual and sequence
  // since we store the trimmed sequence in svabaRead
  SetSequence(std::string());

}

svabaRead::svabaRead(const SeqLib::BamRecord r, const std::string& prefix) {
  
  b = r.shared_pointer(); // copy the BamRecord main read pointer
  //seq = nullptr;
  assert(prefix.length() >= 4);
  p = prefix;
  //memcpy(p, prefix.data(), 4);
  
}

std::string svabaRead::Prefix() const { 
  assert(p.c_str()[0]=='t' || p.c_str()[0] == 'n');
  return p; //std::string(p, 4); 
}

std::string svabaRead::Seq() const {

  if (!seq)
    return Sequence();
  else
    return std::string(seq.get());

}

std::string svabaRead::CorrectedSeq() const {

  if (!seq_corrected)
    return Sequence();
  else
    return std::string(seq_corrected.get());

}

void svabaRead::SetSeq(const std::string& nseq) {
  seq = SeqPointer<char>(strdup(nseq.c_str()));
}

void svabaRead::SetCorrectedSeq(const std::string& nseq) {
  seq_corrected = SeqPointer<char>(strdup(nseq.c_str()));
}

std::string svabaRead::SR() const {
  return(p + "_" + std::to_string(AlignmentFlag()) + "_" + Qname());
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
