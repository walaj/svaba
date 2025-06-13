#include "svabaRead.h"

svabaRead::svabaRead() = default;

void svabaRead::AddR2C(const std::string& contig_name, const r2c& r) {
  //    if (!m_r2c) 
  //    m_r2c = SeqPointer<R2CMap>(new R2CMap());
  
  auto it = m_r2c.find(contig_name);
  if(it != m_r2c.end()) 
    it->second = r;
  else
    m_r2c.insert({contig_name, r});
}

r2c svabaRead::GetR2C(const std::string& contig_name) const {
  //assert(m_r2c);
  R2CMap::const_iterator ff = m_r2c.find(contig_name);
  assert(ff != m_r2c.end());
  return m_r2c.find(contig_name)->second;
}

int svabaRead::CorrectedSeqLength() const {
  if (seq_corrected.length())
    return seq_corrected.length();
  else
    return Length();
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
  if (endpoint != -1 && new_len < Length() &&
      new_len > 0 &&
      new_len - startpoint >= 0 &&
      startpoint + new_len <= Length()) { 
    try { 
      SetCorrectedSeq(Sequence().substr(startpoint, new_len));
    } catch (...) {
      std::cerr << "Subsequence failure with sequence of length "  
		<< Sequence().length() << " and startpoint "
		<< startpoint << " endpoint " << endpoint 
		<< " newlen " << new_len << std::endl;
    }

  } else {
    SetCorrectedSeq(Sequence()); // copies the sequence into private "seq" char
  }

  // remove the HTSlib version of qual and sequence
  // since we store the trimmed sequence in svabaRead
  //SetSequence(std::string());

}

svabaRead::svabaRead(const SeqLib::BamRecord& r, std::string_view prefix)
  : SeqLib::BamRecord()    // base class ctor will init `b` to nullptr
{
  // bam_dup1 will allocate-and-copy a new bam1_t for us:
  bam1_t* dup = bam_dup1(r.raw());
  if (!dup)
    throw std::runtime_error("svabaRead: failed to duplicate BamRecord");

  // wrap it in your shared-ptr with the proper deleter
  b = SeqPointer<bam1_t>(dup, SeqLib::Bam1Deleter());

  // now copy over your svabaRead-specific prefix
  assert(prefix.size() >= 4);
  p.assign(prefix);  
}

std::string svabaRead::Prefix() const { 
  assert(p.c_str()[0]=='t' || p.c_str()[0] == 'n');
  return p; //std::string(p, 4); 
}

std::string svabaRead::CorrectedSeq() const {
  
  if (!seq_corrected.length())
    return Sequence();
  else
    return seq_corrected;

}

void svabaRead::SetCorrectedSeq(std::string_view nseq) {
  seq_corrected = nseq;
}

std::string svabaRead::UniqueName() const {
  return(p + "_" + std::to_string(AlignmentFlag()) + "_" + Qname());
}
