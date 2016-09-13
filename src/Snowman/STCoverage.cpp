#include "STCoverage.h"
#include "SeqLib/SeqLibCommon.h"
#include <stdexcept>
#include <algorithm>


using namespace SeqLib;

  void STCoverage::clear() {
    m_map.clear();
  }

  void STCoverage::settleCoverage() {
    //SeqLib::GRC tmp = m_grc;
    //m_grc.MergeOverlappingIntervals();
  }
  
  STCoverage::STCoverage(const SeqLib::GenomicRegion& gr) {
    //m_gr = gr;
    //v = uint16_sp(new std::vector<uint16_t>(gr.Width(),0));
  }

  uint16_t STCoverage::maxCov() const {
    return (*std::max_element(v->begin(), v->end()));
  }

  /*void STCoverage::addRead2(const BamRecord& r) {

    int p = r.Position();
    int e = r.PositionEnd();

    if (p < 0 || e < 0)
      return;

    if (p < m_gr.pos1 || e > m_gr.pos2)
      return;

    if (r.ChrID() != m_gr.chr)
      return;

    assert(p - m_gr.pos1 < v->size());
    ++v[p - m_gr.pos1];
    }
  */
  void STCoverage::addRead(const BamRecord &r, int buff, bool full_length) {
    
    //m_settled = false;
    //m_grc.add(GenomicRegion(r.ChrID(), r.Position(), r.PositionEnd()));
    
    // out of bounds
    //if (r.Position() < m_gr.pos1 || r.PositionEnd() > m_gr.pos2 || r.ChrID() != m_gr.chr)
    //  return;

    //int p = std::min(r.Position() - m_gr.pos1, m_gr.pos2);
    //int e = std::min(r.PositionEnd() - m_gr.pos1, m_gr.pos2);
    int p = -1; 
    int e = -1;

    if (full_length) {
      Cigar c = r.GetCigar();
      // get beginning
      if (c.size() && c[0].RawType() == BAM_CSOFT_CLIP)
	p = std::max((int32_t)0, r.Position() - (int32_t)c[0].Length()); // get prefixing S
      else
	p = r.Position();
      // get end
      if (c.size() && c.back().RawType() == BAM_CSOFT_CLIP)
	e = r.PositionEnd() + c.back().Length();
      else
	e = r.PositionEnd();
    }
    else {
      p = r.Position() + buff;
      e = r.PositionEnd() - buff;
    }

    if (p < 0 || e < 0)
      return;

    // if we don't have an empty map for this, add
    if (r.ChrID() >= (int)m_map.size()) {
      int k = m_map.size();
      while (k <= r.ChrID()) {
	m_map.push_back(CovMap());
	//m_map.back().reserve(reserve_size);
	++k;
      }
    }

    assert(e - p < 1e6); // limit on read length
    assert(r.ChrID() >= 0);
    assert(r.ChrID() < (int)m_map.size());

    try {
       while (p <= e) {
	//CovMap::iterator iter = m_map.find(p);
	++(m_map[r.ChrID()][p]); // add one to this position
	++p;
	//if (v->at(p) < 60000) // 60000 is roughly int16 lim
	//  v->at(p)++;
	//++p;
	
      }
    } catch (std::out_of_range &oor) {
      std::cerr << "Position " << p << " on tid " << r.ChrID()
		<< " is greater than expected max of " << v->size() << " -- skipping" << std::endl;
      
    }
    
  }
  
  std::ostream& operator<<(std::ostream &out, const STCoverage &c) {
    out << "Region " << c.m_gr << " v.size() " << c.v->size() << std::endl;
    return out;
  }

  void STCoverage::ToBedgraph(std::ofstream * o, const bam_hdr_t * h) const {

    //settleCoverage();

    // unitialized so nothing to do
    if (m_gr.chr == -1 || v->size() == 0)
      return;

    size_t curr_start = 0;
    size_t curr_val = v->at(0);
    for (size_t i = 0; i < v->size(); ++i) {
      if (v->at(i) != curr_val) {
	(*o) << m_gr.ChrName(h) << "\t" << (curr_start + m_gr.pos1) << "\t" << (i+m_gr.pos1) << "\t" << curr_val << std::endl;
	curr_start = i;
	curr_val = v->at(i);
      }
    }
    
    // need to dump last one
    if ( (curr_start+1) != v->size()) 
      (*o) << m_gr.ChrName(h) << "\t" << (curr_start + m_gr.pos1) << "\t" << (v->size()+m_gr.pos1-1) << "\t" << curr_val << std::endl;
  }
  
  int STCoverage::getCoverageAtPosition(int chr, int pos) const {

    //CovMapMap::iterator it = m_map.find(chr);
    //if (it == m_map.end())
    //  return 0;
    if (chr >= (int)m_map.size())
      return 0;
    
    //std::cerr << " MAP " << std::endl;
    //for (auto& i : m_map)
    //  std::cerr << i.first << " " << i.second << std::endl;
    
    //if (!m_settled)
    //  settleCoverage();

    //if (pos < m_gr.pos1 || pos > m_gr.pos2) {
      //std::cerr << "Coverage query out of bounds for location " << m_gr.chr << ":" << pos << std::endl;
    //  return 0;
    //}
    
    //size_t q = pos - m_gr.pos1;
    //if (q >= v->size()) {
    // std::cerr << "Coverage query out of bounds for location " << m_gr.chr << ":" << pos << " with pos-start of " << q << " attempt on v of size " << v->size() << std::endl;
    //  return 0;
    //}
    //return it->second[pos];
    
    CovMap::const_iterator ff = m_map[chr].find(pos);
    if (ff == m_map[chr].end()) {
      return 0;
    }

    return ff->second;

    //return (v->at(q));
    

}
