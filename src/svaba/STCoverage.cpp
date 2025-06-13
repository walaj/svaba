#include "STCoverage.h"
#include "SeqLib/SeqLibCommon.h"
#include "SeqLib/BamRecord.h"
#include <stdexcept>
#include <algorithm>

using SeqLib::BamRecord;
using SeqLib::GenomicRegion;
using SeqLib::Cigar;

  void STCoverage::clear() {
    m_map.clear();
  }
  
void STCoverage::addRead(const BamRecord &r, int buff) { //, bool full_length) {
    
    int p = -1; 
    int e = -1;

    // old code to also cover soft clips, not needed
    /*    if (full_length) {
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
	  else { */
    p = r.Position() + buff;
    e = r.PositionEnd() - buff;
    //    }
    
    if (p < 0 || e < 0)
      return;
    
    assert(e - p < 1e6); // limit on read length
    assert(r.ChrID() >= 0);
    assert(r.ChrID() < (int)m_map.size());
    
    try {
      while (p <= e) {
	++m_map[r.ChrID()][p]; // add one to this position
	++p;
      }
    } catch (std::out_of_range &oor) {
      std::cerr << "Position " << p << " on tid " << r.ChrID()
		<< " is greater than expected max -- skipping" << std::endl;
      
    }
    
  }

int STCoverage::getCoverageAtPosition(int chr, int pos) const {

  auto chr_it = m_map.find(chr);
  if (chr_it == m_map.end()) {
    return 0;
  }
  
  auto pos_it = chr_it->second.find(pos);
  if (pos_it == chr_it->second.end()) {
    return 0;
  }
  
  return pos_it->second;
  
}
