#include "Util.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include <vector>

class AlignedContig {

 public:  
  AlignedContig(const BamTools::BamAlignment align) { m_align.push_back(align); }
  ~AlignedContig() {}

  void addAlignement(const BamTools::BamAlignment align) { 
        m_align.push_back(align); 
        m_num++;
  }

  void addRead(const SeqRecord sr) { m_sr.push_back(sr); }

 private:
  std::vector<BamTools::BamAlignment> m_align;
  std::vector<SeqRecord> m_sr;
  unsigned m_num;
 
}
