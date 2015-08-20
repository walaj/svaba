#include "SimGenome.h"

#include "SnowTools/GenomicRegionCollection.h"



SimGenome::SimGenome(const SnowTools::GenomicRegion& gr, int nbreaks, int ndels, faidx_t * findex, bool scramble) : m_gr(gr) {
  
  std::vector<size_t> index = {0};
  size_t ii = 1;

  // choose a set of random break points 
  SnowTools::GRC grc;
  grc.add(SnowTools::GenomicRegion(m_gr.chr, m_gr.pos1, m_gr.pos1, '+'));
  while(grc.size() < (size_t)(nbreaks+1)) {
  
    //uint32_t r;
    //SnowTools::genRandomValue(r, m_gr.width(), rand() % maxrand); // nu
    uint32_t r = rand() % m_gr.width();
    r += m_gr.pos1; 
    
    char strand = (rand() % 2) ? '+' : '-';
    grc.add(SnowTools::GenomicRegion(0, r, r, strand));
    index.push_back(ii++);

  }

  // setup the scramble probability distribution
  std::vector<size_t> scramble_size;
  for (size_t i = 1; i <= 200; ++i) {
    if (i <= 100)
      scramble_size.push_back(0);
    else 
      scramble_size.push_back(i - 100);
  }

  grc.SortAndStretchRight(m_gr.pos2);

  // shuffle the indices
  if (index.size() > 2)
    std::random_shuffle(index.begin() + 1, index.end() - 1); // seeded by srand

  // get the sequence and add insertions / deletions
  double grcwidth = grc.width();
  size_t id = 0;
  for (auto& i : index) {
    SeqFrag sf(grc.at(i));
    sf.frag_id = id;
    ++id;
    sf.getSeqFromRef(findex);

    int nd = ceil((double)sf.m_seq.length() / grcwidth * (double)ndels);
    sf.addDels(nd, findex);

    // scramble the ends
    if (scramble) {
      size_t left = scramble_size[rand() % 200];
      size_t right = scramble_size[rand() % 200];
      sf.addScrambledEnds(left, right);
    }

    m_sfv.push_back(sf);
    m_indels.insert(m_indels.end(), sf.m_indels.begin(), sf.m_indels.end());


  }
  
  //debug
  size_t mmm = 0;
  for (auto& i : m_sfv)
    mmm += i.m_indels.size();
  std::cerr << "TOTAL NUMBER OF INDELS IS " << mmm << std::endl;
  
}

std::string SimGenome::printBreaks() const {

  std::stringstream ss;

  // MARCIN strand convention here. Left side strand for revComp frag is "-". 
  // right side strand for non-revcomp is "-"

  for (size_t i = 1; i < m_sfv.size(); ++i) {
    ss << m_sfv[i-1].m_gr.chr << "\t" << 
      m_sfv[i-1].getRightSide() << "\t" << (m_sfv[i-1].getStrand() == '+' ? '-' : '+')  << "\t" << m_sfv[i-1].left_scramble << "\t" << 
      m_sfv[i].getLeftSide()    << "\t" << (m_sfv[i].getStrand())                       << "\t" << m_sfv[i].right_scramble << std::endl;
  }
  return ss.str();

}

std::string SimGenome::getSequence() const {

  std::string s;
  for (auto& i : m_sfv)
    s += i.m_seq;
  return s;
}

std::ostream& operator<<(std::ostream& out, const SimGenome& s) {
  
  //for (auto& i : s.m_sfv) {
  //  out << i << std::endl;
  //}
 
  return out;
}
