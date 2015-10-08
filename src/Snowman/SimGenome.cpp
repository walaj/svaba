#include "SimGenome.h"

#include "SnowTools/GenomicRegionCollection.h"

SimGenome::SimGenome(const SnowTools::GenomicRegion& gr, int nbreaks, int ndels, faidx_t * findex, bool scramble, int viral_count) : m_gr(gr) {
  
  std::vector<size_t> index = {0};
  size_t ii = 1;

  double viral_prob = viral_count ? (double)viral_count/((double)(nbreaks+1)) * 2 * 1000: 0;

  // int
  size_t avg_width  = std::max((size_t)(m_gr.width() / nbreaks  ), (size_t)1000);
  size_t max_breaks = std::min((size_t)(m_gr.width() / avg_width), (size_t)nbreaks);
  std::cerr << "avg width: " << avg_width << " max breaks: " << max_breaks << std::endl;
  size_t real_vcount = 0;

  // choose a set of random break points 
  SnowTools::GRC grc;
  grc.add(SnowTools::GenomicRegion(m_gr.chr, m_gr.pos1, m_gr.pos1, '+'));
  size_t exp_break = 0;
  while(grc.size() < (size_t)(nbreaks)) {
  
    exp_break += avg_width;
    uint32_t r = exp_break + (rand() % 400) - 200; //m_gr.width();
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
  grc.sendToBED("segments.bed");
  
  // shuffle the indices
  if (index.size() > 2)
    std::random_shuffle(index.begin() + 1, index.end() - 1); // seeded by srand
  
  //debug
  //for (auto& i : index)
  //  std::cerr << i << "\t";
  //std::cerr << std::endl;
  
  // get the sequence and add insertions / deletions
  double grcwidth = grc.width();
  size_t id = 0;
  for (auto& i : index) {
    SeqFrag sf(grc.at(i), findex);
    sf.frag_id = id;
    ++id;
    sf.getSeqFromRef(findex);

    if (rand() % 1000 < viral_prob && grc.at(i).strand == '+') {
      sf.spikeMicrobe();
      ++real_vcount;
    } else {
      int nd = ceil((double)sf.m_seq.length() / grcwidth * (double)ndels);
      sf.addIndels(nd);
    }

    // scramble the ends
    if (scramble) {
      size_t left = 0, right = 0;
      //if (i != 0 && sf.m_gr.strand != '+') // don't scramble left side of +
      //left = scramble_size[rand() % 200];
      //else if (i != index.size() && sf.m_gr.strand != '-') // dont' scramble right side - 
      if (i != index.size())
	right = scramble_size[rand() % 200];
      sf.addScrambledEnds(0, right);
    }

    m_sfv.push_back(sf);
    m_indels.insert(m_indels.end(), sf.m_indels.begin(), sf.m_indels.end());

  }

  //debug
  std::ofstream sg;
  sg.open("segments.tsv", std::ios::out);
  for (auto& i : m_sfv) {
    sg << i << std::endl;
  }
  sg.close();
  
  //debug
  size_t mmm = 0;
  for (auto& i : m_sfv)
    mmm += i.m_indels.size();
  std::cerr << "TOTAL NUMBER OF INDELS IS " << mmm << " AND VIRAL INTEGRATION IS " << real_vcount << std::endl;

}

std::string SimGenome::printBreaks() const {

  std::stringstream ss;

  // MARCIN strand convention here. Left side strand for revComp frag is "-". 
  // right side strand for non-revcomp is "-"

  for (size_t i = 1; i < m_sfv.size(); ++i) {
    //ss << m_sfv[i-1].m_gr.chr << "\t" << 
      //      m_sfv[i-1].getRightSide() << "\t" << (m_sfv[i-1].getStrand() == '+' ? '-' : '+')  << "\t" << m_sfv[i-1].left_scramble << "\t" << 
      //m_sfv[i].getLeftSide()    << "\t" << (m_sfv[i].getStrand())                       << "\t" << m_sfv[i].right_scramble << std::endl;
    ss << m_sfv[i-1].m_gr.chr << "\t" << 
      m_sfv[i-1].getRightSide() << "\t" << (m_sfv[i-1].getStrand() == '+' ? '-' : '+')  << "\t" << 
      m_sfv[i].getLeftSide()    << "\t" << (m_sfv[i].getStrand())                       << "\t" << m_sfv[i-1].right_scramble << std::endl;
  }

  // loop through the tandem duplications and indels
  for (auto& i : m_indels) {
    //if (i.len >= 50 && i.type == 'D') 
    //  ss << i.gr.chr << "\t" << i.gr.pos1 << "\t-" << i.gr.pos2 << "\t+" << std::endl;
    if (i.len >= 50 && i.type == 'I') // tandem dup
      ss << i.gr.chr << "\t" << i.gr.pos2 << "\t-\t" << i.gr.pos1 << "\t+" << std::endl;      
  }
  return ss.str();

}

std::string SimGenome::printMicrobeSpikes() const {

  std::stringstream ss;
  for (auto& i : m_sfv) {
    if (i.phage_string.length()) {
      ss << i.m_gr.chr << "\t" << (i.m_gr.pos1 + i.phage_site) << "\t" << i.phage_string << std::endl;
    }
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
