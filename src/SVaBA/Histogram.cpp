#include "Histogram.h"
#include "SeqLib/SeqLibUtils.h"
#include <fstream>
#include <cmath>
#include <algorithm>
#include <sstream>

#define BINARY_SEARCH 1

#define DEBUG_HISTOGRAM

using namespace SeqLib;

Histogram::Histogram(const int32_t& start, const int32_t& end, const uint32_t& width)
{
  
  if (end >= start)
    throw std::invalid_argument("Histogram end must be > start");

  Bin bin;
  bin.bounds.first = start;

  int32_t next_end = start + width - 1; // -1 because width=1 is bound.first = bounds.second

  while (next_end < end) 
    {
      // finish this bin
      bin.bounds.second = next_end;
      m_bins.push_back(bin);
      m_ind.push_back(bin.bounds.first); // make the index for the lower bound

      // start a new one
      bin.bounds.first = next_end+1;
      next_end += width;
    }
  
  // finish the last bin
  bin.bounds.second = end;
  m_bins.push_back(bin);
  m_ind.push_back(bin.bounds.first);

  // add a final bin
  //bin.bounds.first = end+1;
  //bin.bounds.second = INT_MAX;
  //m_bins.push_back(bin);
  //m_ind.push_back(bin.bounds.first);

}

void Histogram::toCSV(std::ofstream &fs) {

  for (auto& i : m_bins) 
    fs << i << std::endl;

}
void Histogram::removeElem(const int32_t& elem) {
  --m_bins[retrieveBinID(elem)];
}

void Histogram::addElem(const int32_t& elem) {
  ++m_bins[retrieveBinID(elem)];
}

std::string Histogram::toFileString() const {
  std::stringstream ss;
  for (auto& i : m_bins)
    if (i.m_count)
      ss << i.bounds.first << "_" << i.bounds.second << "_" << i.m_count << ",";
  std::string out = ss.str();
  out.pop_back(); // trim off last comma
  return(out);
  
}

size_t Histogram::retrieveBinID(const int32_t& elem) const {

  if (elem < m_bins[0].bounds.first) 
    {
#ifdef DEBUG_HISTOGRAM
      std::cerr << "removeElem: elem of value " <<  elem << " is below min bin " << m_bins[0] << std::endl;
      exit(1);
#endif
      return 0;
    }

  if (elem > m_bins.back().bounds.second) 
    {
#ifdef DEBUG_HISTOGRAM
      std::cerr << "removeElem: elem of value " <<  elem << " is above max bin " << m_bins.back() << std::endl;
      exit(1);
#endif
      return m_bins.size();
    }


  if (m_bins[0].contains(elem)) 
    return 0;
  if (m_bins.back().contains(elem)) 
    return m_bins.size();

#ifdef BINARY_SEARCH
  // binary search
  std::vector<int32_t>::const_iterator it = std::upper_bound(m_ind.begin(), m_ind.end(), elem);
  size_t i = it - m_ind.begin()-1;
  assert(i < m_ind.size());
  return i;
#else
  for (size_t i = 0; i < m_bins.size(); i++) {
    if (m_bins[i].contains(elem)) {
      return i;
    }
  }
#endif
  std::cerr << "bin not found for element " << elem << std::endl;
  return 0;
}

/*void Histogram::initialize(size_t num_bins, std::vector<int32_t>* pspanv, size_t min_bin_width) {

  // ensure that they spans are sorted
  std::sort(pspanv->begin(), pspanv->end());

  // fill the histogram bins with matrix pairs (pre-sorted by distance)
  Bin bin; 

  // get number of inter-chr events
  size_t intra = 0;
  for (auto& i : *pspanv)
    if (i != INTERCHR)
      intra++;

  int bin_cut = 0;
  try {
    bin_cut = floor((double)intra / (double)num_bins);
    if (bin_cut == 0)
      throw 1;
  } catch(...) {
    std::cerr << "Error in determining bin cut. Not enought events or too many bins?" << std::endl;
    std::cerr << "Events: " << pspanv->size() << " Num Bins " << num_bins << " quantile count (hist height) " << bin_cut << std::endl;
  }

  std::cout << "...Events per bin: " << bin_cut << " num bins " << num_bins << std::endl;

  S last_span = 0;
  size_t tcount = 0; // count events put into bins

  // iterate over spans
  for (auto& span : *pspanv) {
    if (span != INTERCHR) {
      
      ++tcount;
      
      // moved into a new bin? (or done?)
      if (bin.getCount() > bin_cut && span != last_span && (last_span - bin.bounds.first) >= min_bin_width) { 

	// finalize, save old bin
	bin.bounds.second = last_span;
	m_bins.push_back(bin);
	
	// new bin
	bin.bounds.first = last_span+1;
	bin.m_count = 0;
	
      }
      ++bin;
      if (bin.getCount() >= bin_cut) {
	last_span = span;
      }
      
      //update the size of current bin
      bin.bounds.second = span;
    }
  }
  // add the last bin
  bin.bounds.second = INTERCHR-1; 
  m_bins.push_back(bin);

  // add a bin for interchr events
  bin.bounds = {INTERCHR, INTERCHR};
  bin.m_count = pspanv->size() - intra;
  m_bins.push_back(bin);

  // make the indices of lower bound
  for (auto& i : m_bins)
    m_ind.push_back(i.bounds.first);

  if (m_bins.size() != (num_bins+1)) {
    //std::cout << " bin cut " << bin_cut << std::endl;
    //std::cout << " num bins " << num_bins << " bins.size() " << m_bins.size() << std::endl;
    //assert(bins.size() == (num_bins+1));
  }

  }*/

bool Bin::operator < (const Bin& b) const {
  return (bounds.first < b.bounds.first || (bounds.first==b.bounds.first && bounds.second < b.bounds.second));

}

bool Bin::contains(const int32_t& elem) const {

  return (elem >= bounds.first && elem <= bounds.second); 


}

Bin& Bin::operator++()
{
  ++m_count;
  return *this;
}


Bin& Bin::operator--() {
  assert(m_count > 0); 
  --m_count;
  return *this;
}
