#include "Coverage.h"

#include "SnowTools/SnowToolsCommon.h"

Coverage::Coverage(int tid, int tstart, int tend) : id(tid), start(tstart), end(tend) {
  assert(end > start);
  v = uint16_sp(new std::vector<uint16_t>(tend - tstart,0));
}


void Coverage::addRead(Read &r) {

  try {
    size_t p = r->core.pos - start;
    size_t e = bam_endpos(r.get()) - start;
    while (p <= e) {
      if (v->at(p) < 60000)
	v->at(p)++;
      ++p;
	  
    }
  } catch (out_of_range &oor) {
    std::cerr << "Position " << (r->core.pos) << " on tid " << r->core.tid 
	      << " is greater than expected max of " << v->size() << " -- skipping" << std::endl;

  }

}

std::ostream& operator<<(std::ostream &out, const Coverage &c) {
  
  size_t curr_start = 0;
  size_t curr_val = c.v->at(0);
  for (size_t i = 0; i < c.v->size(); ++i)
    if (c.v->at(i) != curr_val) {
      out << c.id << "\t" << (curr_start + c.start) << "\t" << (i-1+c.start) << "\t" << curr_val << std::endl;
      curr_start = i;
      curr_val = c.v->at(i);
    }
  if ( (curr_start+1) != c.v->size()) // need to dump last one
    out << c.id << "\t" << (curr_start + c.start) << "\t" << (c.v->size()+c.start-1) << "\t" << curr_val << std::endl;
    
  return out;
}

uint16_t Coverage::getCoverageAtPosition(size_t pos) const {

  if (pos < start || pos > end) {
    std::cerr << "Coverage query out of bounds for location " << id << ":" << pos << std::endl;
    return 0;
  }
  
  size_t q = pos - start;
  if (q >= v->size()) {
    std::cerr << "Coverage query out of bounds for location " << id << ":" << pos << " with pos-start of " << q << " attempt on v of size " << v->size() << std::endl;
    return 0;
  }

  return (v->at(q));
    

}
