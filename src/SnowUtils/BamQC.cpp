#include "BamQC.h"

template <typename T> void printQCVec(ostream &out, const vector<T> &vec) {
  for (auto it = vec.begin(); it != vec.end(); it++)
    out << "," << *it;
  return;
}

// instantiate a new read group
BamQCReadGroup::BamQCReadGroup() {

  mapq = vector<size_t>(61,0); // 60 elems of 0
  nm   = vector<size_t>(102,0); // 101 elems of 0
  isize= vector<size_t>(2001,0); // (everything above 2000 is inter)
  clip = vector<size_t>(102,0); 
  as   = vector<size_t>(102,0);
  xp   = vector<size_t>(102,0);
  len  = vector<size_t>(102,0);
  phred= vector<size_t>(61, 0);
}

// make the output
std::ostream& operator<<(std::ostream& out, const BamQC& qc) {

  for (auto it = qc.map.begin(); it != qc.map.end(); it++)
    out << "READGROUP:" << it->first << endl << it->second << endl;
  return out;
}

// make the output
std::ostream& operator<<(std::ostream& out, const BamQCReadGroup& rg) {

  out << "total," << rg.num_reads << endl;
  out << "unmap," << rg.unmap << endl;
  out << "qcfail," << rg.qcfail << endl;
  out << "duplicate," << rg.duplicate << endl;
  out << "supplementary," << rg.supp << endl;

  out << "mapq";
  printQCVec<size_t>(out, rg.mapq);
  out << endl;

  out << "nm";
  printQCVec<size_t>(out, rg.nm);
  out << endl;

  out << "isize";
  printQCVec<size_t>(out, rg.isize);
  out << endl;

  out << "as";
  printQCVec<size_t>(out, rg.as);
  out << endl;
  
  out << "xp";
  printQCVec<size_t>(out, rg.xp);
  out << endl;

  out << "clip";
  printQCVec<size_t>(out, rg.clip);
  out << endl;

  out << "len";
  printQCVec<size_t>(out, rg.len);
  out << endl;

  out << "phred";
  printQCVec<size_t>(out, rg.phred);

  return out;
}

