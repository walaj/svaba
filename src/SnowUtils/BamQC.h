#ifndef BAMQC_H
#define BAMQC_H

#include <unordered_map>
#include <vector>
#include <cstdlib>
#include <iostream>

using namespace std;

struct BamQCReadGroup {

  BamQCReadGroup();

  size_t num_reads = 0;
  size_t num_supp = 0;
  size_t unmap = 0;  
  size_t qcfail = 0;
  size_t duplicate = 0;
  size_t supp = 0;
  vector<size_t> mapq;
  vector<size_t> nm;
  vector<size_t> isize;
  vector<size_t> as;
  vector<size_t> xp;
  vector<size_t> clip;
  vector<size_t> len;
  vector<size_t> phred;

  friend ostream& operator<<(std::ostream& out, const BamQCReadGroup& rg);
  
};

struct BamQC {
  
  unordered_map<string, BamQCReadGroup> map;
  friend ostream& operator<<(std::ostream& out, const BamQC& qc);

};


#endif
