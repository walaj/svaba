#ifndef BAMQC_H
#define BAMQC_H

#include <unordered_map>
#include <vector>
#include <cstdlib>
#include <iostream>
#include "api/BamReader.h"
#include "api/BamWriter.h"


using namespace std;

struct BamQCReadGroup {

  BamQCReadGroup();

  size_t num_reads;
  size_t num_supp;
  size_t unmap;  
  size_t qcfail;
  size_t duplicate;
  size_t supp;
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
  bool use = true; // hack so that you can pass qc but not use, in case you want to speed of bam processing without the BuildCharData

  // count an additional read
  void addRead(BamTools::BamAlignment &a);
  
};


#endif
