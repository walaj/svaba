#ifndef VCFRECORD_H
#define VCFRECORD_H

#include <string>
#include <vector>

using namespace std;

struct VCFRecord {
  
  int chr;
  int pos;
  string record;

  VCFRecord() {}
  ~VCFRecord() {}
  VCFRecord(int t_chr, int t_pos, string t_record) : chr(t_chr), pos(t_pos), record(t_record) {}
  
  bool operator < (const VCFRecord &v) const {
    bool out = false;
    if (chr < v.chr)
      out = true;
    else if (chr == v.chr && pos < v.pos)
      out = true;
    return out;
  }

};

typedef vector<VCFRecord> VCFRecordVector;

#endif
