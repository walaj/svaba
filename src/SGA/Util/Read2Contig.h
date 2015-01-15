#ifndef READ2CONTIG_H
#define READ2CONTIG_H

#include <string>
#include <vector>
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include <sstream>

struct Read2Contig {

  // Constructors
  Read2Contig(std::string tcname, std::string trname, unsigned tpos, std::string tseq) : cname(tcname), rname(trname), pos(tpos), seq(tseq) {}
  Read2Contig(std::string tcname, std::string trname, unsigned tpos, BamTools::BamAlignment ta, std::string tseq) : cname(tcname), rname(trname), pos(tpos), a(ta), seq(tseq) {}
  Read2Contig() {}

  ~Read2Contig() {}

  //data
  std::string cname;
  std::string rname;
  int pos;
  BamTools::BamAlignment a;
  std::string seq;
  bool valid = true;
  double sw_score = -1;

  // set the sort order
  bool operator < (const Read2Contig& str) const { return (pos < str.pos); }

  // print it out
  std::string toString() const {
    std::stringstream ss;
    std::string sep = ",";
    ss << cname << sep << rname << sep << pos << sep << seq << sep << sw_score;
    return ss.str();
  }

};

typedef std::vector<Read2Contig> R2CVec;

#endif
