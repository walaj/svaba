#ifndef ALIGNED_CONTIG_H
#define ALIGNED_CONTIG_H

#include "GenomicRegion.h"
#include <algorithm>
#include "BreakPoint.h"
#include "AuxUtils.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/algorithms/Sort.h"
#include <memory>
#include "SnowUtils.h"

using namespace std;

typedef shared_ptr<BamTools::BamAlignment> BamAlignmentUP;
typedef vector<BamAlignmentUP> BamAlignmentUPVector;

typedef vector<BamTools::CigarOp> CigarOpVec;

struct Break {

  int cpos;
  int gpos;
  int refID;
  unordered_map<string, BamAlignmentUP> splits;
  unordered_map<string, BamAlignmentUP> discs;

};
typedef vector<Break> BreakVector;

struct ContigFragment {

  ContigFragment(int b1, int b2, int g1, int g2) {
    cpos = {b1, b2};
    gpos = {g1, g2};

    assert(g1 < g2);
    assert(b1 < b2);
  }

  bool operator < (const ContigFragment& b) const {
    return b.cpos.first < cpos.first;
  }

  pair<int,int> cpos;
  pair<int,int> gpos;

};

struct AlignedContig2 {

  BreakVector brv;
  BamAlignmentVector aln; // vector of CONTIG alignments
  string name;
  string seq;
  BamAlignmentUPVector reads; // vector or READS to contigs
  unordered_map<string, vector<size_t> > cov; 

  AlignedContig2(const string& sam, const BamTools::BamReader *treader);

  // convert a cigar string to a CigarOp
  CigarOpVec stringToCigar(const string& val);
  
  // parse tag data from a sam record, add to BamAlignment
  void parseTags(const string& val, BamAlignment &a);
  
  // perform SW alignment of reads to contig
  void alignReadsToContigs(BamAlignmentUPVector &bav);
  
  // convert cigars of alignment to generate break points
  void setBreaks();

  CigarOpVec orientCigar(const BamTools::BamAlignment& align);
};


#endif
