#ifndef VARIANT_BAM_READER_H
#define VARIANT_BAM_READER_H

#include "MiniRules.h"
#include "GenomicRegion.h"
#include "BamQC.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include <unordered_map>
#include <memory>

typedef shared_ptr<BamAlignment> BamAlignmentUP;
typedef vector<BamAlignmentUP> BamAlignmentUPVector;

using namespace std;
using namespace BamTools;

// Phred score transformations
inline int char2phred(char b) {
  uint8_t v = b;
  assert(v >= 33);
  return v - 33;
}

/////////////// 
// Hold read counts
//////////////
struct ReadCount {

  int keep = 0;
  int total = 0;
  
  int percent () const {
    int perc  = SnowUtils::percentCalc<int>(keep, total); 
    return perc;
  }

  string totalString() const {
    return SnowUtils::AddCommas<int>(total);
  }

  string keepString() const {
    return SnowUtils::AddCommas<int>(keep);
  }

};

//////////
// read a variant bam and write it to disc or store in memory
/////////
class VariantBamReader {

 public:
  VariantBamReader() {}
  ~VariantBamReader() {
    delete m_writer;
    delete m_reader;
  }

  void setPrefix(string prefix) { m_prefix = prefix; } 
  
  VariantBamReader(string inbam, string outbam, MiniRulesCollection* mr, int verbose);

  // remove duplicate reads by name and alignment flag
  static void deduplicateReads(const BamAlignmentVector &inbav, BamAlignmentVector &outbav);

  static unsigned getClipCount(BamAlignment &a);
  static void qualityTrimRead(int qualTrim, string &seq, string &qual);

  //bool writeVariantBam(BamQC &qc, bool qc_only);
  //bool writeVariantBam(BamQC &qc, BamAlignmentVectorP &bav);
  
  // set which part of the bam to read
  bool setBamRegion(GenomicRegion gp);

  // create the index file for the output bam
  void MakeIndex();

  // print to stdout
  void printMessage(const ReadCount &rc_main, const BamAlignment * a) const;

 private:
  
  string m_bam;
  string m_out;
  BamReader * m_reader;
  BamWriter * m_writer;
  GenomicRegion m_region;
  MiniRulesCollection * m_mr;
  int m_verbose;

  string m_prefix = "";
};

typedef unordered_map<string, VariantBamReader*> VariantBamReaderMap;


#endif 
