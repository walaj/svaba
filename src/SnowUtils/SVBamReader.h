#ifndef SVBAMREADER_H
#define SVBAMREADER_H

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include <vector>
#include <string>
#include "GenomicRegion.h"
#include <unordered_map>

#include "api/api_global.h"
#include "api/BamAlignment.h"
#include "api/BamReader.h"
#include "api/BamMultiReader.h"
#include "api/algorithms/Sort.h"

using namespace std;
using namespace BamTools;

typedef vector<BamTools::BamAlignment> BamAlignmentVector;

class SVBamReader {

 public:

  bool clipOnly = false;
  
  SVBamReader() {}

  SVBamReader(const string filename, const string prefix, const int isize = 2000, 
	      const unsigned mapq = 30, const int qualthresh = 4, const int mol = 40, 
	      const bool skip_supp = false, const bool skip_r2 = false, 
	      const int minclip = 10, const unsigned verbose = 1) 
    : m_bam(filename), m_prefix(prefix), m_isize(isize), m_mapq(mapq), m_qualthresh(qualthresh), 
    m_mol(mol), m_skip_supp(skip_supp), m_skip_r2(skip_r2), m_minclip(minclip), m_verbose(verbose) {}

  ~SVBamReader() { }

  SVBamReader(const string filename, const string prefix, const int isize = 2000, 
	      const unsigned mapq = 30, const unsigned verbose = 1)
    : m_bam(filename), m_prefix(prefix), m_isize(isize), 
    m_mapq(mapq), m_verbose(verbose) {}
  
    // appropriate constructor for reading full bam (no filtering)
    SVBamReader(const string filename) : m_bam(filename) {}

    bool disc_cluster_only = false;

  static bool IsTumorRead(const BamAlignment &a);

  bool setBamRegion(int refid, int pos1, int pos2);
  bool setBamRegion(GenomicRegion gr);

  bool findBamIndex();

  bool contigBamToBAVec(BamAlignmentVector &bav);

  void softClip(int qualTrim, string &seq, string const &qual, unsigned clipnum, bool &tooshort, bool &not_real_clip, int &startpoint, int &endpoint);

  bool updateCounters(BamAlignment a, unsigned clipnum);

  bool passNtest(string str) const;

  bool GetR2Read(BamAlignment const &a, string const &r2, string &q2, string &trim_seq);

  bool bamToBAVec(BamAlignmentVector &bav);

  void setReadLimit(unsigned lim) { limit = lim; };
  void setNMLimit(unsigned nmlim) { m_nmlim = nmlim; };

  static string findBamIndex(string bam);

  static void getRefVector(string bamfile, RefVector &ref);

  static void getSamHeader(string bamfile, SamHeader &sam);

  // deduplicate the reads by name
  static void deduplicateReads(const BamAlignmentVector &inbav, BamAlignmentVector &outbav);
  
  // clean out the tags
  static void clearFinalTags(BamAlignment * align);
  
  // deduplicate the reads by position (deals with reads which are not correctly marked as duplicates)
  static void deduplicateReadsPos(const BamAlignmentVector &inbav, BamAlignmentVector &outbav);
  
  static unsigned getClipCount(BamTools::BamAlignment a);
  
  bool bamToBAVecCluster(BamAlignmentVector &bav);
  
  string getBamFilename() { return m_bam;}
  
  bool R2CbamToBAVec(BamAlignmentVector &bav);
  
  void setSkipR2(bool val) { m_skip_r2 = val; }
  
  size_t discordantCount(GenomicRegion const &gr1, GenomicRegion const &gr2, int span);

  bool preprocessBam(BamTools::BamWriter &writer);

  int perclimit = 30;

 private:

  string m_bam;
  string m_prefix;
  int m_isize;
  unsigned m_mapq;
  unsigned m_qualthresh;
  unsigned m_mol; // min overlap. Used to determine if to throw seq away
  unsigned m_skip_supp = false;
  unsigned m_skip_r2 = false;
  string m_bai;
  string m_bam_bai;
  BamTools::BamReader m_reader;
  //unsigned m_region[4];
  GenomicRegion m_region;
  size_t m_minclip = 10;
  size_t m_verbose;
  
  //keep track of what is happening
  unsigned disc_num = 0;
  unsigned unmap_num = 0;
  unsigned clip_num = 0;
  unsigned disc_clip_num = 0;
  unsigned mapq_num = 0;
  unsigned used_num = 0;
  unsigned total_reads = 0;
  unsigned low_qual = 0;
  unsigned non_num = 0;

  // maximimum number of reads before aborting
  unsigned limit = 1000000000;
  unsigned m_nmlim = 3;

  

};

// define a sorter to sort by the Mate Position
typedef std::binary_function<BamAlignment, BamAlignment, bool> AlignmentSortBase;
struct ByMatePosition : public AlignmentSortBase {
  
  // ctor
  ByMatePosition(const Algorithms::Sort::Order& order = Algorithms::Sort::AscendingOrder)
    : m_order(order) { }
  
  // comparison function
  bool operator()(const BamTools::BamAlignment& lhs, const BamTools::BamAlignment& rhs) {
    
    // force unmapped aligmnents to end
    if ( lhs.MateRefID == -1 ) return false;
    if ( rhs.MateRefID == -1 ) return true;
    
    // if on same reference, sort on position
    if ( lhs.MateRefID == rhs.MateRefID )
      return Algorithms::Sort::sort_helper(m_order, lhs.MatePosition, rhs.MatePosition);
    
    // otherwise sort on reference ID
    return Algorithms::Sort::sort_helper(m_order, lhs.MateRefID, rhs.MateRefID);
  }
  // used by BamMultiReader internals
  static inline bool UsesCharData(void) { return false; }
  
  // data members
  private:
  const Algorithms::Sort::Order m_order;
};


#endif
