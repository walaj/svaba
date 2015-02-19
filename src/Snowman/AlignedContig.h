#ifndef ALIGNED_CONTIG
#define ALIGNED_CONTIG

#include "Util.h" //TODO
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
using namespace BamTools;

typedef shared_ptr<BamAlignment> BamAlignmentUP;
typedef vector<BamAlignmentUP> BamAlignmentUPVector;

typedef vector<BamTools::BamAlignment> BAVec;
typedef vector<BamTools::CigarOp> CigarOpVec;

typedef std::vector<std::string> StringVec;

class AlignedContig;
typedef unordered_map<string, AlignedContig> ContigMap;
typedef vector<AlignedContig> AlignedContigVec;

struct CAlignment {

  BamTools::BamAlignment align;
  CigarOpVec cigar;
  string cigstring; 
  int break1 = -1;
  int break2 = -1;
  int gbreak1 = -1;
  int gbreak2 = -1;

  unsigned start;

  BamAlignmentUPVector reads_b1;
  BamAlignmentUPVector reads_b2;

  unsigned nsplit1 = 0;
  unsigned tsplit1 = 0;
  unsigned nsplit2 = 0;
  unsigned tsplit2 = 0;

  size_t ndisc = 0;
  
  bool ca_local = false; // is an alignmend to anchor window

  // converts the BamTools cigar format to string
  string cigarToString(CigarOpVec cig);

  // CONSTRUCTOR
  CAlignment(BamTools::BamAlignment talign, CigarOpVec tcigar);

  // define how to sort these 
  bool operator < (const CAlignment& str) const { return (start < str.start); }

  // define how to sort these 
  friend ostream& operator<<(ostream &out, const CAlignment& c); 

};

// define a way to order the contigs by start
struct AlignmentOrdering {
  inline bool operator() (const CAlignment& struct1, const CAlignment& struct2) {
    return (struct1.start < struct2.start);
  }
};

typedef vector<CAlignment> AlignVec;

class AlignedContig {

 public:  

  // constructor taking in a BAM record from the contig BAM
  AlignedContig(const BamTools::BamAlignment align); 
  AlignedContig() {}
  ~AlignedContig() {}
  AlignedContig(std::string sam, const BamTools::BamReader * reader, GenomicRegion &twindow);

  bool hasvariant = false;
  BPVec m_breaks;  
  BreakPoint m_farbreak;
  BreakPoint m_farbreak_filt; // global breakpoint, but remove bad fragments (e.g. 60, 60, 0 mapqs, keep 60-60 connection)
  //vector<BamAlignment> m_bamreads;
  BamAlignmentUPVector m_bamreads;
  AlignVec m_align;
  string samrecord;
  GenomicRegion window;

  // new
  BamAlignmentVector aln;

  void parseTags(const string& val, BamTools::BamAlignment &a);
  CigarOpVec stringToCigar(const string& val);

  void splitIndelCoverage();

  bool hasVariant() const { return hasvariant; }

  int maxSplit() const;

  // add a new contig alignment
  void addAlignment(const BamTools::BamAlignment align);

  void fillExtraReads();

  void addDiscordantCluster(DiscordantCluster dc) { m_dc.push_back(dc); } 

  int numDiscordantClusters() const { return m_dc.size(); }

  string printDiscordantClusters() const;

  void setIndelBreaks(CAlignment &align);

  // return the name of the contigs
  string getContigName() const { return m_align[0].align.Name; }
 
  // return the number of alignments
  int getNumPrimaryAlign() const { return m_align.size(); }
  int getNumSecondaryAlign() const { return m_align_second.size(); }

  int getMaxMapq() const { 
    int m = -1;
    for (auto& i : m_align)
      if (i.align.MapQuality > m)
	m = i.align.MapQuality;
    return m;
  }

  int getMinMapq() const { 
    int m = 1000;
    for (auto& i : m_align)
      if (i.align.MapQuality < m)
	m = i.align.MapQuality;
    return m;
  }

  //bool hasLocal() const { return m_local; }

  /*int getNumReads(const char type) const { 
    int countr = 0;
    for (R2CVec::const_iterator it = m_reads.begin(); it != m_reads.end(); it++) 
      if (it->rname.at(0) == type)
	countr++;
    return countr; 
    };*/

  // sort the read2contigs
  void sortReads();

  // make work function for getting PER-ALIGNMENT breaks from BWA-MEM
  void setBreaks(CAlignment &align);

  void splitCoverage();

  string printForR() const;

  // set whether this alignment has a somatic breakpoint
  //void setSomatic();

  string printAlignments() const;
  void printContigFasta(ofstream &ostream) const;

  bool isGermline() const { return m_germline; };
  bool isSomatic() const { return m_somatic; };

  bool intersect(const AlignedContig *al) const; 

  // flips the cigar if the contig is aligned to the opposite strand
  CigarOpVec orientCigar(const BamTools::BamAlignment align);

  // parses the contig file name to determine where the anchor window was
  //void setWindow(const string s);

  // 
  //Window getWindow() const { return m_window; }

  // find the breakpoint pairs by looping through ALL the alignments
  void getBreakPairs();
   
  // get the break pairs, but update incase things have changed
  BPVec getBreaks() const { 
    return m_breaks; 
  }

  // add masked seq
  void addMaskedSeq(string seq) { m_masked_seq  = seq; }

  // add repeat masker
  void addRepeatMasker(RepeatMasker rp) { 
    m_rep_vec.push_back(rp); 
    sort(m_rep_vec.begin(), m_rep_vec.end());
  }

  // add SBlat alignment
  void addSBlat(SBlat b) { 
    m_blat_vec.push_back(b); 
    sort(m_blat_vec.begin(), m_blat_vec.end());
  }
  
  // run all the algorithms for updating the breakpoints, in case other parts changed
  /*void updateBreakpointData(bool skip_realign, bool no_r2c_matched) {

    // realign the reads to the contigs

    if (!skip_realign) 
      realignReads();
    else
      readR2Creads();

    // sort the reads for better visualization
    if (!no_r2c_matched)
      sortReads(); 

    // get the split coverage
    splitCoverage();

    // get the break pairs
    getBreakPairs();

    //debug
    if (m_breaks.size() == 0)
      cerr << "m_align size: "<< m_align.size() << endl;

    // if we don't need to write the output, clear most of the read data
    if (no_r2c_matched) {
      for (BamAlignmentVector::iterator it = m_bamreads.begin(); it != m_bamreads.end(); it++) {
	it->QueryBases = "";
	it->Qualities = "";
	it->RemoveTag("JW");
	it->RemoveTag("TS");
      }
    }

  }
  */
  void readR2Creads();

  BreakPoint getGlobalBreak() const { assert(m_farbreak.gr1.mapq <= 60 && m_farbreak.gr2.mapq <= 60); return m_farbreak; }

  // doing a more stringent alignment of reads to the contig
  // this is to remove normals that ruin somatic calls
  void realignReads();

  string getName() const { return m_align[0].align.Name; }

  void alignReadsToContigs(BamAlignmentUPVector &bav);

  AlignVec getAlignments() const { return m_align; }

  void setSomatic(const bool somatic) { m_somatic = somatic; }

  void setGermline(const bool germline) { m_germline = germline; }

  void addBams(string tum, string norm, string pan, string r2c) { 
    tbam = tum;
    nbam = norm;
    pbam = pan;
    rbam = r2c;
  }

  //vector<BamAlignment> getBamReads() const { return m_bamreads; }

  void settleContigs();

  string getSequence() const { return m_seq; }
  

  friend ostream& operator<<(ostream &out, const AlignedContig &ac);
 private:

  AlignVec m_align_second;

  bool m_somatic = false;
  bool m_germline = false;
  string m_seq = "";
  string m_masked_seq = "";
  RepeatMaskerVec m_rep_vec;
  SBlatVec m_blat_vec;
  int mapq_threshold = 60;
  // store the raw alignments. Move into AlignVec if keeping the contig

  // the bam files that
  string tbam, nbam, pbam, rbam;
  //SVBamReader treader, nreader, preader;

  vector<DiscordantCluster> m_dc;

};

// define a sorter to sort by the Mate Position
typedef std::binary_function<BamAlignmentUP, BamAlignmentUP, bool> AlignmentSortBase;
struct ByALTag : public AlignmentSortBase {
  
  // ctor
  ByALTag(const Algorithms::Sort::Order& order = Algorithms::Sort::AscendingOrder)
    : m_order(order) { }

  // comparison function
  bool operator()(const BamAlignmentUP& lhs, const BamAlignmentUP& rhs) {

    int alL = SnowUtils::GetIntTag(lhs, "AL").back();
    int alR = SnowUtils::GetIntTag(rhs, "AL").back();

    // otherwise sort on reference ID
    return Algorithms::Sort::sort_helper(m_order, alL, alR);
  }
  // used by BamMultiReader internals
  static inline bool UsesCharData(void) { return false; }
  
  // data members
  private:
  const Algorithms::Sort::Order m_order;
};


#endif


