#ifndef SNOW_BREAKPOINT_H
#define SNOW_BREAKPOINT_H

#include <cstdlib>
#include <string>
#include <unordered_map>
#include "GenomicRegion.h"
#include "DiscordantCluster.h"

#include "reads.h"

using namespace std;

struct BreakPoint;
typedef vector<BreakPoint> BPVec;
typedef unordered_map<string, size_t> PON;
//typedef unordered_map<string, BreakPoint> BPMap;
//typedef shared_ptr<BamAlignment> BamAlignmentUP;
//typedef vector<BamAlignmentUP> BamAlignmentUPVector;

/**
 *
 */
void runRefilterBreakpoints(int argc, char** argv);

/**
 *
 */
void parseBreakOptions(int argc, char** argv);

struct BreakPoint {

  static string header() { return "chr1\tpos1\tstrand1\tchr2\tpos2\tstrand2\tspan\tmapq1\tmapq2\tnsplit\ttsplit\tndisc\ttdisc\tncigar\ttcigar\thomology\tinsertion\tcontig\tnumalign\tconfidence\tevidence\treads\tpon_samples\trepeat_seq"; }

  // Discovar information
  size_t disco_tum = 0;
  size_t disco_norm = 0;
  bool discovar = false;
  
  // reads spanning this breakpoint
  ReadVec reads;

  // discordant reads supporting this aseembly bp
  DiscordantCluster dc;

  // breakpoints on the reference
  GenomicRegion gr1;
  GenomicRegion gr2;

  string read_names;

  int cpos1 = 0;  
  int cpos2 = 0;

  string seq;

  string cname;

  string insertion = "";
  string homology = "";

  string repeat_seq = "";

  string id1;
  string id2;
  int matchlen1 = 0;
  int matchlen2 = 0;

  size_t tcigar = 0;
  size_t ncigar = 0;

  size_t min_end_align_length = 0; // minimum length of alignment on end. real short are not to be trusted
  
  //char strand1;
  //char strand2;

  bool isSomatic = false;
  bool isGermline = false;

  size_t pon = 0;

  //unsigned mapq1; 
  //unsigned mapq2; 

  unsigned tsplit1 = 0;
  unsigned tsplit2 = 0;

  unsigned nsplit1 = 0;
  unsigned nsplit2 = 0;

  size_t nsplit = 0;
  size_t tsplit = 0;

  unsigned tall = 0;
  unsigned nall = 0; 

  //int nm1 = 0;
  //int nm2 = 0;

  unsigned num_dups = 0;
   
  //Window window;
  GenomicRegion window;

  //int span;

  unsigned num_align = 0;

  bool part_of_local = false;

  bool local1 = false;
  bool local2 = false;

  string evidence = "";
  string confidence = "";

  bool isindel = false;

  BreakPoint(DiscordantCluster tdc);
  BreakPoint() {
    gr1.pos1 = 0;
    gr2.pos1 = 0;
  }

  BreakPoint(string &line);

  
  static void readPON(string &file, unique_ptr<PON> &pmap);
  
  /*! @function determine if the breakpoint has split read support
   * @param reference to a vector of read smart pointers that have been aligned to a contig
   * @discussion Note: will cause an error if the AL tag not filled in for the reads. 
   * The AL tag is filled in by AlignedContig::alignReadsToContigs.
   */
  void splitCoverage(ReadVec &bav);

  /*! @function get the span of the breakpoints (in bp). -1 for interchrom
   * @return int distance between breakpoints
   */
  int getSpan() const { 
    if (isindel && insertion == "" )// deletion
      return (abs(gr1.pos1 - gr2.pos1) - 1);
    if (isindel)
      return (insertion.length()); // insertion
    if (gr1.chr == gr2.chr)
      return abs(gr1.pos1-gr2.pos1);
    else
      return -1;
  }

  /*! @function check the breakpoint against a panel of normals
   * @param Panel of normals hash
   * @return number of normal samples with this variant
   */
  int checkPon(unique_ptr<PON> &p);
  
  /*! @function get a unique string representation of this breakpoint.
   * Format for indel is chr_breakpos_type (eg. 0_134134_I)
   * @return string with breakpoint info
   */
  string getHashString() const;

  bool hasMinimal() const;

  string toString() const; 
 
  bool sameBreak(BreakPoint &bp) const;

  void order();

  bool isEmpty() const { return (gr1.pos1 == 0 && gr2.pos1 == 0); }

  // return whether a bp is good to move on
  bool isGoodSomatic(int mapq, size_t tsplit_cutoff, size_t nsplit_cutoff) const;

  string toFileString(bool noreads = false);
  
  bool hasDiscordant() const;

  // return whether a bp is good to move on
  bool isGoodGermline(int mapq, size_t allsplit) const;

  bool operator==(const BreakPoint& bp) const;

  // define how to sort these 
  bool operator < (const BreakPoint& bp) const { 
    return (gr1 < bp.gr1); // || (gr1 == gr2 && nsplit1 > bp.nsplit1)
    //(bp.gr1.ref == refID1 && bp.pos1 > pos1) || // low pos is first
      //  (bp.refID1 == refID1 && bp.pos1 == pos1 && bp.pos2 == pos2 && nsplit1 > bp.nsplit1) || // if same, check nsplit
      // (bp.refID1 == refID1 && bp.pos1 == pos1 && bp.pos2 == pos2 && nsplit1 == bp.nsplit1 && tsplit1 > bp.tsplit1); // if also same, check tsplit
  }
  friend ostream& operator<<(std::ostream& out, const BreakPoint& bp) { out << bp.toString(); return out; }

  // print to file
  void printToFile(ofstream &of, const BamAlignmentVector &bamreads);

};

#endif
