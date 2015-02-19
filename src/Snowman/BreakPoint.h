#ifndef SNOW_BREAKPOINT_H
#define SNOW_BREAKPOINT_H

#include <cstdlib>
#include <string>
#include <unordered_map>
#include "GenomicRegion.h"
#include "DiscordantCluster.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include <memory>

using namespace std;

struct BreakPoint;
typedef vector<BreakPoint> BPVec;
//typedef unordered_map<string, BreakPoint> BPMap;
//typedef shared_ptr<BamAlignment> BamAlignmentUP;
//typedef vector<BamAlignmentUP> BamAlignmentUPVector;

struct BreakPoint {

  static string header() { return "chr1\tpos1\tstrand1\tchr2\tpos2\tstrand2\tspan\tmapq1\tmapq2\tnsplit\ttsplit\tndisc\ttdisc\thomology\tinsertion\tcontig\tnumalign\tconfidence\tevidence\treads"; }

  // Discovar information
  size_t disco_tum = 0;
  size_t disco_norm = 0;
  bool discovar = false;
  
  // reads spanning this breakpoint
  BamAlignmentUPVector reads;

  // discordant reads supporting this aseembly bp
  DiscordantCluster dc;

  // breakpoints on the reference
  GenomicRegion gr1;
  GenomicRegion gr2;


  //unsigned pos1 = 0;
  //unsigned pos2 = 0;

  unsigned cpos1 = 0;  
  unsigned cpos2 = 0;

  //unsigned refID1 = 0;
  //unsigned refID2 = 0;

  string seq;

  string cname;

  string insertion;
  string homology;

  string id1;
  string id2;
  int matchlen1 = 0;
  int matchlen2 = 0;
  
  //char strand1;
  //char strand2;

  bool isSomatic = false;
  bool isGermline = false;

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

  int span;

  unsigned num_align = 0;

  bool part_of_local = false;

  bool local1 = false;
  bool local2 = false;

  string evidence = "";
  string confidence = "";

  BreakPoint(DiscordantCluster tdc);
  BreakPoint() {
    gr1.pos1 = 0;
    gr2.pos1 = 0;
  }

  static string BreakPointHeader();

  string toString() const; 
 
  bool sameBreak(BreakPoint &bp) const;

  void order();

  bool isEmpty() const { return (gr1.pos1 == 0 && gr2.pos1 == 0); }

  // return whether a bp is good to move on
  bool isGoodSomatic(int mapq, size_t tsplit_cutoff, size_t nsplit_cutoff) const;

  string toFileString();
  
  bool hasDiscordant() const;

  // return whether a bp is good to move on
  bool isGoodGermline(int mapq, size_t allsplit) const;

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
