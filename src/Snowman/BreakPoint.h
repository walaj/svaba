#ifndef SNOW_BREAKPOINT_H
#define SNOW_BREAKPOINT_H

#include <cstdlib>
#include <string>
#include <unordered_map>
#include "GenomicRegion.h"

struct BreakPoint;
typedef vector<BreakPoint> BPVec;
typedef unordered_map<string, BreakPoint> BPMap;

struct BreakPoint {

  bool isBest = false; // marked for best breakpoint when mulitple are redundant
 
  size_t disco_tum = 0;
  size_t disco_norm = 0;

  //string idcommon = "";
  //string pairid = "";
  
  bool discovar = false;

  DiscordantCluster dc;

  unsigned pos1;
  unsigned pos2;

  unsigned cpos1;  
  unsigned cpos2;

  unsigned refID1;
  unsigned refID2;

  string seq;

  string cname;

  string insertion;
  string homology;

  string id1;
  string id2;
  int matchlen1 = 0;
  int matchlen2 = 0;
  
  char strand1;
  char strand2;

  bool isSomatic = false;
  bool isGermline = false;

  unsigned mapq1; 
  unsigned mapq2; 

  unsigned tsplit1 = 0;
  unsigned tsplit2 = 0;

  unsigned nsplit1 = 0;
  unsigned nsplit2 = 0;

  size_t nsplit = 0;
  size_t tsplit = 0;

  unsigned tall = 0;
  unsigned nall = 0; 

  int nm1 = 0;
  int nm2 = 0;

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
  BreakPoint() {}

  static string BreakPointHeader();

  string toString() const; 
 
  bool sameBreak(BreakPoint &bp) const;

  void order();

  // return whether a bp is good to move on
  bool isGoodSomatic(int mapq, size_t tsplit_cutoff, size_t nsplit_cutoff) const;

  bool hasDiscordant() const;

  // return whether a bp is good to move on
  bool isGoodGermline(int mapq, size_t allsplit) const;

  // define how to sort these 
  bool operator < (const BreakPoint& bp) const { 
    return bp.refID1 > refID1 || 
       (bp.refID1 == refID1 && bp.pos1 > pos1) || // low pos is first
       (bp.refID1 == refID1 && bp.pos1 == pos1 && bp.pos2 == pos2 && nsplit1 > bp.nsplit1) || // if same, check nsplit
       (bp.refID1 == refID1 && bp.pos1 == pos1 && bp.pos2 == pos2 && nsplit1 == bp.nsplit1 && tsplit1 > bp.tsplit1); // if also same, check tsplit
  }
  friend ostream& operator<<(std::ostream& out, const BreakPoint& bp) { out << bp.toString(); return out; }

  // print to file
  void printToFile(ofstream &of, const BamAlignmentVector &bamreads);

};

#endif
