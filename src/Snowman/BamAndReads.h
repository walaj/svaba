#ifndef BAM_AND_READS_H
#define BAM_AND_READS_H

#include <vector>
#include "GenomicRegion.h"
#include "VariantBamReader.h"
#include "MiniRules.h"
#include <time.h>
#include <memory>
#include "MinRead.h"

using namespace std;
using namespace BamTools;

typedef shared_ptr<BamAlignment> BamAlignmentUP;
typedef vector<BamAlignmentUP> BamAlignmentUPVector;

// a smaller region to assemble
struct AssemblyRegion {

  AssemblyRegion() {}
  AssemblyRegion(GenomicRegion gr) { region = gr; }
  ~AssemblyRegion() {
    //for (vector<BamAlignment*>::iterator it = reads.begin() ; it != reads.end(); ++it) 
    //  delete (*it);
    //reads.clear();
  }
  
  BamAlignmentUPVector reads;
  MinReadSPVector mrv;
  GenomicRegion region;

  GenomicRegionVector partner_windows;
  GenomicIntervalTreeMap tree_pw;
};

typedef vector<AssemblyRegion> AssemblyRegionVector;
typedef unique_ptr<AssemblyRegion> AssemblyRegionUP;
typedef vector<AssemblyRegionUP> AssemblyRegionUPVector;


// store a reader of a bam file and vectors of pointers
// to reads which live on the heap, for both anchor and partner
struct BamAndReads {

  BamAndReads() {}
  ~BamAndReads() { 
    delete reader; 
    arvec.clear();
    //for (AssemblyRegionVector::iterator it = arvec.begin(); it != arvec.end(); it++)
    //  delete (*it);
    //arvec.clear();
  }
  BamAndReads(GenomicRegion gr, MiniRulesCollection *tmr, int tverb, string tbam, string tprefix);

  BamReader * reader;

  int read_time = 0; // timer in seconds for reading BAM
  int unique_reads = 0; // total number of unique reads to be assembled
  int reads = 0; // total number of reads to be assembled (allows doubles, etc)

  int mate_read_time = 0; // timer in seconds for reading BAM
  int mate_unique_reads = 0; // total number of unique reads to be assembled
  int mate_reads = 0; // total number of reads to be assembled (allows doubles, etc)

  // what interval is this defined on
  GenomicRegion interval;
  GenomicIntervalTreeMap tree;

  GenomicRegionVector mate_regions;
  GenomicIntervalTreeMap tree_with_mate;

  // define how to read the bam
  MiniRulesCollection * mr;

  int verbose = 1;
  string bam;
  string prefix;

  // store all of the mini regions to run
  AssemblyRegionUPVector arvec;

  // store all of the discordant reads that map outside of BIGCHUNK
  // we'll need these to lookup extra reads for assembly
  BamAlignmentUPVector disc;

  // break up the interval into small pieces
  void divideIntoRegions();

  // read the bam
  void readBam();

  // read the bam
  void readMateBam();

  // get mate regions
  void calculateMateRegions();

  void _read_bam(BamAlignmentVector &reads, int limit);

  // use the interval tree to see if a read should be added to a group, based on itself or its mate
  void addRead(BamAlignment &a);
  void addMateRead(BamAlignment &a);

  friend ostream& operator<<(ostream &out, BamAndReads &bar);
  
};

#endif
