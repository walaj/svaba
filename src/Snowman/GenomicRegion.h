#ifndef GENOMIC_REGION_H
#define GENOMIC_REGION_H

#include <vector>
#include <regex>
#include <iostream>
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include <unordered_map>
#include "IntervalTree.h"

using namespace std;
using namespace BamTools;

typedef signed long long slong;

// forward declare 
struct GenomicRegion;
typedef vector<GenomicRegion> GenomicRegionVector;

static const vector<string> CHR_NAME {"1", "2", "3", "4", "5", "6", "7", "8", "9",
				  "10", "11", "12", "13", "14", "15", "16", "17", 
                                  "18", "19", "20", "21", "22", "X", "Y", "M"};
static const vector<string> CHR_NAME_NUM {"1", "2", "3", "4", "5", "6", "7", "8", "9",
				  "10", "11", "12", "13", "14", "15", "16", "17", 
                                  "18", "19", "20", "21", "22", "23", "24"};

string const REFHG19 = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";

static const int CHR_LEN [25] = {249250621, 243199373, 198022430, 191154276, //1-4
				 180915260, 171115067, //5-6
				 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, //7-12
				 115169878, 107349540, 102531392, 90354753,  81195210,  78077248, //13-18
				 59128983,  63025520,  48129895,  51304566,  155270560, 59373566, //19-24
                                 16571}; //25

static const slong CHR_CLEN [25] = {0, 249250621, 492449994,  690472424, 881626700, 1062541960, 1233657027,
				       1392795690,1539159712,1680373143,1815907890,1950914406,2084766301,
				       2199936179, 2307285719, 2409817111, 2500171864, 2581367074, 2659444322,                                                                                                                   
				       2718573305, 2781598825, 2829728720, 2881033286, 3036303846, 3095677412};


static const int NONCENT_CHR [44] = {1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,
    11,11,12,12,13,14,15,16,16,17,17,18,18,19,19,20,20,21,21,22,23,23,24,24};

struct GenomicRegion {

  GenomicRegion() {};
  GenomicRegion(int t_chr, int t_pos1, int t_pos2);
  GenomicRegion(int t_chr, int t_pos1);
  GenomicRegion(string t_chr, string t_pos1, string t_pos2);

  // members
  int chr = 0;
  int pos1 = 0;
  int pos2 = 0;
  int tcount = 0;
  int ncount = 0;
  string cluster;

  //unordered_map<string, BamAlignment> reads;
  unordered_map<string, size_t> rname;
  vector<int> mapq;
  char strand = '*';
  slong abspos1;
  slong abspos2;
  
  static GenomicRegionVector non_centromeres;

  static slong convertPos(unsigned refid, unsigned pos, bool revstrand = false);
  static int chrToNumber(string ref);
  static string chrToString(int ref);
  static slong RPtoNum(string rp);

  // intersect region with a vector
  GenomicRegionVector intersection(GenomicRegionVector &grv);

  // checks whether a GenomicRegion is empty
  bool isEmpty() const;

  // define how these are to be sorted
  bool operator < (const GenomicRegion& b) const;

  // return regions containing the whole genome
  static GenomicRegionVector getWholeGenome();

  // determine if something overlaps with centromere 
  int centromereOverlap() const;

  // determine if something overlaps with blacklist regions
  int blacklistOverlap() const;

  // check if there is an overlap
  int getOverlap(const GenomicRegion gr) const;

  friend ostream& operator<<(std::ostream& out, const GenomicRegion& gr);
  string toString() const;
  string toStringOffset() const;
  void pad(int pad);
  
  int width() const;

  // parse a contig name to get window
  // formatted as contig_8:12123-123123_whatever
  GenomicRegion (const string s);
    
};


// Simplified BED format that only uses the first 4 columns.
// Copyright (c) 2013-2014 Kamil Slowikowski
// See LICENSE for GPLv3 license.
class BEDRow {
 public:
  std::string name;
  GenomicRegion gr;
  void readNextRow(std::istream & stream);
};

/*static std::istream & operator>>(std::istream & stream, BEDRow & a) {
  a.readNextRow(stream);
  return stream;
}*/

// define the interval trees
typedef IntervalT<GenomicRegion> GenomicInterval;
typedef vector<GenomicInterval> GenomicIntervalVector;
typedef IntervalTTree<GenomicRegion> GenomicTree;


// define a structure for holding tumor and normal counts for
// Discovar contigs
struct ATuple {
  size_t tum = 0;
  size_t norm= 0;

  // constructor from strings
  ATuple(string ttum, string nnorm) { 
    tum = static_cast<size_t>(stoi(ttum));
    norm = static_cast<size_t>(stoi(nnorm));
  }

  // constructor from ints
  ATuple(int ttum, int nnorm) { 
    tum = static_cast<size_t>(ttum);
    norm = static_cast<size_t>(nnorm);
  }


};

// define a structure to hold discordant clusters
struct DiscordantCluster {

  string cluster;
  size_t tcount = 0;
  size_t ncount = 0; 
  unordered_map<string, bool> qnames;
  vector<int> mapq;
  string contig = "";

  GenomicRegion reg1;
  GenomicRegion reg2;

  DiscordantCluster() {}
  ~DiscordantCluster() {}
  DiscordantCluster(string tcluster);

  // return the mean mapping quality for this cluster
  double getMeanMapq() const;

  // add the read names supporting this cluster
  void addRead(string name);
  
  // define how to print this to stdout
  friend ostream& operator<<(std::ostream& out, const DiscordantCluster& dc);

  // define how to print to file
  string toFileString() const;

  // define how these are to be sorted
  bool operator < (const DiscordantCluster& b) const;


};

typedef unordered_map<string, DiscordantCluster> DMap;
typedef vector<DiscordantCluster> DVec;

/////////////////////////
/////////////////////////
#endif
