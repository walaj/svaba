#ifndef GENOMIC_REGION_H
#define GENOMIC_REGION_H

#include <vector>
#include <regex>
#include <iostream>
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include <unordered_map>
#include "SnowUtils.h"
#include "IntervalTree.h"

#include <utility>
#include <list>

using namespace std;
using namespace BamTools;

typedef signed long long slong;

// forward declare 
struct GenomicRegion;
typedef vector<GenomicRegion> GenomicRegionVector;
typedef list<GenomicRegion> GenomicRegionList;
typedef EInterval<GenomicRegion> GenomicInterval;
typedef unordered_map<int, vector<GenomicInterval> > GenomicIntervalMap;
typedef EIntervalTree<GenomicRegion> GenomicIntervalTree;
typedef unordered_map<int, GenomicIntervalTree> GenomicIntervalTreeMap;
typedef vector<GenomicInterval> GenomicIntervalVector;

// declare stuff for 2D genome (for cluster pairs
typedef unsigned long long uslong;
typedef EInterval<uslong> Genomic2DInterval;
typedef vector<Genomic2DInterval> Genomic2DIntervalVector;

static const vector<string> CHR_NAME {"1", "2", "3", "4", "5", "6", "7", "8", "9",
				  "10", "11", "12", "13", "14", "15", "16", "17", 
                                  "18", "19", "20", "21", "22", "X", "Y", "M"};
static const vector<string> CHR_NAME_NUM {"1", "2", "3", "4", "5", "6", "7", "8", "9",
				  "10", "11", "12", "13", "14", "15", "16", "17", 
                                  "18", "19", "20", "21", "22", "23", "24"};

static const int CHR_LEN [25] = {249250621, 243199373, 198022430, 191154276, //1-4
				 180915260, 171115067, //5-6
				 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, //7-12
				 115169878, 107349540, 102531392, 90354753,  81195210,  78077248, //13-18
				 59128983,  63025520,  48129895,  51304566,  155270560, 59373566, //19-24
                                 16571}; //25

static const uslong CHR_CLEN [25] = {0, 249250621, 492449994,  690472424, 881626700, 1062541960, 1233657027,
				    1392795690,1539159712,1680373143,1815907890,1950914406,2084766301,
				    2199936179, 2307285719, 2409817111, 2500171864, 2581367074, 2659444322,                                                                                                                   
				    2718573305, 2781598825, 2829728720, 2881033286, 3036303846, 3095677412};

static string const REFHG19 = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";

static const int NONCENT_CHR [44] = {1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,
    11,11,12,12,13,14,15,16,16,17,17,18,18,19,19,20,20,21,21,22,23,23,24,24};

struct GenomicRegion {

  GenomicRegion() {};
  GenomicRegion(int t_chr, int t_pos1, int t_pos2);
  GenomicRegion(int t_chr, int t_pos1, int t_pos2, char t_strand, int t_mapq) : 
    chr(t_chr), pos1(t_pos1), pos2(t_pos2), strand(t_strand), mapq(t_mapq) {}
  GenomicRegion(string t_chr, string t_pos1, string t_pos2);

  // members
  int chr = 0;
  int pos1 = 0;
  int pos2 = 0;
  int ncount = 0;
  int tcount = 0;
  char strand = '*';
  int mapq = 0;
  string cluster;
  string id;
  int nm;
  

  static GenomicRegionVector non_centromeres;

  static int chrToNumber(string ref);
  static string chrToString(int ref);
  static GenomicRegionVector regionFileToGRV(string file, int pad);
  static GenomicRegionVector mergeOverlappingIntervals(const GenomicRegionVector &grv);
  static GenomicIntervalTreeMap createTreeMap(const GenomicRegionVector &grv);
  static void sendToBed(const GenomicRegionVector &grv, const string file);

  // divide a region into pieces of width and overlaps
  GenomicRegionVector divideWithOverlaps(int width, int ovlpo);

  static uslong posToBigPos(int refid, int pos);

  // checks whether a GenomicRegion is empty
  bool isEmpty() const;

  // define how these are to be sorted
  bool operator < (const GenomicRegion& b) const;
  bool operator==(const GenomicRegion& b) const;
  bool operator<=(const GenomicRegion &b) const;

  static bool GenomicIntervalLessThan (const GenomicInterval& a, const GenomicInterval& b);
  
  friend ostream& operator<<(std::ostream& out, const GenomicInterval& gi);

  // return regions containing the whole genome
  static GenomicRegionVector getWholeGenome();

  // determine if something overlaps with centromere 
  int centromereOverlap() const;

  // check if there is an overlap
  int getOverlap(const GenomicRegion gr) const;

  friend ostream& operator<<(std::ostream& out, const GenomicRegion& gr);
  string toString() const;
  void pad(int pad);
  
  int width() const;

};

#endif
