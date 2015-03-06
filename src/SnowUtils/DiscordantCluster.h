#ifndef DISCORDANT_CLUSTER_H
#define DISCORDANT_CLUSTER_H

#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include "GenomicRegion.h"
#include "VariantBamReader.h"

#include "reads.h"

using namespace std;

// define a structure to hold discordant clusters
struct DiscordantCluster {

  string cluster;
  string id;
  size_t tcount = 0;
  size_t ncount = 0; 
  unordered_map<string, bool> qnames; // TODO get rid of it
  unordered_map<string, Read> reads;
  unordered_map<string, Read> mates;

  double reads_mapq; 
  double mates_mapq;
  vector<int> mapq; // TODO remoe
  string contig = "";

  static string header() { return "chr1\tpos1\tstrand1\tchr2\tpos2\tstrand2\ttcount\tncount\tmapq1\tmapq2\treads"; }

  GenomicRegion reg1;
  GenomicRegion reg2;

  DiscordantCluster() {}
  ~DiscordantCluster() {}
  DiscordantCluster(ReadVec &this_reads, ReadVec &all_reads);
  DiscordantCluster(string tcluster);
  DiscordantCluster(string tcluster, GenomicRegion gr1, GenomicRegion gr2) {
    cluster = tcluster;
    reg1 = gr1;
    reg2 = gr2;
  }

  void addMateReads(ReadVec &bav);

  // return the mean mapping quality for this cluster
  double getMeanMapq(bool mate) const;
  
  double getMeanMapq() const;
  string toRegionString() const;

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

#endif
