#ifndef DISCORDANT_CLUSTER_H
#define DISCORDANT_CLUSTER_H

#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include "GenomicRegion.h"
#include "VariantBamReader.h"

using namespace std;

// define a structure to hold discordant clusters
struct DiscordantCluster {

  string cluster;
  string id;
  size_t tcount = 0;
  size_t ncount = 0; 
  unordered_map<string, bool> qnames; // TODO get rid of it
  unordered_map<string, BamAlignmentUP> reads;
  unordered_map<string, BamAlignmentUP> mates;
  vector<int> mapq;
  string contig = "";

  GenomicRegion reg1;
  GenomicRegion reg2;

  DiscordantCluster() {}
  ~DiscordantCluster() {}
  DiscordantCluster(BamAlignmentUPVector &this_reads, BamAlignmentUPVector &all_reads);
  DiscordantCluster(string tcluster);
  DiscordantCluster(string tcluster, GenomicRegion gr1, GenomicRegion gr2) {
    cluster = tcluster;
    reg1 = gr1;
    reg2 = gr2;
  }

  void addMateReads(BamAlignmentUPVector &bav);

  // return the mean mapping quality for this cluster
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
