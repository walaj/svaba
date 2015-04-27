#ifndef DISCORDANT_CLUSTER_H
#define DISCORDANT_CLUSTER_H

#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>

#include "SnowTools/GenomicRegion.h"
#include "SnowTools/HTSTools.h"


// define a structure to hold discordant clusters
struct DiscordantCluster {

  std::string cluster;
  std::string id;
  size_t tcount = 0;
  size_t ncount = 0; 
  std::unordered_map<std::string, bool> qnames; // TODO get rid of it
  std::unordered_map<std::string, Read> reads;
  std::unordered_map<std::string, Read> mates;

  double reads_mapq; 
  double mates_mapq;
  std::vector<int> mapq; // TODO remoe
  std::string contig = "";

  static std::string header() { return "chr1\tpos1\tstrand1\tchr2\tpos2\tstrand2\ttcount\tncount\tmapq1\tmapq2\treads"; }

  SnowTools::GenomicRegion reg1;
  SnowTools::GenomicRegion reg2;

  DiscordantCluster() {}
  ~DiscordantCluster() {}
  DiscordantCluster(ReadVec &this_reads, ReadVec &all_reads);
  DiscordantCluster(std::string tcluster);
  DiscordantCluster(std::string tcluster, SnowTools::GenomicRegion gr1, SnowTools::GenomicRegion gr2) {
    cluster = tcluster;
    reg1 = gr1;
    reg2 = gr2;
  }

  void addMateReads(ReadVec &bav);

  // return the mean mapping quality for this cluster
  double getMeanMapq(bool mate) const;
  
  double getMeanMapq() const;
  std::string toRegionString() const;

  // add the read names supporting this cluster
  void addRead(std::string name);
  
  // define how to print this to stdout
  friend std::ostream& operator<<(std::ostream& out, const DiscordantCluster& dc);

  // define how to print to file
  std::string toFileString(bool with_read_names = false) const;

  // define how these are to be sorted
  bool operator < (const DiscordantCluster& b) const;


};

typedef std::unordered_map<std::string, DiscordantCluster> DMap;
typedef std::vector<DiscordantCluster> DVec;

#endif
