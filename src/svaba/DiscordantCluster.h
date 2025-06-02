#pragma once

#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include <unordered_map>

#include "svabaRead.h"
#include "SvabaSharedConfig.h"

typedef std::vector<svabaReadVector> svabaReadClusterVector;

/** Class to hold clusters of discordant reads */
class DiscordantCluster 
{
  
  friend struct BreakPoint;
  
public:
  
  /** Create an empty cluster */
  DiscordantCluster() { 
    m_reg1 = SeqLib::GenomicRegion(); 
    m_reg2 = SeqLib::GenomicRegion(); 
    assert(m_reg1.IsEmpty()); 
    ncount = 0; tcount = 0;
    mapq1 = -1; mapq2 = -1;
  }
  
  /** Make a cluster from a set of reads (pre-clustered) and look up a larger set to find 
   * their mates 
   * @param this_reads Pre-clustered set of discordant reads (but not their mates)
   * @param all_reads A pile of reads to search for mates
   */
  DiscordantCluster(const svabaReadVector& this_reads,
		    const svabaReadVector& all_reads,
		    int max_mapq_possible);
  
  /** Is this discordant cluster empty? */
  bool isEmpty() const;
  
  /** Return a string representing the output file header */
  static std::string header() { 
    return "chr1\tpos1\tstrand1\tchr2\tpos2\tstrand2\ttcount\tncount\ttcount_hq\tncount_hq\t\tmapq1\tmapq2\tcname\tregion_string\treads\tcompeting_id"; 
  }
  
  bool hasAssociatedAssemblyContig() const { return m_contig.length(); }
  
  void addMateReads(const svabaReadVector& bav);
  
  /** Return the discordant cluster as a string with just coordinates */
  std::string toRegionString(const SeqLib::BamHeader& h) const;
  
  /** Return the ID associated with this cluster */
  std::string ID() const { return m_id; } 
  
  /** Print this with region string and read counts and mapq */
  //friend std::ostream& operator<<(std::ostream& out, const DiscordantCluster& dc);
  
  /** Print this with region string and read counts and mapq */
  std::string print(const SeqLib::BamHeader& h) const;
  
  /** Return as a string for writing to a file */
  std::string toFileString(const SeqLib::BamHeader& h, bool with_read_names) const;
  
  /** Sort by coordinate */
  bool operator < (const DiscordantCluster& b) const;
  
  /** Is this a valid cluster? */
  bool valid() const;
  
  static void __remove_singletons(svabaReadClusterVector& b);
  
  static std::unordered_map<std::string, DiscordantCluster>
  clusterReads(svabaReadVector& bav,
	       const SeqLib::GenomicRegion& interval,
	       int max_mapq_possible); 
  
  static bool __add_read_to_cluster(svabaReadClusterVector &cvec, svabaReadVector &clust, const svabaRead &a, bool mate);
  
  static void __cluster_reads(const svabaReadVector& brv, svabaReadClusterVector& fwd, svabaReadClusterVector& rev, int orientation);
  
  static void __cluster_mate_reads(svabaReadClusterVector& brcv, svabaReadClusterVector& fwd, svabaReadClusterVector& rev);
  
  static void __convertToDiscordantCluster(std::unordered_map<std::string, DiscordantCluster> &dd, const svabaReadClusterVector& cvec, const svabaReadVector& bav, int max_mapq_possible);
  
  /** Query an interval against the two regions of the cluster. If the region overlaps
   * with one region, return the other region. This is useful for finding the partner 
   * region give a query region */
  SeqLib::GenomicRegion GetMateRegionOfOverlap(const SeqLib::GenomicRegion& gr) const; 
  
  int tcount = 0;
  int ncount = 0; 
  
  int tcount_hq = 0;
  int ncount_hq = 0;
  
  int max_possible_mapq = 0;
  
  std::unordered_map<std::string, int> counts; // supporting read counts per sample (e.g. t001 - 4, n001 - 6)
  
  std::unordered_map<std::string, svabaRead> reads;
  std::unordered_map<std::string, svabaRead> mates;
  
  std::string m_contig = "";
  
  double read_score = 0;
  double mate_score = 0;
  
  //int rp_orientation = -1; // FR, FF,  RR, RF
  
  SeqLib::GenomicRegion m_reg1;
  SeqLib::GenomicRegion m_reg2;
  
  int mapq1;
  int mapq2;
  
  std::string m_id_competing; // id of discordant cluster with same span, different strands
  
private:    
  std::string m_id;
  
  // return the mean mapping quality for this cluster
  double __getMeanMapq(bool mate = false) const;
};

//! vector of AlignmentFragment objects
typedef std::vector<DiscordantCluster> DiscordantClusterVector;

//! Store a set of DiscordantCluster objects, indexed by the "id" field
typedef std::unordered_map<std::string, DiscordantCluster> DiscordantClusterMap;


