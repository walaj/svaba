#pragma once

#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include <unordered_map>

#include "svabaRead.h"
#include "SvabaSharedConfig.h"

class DiscordantCluster;

typedef std::unordered_map<std::string, svabaReadPtr> DiscordantReadMap;
typedef std::vector<svabaReadPtrVector> svabaReadClusterVector;
typedef std::vector<DiscordantCluster> DiscordantClusterVector;
typedef std::unordered_map<std::string, DiscordantCluster> DiscordantClusterMap;

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
  DiscordantCluster(const svabaReadPtrVector& this_reads,
		    const svabaReadPtrVector& all_reads,		    
		    const SeqLib::BamHeader& header);
  
  /** Is this discordant cluster empty? */
  bool isEmpty() const;
  
  /** Return a string representing the output file header */
  static std::string header() { 
    return "chr1\tpos1\tstrand1\tchr2\tpos2\tstrand2\ttcount\tncount\t\tmapq1\tmapq2\tcname\tid";
  }
  
  bool hasAssociatedAssemblyContig() const { return m_contig.length(); }
  
  void addMateReads(const svabaReadPtrVector& bav); 
  
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
  
  static std::unordered_map<std::string, DiscordantCluster>
  clusterReads(svabaReadPtrVector& bav,
	       const SeqLib::GenomicRegion& interval,
	       const SeqLib::BamHeader& header);

  void labelReads();
  
  static bool __add_read_to_cluster(svabaReadClusterVector &cvec,
				    svabaReadPtrVector &clust,
				    svabaReadPtr&, bool ismate);
  
  static void __cluster_reads(svabaReadPtrVector& brv,
			      svabaReadClusterVector& fwd,
			      svabaReadClusterVector& rev, int orientation);
  
  static void __cluster_mate_reads(svabaReadClusterVector& brcv,
				   svabaReadClusterVector& fwd,
				   svabaReadClusterVector& rev);
  
  static void __convertToDiscordantCluster(DiscordantClusterMap& dd,
					   const svabaReadClusterVector& cvec,
					   const svabaReadPtrVector& bav,
					   const SeqLib::BamHeader& header);

  static bool __valid_cluster(svabaReadPtrVector& clust, bool ismate);
  
  /** Query an interval against the two regions of the cluster. If the region overlaps
   * with one region, return the other region. This is useful for finding the partner 
   * region give a query region */
  SeqLib::GenomicRegion GetMateRegionOfOverlap(const SeqLib::GenomicRegion& gr) const; 
  
  int tcount = 0;
  int ncount = 0; 

  // supporting read counts per sample (e.g. t001 - 4, n001 - 6)  
  std::unordered_map<std::string, int> counts;

  DiscordantReadMap reads, mates;
  
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
  double __getMeanMapq(const DiscordantReadMap& m) const;
};

//! vector of AlignmentFragment objects
typedef std::vector<DiscordantCluster> DiscordantClusterVector;

//! Store a set of DiscordantCluster objects, indexed by the "id" field
typedef std::unordered_map<std::string, DiscordantCluster> DiscordantClusterMap;


