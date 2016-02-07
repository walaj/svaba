#ifndef DISCORDANT_CLUSTER_H__
#define DISCORDANT_CLUSTER_H__

#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <cassert>

#include "SnowTools/GenomicRegion.h"
#include "SnowTools/BamRead.h"

namespace SnowTools 
{

  //class DiscordantCluster;
  
  /** Class to hold clusters of discordant reads */
  class DiscordantCluster 
  {

    friend class BreakPoint;
    
  public:

    /** Create an empty cluster */
    DiscordantCluster() { 
      m_reg1 = GenomicRegion(); 
      m_reg2 = GenomicRegion(); 
      assert(m_reg1.isEmpty()); 
      ncount = 0; tcount = 0;
      mapq1 = -1; mapq2 = -1;
    }

    /** Make a cluster from a set of reads (pre-clustered) and look up a larger set to find 
     * their mates 
     * @param this_reads Pre-clustered set of discordant reads (but not their mates)
     * @param all_reads A pile of reads to search for mates
     */
    DiscordantCluster(const BamReadVector& this_reads, const BamReadVector& all_reads);
    
    /** Is this discordant cluster empty? */
    bool isEmpty() const;

    /** Return a string representing the output file header */
    static std::string header() { 
      return "chr1\tpos1\tstrand1\tchr2\tpos2\tstrand2\ttcount\tncount\tmapq1\tmapq2\treads"; 
    }
    
    bool hasAssociatedAssemblyContig() const { return m_contig.length(); }

    void addMateReads(const BamReadVector& bav);
    
    /** Return the discordant cluster as a string with just coordinates */
    std::string toRegionString() const;
    
    /** Add the read names supporting this cluster */
    void addRead(std::string name);
    
    /** Print this with region string and read counts and mapq */
    friend std::ostream& operator<<(std::ostream& out, const DiscordantCluster& dc);
    
    /** Return as a string for writing to a file */
    std::string toFileString(bool with_read_names = false) const;
    
    /** Sort by coordinate */
    bool operator < (const DiscordantCluster& b) const;

    static std::unordered_map<std::string, DiscordantCluster> clusterReads(const BamReadVector& bav, const GenomicRegion& interval);

    static bool __add_read_to_cluster(BamReadClusterVector &cvec, BamReadVector &clust, const BamRead &a, bool mate);

    static void __cluster_reads(const BamReadVector& brv, BamReadClusterVector& fwd, BamReadClusterVector& rev);

    static void __cluster_mate_reads(BamReadClusterVector& brcv, BamReadClusterVector& fwd, BamReadClusterVector& rev);

    static void __convertToDiscordantCluster(std::unordered_map<std::string, DiscordantCluster> &dd, const BamReadClusterVector& cvec, const BamReadVector& bav);

    /** Query an interval against the two regions of the cluster. If the region overlaps
     * with one region, return the other region. This is useful for finding the partner 
     * region give a query region */
    GenomicRegion GetMateRegionOfOverlap(const GenomicRegion& gr) const; 

    int tcount = 0;
    int ncount = 0; 
    std::unordered_map<std::string, int> counts;

    std::unordered_map<std::string, BamRead> reads;
    std::unordered_map<std::string, BamRead> mates;

    std::string m_contig = "";

    double read_score = 0;
    double mate_score = 0;
    
    GenomicRegion m_reg1;
    GenomicRegion m_reg2;

    int mapq1;
    int mapq2;

  private:    
    std::string m_id;
    std::unordered_map<std::string, bool> qnames; // TODO get rid of it

    // return the mean mapping quality for this cluster
    double __getMeanMapq(bool mate = false) const;
  };
  
  //! vector of AlignmentFragment objects
  typedef std::vector<DiscordantCluster> DiscordantClusterVector;
  
  //! Store a set of DiscordantCluster objects, indexed by the "id" field
  typedef std::unordered_map<std::string, DiscordantCluster> DiscordantClusterMap;

  
}

#endif
