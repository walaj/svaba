#ifndef SEQLIB_BAMSTATS_H__
#define SEQLIB_BAMSTATS_H__

#include <unordered_map>
#include <cstdint>
#include <iostream>

#include "Histogram.h"
#include "SeqLib/BamRecord.h"

  /** Store information pertaining to a given read group *
   *
   * This class will collect statistics on number of: read, supplementary reads, unmapped reads, qcfail reads, duplicate reads.
   * It will also create Histogram objects to store counts of: mapq, nm, isize, clip, mean phred score, length
   */
class BamReadGroup {

  friend class BamStats;
  
 public:

  /** Construct an empty BamReadGroup */
  BamReadGroup() {}

  /** Construct an empty BamReadGroup for the specified read group
   * @param name Name of the read group
   */
  BamReadGroup(const std::string& name);

  /** Display basic information about this read group
   */
  friend std::ostream& operator<<(std::ostream& out, const BamReadGroup& rg);

  /** Add a BamRecord to this read group */
  void addRead(SeqLib::BamRecord &r);

 private:

  // count number of reads 
  size_t reads;
  size_t supp;
  size_t unmap;  
  size_t qcfail;
  size_t duplicate;
  size_t mate_unmap;

  int mapq_max = 0;
  int realen_max = 0;
  std::vector<int> isize_vec;
  
  std::string m_name;

};

/** Class to store statistics on a BAM file.
 *
 * BamStats currently stores a map of BamReadGroup objects. Bam statistics
 * are collected then on a read-group basis, but can be output in aggregate. See
 * BamReadGroup for description of relevant BAM statistics.
 */
class BamStats
{

 public:
  
  /** Loop through the BamReadGroup objections and print them */
  friend std::ostream& operator<<(std::ostream& out, const BamStats& qc);

  /** Add a read by finding which read group it belongs to and calling the 
   * addRead function for that BamReadGroup.
   */
  void addRead(SeqLib::BamRecord &r);

 private:
  
  std::unordered_map<std::string, BamReadGroup> m_group_map;

};

#endif
