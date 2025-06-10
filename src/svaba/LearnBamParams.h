// LearnBamParams.h
#pragma once

#include <string>
#include <map>
#include <ostream>
#include <vector>
#include <unordered_map>
#include <memory>

class SvabaSharedConfig;
namespace SeqLib {
  class BamRecord;
  class BamReader;
}

/** Store information pertaining to a given read group *
 *
 * This class will collect statistics on number of: read, supplementary reads, unmapped reads, qcfail reads, duplicate reads.
 * It will also create Histogram objects to store counts of: mapq, nm, isize, clip, mean phred score, length
 */
class BamReadGroup {
  
 public:

  /** Construct an empty BamReadGroup */
  BamReadGroup() {}

  /** Add a BamRecord to this read group */
  void addRead(const SeqLib::BamRecord &r);

  void computeStats();

  friend std::ostream& operator<<(std::ostream& os, const BamReadGroup& bg);
  
  // count number of reads 
  size_t reads = 0;
  size_t supp = 0;
  size_t unmap = 0;  
  size_t qcfail = 0;
  size_t duplicate = 0;
  size_t mate_unmap = 0;

  int mapq_max = 0;
  int readlen_max = 0;
  std::vector<uint32_t> isize_vec;

  // statistics
  double isize_mean = 0;
  double sd_isize = 0;
  
};


/**
 * LearnBamParams encapsulates logic to learn insert-size and related statistics
 * from a single BAM file or across multiple BAM files. It:
 *   - Opens and iterates through each read in the BAM (via SeqLib::BamReader)
 *   - Groups reads by read-group and collects metrics (read length, MAPQ, insert-size distribution)
 *   - Computes per-read-group statistics (mean, median, SD of insert size, coverage, clip fraction, etc.)
 *   - Provides a static helper to run this process on a set of BAM files and
 *     produce both per-file (per RG) and global summary statistics (max read length,
 *     max mapQ, and a discordant read size cutoff).
 */
class LearnBamParams {

public:

  /// Learn parameters from a single BAM file
  LearnBamParams(SvabaSharedConfig& sc_,
		 const std::string& bamPath); // opens the BAM

  // map of RG : params, for a single bam
  void learnParams();  // learn the RGs

  // this is map of RG-name : BamReadGroup
  std::unordered_map<std::string, BamReadGroup> bam_read_groups;

  // per BAM readlen / mapq max
  int readlen_max = 0;
  int mapq_max = 0;
  
private:
 
 std::string       bam_; // file path of the BAM
 SvabaSharedConfig&  sc; 
 std::shared_ptr<SeqLib::BamReader> reader_;

};
