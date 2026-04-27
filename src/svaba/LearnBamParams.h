// LearnBamParams.h
#pragma once

#include <string>
#include <map>
#include <ostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
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
 *
 * Sampling strategy (SvABA2.0): coordinate-sorted BAMs cluster read groups
 * by genomic position, so scanning from the start only sees RGs present in
 * early chromosomes. To ensure all RGs are represented:
 *   1. Build sampling windows at the midpoint of each reference contig
 *   2. For each window, read up to `reads_per_window` reads (100k default)
 *   3. Track per-RG counts; skip reads from already-saturated RGs
 *   4. Stop early once all header-declared RGs are satisfied
 * This replaces the old sequential-from-start scan that missed 20/21 RGs.
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
  double isize_max = 0;

  /// Write per-RG isize distributions to a TSV for R plotting.
  /// Called after learnParams() and before computeStats() clears the vectors.
  /// Output: ${prefix}.learn.tsv.gz with columns: bam, rg, isize
  void dumpLearnData(const std::string& prefix) const;

private:

  /// Process reads from the current reader position, updating rg_count
  /// and bam_read_groups. Returns number of reads consumed.
  size_t consumeReads(size_t max_reads,
		      std::unordered_map<std::string, size_t>& rg_count,
		      std::unordered_set<std::string>& satisfied_rgs,
		      const std::vector<std::string>& groups);

 std::string       bam_; // file path of the BAM
 SvabaSharedConfig&  sc;
 std::shared_ptr<SeqLib::BamReader> reader_;

};
