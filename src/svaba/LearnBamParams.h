// LearnBamParams.h
#pragma once

#include <string>
#include <map>
#include <ostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <array>

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

  /// Detect a recurrent 5' soft-clip "tag" (untrimmed UMI / inline barcode /
  /// adapter or library spacer) from the accumulated clip5 histogram. Sets
  /// tag_trim_5p / tag_frac / tag_seq. Cheap; uses only the clip5_* counters.
  void detectTag();

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

  // statistics (set by computeStats)
  size_t n_isize_pairs = 0;  // number of FR pairs used (after 98% trim)
  double isize_median = 0;
  double sd_isize = 0;      // SD computed around the median

  // --- 5' soft-clip tag detection (set by addRead + detectTag) ---
  // Histogram of 5' soft-clip length (read's 5' end: leading clip for forward
  // reads, trailing clip for reverse, since SEQ is stored reference-forward).
  std::array<size_t, 16> clip5_hist{};
  size_t clip5_total = 0;    // mapped reads with a CIGAR considered (frac denom)
  std::unordered_map<std::string, size_t> clip5_seq;  // forward-read clipped bases (fixed vs random)
  int    tag_trim_5p = 0;    // detected tag length, bp (0 = none); the trim cutoff
  double tag_frac    = 0.0;  // fraction of reads carrying the modal 5' clip
  std::string tag_seq;       // dominant clipped sequence if fixed; "" if variable (UMI)

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
  int tag_trim_5p = 0;  // max detected 5' tag length across this BAM's RGs (0 = none)

  /// Write per-RG isize distributions to a TSV for R plotting.
  /// Called after learnParams() and before computeStats() clears the vectors.
  /// Output: ${prefix}.learn.tsv.gz with columns: bam, rg, isize
  void dumpLearnData(const std::string& prefix) const;

  /// Path of the BAM this object learned from.
  const std::string& bamPath() const { return bam_; }

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

// --- insert-size param cache (cloud scatter: learn once, reuse per shard) ------
// writeBamParams: dump every learned BAM's per-RG isize_median/sd_isize (+ the
//   per-bam readlen/mapq/isize maxima) to a small TSV.
// loadBamParams: restore them into sc.bamStats, matched to the current run's BAMs
//   by PATH (so the prefix scheme can differ), letting svaba skip the genome-wide
//   insert-size sweep. Returns true only if EVERY BAM in the current run is covered.
void writeBamParams(const SvabaSharedConfig& sc, const std::string& file);
bool loadBamParams(SvabaSharedConfig& sc, const std::string& file);
