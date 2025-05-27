// LearnBamParams.h
#pragma once

#include <string>
#include <map>
#include <ostream>
#include "SvabaSharedConfig.h"

/**
 * BamParams collects and stores basic summary statistics for a single read group:
 *   - readlen:       the maximum read length observed
 *   - max_mapq:      the maximum mapping quality observed
 *   - mean_isize:    the mean insert size of properly paired reads
 *   - sd_isize:      the standard deviation of the insert size distribution
 *
 * After reading through a BAM file, LearnBamParams::collectStats() fills
 * these fields and they are used downstream to set thresholds for
 * SV/indel calling (e.g., discordant read size cutoffs).
 */
struct BamParams {
  int readlen      = 0;
  int max_mapq     = 0;
  double mean_isize= 0.0, sd_isize = 0.0;

  void collectStats();
  
  friend std::ostream& operator<<(std::ostream& out, const BamParams& p);
};
using BamParamsMap = std::map<std::string,BamParams>;


/**
 * BamLearningResult aggregates the insert-size learning results for all input BAM files:
 *   - perFile:                   maps each sample name to its BamParamsMap (per-read group statistics)
 *   - globalReadLen:             the maximum read length observed across all files
 *   - globalMaxMapQ:             the maximum mapping quality observed across all files
 *   - globalMinDiscordantSize:   the highest discordant-read size cutoff computed 
 *                               (mean+SD*cutoff) across all read groups in all files
 */
struct BamLearningResult {
  std::map<std::string, BamParamsMap> perFile;
  int globalReadLen             = 0;
  int globalMaxMapQ             = 0;
  int globalMinDiscordantSize   = 0;
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
  LearnBamParams(const SvabaSharedConfig& sc,
		 const std::string& bamPath);

  // map of RG : params, for a single bam
  BamParamsMap learnParams();  // scan this->bam_, return per RG stats

  /// Scan *all* the given BAMs and return both per-file and global summaries
  static BamLearningResult learnAll(const SvabaSharedConfig& sc)
    
private:
  std::string       bam_;
  const SvabaSharedConfig&  sc_; 
  SeqLib::BamReader reader_;
};
