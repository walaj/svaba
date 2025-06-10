#pragma once

#include <map>
#include <memory>
#include <vector>
#include <mutex>

#include "svabaBamWalker.h"
#include "svabaLogger.h"
#include "svabaOptions.h"
#include "AlignedContig.h"
#include "DiscordantCluster.h"
#include "BreakPoint.h"
#include "SeqLib/RefGenome.h"

#include "SvabaSharedConfig.h"
#include "SeqLib/BWAAligner.h"

using WalkerMap = std::map<std::string, svabaBamWalker>;

class svabaThreadUnit {
  
public:
  
  //svabaThreadUnit() = default;
  ~svabaThreadUnit() = default;
  
  svabaThreadUnit(SvabaSharedConfig& sc_);
  
  // local version of aligner class, but will hold shared memory index
  std::shared_ptr<SeqLib::BWAAligner> bwa_aligner; //(sc.bwa_idx);

  size_t processed_count = 0;
  size_t total_count = 0; // total to process
  size_t processed_since_memory_dump = 0;
  
  // results
  std::vector<AlignedContig>                 master_alc;
  SeqLib::BamRecordVector                    master_contigs;
  BPVec                                      m_bps;
  DiscordantClusterMap                       m_disc;
  //size_t                                     m_bamreads_count = 0;
  //size_t                                     m_disc_reads     = 0;
  SeqLib::GRC                                badd; // bad regions
  int                                        threadId;

  // very verbose outpout
  svabaReadVector                            all_weird_reads;
  svabaReadVector                            all_discordant_reads;
  SeqLib::BamRecordVector                            all_corrected_reads;    
  
  // store the BAM .bai indicies for for this thread
  WalkerMap                            walkers;
  
  // non-copyable, movable
  svabaThreadUnit(const svabaThreadUnit&) = delete;
  svabaThreadUnit& operator=(const svabaThreadUnit&) = delete;
  svabaThreadUnit(svabaThreadUnit&&) noexcept = default;
  svabaThreadUnit& operator=(svabaThreadUnit&&) noexcept = default;

  void clear();

  bool MemoryLimit(size_t readLimit, size_t contLimit) const;

  // store the faidx index for this thread
  std::unique_ptr<SeqLib::RefGenome>   ref_genome;

  SvabaSharedConfig& sc;

private:

};
