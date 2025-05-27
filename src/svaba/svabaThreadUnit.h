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

#include "SeqLib/BWAIndex.h"
#include "SeqLib/BWAAligner.h"

using WalkerMap = std::map<std::string, svabaBamWalker>;

class svabaThreadUnit {
  
public:
  
  //svabaThreadUnit() = default;
  ~svabaThreadUnit() = default;
  
  svabaThreadUnit(
		  BWAIndexPtr idx,
		  const std::string& reference,
		  const std::map<std::string,std::string>& bamFiles,
		  const SvabaOptions& opts
		  ) {

    // load the samtools faidx
    ref_genome_ = std::make_unique<SeqLib::RefGenome>();
    ref_genome_->LoadIndex(reference);

    // set the *shared* memory BWAIndex
    bwa_aligner = SeqLib::BWAAligner(idx);

    // load the .bai files for each bam
    for (const auto&p : bamFiles) {
      auto& walker = walkers_[p.first];
      walker.Open(p.second);
      walker.bwa_aligner = bwa_aligner;
      
      
    }

    // set the parameters
    
    
  }

  // local version of aligner class, but will hold shared memory index
  SeqLib::BWAAligner bwa_aligner;
  
  // results
  std::vector<AlignedContig>                 m_alc;
  SeqLib::BamRecordVector                    m_contigs;
  BPVec                                      m_bps;
  DiscordantClusterMap                       m_disc;
  size_t                                     m_bamreads_count = 0;
  size_t                                     m_disc_reads     = 0;
  SeqLib::GRC                                badd; // bad regions

  // non-copyable, movable
  svabaThreadUnit(const svabaThreadUnit&) = delete;
  svabaThreadUnit& operator=(const svabaThreadUnit&) = delete;
  svabaThreadUnit(svabaThreadUnit&&) noexcept = default;
  svabaThreadUnit& operator=(svabaThreadUnit&&) noexcept = default;

  void clear() {
    m_alc.clear();
    m_contigs.clear();
    m_bps.clear();
    m_disc.clear();
    m_bamreads_count = 0;
    m_disc_reads     = 0;
    badd.clear();
  }

  bool MemoryLimit(size_t readLimit, size_t contLimit) const {
    return m_bamreads_count > readLimit
        || m_contigs.size()   > contLimit
        || m_disc_reads       > readLimit;
  }

private:

  // store the BAM .bai indicies for for this thread
  WalkerMap                            walkers_;

  // store the faidx index for this thread
  std::unique_ptr<SeqLib::RefGenome>   ref_genome_;  
};
