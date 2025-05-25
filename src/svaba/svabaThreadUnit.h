#pragma once

#include <map>
#include <memory>
#include <vector>
#include <mutex>

#include "svabaBamWalker.h"
#include "AlignedContig.h"
#include "DiscordantCluster.h"
#include "BreakPoint.h"
#include "SeqLib/RefGenome.h"

using WalkerMap = std::map<std::string, svabaBamWalker>;

struct svabaThreadUnit;
void WriteFilesOut(svabaThreadUnit&);

struct svabaThreadUnit {

  // per-thread BAM walkers and reference genomes
  WalkerMap                                  walkers;
  std::unique_ptr<SeqLib::RefGenome>         ref_genome;

  // results
  std::vector<AlignedContig>                 m_alc;
  SeqLib::BamRecordVector                    m_contigs;
  BPVec                                      m_bps;
  DiscordantClusterMap                       m_disc;
  size_t                                     m_bamreads_count = 0;
  size_t                                     m_disc_reads     = 0;
  SeqLib::GRC                                badd; // bad regions

  svabaThreadUnit() = default;
  ~svabaThreadUnit() = default;

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
};
