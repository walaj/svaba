#ifndef SNOWMAN_WORK_UNIT_H__
#define SNOWMAN_WORK_UNIT_H__

// hold outputs for a single (or set on the same thread) of local assemblies

#include <map>
#include "SnowmanBamWalker.h"
#include "AlignedContig.h"
#include "DiscordantCluster.h"
#include "BreakPoint.h"
#include "DiscordantCluster.h"
#include "SeqLib/RefGenome.h"

typedef std::map<std::string, SnowmanBamWalker> WalkerMap;

struct SnowmanWorkUnit {
  
  // its own thread-safe versions of readers and genomes
  WalkerMap walkers;
  SeqLib::RefGenome * ref_genome = nullptr;
  SeqLib::RefGenome * vir_genome = nullptr;
  SeqLib::GRC m_bad_regions;// bad region tracker for this thread
  
  // other structures to hold results
  std::vector<AlignedContig> m_alc;
  SeqLib::BamRecordVector m_contigs, m_vir_contigs;
  BPVec m_bps;
  DiscordantClusterMap m_disc;
  size_t m_bamreads_count = 0;
  size_t m_disc_reads = 0;
  SeqLib::GRC badd;

  void clear() {
    m_alc.clear();
    m_contigs.clear();
    m_vir_contigs.clear();
    m_bps.clear();
    m_disc.clear();
    m_bamreads_count = 0;
  }
  
  bool MemoryLimit(size_t read, size_t cont) const {
    const size_t readlim = read;
    const size_t contlim = cont;
    return m_bamreads_count > readlim || m_contigs.size() > contlim || m_vir_contigs.size() > contlim || m_disc_reads > readlim;
  }
  
  ~SnowmanWorkUnit() {
    clear();
    if (ref_genome)
      delete ref_genome;
    if (vir_genome)
      delete vir_genome;
  }
  
};


#endif
