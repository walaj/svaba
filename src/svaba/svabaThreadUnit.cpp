#include "svabaThreadUnit.h"

svabaThreadUnit::svabaThreadUnit(SvabaSharedConfig& sc_) : sc(sc_) {
  
  // load the samtools faidx
  ref_genome = std::make_unique<SeqLib::RefGenome>();
  ref_genome->LoadIndex(sc.opts.refGenome);
  
  // set the *shared* memory BWAIndex
  bwa_aligner = std::make_shared<SeqLib::BWAAligner>(sc.bwa_idx);
  
  // load the .bai files for each bam and set parameters
  for (const auto&p : sc.opts.bams) {
    
    // instantiate a new BAM walker
    auto [it, inserted] = walkers.try_emplace(p.first, sc);
    svabaBamWalker& walker = it->second;
    walker.Open(p.second);
  }
  
}

void svabaThreadUnit::clear() {
  master_alc.clear();
  master_contigs.clear();
  m_bps.clear();
  m_disc.clear();
  m_bamreads_count = 0;
  m_disc_reads     = 0;
  badd.clear();
}

bool svabaThreadUnit::MemoryLimit(size_t readLimit, size_t contLimit) const {
  return m_bamreads_count > readLimit
    || master_contigs.size()   > contLimit
    || m_disc_reads       > readLimit;
}

