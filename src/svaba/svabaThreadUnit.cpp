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
  //m_bamreads_count = 0;
  // m_disc_reads     = 0;
  badd.clear();

  all_weird_reads.clear();
  all_discordant_reads.clear();
  all_corrected_reads.clear();  
}

bool svabaThreadUnit::MemoryLimit(size_t readLimit, size_t contLimit) const {

  size_t stored_reads = all_weird_reads.size() +
    all_discordant_reads.size() +
    all_corrected_reads.size();
  
  bool mem_exceeded = stored_reads > readLimit ||
                      master_contigs.size()   > contLimit;
  
  if (mem_exceeded) {
    sc.logger.log(true, true, "...writing files on thread ",
		  threadId, " with limit hit of ",
		  stored_reads, " reads and ", master_contigs.size(),
		  " contigs");
  }
  
  return mem_exceeded;
}

