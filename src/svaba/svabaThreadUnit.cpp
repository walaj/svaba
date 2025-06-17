#include "svabaThreadUnit.h"
#include "svabaOutputWriter.h"

#if defined(__GLIBC__)
  #include <malloc.h>
#endif

svabaThreadUnit::svabaThreadUnit(SvabaSharedConfig& sc_) : sc(sc_) {
  
  // load the samtools faidx
  ref_genome = std::make_unique<SeqLib::RefGenome>();
  ref_genome->LoadIndex(sc.opts.refGenome);
  
  // set the *shared* memory BWAIndex
  bwa_aligner = std::make_shared<SeqLib::BWAAligner>(sc.bwa_idx);


  
  // load the .bai files for each bam and set parameters
  for(const auto& [pref, path] : sc.opts.bams)  {
    
    // instantiate a new BAM walker
    auto [it, inserted] = walkers.try_emplace(pref, std::make_shared<svabaBamWalker>(sc));
    it->second->Open(path);
    it->second->SetPrefix(pref);
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
  
#if defined(__GLIBC__)
  malloc_trim(0);
#endif
}

void svabaThreadUnit::flush() {
  
  if (MemoryLimit(THREAD_READ_LIMIT, THREAD_CONTIG_LIMIT)) {
    sc.writer.writeUnit(*this, sc); // mutex and flush are inside this call
    clear();
  }
  
  // clear out the reads and reset the walkers
  for (auto& [_, walker] : walkers) {
    walker->clear();
  }
#if defined(__GLIBC__)
  malloc_trim(0);
#endif
}

bool svabaThreadUnit::MemoryLimit(size_t readLimit, size_t contLimit) const {
  
  size_t stored_reads = all_weird_reads.size() +
    all_discordant_reads.size() +
    all_corrected_reads.size();

  bool mem_exceeded = stored_reads > readLimit ||
                      master_contigs.size()   > contLimit;
  
  if (mem_exceeded) {
    sc.logger.log(sc.opts.verbose > 1, false, "...writing files on thread ",
		  threadId, " with limit hit of ",
		  stored_reads, " reads and ", master_contigs.size(),
		  " contigs=========================================================");
    
    size_t all_walker_reads = 0;
    for (const auto& [_, walker] : walkers) {
      all_walker_reads += walker->reads.size();
    }    
    
    // more memory mapping
    /*    sc.logger.log(true, true, " bad regions: ", badd.size(),
	  " master_alc.size() ", master_alc.size(),
	  " m_disc.size() ", m_disc.size(),
	  " walkers.reads.size() ", all_walker_reads);
    */
  } 
  
  return mem_exceeded;
}

