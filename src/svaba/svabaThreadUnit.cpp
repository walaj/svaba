#include "svabaThreadUnit.h"
#include "svabaOutputWriter.h"

#if defined(__GLIBC__)
  #include <malloc.h>
#endif

// Build a new BamHeader that prepends an @HD line (with the given sort order)
// to the @SQ block of `src`. Avoids SeqLib::BamHeader::AsString(), which
// blindly dereferences h->text; in modern htslib that field is lazy and is
// often NULL until sam_hdr_str() materializes it, causing a segfault.
static SeqLib::BamHeader make_header_with_hd(const SeqLib::BamHeader& src,
                                             const std::string& so) {
  std::string txt = "@HD\tVN:1.6\tSO:" + so + "\n";

  const bam_hdr_t* h = src.get();
  if (h) {
    // sam_hdr_str forces htslib to regenerate the textual form from its
    // internal representation and caches it in h->text.
    const char* cur = sam_hdr_str(const_cast<bam_hdr_t*>(h));
    if (cur) txt += cur;
  }

  return SeqLib::BamHeader(txt);
}

svabaThreadUnit::svabaThreadUnit(SvabaSharedConfig& sc_,
				 int thread) : sc(sc_), threadId(thread) {
  
  // load the samtools faidx
  ref_genome = std::make_unique<SeqLib::RefGenome>();
  ref_genome->LoadIndex(sc.opts.refGenome);
  
  // set the *shared* memory BWAIndex
  bwa_aligner = std::make_shared<SeqLib::BWAAligner>(sc.bwa_idx);
  //bwa_aligner->allocBuffer(4096); // longest contig

  // load the .bai files for each bam and set parameters
  for(const auto& [pref, path] : sc.opts.bams)  {
    
    // instantiate a new BAM walker
    auto [it, inserted] = walkers.try_emplace(pref, std::make_shared<svabaBamWalker>(sc));
    it->second->Open(path);
    it->second->SetPrefix(pref);
  }

  // make the header in case dumping reads for debugging
  auto hdr_unsorted = make_header_with_hd(sc.header, "unsorted");
  
  // setup the weird read writers
  if (sc.opts.dump_weird_reads) {
    
    // instantiate a new BAM writer
    std::string bamname =sc.opts.analysisId +
      ".thread" + std::to_string(threadId) + ".weird.bam";
    
    auto [it, inserted] = writers.try_emplace("w", std::make_shared<SeqLib::BamWriter>(SeqLib::BAM));
    it->second->SetHeader(hdr_unsorted);
    if (!it->second->Open(bamname)) {
      std::cerr << "ERROR: could not open weird read writer for thread " <<
	threadId << " at path " << bamname << "/n";
      exit(EXIT_FAILURE);
    }
    it->second->WriteHeader();
  }
  
  // setup the corrected read writers
  if (sc.opts.dump_corrected_reads) {
    
    // instantiate a new BAM writer
    std::string bamname =sc.opts.analysisId +
      ".thread" + std::to_string(threadId) + ".corrected.bam";
    
    auto [it, inserted] = writers.try_emplace("c", std::make_shared<SeqLib::BamWriter>(SeqLib::BAM));
    it->second->SetHeader(hdr_unsorted);
    if (!it->second->Open(bamname)) {
      std::cerr << "ERROR: could not open corrected read writer for thread " <<
	threadId << " at path " << bamname << "/n";
      exit(EXIT_FAILURE);
    }
    it->second->WriteHeader();    
  }

  // setup the discordant read writers
  if (sc.opts.dump_discordant_reads) {
    
    // instantiate a new BAM writer
    std::string bamname =sc.opts.analysisId +
      ".thread" + std::to_string(threadId) + ".discordant.bam";
    
    auto [it, inserted] = writers.try_emplace("d", std::make_shared<SeqLib::BamWriter>(SeqLib::BAM));
    it->second->SetHeader(hdr_unsorted);
    if (!it->second->Open(bamname)) {
      std::cerr << "ERROR: could not open discordant read writer for thread " <<
	threadId << " at path " << bamname << "/n";
      exit(EXIT_FAILURE);
    }
    it->second->WriteHeader();
  }

  
}

void svabaThreadUnit::clear() {
  
  master_alc.clear();
  master_contigs.clear();
  
  m_bps.clear();
  m_disc.clear();
  //m_bamreads_count = 0;
  // m_disc_reads     = 0;
  //badd.clear();

  all_weird_reads.clear();
  all_corrected_reads.clear();
  
  ss.str("");       // Clear the string content
  ss.clear();
  
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

svabaThreadUnit::~svabaThreadUnit() {

  //
  if (sc.opts.dump_weird_reads) {
    auto it = writers.find("w");
    it->second->Close();
  }

  //
  if (sc.opts.dump_corrected_reads) {
    auto it = writers.find("c");
    it->second->Close();
  }

}
