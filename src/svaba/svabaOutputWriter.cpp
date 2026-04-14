// svabaOutputWriter.cpp

#include "SvabaOutputWriter.h"
#include "SvabaThreadUnit.h"
#include "SvabaUtils.h"
#include "BreakPoint.h"
#include "DiscordantCluster.h"
#include "SvabaSharedConfig.h"

#include <htslib/sam.h>

#include <mutex>

using namespace std;

namespace {

  // one mutex to serialize ALL threads writes
  static mutex writeMutex_;

}

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

SvabaOutputWriter::SvabaOutputWriter(SvabaLogger& logger_, SvabaOptions& opts_)
  : logger(logger_), opts(opts_)
{ }

// init(): call once, from run_svaba.cpp immediately after you've parsed
// your opts and before you shoot off any threads.
void SvabaOutputWriter::init(const string& analysis_id,
			     const SeqLib::BamHeader& b_header) {

  // keep a local reference to the header
  bam_header_ = b_header;

  // open our gzipped text outputs
  svabaUtils::fopen(analysis_id + ".alignments.txt.gz", all_align_);
  svabaUtils::fopen(analysis_id + ".bps.txt.gz",        os_allbps_);
  svabaUtils::fopen(analysis_id + ".discordant.txt.gz",  os_discordant_);
  svabaUtils::fopen(analysis_id + ".runtime.txt",  os_runtime_);  

  // write the header line for breakpoints file:
  os_allbps_ << BreakPoint::header();
  for (auto &p : opts.bams) 
    os_allbps_ << "\t" << p.first << "_" << p.second;
  os_allbps_ << "\n";

  // write the header line for runtime file
  os_runtime_ <<
    svabaUtils::svabaTimer::header << "\n";
  
  // write the header line for discordant clusters:
  os_discordant_ << DiscordantCluster::header() << "\n";

  // BAM output for aligned contigs
  std::string aligned_contigs_bam_path = analysis_id;
  aligned_contigs_bam_path.append(".contigs.bam");
  b_contig_writer_ = SeqLib::BamWriter(SeqLib::BAM);
  auto hdr_unsorted = make_header_with_hd(b_header, "unsorted");
  b_contig_writer_.SetHeader(hdr_unsorted);  
  if (!b_contig_writer_.Open(aligned_contigs_bam_path)) {    
    std::cerr << "ERROR: could not open aligned contig writer " <<
      aligned_contigs_bam_path << std::endl;
    exit(EXIT_FAILURE);
  }
  b_contig_writer_.WriteHeader();
  
  
}

void SvabaOutputWriter::writeUnit(svabaThreadUnit& unit,
				  SvabaSharedConfig& sc) {

  // write the weird reads
  if (sc.opts.dump_weird_reads) {
    
    auto it = unit.writers.find("w");
    if (it == unit.writers.end()) {
      std::cerr << " BAM writer for weird reads not found" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    for (const auto& r : unit.all_weird_reads) {

      // switch out the qname so it indicates also which BAM it was from
      r->SetQname(r->UniqueName());

      r->RemoveTag("RG");
      r->RemoveTag("PG");
      bool ok = it->second->WriteRecord(*r);
      if (!ok)
	std::cerr << "...unable to write weird read record" << std::endl;
    }
  }

  // write the discordant reads
  if (sc.opts.dump_discordant_reads) {
    auto it = unit.writers.find("d");
    if (it == unit.writers.end()) {
      std::cerr << " BAM writer for discordant reads not found" << std::endl;
    }
    
    for (const auto& r : unit.all_weird_reads) {
      if (r->dd > 0) {

	// switch out the qname so it indicates also which BAM it was from
	r->SetQname(r->UniqueName());

	// write it
	r->RemoveTag("RG");
	r->RemoveTag("PG");
	bool ok = it->second->WriteRecord(*r);
	if (!ok)
	  std::cerr << "...unable to write discordant read record" << std::endl;
      }
    }

  }
  
  // write the corrected reads
  if (sc.opts.dump_corrected_reads) {
    auto it = unit.writers.find("c");
    if (it == unit.writers.end()) {
      std::cerr << " read found with prefix " << " c " <<
	" that is not in the weird read writers\n";
      exit(EXIT_FAILURE);
    }
    
    for (const auto& r : unit.all_corrected_reads) {
      r->RemoveTag("RG");
      r->RemoveTag("PG");
      bool ok = it->second->WriteRecord(*r);
      if (!ok)
	std::cerr << "...unable to write corrected read record" << std::endl;
    }
  }    
  lock_guard<mutex> guard(writeMutex_); // lock the writers
  
  sc.total_regions_done += unit.processed_since_memory_dump;
  unit.processed_since_memory_dump = 0;

  // alignment plot lines
  for (const auto& alc : unit.master_alc) {
    if (alc.hasVariant()) 
      all_align_ << alc.printToAlignmentsFile(bam_header_) << "\n";
  }

  // discordant clusters
  // hardcoding "false" for readtracking for simplicity
  for (const auto& kv : unit.m_disc) {
    const auto& dc = kv.second;
    if (dc.valid())
      os_discordant_ << dc.toFileString(bam_header_, false) << "\n";
  }

  // write contig alignments to BAM
  for (auto& i : unit.master_contigs) {
    bool ok = b_contig_writer_.WriteRecord(*i);
    assert(ok);
  }

  // breakpoints
  for (auto& bp : unit.m_bps) {
    
    if ( bp->hasMinimal() )
      //	 && (bp.confidence != "NOLOCAL" || bp.complex_local) )
      {
	os_allbps_ << bp->toFileString(sc.header) << "\n";
      }
  }
  
  // discordant reads
  //std::cerr << " DC " << unit.all_discordant_reads.size() << std::endl;    
  // if (opts.dump_discordant_reads) {
  //   for (const auto& r : unit.all_discordant_reads) {
  //     bool ok = b_discordant_read_writer_.WriteRecord(*r);
  //     assert(ok);
  //   }
  // }

  // weird reads
  // if (opts.dump_weird_reads) {
  //   for (const auto& r : unit.all_weird_reads) {
  //     bool ok = b_weird_read_writer_.WriteRecord(*r);
  //     assert(ok);
  //   }
  // }

  // // corrected reads
  // if (opts.dump_corrected_reads) {
  //   for (const auto& r : unit.all_corrected_reads) {
  //     bool ok = b_corrected_read_writer_.WriteRecord(*r);
  //     assert(ok);
  //   }
  // }

  // runtime
  os_runtime_ << unit.ss.str();
}

void SvabaOutputWriter::close() {

  all_align_.close();
  os_allbps_.close();
  os_discordant_.close();
  
  // if (opts.dump_discordant_reads) {
  //   if (!b_discordant_read_writer_.Close()) {
  //     std::cerr << "Unable to close discordant read writer" << std::endl;
  //   }
  // }
  
  // if (opts.dump_weird_reads) {
  //   if (!b_weird_read_writer_.Close()) {
  //     std::cerr << "Unable to close weird read writer" << std::endl;
  //   }
  // }
  
  // if (opts.dump_corrected_reads) {
  //   if (!b_corrected_read_writer_.Close()) {
  //     std::cerr << "Unable to close corrected read writer" << std::endl;
  //   }
  // }
  
  if (!b_contig_writer_.Close()) {
    std::cerr << "Unable to close contigs bam writer" << std::endl;
  }
  
}
