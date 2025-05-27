// svabaOutputWriter.cpp

#include "svabaOutputWriter.h"
#include "svabaThreadUnit.h"
#include "svabaUtils.h"
#include "BreakPoint.h"
#include "DiscordantCluster.h"

#include <mutex>

using namespace std;

namespace {
  
  // one mutex to serialize ALL threads writes
  static mutex writeMutex_;

}

SvabaOutputWriter::SvabaOutputWriter(SvabaLogger& logger, const SvabaOptions& opts)
  : logger_(logger), opts_(opts)
{ }

// init(): call once, from run_svaba.cpp immediately after you've parsed
// your opts and before you shoot off any threads.
void SvabaOutputWriter::init(const string& analysis_id,
                             const map<string,string>& bamFiles,
			     const SeqLib::BamHeader& b_header) {
  
  // remember these flags for writeUnit:
  read_tracking_      = false; 

  // make a copy of the header for use by this writer
  bam_header_ = b_header;

  // open our gzipped text outputs
  svabaUtils::fopen(analysis_id + ".alignments.txt.gz", all_align_);
  svabaUtils::fopen(analysis_id + ".bps.txt.gz",        os_allbps_);
  svabaUtils::fopen(analysis_id + ".discordant.txt.gz",  os_discordant_);

  // write the header line for breakpoints file:
  os_allbps_ << BreakPoint::header();
  for (auto &p : bamFiles) 
    os_allbps_ << "\t" << p.first << "_" << p.second;
  os_allbps_ << "\n";

  // write the header line for discordant clusters:
  os_discordant_ << DiscordantCluster::header() << "\n";

  // open the contig bam for writing contig alignments
  string wname = analysis_id + ".contigs.bam";
  b_contig_writer_.SetHeader(bam_header_);
  if (!b_contig_writer_.Open(wname)) {
    std::cerr << "Unable to open " << wname << std::endl;
    assert(false);
  }
  b_contig_writer_.WriteHeader();
  
}

void SvabaOutputWriter::writeUnit(svabaThreadUnit& unit) {
  lock_guard<mutex> guard(writeMutex_); // lock the writers

  // alignment plot lines
  for (const auto& alc : unit.m_alc) {
    if (alc.hasVariant()) 
      all_align_ << alc.print(bam_header_) << "\n";
  }

  // discordant clusters
  // hardcoding "false" for readtracking for simplicity
  for (const auto& kv : unit.m_disc) {
    const auto& dc = kv.second;
    if (dc.valid())
      os_discordant_ << dc.toFileString(bam_header_, read_tracking_) << "\n";
  }

  // write contig alignments to BAM
  for (auto& i : unit.m_contigs) {
    i.RemoveTag("MC");
    b_contig_writer_.WriteRecord(i);
  }

  // breakpoints
  for (auto& bp : unit.m_bps) {
    if ( bp.hasMinimal() 
	 && (bp.confidence != "NOLOCAL" || bp.complex_local) )
      {
	os_allbps_ << bp.toFileString(!read_tracking_) << "\n";
      }
  }
  
  // clear them out so this unit can reaccumulate
  unit.clear();
}
