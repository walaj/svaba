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

  // write the header line for breakpoints file:
  os_allbps_ << BreakPoint::header();
  for (auto &p : opts.bams) 
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

  // BAM output for aligned contigs
  std::string aligned_contigs_bam_path = analysis_id;
  aligned_contigs_bam_path.append(".discordant.bam");
  b_contig_writer_ = SeqLib::BamWriter(SeqLib::BAM);
  b_contig_writer_.SetHeader(b_header);
  if (!b_contig_writer_.Open(aligned_contigs_bam_path)) {    
    std::cerr << "ERROR: could not open aligned contig writer " <<
      aligned_contigs_bam_path << std::endl;
    exit(EXIT_FAILURE);
  }
  b_contig_writer_.WriteHeader();
  
  // BAM output for weird reads
  if (opts.dump_weird_reads) {
    std::string weird_read_bam_path = analysis_id;
    weird_read_bam_path.append(".weird.bam");
    b_weird_read_writer_ = SeqLib::BamWriter(SeqLib::BAM);
    b_weird_read_writer_.SetHeader(b_header);
    if (!b_weird_read_writer_.Open(weird_read_bam_path)) {
      std::cerr << "ERROR: could not open output weird read writer " <<
	weird_read_bam_path << std::endl;
      exit(EXIT_FAILURE);
    }
    b_weird_read_writer_.WriteHeader();
  }
  
  // BAM output for discordant reads
  if (opts.dump_discordant_reads) {
    std::string discordant_read_bam_path = analysis_id;
    discordant_read_bam_path.append(".discordant.bam");
    b_discordant_read_writer_ = SeqLib::BamWriter(SeqLib::BAM);
    b_discordant_read_writer_.SetHeader(b_header);
    if (!b_discordant_read_writer_.Open(discordant_read_bam_path)) {    
      std::cerr << "ERROR: could not open output discordant read writer " <<
	discordant_read_bam_path << std::endl;
      exit(EXIT_FAILURE);
    }
    b_discordant_read_writer_.WriteHeader();
  }

  // BAM output for discordant reads
  if (opts.dump_corrected_reads) {
    std::string corrected_read_bam_path = analysis_id;
    corrected_read_bam_path.append(".corrected.bam");
    b_corrected_read_writer_ = SeqLib::BamWriter(SeqLib::BAM);
    b_corrected_read_writer_.SetHeader(b_header);
    if (!b_corrected_read_writer_.Open(corrected_read_bam_path)) {    
      std::cerr << "ERROR: could not open output corrected read writer " <<
	corrected_read_bam_path << std::endl;
      exit(EXIT_FAILURE);
    }
    b_corrected_read_writer_.WriteHeader();
  }

  
}

void SvabaOutputWriter::writeUnit(svabaThreadUnit& unit) {
  lock_guard<mutex> guard(writeMutex_); // lock the writers

  // alignment plot lines
  for (const auto& alc : unit.master_alc) {
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
  for (auto& i : unit.master_contigs) {
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

  // discordant reads
  if (opts.dump_discordant_reads) {
    for (const auto& r : unit.all_discordant_reads)
      b_discordant_read_writer_.WriteRecord(r);
  }

  // weird reads
  if (opts.dump_weird_reads) {
    for (const auto& r : unit.all_weird_reads)
      b_weird_read_writer_.WriteRecord(r);
  }

  // corrected reads
  if (opts.dump_corrected_reads) {
    for (const auto& r : unit.all_corrected_reads)
      b_corrected_read_writer_.WriteRecord(r);
  }

  // aligned contigs
  for (const auto& r : unit.master_contigs)
    b_contig_writer_.WriteRecord(r);
  
  // clear them out so this unit can reaccumulate
  unit.clear();
}

void SvabaOutputWriter::close() {

  all_align_.close();
  os_allbps_.close();
  os_discordant_.close();
  
  if (opts.dump_discordant_reads) {
    if (!b_discordant_read_writer_.Close()) {
      std::cerr << "Unable to close discordant read writer" << std::endl;
    }
  }
  
  if (opts.dump_weird_reads) {
    if (!b_weird_read_writer_.Close()) {
      std::cerr << "Unable to close weird read writer" << std::endl;
    }
  }

  if (opts.dump_corrected_reads) {
    if (!b_corrected_read_writer_.Close()) {
      std::cerr << "Unable to close corrected read writer" << std::endl;
    }
  }

  b_contig_writer_.Close();
  
}
