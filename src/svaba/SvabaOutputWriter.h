// svabaOutputWriter.h
#pragma once

#include <map>
#include <string>
#include <fstream>
#include "SvabaThreadUnit.h"
#include "SeqLib/BamWriter.h"
#include "gzstream.h"
#include "SvabaLogger.h"
#include "SvabaOptions.h"


using std::ofstream;
struct SvabaSharedConfig;

struct SvabaOutputWriter {

public:
  
  SvabaOutputWriter(SvabaLogger& logger_, SvabaOptions& opts_);
  
  // call this once at startup:
  void init(const std::string& analysis_id,
	    const SeqLib::BamHeader& b_header);

  // call this from any thread under a lock:
  void writeUnit(svabaThreadUnit& unit,
		 SvabaSharedConfig& sc);

  void close();
  
private:

  SvabaLogger&        logger;
  SvabaOptions&       opts;
  
  ogzstream           all_align_;
  // SvABA2.0: structured (TSV) counterpart of all_align_. Same info, but
  // serialized as one row per contig / one row per r2c-aligned read so
  // downstream tools can re-plot on demand. See AlignedContig::printToR2CTsv.
  ogzstream           os_r2c_;
  ogzstream           os_allbps_;
  ogzstream           os_discordant_;
  ofstream            os_runtime_;

  SeqLib::BamWriter   b_contig_writer_;
  
  SeqLib::BamWriter   b_weird_read_writer_;
  SeqLib::BamWriter   b_discordant_read_writer_;
  SeqLib::BamWriter   b_corrected_read_writer_;  
  
  SeqLib::BamHeader bam_header_;
  bool                disc_cluster_only_ = false;
};

