// svabaOutputWriter.h
#pragma once

#include <map>
#include <string>
#include "svabaThreadUnit.h"
#include "SeqLib/BamWriter.h"
#include "gzstream.h"
#include "svabaLogger.h"
#include "svabaOptions.h"

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
  ogzstream           os_allbps_;
  ogzstream           os_discordant_;

  SeqLib::BamWriter   b_contig_writer_;
  
  SeqLib::BamWriter   b_weird_read_writer_;
  SeqLib::BamWriter   b_discordant_read_writer_;
  SeqLib::BamWriter   b_corrected_read_writer_;  
  
  SeqLib::BamHeader bam_header_;
  bool                read_tracking_ = false;
  bool                disc_cluster_only_ = false;
};

