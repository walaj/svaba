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
  
  // SvABA2.0: alignments.txt.gz (pre-rendered ASCII) and the shared
  // r2c.txt.gz stream used to live here. Both are gone:
  //   - alignments.txt.gz is removed entirely; its re-plot-able successor
  //     is r2c.txt.gz (AlignedContig::printToR2CTsv + r2cTsvHeader).
  //   - r2c.txt.gz is now emitted per-thread by svabaThreadUnit::r2c_out_
  //     and merged at postprocess time (gzip is concatenation-safe).
  // The remaining shared streams below are small enough that mutex
  // contention isn't measurable (bps + cluster-level discordant +
  // runtime are KB-MB per run).
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

