// svabaOutputWriter.h
#pragma once

#include <map>
#include <string>
#include "svabaThreadUnit.h"
#include "SeqLib/BamWriter.h"
#include "gzstream.h"

struct SvabaOutputWriter {
  
  // call this once at startup:
  void init(const std::string& analysis_id,
            const std::map<std::string,std::string>& bamFiles,
	    const SeqLib::BamHeader& b_header);

  // call this from any thread under a lock:
  void writeUnit(svabaThreadUnit& unit);
  
private:
  ogzstream           all_align_;
  ogzstream           os_allbps_;
  ogzstream           os_discordant_;

  SeqLib::BamWriter   b_contig_writer_;

  SeqLib::BamHeader bam_header_;
  bool                read_tracking_ = false;
  bool                disc_cluster_only_ = false;
};

