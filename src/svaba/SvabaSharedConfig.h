#pragma once

#include <set>
//#include <atomic>

#include "SeqLib/BWAIndex.h"
#include "SeqLib/BamHeader.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/ReadFilter.h"
#include "LearnBamParams.h"
#include <ctime>

class SvabaOutputWriter;
class SvabaLogger;
class SvabaOptions;
class SvabaOutputWriter;
class DBSnpFilter;

class SvabaSharedConfig {
  
 public:

 SvabaSharedConfig(SvabaLogger&        _logger,
		   SvabaOptions&       _opts,
		   SvabaOutputWriter&  _writer)
   : logger(_logger),
    opts(_opts),
    writer(_writer) {}

  // Delete copy constructor and copy assignment
  SvabaSharedConfig(const SvabaSharedConfig&) = delete;
  SvabaSharedConfig& operator=(const SvabaSharedConfig&) = delete;
  
  // Delete move constructor and move assignment
  SvabaSharedConfig(SvabaSharedConfig&&) = delete;
  SvabaSharedConfig& operator=(SvabaSharedConfig&&) = delete;
  
  size_t total_regions_to_process = 0;
  size_t total_regions_done = 0;

  struct timespec start;

  // store which readgroups we already warned about
  std::unordered_set<std::string> warned; 
  
  SvabaLogger&          logger;

  SvabaOptions&         opts;
  
  SvabaOutputWriter&     writer; 
  
  SeqLib::BWAIndexPtr   bwa_idx;
  
  SeqLib::BamHeader     header;

  SeqLib::GRC           blacklist;

  SeqLib::GRC           file_regions;  

  SeqLib::GRC           germline_svs;

  int readlen = -1;
  double insertsize = -1;
  
  std::shared_ptr<DBSnpFilter>           dbsnp_filter;

  SeqLib::Filter::ReadFilterCollection  mr;

  std::unordered_map<std::string, LearnBamParams> bamStats;

  // needed for aligned contig
  std::set<std::string> prefixes;

  std::string args = "svaba";
};
