#pragma once

#include "svabaLogger.h"
#include "svabaOptions.h"
#include "SeqLib/BWAIndex.h"
#include "svabaOutputWriter.h"
#include "SeqLib/BamHeader.h"
#include "SeqLib/GenomicRegionCollection.h"

struct SvabaSharedConfig {

  SvabaLogger&          logger;

  SvabaOptions&         opts;
  
  SeqLib::BWAIndexPtr   bwa_idx;

  SvabaOutputWriter     writer; 

  SeqLib::BamHeader     header;

  SeqLib::GRC           blacklist;
  
}
