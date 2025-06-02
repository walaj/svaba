#pragma once

#include <cstddef>                            // for size_t
#include "svabaThreadUnit.h"                  // for svabaThreadUnit

class SvabaLogger;
class SvabaOptions;
class SvabaOutputWriter;

/// Encapsulates everything needed to process one genomic region:
///  - shared logger
///  - shared run-time options
///  - shared output writer
/// Its `process()` method is exactly where your old runWorkItem logic goes.

namespace SeqLib {
  class GenomicRegion;
  class BamHeader;
}

class SvabaRegionProcessor {
public:
  SvabaRegionProcessor(SvabaSharedConfig& sh_cf); 

  /// Run the SVABA pipeline on `region` using per-thread state in `unit`.
  /// Returns true on success, false on error.
  bool process(const SeqLib::GenomicRegion& region,
               svabaThreadUnit&             unit,
               size_t                       threadId);

  void runMateCollectionLoop(const SeqLib::GenomicRegion& region,
			     svabaThreadUnit& stu);
    

private:
  SvabaSharedConfig& sc;
};
