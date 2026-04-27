#pragma once

#include <cstddef>                            // for size_t
#include "SvabaThreadUnit.h"                  // for svabaThreadUnit
#include "SvabaUtils.h"
#include "SeqLib/GenomicRegionCollection.h"   // for GRC return type

class SvabaLogger;
class SvabaOptions;
class SvabaOutputWriter;

/// Encapsulates everything needed to process one genomic region:
///  - shared logger
///  - shared run-time options
///  - shared output writer
/// Its `process()` method is exactly where your old runWorkItem logic goes.

namespace SeqLib {
  class BamHeader;
}

using svabaUtils::svabaTimer;

class SvabaRegionProcessor {
public:
  SvabaRegionProcessor(SvabaSharedConfig& sh_cf); 

  /// Run the SVABA pipeline on `region` using per-thread state in `unit`.
  /// Returns true on success, false on error.
  bool process(const SeqLib::GenomicRegion& region,
               svabaThreadUnit&             unit,
               size_t                       threadId);

  /// Run mate-region collection. Returns the merged somatic mate regions
  /// (as a GRC) so the caller can include their reference in a local BWA
  /// index for native realignment.
  SeqLib::GRC runMateCollectionLoop(const SeqLib::GenomicRegion& region,
				    svabaThreadUnit& stu);
    

private:
  SvabaSharedConfig& sc;
};
