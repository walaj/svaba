#pragma once

#include <cstddef>                            // for size_t
#include "SeqLib/GenomicRegion.h"             // for SeqLib::GenomicRegion
#include "svabaThreadUnit.h"                  // for svabaThreadUnit
#include "SvabaLogger.h"                      // your logging class
#include "SvabaOptions.h"                     // your options struct/class
#include "svabaOutputWriter.h"                // your output writer

#include "SeqLib/BamHeader.h"

/// Encapsulates everything needed to process one genomic region:
///  - shared logger
///  - shared run-time options
///  - shared output writer
/// Its `process()` method is exactly where your old runWorkItem logic goes.
class SvabaRegionProcessor {
public:
  SvabaRegionProcessor(SvabaLogger&        logger,
                       const SvabaOptions& opts,
		       const SeqLib::BamHeader& header,
                       SvabaOutputWriter&  writer);

  /// Run the SVABA pipeline on `region` using per-thread state in `unit`.
  /// Returns true on success, false on error.
  bool process(const SeqLib::GenomicRegion& region,
               svabaThreadUnit&             unit,
               size_t                       threadId);

private:
  SvabaLogger&        logger_;
  const SvabaOptions& opts_;
  const SeqLib::BamHeader& header_;
  SvabaOutputWriter&  writer_;
};
