#include "SvabaRegionProcessor.h"

SvabaRegionProcessor::SvabaRegionProcessor(SvabaLogger&        logger,
                                           const SvabaOptions& opts,
					   const SeqLib::BamHeader& header,
                                           SvabaOutputWriter&  writer)
  : logger_(logger)
  , opts_(opts)
  , header_(header), 
  , writer_(writer)
{ }

bool SvabaRegionProcessor::process(const SeqLib::GenomicRegion& region,
                                   svabaThreadUnit&             unit,
                                   size_t                       threadId)
{

  logger_.log(opts_.verbose > 1, false/*TODO*/, "===Running region ", region.ToString(header_),
	      " on thread ", threadID);

  // setup the params
  for (auto& w : stu.walkers)
    setWalkerParams(w.second);
  
  
  // === 1) your per-region setup using opts_ ===
  // e.g. setWalkerParams(unit.walkers[...], opts_.someThreshold);

  // === 2) your main SVABA logic (formerly in runWorkItem) ===
  //   - extraction
  //   - assembly
  //   - realignment
  //   - breakpoint detection
  //   - etc.
  //
  // You can use:
  //    logger_.log(/*toErr*/false, /*toLog*/true, "some message");
  // when you had WRITELOG before, and
  //    writer_.writeUnit(unit);
  // instead of the old ::WriteFilesOut.

  bool success = true;
  try {
    //  paste in your old runWorkItem body here, replacing:
    //    WRITELOG(...)       logger_.log(...)
    //    WriteFilesOut(...)  writer_.writeUnit(unit)
    //
    // e.g.
    // logger_.log(false,true,
    //   "Thread ", threadId,
    //   ": processing region ", region.str());
    //
    // runAssembly(...);
    // collectBreakpoints(...);
    //
    // at end:
    writer_.writeUnit(unit);
  }
  catch (const std::exception& e) {
    logger_.log(true, true,
                "Error in thread ", threadId,
                " on region ", region.str(),
                ": ", e.what());
    success = false;
  }
  return success;
}
