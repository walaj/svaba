// SvabaFileLoader.h
#pragma once

#include <memory>
#include <string>

#include "SeqLib/GenomicRegionCollection.h"

class SvabaLogger;
class SvabaOptions;
class SvabaSharedConfig;
namespace SeqLib {
  class RefGenome;
}

/// A small helper for all the open this file, bail out if it fails logic.
/// Carries around references to your single SvabaLogger and SvabaOptions.
class SvabaFileLoader {
public:
  /// You must give it your one-and-only logger and options.
  SvabaFileLoader(SvabaSharedConfig& sc_);

  /// Loads a BED style file into a GenomicRegionCollection, or returns empty if no path.
  void loadBedRegions(const std::string& path,
		      SeqLib::GRC& gr);

  std::unique_ptr<SeqLib::RefGenome> loadReference();
  
  /// load the regions to run file and parese out jobs
  void countJobs(SeqLib::GRC& runRegions);
  
private:
  SvabaSharedConfig& sc;
};
