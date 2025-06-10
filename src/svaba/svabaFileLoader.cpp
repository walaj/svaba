// SvabaFileLoader.cpp
#include <filesystem> 
#include <memory>                           // for std::unique_ptr
#include <string>                           // for std::string
#include <vector>                           // if you use any std::vector in this file
#include "svabaFileLoader.h"                // your own header

#include "SvabaSharedConfig.h"
#include "svabaLogger.h"
#include "svabaOptions.h"
#include "SeqLib/RefGenome.h"               // defines SeqLib::RefGenome
#include "SeqLib/GenomicRegionCollection.h" // defines SeqLib::GRC
#include "SeqLib/BamHeader.h"               // if you use SeqLib::BamHeader anywhere
#include "SeqLib/BamReader.h"               // if you invoke any BamReader methods

namespace fs = std::filesystem;
using std::string;

SvabaFileLoader::SvabaFileLoader(SvabaSharedConfig& sc_) : sc(sc_) {}

std::unique_ptr<SeqLib::RefGenome> SvabaFileLoader::loadReference() {
  auto rg = std::make_unique<SeqLib::RefGenome>();
  try {
    rg->LoadIndex(sc.opts.refGenome);
  } catch (const std::exception& e) {
    sc.logger.log(true, true,
             "ERROR: Unable to load reference genome '", sc.opts.refGenome,
             "': ", e.what());
    throw;
  }
  sc.logger.log(false, true, "Loaded reference genome: ", sc.opts.refGenome);
  return rg;
}

void SvabaFileLoader::loadBedRegions(const std::string& path,
				     SeqLib::GRC& gr) {
  if (path.empty()) return;
  
  try {
    gr  = SeqLib::GRC(path, sc.header);
    gr.CreateTreeMap();    
    sc.logger.log(true, true,
             "Loaded BED file '", path, "' (", gr.size(), " regions)");
  } catch (const std::exception& e) {
    sc.logger.log(true, true,
             "ERROR: cannot read BED file '", path, "': ", e.what());
    throw;
  }
}

void SvabaFileLoader::countJobs(
    SeqLib::GRC& runRegions) {

  SeqLib::GRC fileRegions;
  const auto& rf = sc.opts.regionFile;
  // 1) If its a real BED/region file on disk
  if (!rf.empty() && fs::exists(rf)) {
    try {
      fileRegions = SeqLib::GRC(rf, sc.header);
    } catch (const std::exception& ex) {
      sc.logger.log(
		  /*toErr=*/true, /*toLog=*/true,
		  "Could not parse region file '", rf, "': ", ex.what(),
		  ".  Skipping region file.");
    }
    
    // 2) If it looks like chr:start-end
  } else if (rf.find(':') != string::npos && rf.find('-') != string::npos) {
    fileRegions.add(SeqLib::GenomicRegion(rf, sc.header));
    
    // 3) If its just a single chromosome name
  } else if (!rf.empty()) {
    // construct region from chr=rf, pos1=1..end
    SeqLib::GenomicRegion gr(rf, "1", "1", sc.header);
    if (gr.chr < 0 || gr.chr >= sc.header.NumSequences()) {
      throw std::runtime_error("Region file '" + rf +
			       "' failed to match any chromosome in BAM header");
    }
    gr.pos2 = sc.header.GetSequenceLength(gr.chr);
    fileRegions.add(gr);
    
    // 4) else: whole genome
  } else {
    for (int i = 0; i < sc.header.NumSequences(); ++i) {
      fileRegions.add(
		      SeqLib::GenomicRegion(i, 1, sc.header.GetSequenceLength(i))
		      );
    }
  }
  
  if (fileRegions.size() == 0) {
    throw std::runtime_error("No regions found"
			     + (rf.empty() ? string("") : " in '" + rf + "'"));
  }
  
  // 5) Chunk them up
  if (sc.opts.chunkSize > 0) {
    for (auto& r : fileRegions) {
      // this SeqLib constructor splits a region into windows of size chunkSize with pad
      SeqLib::GRC pieces(sc.opts.chunkSize, sc.opts.windowpad, r);
      runRegions.Concat(pieces);
    }
  } else {
    // if chunkSize <= 0, we just do one big pass:
    runRegions = fileRegions;
  }
  
  // if no explicit regionFile, clear fileRegions to signal whole genome
  if (rf.empty()) {
    fileRegions.clear();
  }
    
  return;
}
