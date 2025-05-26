#pragma once

// --- minimal STL includes ---
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <cstddef>   // for size_t

// --- forward declare only what's needed from SeqLib ---
namespace SeqLib {
  class BamReader;
  class BamRecord;
  using BamRecordVector = std::vector<BamRecord>;

  class BWAWrapper;
  class RefGenome;
  class CigarMap;

  class SharedHTSFile;
  
  template<typename T> class GenomicRegionCollection;
  using GRC = GenomicRegionCollection<GenomicRegion>;
}

// --- forward declare your own types ---
struct svabaThreadUnit;
struct svabaBamWalker;
struct AlignedContig;
struct DiscordantCluster;
struct svabaRead;  // if you have a svabaRead type



