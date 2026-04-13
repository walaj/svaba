#pragma once

// vcf.h
//
// Read a svaba `bps.txt.gz` file back in, dedupe the breakpoint pairs,
// and write out somatic/germline SV and indel VCFs.
//
// Design notes:
// * `BreakPoint` is non-copyable and non-movable (see BreakPoint.h), so we
//   store every breakpoint as a `std::shared_ptr<BreakPoint>` — both the
//   `VCFEntry` views and the `VCFEntryPair` aggregate share ownership of
//   the same underlying object.
// * `BreakPoint::SampleInfo` is a private nested type. vcf.cpp therefore
//   iterates `bp->allele` with `auto` and never names the type directly.
// * The bps.txt.gz file is read one line at a time (memory-efficient); only
//   parsed records that survive the per-contig secondary cap are retained
//   in memory.
// * The external API exposed to `run_svaba.cpp` is preserved verbatim:
//   `VCFFile(file, id, sc, header, nopass, verbose)` +
//   `writeIndels / writeSVs(basename, zip, onefile, bam_header)` +
//   `include_nonpass`.

#include <cstdint>
#include <memory>
#include <ostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "SeqLib/BamHeader.h"
#include "SeqLib/GenomicRegion.h"

#include "BreakPoint.h"
#include "SvabaSharedConfig.h"

// ---------------------------------------------------------------------------
// Header-field dictionaries.
// ---------------------------------------------------------------------------
using InfoMap        = std::unordered_map<std::string, std::string>;
using FilterMap      = std::unordered_map<std::string, std::string>;
using FormatMap      = std::unordered_map<std::string, std::string>;
using SampleMap      = std::unordered_map<std::string, std::string>;
using ContigFieldMap = std::unordered_map<std::string, std::string>;

// ---------------------------------------------------------------------------
// VCFHeader — the ##-lines that precede a VCF body.
// ---------------------------------------------------------------------------
struct VCFHeader {

  VCFHeader()  = default;
  ~VCFHeader() = default;

  std::string fileformat = "VCFv4.2";
  std::string filedate;
  std::string source;
  std::string reference = "hg19";
  std::string colnames  = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

  InfoMap        infomap;
  FilterMap      filtermap;
  FormatMap      formatmap;
  SampleMap      samplemap;
  ContigFieldMap contigfieldmap;

  friend std::ostream& operator<<(std::ostream& out, const VCFHeader& v);

  void addInfoField  (std::string field, std::string number,
                      std::string type,  std::string description);
  void addFilterField(std::string field, std::string description);
  void addFormatField(std::string field, std::string number,
                      std::string type,  std::string description);
  void addSampleField(std::string field);
  void addContigField(std::string id, int len);
};

// ---------------------------------------------------------------------------
// VCFEntry — a single breakend row (one side of a BreakPoint, or an indel).
// ---------------------------------------------------------------------------
struct VCFEntry {

  VCFEntry()  = default;
  ~VCFEntry() = default;

  // Shared ownership with the owning VCFEntryPair.
  std::shared_ptr<BreakPoint> bp;

  // Unique id across all VCFEntryPairs in a VCFFile. `id_num` is 1 for
  // the left breakend (b1) and 2 for the right (b2); indels always use 1.
  uint32_t id      = 0;
  uint8_t  id_num  = 0;

  std::string getRefString() const;
  std::string getAltString(const SeqLib::BamHeader& header) const;
  std::string getIdString () const;

  // Build the full INFO map (keys → values) for this entry.
  std::unordered_map<std::string, std::string> fillInfoFields() const;

  // Serialize the entry as one VCF body line (no trailing newline).
  std::string toFileString(const SeqLib::BamHeader& header) const;

  // Sort primarily by genomic position of the represented breakend.
  bool operator< (const VCFEntry& v) const;
  bool operator==(const VCFEntry& v) const;
};

using VCFEntryVec = std::vector<VCFEntry>;

// ---------------------------------------------------------------------------
// VCFEntryPair — holds both sides (e1, e2) of a rearrangement, plus the
// shared BreakPoint they are built from. Indels still use a pair but only
// populate e1 as a real record.
// ---------------------------------------------------------------------------
struct VCFEntryPair {

  VCFEntryPair() = default;
  ~VCFEntryPair() = default;

  // Build both sides from a shared BreakPoint. Assigns a collision-avoided
  // hash id to e1.id / e2.id.
  explicit VCFEntryPair(std::shared_ptr<BreakPoint>& b);

  VCFEntry                    e1, e2;
  std::shared_ptr<BreakPoint> bp;

  std::string toFileString(const SeqLib::BamHeader& header) const;
};

using VCFEntryPairMap = std::unordered_map<int, std::shared_ptr<VCFEntryPair>>;

// ---------------------------------------------------------------------------
// VCFFile — the top-level object, owning all entry pairs and the two headers.
// ---------------------------------------------------------------------------
struct VCFFile {

  VCFFile()  = default;
  ~VCFFile() = default;

  // Read a svaba bps.txt.gz file, parse line-by-line into BreakPoints,
  // build VCFEntryPairs, and deduplicate. `nopass` controls whether
  // initially non-PASS records are retained (they can still be filtered
  // later by writeIndels/writeSVs via `include_nonpass`).
  VCFFile(std::string              file,
          std::string              id,
          const SvabaSharedConfig& sc,
          const VCFHeader&         vheader,
          bool                     nopass,
          bool                     verbose);

  std::string filename;
  std::string analysis_id;

  bool verbose         = false;
  bool include_nonpass = false;

  VCFHeader indel_header;
  VCFHeader sv_header;

  VCFEntryPairMap entry_pairs; // SVs
  VCFEntryPairMap indels;      // indels
  std::unordered_set<int> dups; // keys in entry_pairs marked as duplicates

  // Dedupe SVs via paired breakend interval-tree matching; dedupe indels
  // via (chr,pos,ref,alt) hash set.
  void deduplicate(const SeqLib::BamHeader& header);

  // Write the body VCFs to disk. `zip` is currently ignored (always writes
  // plain .vcf); `onefile` merges somatic/germline into one emission.
  void writeIndels(std::string              basename,
                   bool                     zip,
                   bool                     onefile,
                   const SeqLib::BamHeader& header) const;

  void writeSVs(std::string              basename,
                bool                     zip,
                bool                     onefile,
                const SeqLib::BamHeader& header) const;
};
