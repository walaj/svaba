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
// Knobs for the `svaba tovcf` subcommand.
// ---------------------------------------------------------------------------
//
// QualMode — what to put in the VCF QUAL column.
//   SUM_LO_PHRED (legacy): -10 log10 P(all alt reads are errors), taken from
//                          BreakPoint::quality (col 34 of bps.txt). Keeps
//                          behavior compatible with `svaba run` / refilter.
//   MAXLOD_PHRED:          round(10 * bp->max_lod), cap 99. One-to-one with
//                          INFO/MAXLOD, meaningful as a VCF QUAL field.
//   MISSING:               emit '.' (VCF spec-valid missing value). Forces
//                          users to look at FILTER / INFO for confidence,
//                          which is usually what you want when somlod and
//                          maxlod are the canonical scores.
enum class QualMode { SUM_LO_PHRED, MAXLOD_PHRED, MISSING };

// SvFormat — how to serialize an SV call.
//   BND_ALWAYS (legacy):          every SV emits as two paired BND records
//                                 with mate-bracket ALT notation. Matches
//                                 the current `svaba run` output exactly.
//   SYMBOLIC_WHEN_OBVIOUS:        intrachromosomal events with unambiguous
//                                 orientation emit as a single symbolic
//                                 record (<DEL> / <DUP> / <INV>). Inter-
//                                 chromosomal and complex events stay BND.
//                                 Nicer for IGV / samplot / bcftools
//                                 pipelines that render symbolic alleles
//                                 natively.
enum class SvFormat { BND_ALWAYS, SYMBOLIC_WHEN_OBVIOUS };

// ---------------------------------------------------------------------------
// VCFHeader — the ##-lines that precede a VCF body.
// ---------------------------------------------------------------------------
struct VCFHeader {

  VCFHeader()  = default;
  ~VCFHeader() = default;

  // VCF spec version this header declares. Defaulted to 4.5 (the most
  // recent formal spec as of 2024) — backwards-compatible at the record
  // level with anything accepting 4.2+, and enables the clean symbolic-
  // SV / EVENT / SVCLAIM idioms the tovcf path uses.
  std::string fileformat = "VCFv4.5";
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

  // Set by VCFFile::writeSvsSingleFile when the SvFormat::
  // SYMBOLIC_WHEN_OBVIOUS mode decides this pair should be emitted as
  // a single `<DEL>`/`<DUP>`/`<INV>` record instead of two paired BND
  // records. When true, only the e1 entry of the pair gets serialized
  // (the caller skips e2 for that pair). Kind is stashed here so
  // getAltString / fillInfoFields can key off it without reclassifying.
  bool        symbolic_rep  = false;
  std::string symbolic_kind;  // "DEL" / "DUP" / "INV" when symbolic_rep is true

  // Somlod cutoff for stamping the INFO/SOMATIC flag. A record gets the
  // flag iff bp->LO_s >= somatic_threshold AND bp->somatic != FAILED.
  // Writers set this per-entry from VCFFile::somatic_threshold before
  // calling toFileString(); default 1.0 matches the tovcf default.
  double      somatic_threshold = 1.0;

  std::string getRefString() const;
  std::string getAltString(const SeqLib::BamHeader& header) const;
  std::string getIdString () const;

  // Build the full INFO map (keys → values) for this entry.
  // `somatic_flag` = true adds `SOMATIC` (a flag INFO field).
  std::unordered_map<std::string, std::string> fillInfoFields() const;

  // Serialize the entry as one VCF body line (no trailing newline).
  //   qual_mode   — what goes in the QUAL column (see QualMode).
  //   NB: the symbolic/BND choice is read off this->symbolic_rep,
  //       which the writer sets before calling.
  std::string toFileString(const SeqLib::BamHeader& header,
                           QualMode qual_mode = QualMode::SUM_LO_PHRED) const;

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

  std::string toFileString(const SeqLib::BamHeader& header,
                           QualMode qual_mode = QualMode::SUM_LO_PHRED) const;
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
  //
  // `skip_dedup_in_ctor` = true has the ctor set `this->skip_dedup`
  // BEFORE the internal `deduplicate()` call fires, so no dedup runs
  // on input that was already deduped upstream (the `svaba tovcf` path).
  // Default false keeps backward-compat with run_svaba and refilter.
  VCFFile(std::string              file,
          std::string              id,
          const SvabaSharedConfig& sc,
          const VCFHeader&         vheader,
          bool                     nopass,
          bool                     verbose,
          bool                     skip_dedup_in_ctor = false);

  std::string filename;
  std::string analysis_id;

  bool verbose         = false;
  bool include_nonpass = false;

  // SvABA2.0 tovcf-path knobs. All default to legacy-compatible values so
  // run_svaba.cpp and refilter.cpp pick up no behavior change. The `svaba
  // tovcf` entry point flips them to the new defaults (MISSING qual,
  // SYMBOLIC_WHEN_OBVIOUS sv format, skip_dedup=true).
  bool     skip_dedup = false;  // deduplicate() returns immediately when set
  QualMode qual_mode  = QualMode::SUM_LO_PHRED;
  SvFormat sv_format  = SvFormat::BND_ALWAYS;

  // Minimum somlod (LO_s) to mark a record with the INFO/SOMATIC flag.
  // Writers stamp this onto each VCFEntry before emission. Default 1.0
  // is the recommended somatic-call gate; lower it to be permissive,
  // raise it to be stricter.
  double   somatic_threshold = 1.0;

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

  // SvABA2.0: single-file writers used by `svaba tovcf`. Unlike
  // writeSVs / writeIndels above, these write EVERY record (somatic and
  // germline together) to one path, with the SOMATIC INFO flag set on
  // somatic rows for downstream filtering. Honors sv_format, qual_mode,
  // and include_nonpass set on the VCFFile. `gzip=true` uses ogzstream;
  // otherwise a plain ofstream is used.
  void writeSvsSingleFile   (const std::string&       path,
                             bool                     gzip,
                             const SeqLib::BamHeader& header) const;
  void writeIndelsSingleFile(const std::string&       path,
                             bool                     gzip,
                             const SeqLib::BamHeader& header) const;
};
