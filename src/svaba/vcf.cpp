// vcf.cpp
//
// Read svaba breakpoints back from `bps.txt.gz`, deduplicate them, and emit
// VCFs. See vcf.h for the external API. Only ~three things change vs the
// pre-overhaul vcf.cpp:
//
//   1. BreakPoint is non-copyable/non-movable, so everything flows through
//      std::shared_ptr<BreakPoint>.
//   2. BreakPoint::SampleInfo is private; we iterate `bp->allele` via
//      structured bindings (auto) and only use its public data members.
//   3. The `BreakPoint(line, &sc)` ctor now takes a SvabaSharedConfig*, and
//      the parsed bps row schema (39 fixed cols + per-sample tail) is
//      described in BreakPoint.cpp at the parse loop.

#include "vcf.h"

#include <algorithm>
#include <cassert>
#include <cmath>       // std::lround (QUAL phred-rescaling)
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "gzstream.h"

#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"

#include "BreakPoint.h"
#include "SvabaOptions.h"

// ---------------------------------------------------------------------------
// File-local constants and helpers.
// ---------------------------------------------------------------------------
namespace {

// Maximum overlaps one breakend can have before we blacklist the whole
// region as an SV-pileup and drop it from the output.
constexpr int HIGH_OVERLAP_LIMIT = 500;

// How many bp of padding to apply when we blacklist a high-overlap region.
constexpr int BAD_REGION_PAD = 200;

// FORMAT column strings — must stay in sync with BreakPoint::SampleInfo::
// toFileString (see BreakPoint.cpp ~line 1710) and with the ##FORMAT lines
// we emit below.
const std::string kSvFormat    = "GT:AD:DP:SR:DR:GQ:PL:LO:LO_n";
const std::string kIndelFormat = "GT:AD:DP:SR:CR:GQ:PL:LO:LO_n";

// INFO fields that were declared `Type=Flag` in the header — flags are
// printed without an `=value` tail.
std::unordered_set<std::string> g_flag_fields;

// Track hashed VCFEntryPair ids we've already handed out so we can retry on
// collision. File-scope because the constructor consumes it across all pairs.
std::unordered_set<uint32_t> g_hash_avoid;

// Stable 32-bit string hash helper; we avoid pulling in htslib/khash.h just
// for this and use std::hash truncated to 32 bits.
inline uint32_t hash_string32(const std::string& s) {
  const size_t h = std::hash<std::string>{}(s);
  return static_cast<uint32_t>(h ^ (h >> 32));
}

// Sort INFO key/value pairs so that READ_ID (if present) always comes last,
// and everything else is sorted alphabetically. This matches legacy output.
bool compareInfoFields(const std::pair<std::string, std::string>& lhs,
                       const std::pair<std::string, std::string>& rhs) {
  const bool lhs_rid = (lhs.first == "READ_ID");
  const bool rhs_rid = (rhs.first == "READ_ID");
  if (rhs_rid && !lhs_rid) return true;
  if (lhs_rid && !rhs_rid) return false;
  return lhs.first < rhs.first;
}

// For the contig-by-length sort in the header emitter.
bool pairCompareDesc(
    const std::pair<int, std::pair<std::string, std::string>>& a,
    const std::pair<int, std::pair<std::string, std::string>>& b) {
  return a.first > b.first;
}

// SVType → short string used in the EVDNC INFO field.
std::string svtype_to_string(SVType t) {
  switch (t) {
    case SVType::NOTSET:     return "NOTSET";
    case SVType::TSI_LOCAL:  return "TSI_LOCAL";
    case SVType::TSI_GLOBAL: return "TSI_GLOBAL";
    case SVType::ASSMB:      return "ASSMB";
    case SVType::ASDIS:      return "ASDIS";
    case SVType::DSCRD:      return "DSCRD";
    case SVType::INDEL:      return "INDEL";
  }
  return "UNKNOWN_SVTYPE";
}

// Return the maximum per-sample LO across all alleles on a BreakPoint.
// Used as the legacy "LOD" INFO field for indels.
double max_allele_lo(const BreakPoint& bp) {
  double m = 0.0;
  for (const auto& [_, al] : bp.allele)
    m = std::max(m, al.LO);
  return m;
}

} // namespace

// ===========================================================================
// VCFHeader
// ===========================================================================

std::ostream& operator<<(std::ostream& out, const VCFHeader& v) {

  out << "##fileformat=" << v.fileformat << '\n';
  out << "##fileDate="   << v.filedate   << '\n';
  out << "##source="     << v.source     << '\n';
  out << "##reference="  << v.reference  << '\n';

  // Emit contigs sorted by length, largest first (matches legacy output).
  using CC = std::pair<std::string, std::string>;
  using CI = std::pair<int, CC>;
  std::vector<CI> contig_vec;
  contig_vec.reserve(v.contigfieldmap.size());
  for (const auto& kv : v.contigfieldmap)
    contig_vec.emplace_back(std::stoi(kv.second), CC(kv.first, kv.second));
  std::sort(contig_vec.begin(), contig_vec.end(), pairCompareDesc);

  for (const auto& c : contig_vec)
    out << "##contig=<ID=" << c.second.first
        << ",length=" << c.second.second << ">\n";

  for (const auto& kv : v.infomap)
    out << "##INFO=<ID="   << kv.first << ',' << kv.second << ">\n";
  for (const auto& kv : v.filtermap)
    out << "##FILTER=<ID=" << kv.first << ',' << kv.second << ">\n";
  for (const auto& kv : v.formatmap)
    out << "##FORMAT=<ID=" << kv.first << ',' << kv.second << ">\n";
  for (const auto& kv : v.samplemap)
    out << "##SAMPLE=<ID=" << kv.first << ">\n";

  out << v.colnames;
  return out;
}

void VCFHeader::addInfoField(std::string field, std::string number,
                             std::string type,  std::string description) {
  if (infomap.find(field) != infomap.end()) {
    std::cerr << "Warning: Info field already exists: " << field << '\n';
    return;
  }
  if (type == "Flag")
    g_flag_fields.insert(field);

  infomap[field] = "Number=" + number +
                   ",Type=" + type +
                   ",Description=\"" + description + "\"";
}

void VCFHeader::addFilterField(std::string field, std::string description) {
  if (filtermap.find(field) != filtermap.end()) {
    std::cerr << "Warning: Filter field already exists: " << field << '\n';
    return;
  }
  filtermap[field] = "Description=\"" + description + "\"";
}

void VCFHeader::addFormatField(std::string field, std::string number,
                               std::string type,  std::string description) {
  if (formatmap.find(field) != formatmap.end()) {
    std::cerr << "Warning: Format field already exists: " << field << '\n';
    return;
  }
  formatmap[field] = "Number=" + number +
                     ",Type=" + type +
                     ",Description=\"" + description + "\"";
}

void VCFHeader::addSampleField(std::string field) {
  if (samplemap.find(field) != samplemap.end()) {
    std::cerr << "Warning: Sample field already exists: " << field << '\n';
    return;
  }
  samplemap[field] = field;
}

void VCFHeader::addContigField(std::string id, int len) {
  contigfieldmap[id] = std::to_string(len);
}

// ===========================================================================
// VCFEntry
// ===========================================================================

std::string VCFEntry::getRefString() const {
  std::string p = (bp->svtype == SVType::INDEL || id_num == 1) ? bp->ref : bp->alt;
  if (p.empty()) {
    std::cerr << "WARNING: Empty ref/alt field for bp\n";
    return "N";
  }
  return p;
}

std::string VCFEntry::getAltString(const SeqLib::BamHeader& header) const {

  // Indels just use bp->alt as-is (VCF-style insert/delete REF/ALT).
  if (bp->svtype == SVType::INDEL) {
    if (bp->alt.empty()) {
      std::cerr << "WARNING: Empty alt field for indel bp\n";
      return "N";
    }
    return bp->alt;
  }

  // SvABA2.0: symbolic SV emission for SvFormat::SYMBOLIC_WHEN_OBVIOUS.
  // The writer classifies each pair before calling; when symbolic_rep
  // is set we emit `<DEL>` / `<DUP>` / `<INV>` as the single-record ALT
  // and rely on END / SVLEN / SVTYPE INFO to describe the event.
  if (symbolic_rep)
    return "<" + symbolic_kind + ">";

  // SV breakends get the BND-style `N]chr:pos]` mate notation.
  const std::string ref = getRefString();

  std::stringstream ptag;
  if (id_num == 1) {
    ptag << bp->b2.gr.ChrName(header) << ':' << bp->b2.gr.pos1;
  } else {
    ptag << bp->b1.gr.ChrName(header) << ':' << bp->b1.gr.pos1;
  }

  std::stringstream alt;
  const char s1 = bp->b1.gr.strand;
  const char s2 = bp->b2.gr.strand;

  if (s1 == '+' && s2 == '+') {
    alt << ref << ']' << ptag.str() << ']';
  } else if (s1 == '+' && s2 == '-') {
    if (id_num == 1) alt << ref << '[' << ptag.str() << '[';
    else             alt << ']' << ptag.str() << ']' << ref;
  } else if (s1 == '-' && s2 == '+') {
    if (id_num == 1) alt << ']' << ptag.str() << ']' << ref;
    else             alt << ref << '[' << ptag.str() << '[';
  } else {
    alt << '[' << ptag.str() << '[' << ref;
  }

  return alt.str();
}

std::string VCFEntry::getIdString() const {
  if (bp->svtype != SVType::INDEL)
    return std::to_string(id) + ':' + std::to_string(id_num);
  return std::to_string(id);
}

std::unordered_map<std::string, std::string> VCFEntry::fillInfoFields() const {

  std::unordered_map<std::string, std::string> info;

  info["SPAN"] = std::to_string(bp->getSpan());
  info["SCTG"] = bp->cname;

  // SvABA2.0: common per-record metadata, populated for both SVs and
  // indels so downstream filters don't need to condition on shape.
  //
  // SOMATIC flag: gated on the somlod threshold carried per-entry
  // (`somatic_threshold`, defaulted 1.0, overridable via `svaba tovcf
  // --somlod N`). The svaba-run-time SomaticState::FAILED state is
  // honored as a hard veto — something earlier in the pipeline
  // explicitly said "not somatic", so we don't mark it regardless
  // of what LO_s looks like.
  if (bp->LO_s >= somatic_threshold &&
      bp->somatic != SomaticState::FAILED)
    info["SOMATIC"] = "";   // flag

  if (!bp->id.empty())
    info["EVENT"] = bp->id; // v3 bp_id (col 52 of bps.txt)

  {
    std::stringstream ss; ss << std::setprecision(4) << bp->max_lod;
    info["MAXLOD"] = ss.str();
  }
  {
    std::stringstream ss; ss << std::setprecision(4) << bp->LO_s;
    info["SOMLOD"] = ss.str();
  }

  if (bp->svtype != SVType::INDEL) {
    info["EVDNC"]  = svtype_to_string(bp->svtype);

    // SvABA2.0: SVTYPE is either BND (paired-breakend notation) or a
    // symbolic type (DEL/DUP/INV) when the writer has opted into
    // symbolic alleles and this record qualifies. symbolic_rep /
    // symbolic_kind are set by writeSvsSingleFile before toFileString
    // runs; when they're unset we fall back to legacy BND.
    if (symbolic_rep) {
      info["SVTYPE"] = symbolic_kind;  // "DEL" / "DUP" / "INV"
      // For symbolic alleles we must set END (SV end position) and
      // SVLEN (positive for INS/DUP, negative for DEL; absolute length
      // for INV because the inverted segment length is what matters).
      const int pos1_1idx = bp->b1.gr.pos1 + 1;
      const int pos2_1idx = bp->b2.gr.pos1 + 1;
      const int end_pos   = std::max(pos1_1idx, pos2_1idx);
      const int length    = std::abs(pos2_1idx - pos1_1idx);
      info["END"] = std::to_string(end_pos);
      if      (symbolic_kind == "DEL") info["SVLEN"] = std::to_string(-length);
      else                              info["SVLEN"] = std::to_string(length);
      // svaba uses primarily junction evidence (+ optional depth via
      // discordant pairs). "J" is the honest claim for assembly-only,
      // "DJ" for ASDIS (assembly+discordant). DSCRD-only is rare for
      // symbolic-eligible intrachrom events but we'd emit "D" there.
      if      (bp->svtype == SVType::ASSMB) info["SVCLAIM"] = "J";
      else if (bp->svtype == SVType::ASDIS) info["SVCLAIM"] = "DJ";
      else if (bp->svtype == SVType::DSCRD) info["SVCLAIM"] = "D";
      else                                   info["SVCLAIM"] = "J";
    } else {
      info["SVTYPE"]  = "BND";
      info["SVCLAIM"] = "J";
    }
  }

  if (!bp->repeat_seq.empty())
    info["REPSEQ"] = bp->repeat_seq;

  // NM / MATENM / MATEID: only for multi-alignment (BND-style) breakpoints.
  if (bp->num_align != 1) {
    info["MATEID"] = std::to_string(id) + ':' +
                     std::to_string(id_num == 1 ? 2 : 1);
    if (id_num == 1) {
      info["NM"]     = std::to_string(bp->b1.nm);
      info["MATENM"] = std::to_string(bp->b2.nm);
    } else {
      info["NM"]     = std::to_string(bp->b2.nm);
      info["MATENM"] = std::to_string(bp->b1.nm);
    }
  } else {
    info["NM"] = std::to_string(bp->b1.nm);
  }

  // MAPQ of the breakend this entry represents.
  info["MAPQ"] = std::to_string(id_num == 1 ? bp->b1.mapq : bp->b2.mapq);

  if (bp->num_align != 1) {
    // --- SV-only fields ---
    if (id_num == 1) {
      if (bp->b1.sub)      info["SUBN"] = std::to_string(bp->b1.sub);
      else if (bp->b2.sub) info["SUBN"] = std::to_string(bp->b2.sub);
    }

    if (!bp->homology.empty())  info["HOMSEQ"]    = bp->homology;
    if (!bp->insertion.empty()) info["INSERTION"] = bp->insertion;

    info["NUMPARTS"] = std::to_string(bp->num_align);

    if (bp->imprecise > 0) info["IMPRECISE"] = "";
    if (bp->secondary > 0) info["SECONDARY"] = "";

    if (info["EVDNC"] != "ASSMB") {
      info["DISC_MAPQ"] = std::to_string(
          id_num == 1 ? bp->dc.mapq1 : bp->dc.mapq2);
    }

  } else {
    // --- indel-only fields ---
    std::stringstream lod_ss;
    lod_ss << std::setprecision(4) << max_allele_lo(*bp);
    info["LOD"] = lod_ss.str();

    if (!bp->rs.empty())
      info["DBSNP"] = bp->rs;
  }

  return info;
}

std::string VCFEntry::toFileString(const SeqLib::BamHeader& header,
                                   QualMode qual_mode) const {

  std::stringstream out;

  // Sort INFO fields for stable, legacy-compatible output.
  auto info_map = fillInfoFields();
  std::vector<std::pair<std::string, std::string>> info_vec(
      info_map.begin(), info_map.end());
  std::sort(info_vec.begin(), info_vec.end(), compareInfoFields);

  std::string info;
  for (const auto& kv : info_vec) {
    // Don't emit HOMSEQ/HOMLEN/INSERTION when imprecise (legacy behavior).
    if (bp->imprecise && (kv.first == "HOMSEQ" ||
                          kv.first == "HOMLEN" ||
                          kv.first == "INSERTION"))
      continue;

    const bool is_flag = (g_flag_fields.count(kv.first) > 0);
    info += kv.first + (is_flag ? "" : "=") + kv.second + ';';
  }
  if (!info.empty())
    info.pop_back(); // trim trailing ';'

  // SvABA2.0: pick QUAL representation. Legacy (SUM_LO_PHRED) keeps
  // exactly the pre-2.0 behavior by passing bp->quality through; the
  // two new modes are for the `svaba tovcf` path, where the canonical
  // confidence lives in INFO/MAXLOD / INFO/SOMLOD and writing the sum-
  // of-LOs Phred into QUAL just invites confused downstream filters.
  std::string qual_str;
  switch (qual_mode) {
    case QualMode::MISSING:
      qual_str = ".";
      break;
    case QualMode::MAXLOD_PHRED: {
      // maxlod is already log10 — scale by 10 for Phred and cap at 99
      // to stay inside the conventional VCF QUAL range.
      const double v = std::max(0.0, bp->max_lod);
      const long   p = std::lround(10.0 * v);
      qual_str = std::to_string(std::min<long>(p, 99));
      break;
    }
    case QualMode::SUM_LO_PHRED:
    default:
      qual_str = std::to_string(bp->quality);
      break;
  }

  // For symbolic SV records, POS is the lower of the two breakend
  // positions (the event's 5' anchor); the INFO/END field carries the
  // other boundary. For BND records, POS is just the breakend this
  // entry represents (legacy behavior).
  int pos;
  std::string chr_name;
  if (symbolic_rep) {
    const int p1 = bp->b1.gr.pos1;
    const int p2 = bp->b2.gr.pos1;
    pos      = std::min(p1, p2);
    // Both breakends live on the same chrom when symbolic; use b1's.
    chr_name = bp->b1.gr.ChrName(header);
  } else {
    const BreakEnd* be = (id_num == 1) ? &bp->b1 : &bp->b2;
    pos      = be->gr.pos1;
    chr_name = be->gr.ChrName(header);
  }

  const char sep = '\t';
  out << chr_name               << sep
      << pos                     << sep
      << getIdString()           << sep
      << getRefString()          << sep
      << getAltString(header)    << sep
      << qual_str                << sep
      << bp->confidence          << sep
      << info                    << sep
      << (bp->svtype == SVType::INDEL ? kIndelFormat : kSvFormat);

  // SvABA2.0: per-sample columns (VCF column 10+). Iterating bp->allele
  // in-order is deterministic because SampleInfo is stored in a
  // std::map keyed by the 4-char bam prefix (n001, t001, ...). The
  // header's `colnames` must list samples in the same order for the
  // VCF to be well-formed — callers that build the header from the
  // bps.txt.gz row-0 sample list implicitly match this ordering.
  for (const auto& [_, al] : bp->allele) {
    out << sep << al.toFileString(bp->svtype);
  }

  return out.str();
}

bool VCFEntry::operator<(const VCFEntry& v) const {
  const BreakEnd* a = (id_num   == 1) ? &bp->b1   : &bp->b2;
  const BreakEnd* b = (v.id_num == 1) ? &v.bp->b1 : &v.bp->b2;
  return a->gr < b->gr;
}

bool VCFEntry::operator==(const VCFEntry& v) const {
  const BreakEnd* a = (id_num   == 1) ? &bp->b1   : &bp->b2;
  const BreakEnd* b = (v.id_num == 1) ? &v.bp->b1 : &v.bp->b2;
  return a->gr == b->gr;
}

// ===========================================================================
// VCFEntryPair
// ===========================================================================

VCFEntryPair::VCFEntryPair(std::shared_ptr<BreakPoint>& b) {
  bp     = b;
  e1.bp  = bp;
  e2.bp  = bp;
  e1.id_num = 1;
  e2.id_num = 2;

  // Compose a stable string describing this pair and hash it. Retry on
  // collision by salting with an increment counter, so ids remain unique
  // even if two genuinely-distinct pairs happen to collide.
  uint32_t hashed = 0;
  int salt = 0;
  do {
    const std::string s =
        std::to_string(bp->b1.gr.pos1) + '-' +
        std::to_string(bp->b2.gr.pos1) + '-' +
        std::to_string(bp->b1.gr.chr)  + '_' +
        std::to_string(bp->b2.gr.chr)  +
        bp->cname                      +
        std::to_string(salt);
    hashed = hash_string32(s);
    ++salt;
  } while (g_hash_avoid.find(hashed) != g_hash_avoid.end());

  g_hash_avoid.insert(hashed);
  e1.id = hashed;
  e2.id = hashed;
}

std::string VCFEntryPair::toFileString(const SeqLib::BamHeader& header,
                                       QualMode qual_mode) const {
  std::stringstream out;
  // If e1 is marked symbolic_rep, only emit it (single-record SV event).
  // Otherwise emit both breakends as a paired BND event.
  out << e1.toFileString(header, qual_mode) << '\n';
  if (!e1.symbolic_rep)
    out << e2.toFileString(header, qual_mode) << '\n';
  return out.str();
}

// ===========================================================================
// VCFFile::VCFFile — read bps.txt.gz, build entries, deduplicate.
// ===========================================================================

namespace {

// All of the boilerplate ##INFO/##FILTER/##FORMAT lines for SV and indel
// headers. Kept here rather than inline in the ctor so the ctor stays
// readable.
void populate_sv_header(VCFHeader& h) {

  // FILTERs
  h.addFilterField("NOLOCAL",        "Contig realigned to region outside of local assembly region, and no disc support.");
  h.addFilterField("LOCALMATCH",     "Contig realigned to assembly region without clipping");
  h.addFilterField("HIGHHOMOLOGY",   "Contig realigns with > 25% of readlength of homology. High probaility of assembly/mapping artifact");
  h.addFilterField("DUPREADS",       "Contig built from what appear to be duplicate reads (split reads all same contig cov))");
  h.addFilterField("NODISC",         "Rearrangement was not detected independently by assembly");
  h.addFilterField("LOWSUPPORT",     "Fewer than 2 split reads or < 4 total alt reads for ASDISC");
  h.addFilterField("COMPETEDISC",    "Discordant cluster found with nearly same breakpoints, but different strands for DSCRD event");
  h.addFilterField("LOWSPAN",        "Discordant read cluster (no split read support), and less than 10kb span and < 12 reads");
  h.addFilterField("LOWMAPQ",        "Assembly contig has non 60/60 mapq and no discordant support");
  h.addFilterField("LOWMATCHLEN",    "Assembly contig has fewer than 40 bases mapping uniquely to a reference locus (<100 if complex mapping or )");
  h.addFilterField("SINGLEBX",       "Variant is supported by only a single BX tag (if run with 10X Genomics data)");
  h.addFilterField("LOWQINVERSION",  "Assembly-only inversion of span < 300 and < 6 split reads. Common artifact in Illumina data");
  h.addFilterField("LOWMAPQDISC",    "Both clusters of reads failed to achieve mean mapq of > 30 for DSCRD");
  h.addFilterField("LOWSPLITSMALL",  "Fewer than 4 split reads for small events ( < 1500 bp)");
  h.addFilterField("LOWICSUPPORT",   "Less than 60bp of contig match on one end of an inter-chromosomal break");
  h.addFilterField("LOWAS",          "Alignment score of one end is less than 80% of contig length, or number of mismatch bases (NM) on one end is >= 10");
  h.addFilterField("WEAKSUPPORTHIREP","Fewer then 7 split reads for variant with >= 10 bases of repeat sequence (need to be more strict)");
  h.addFilterField("WEAKDISC",       "Fewer than 7 supporting discordant reads and no assembly support");
  h.addFilterField("TOOSHORT",       "Contig alignment for part of this rearrangement has <= 25bp match to reference");
  h.addFilterField("PASS",           "Strong assembly support, strong discordant support, or combined support. Strong MAPQ");
  h.addFilterField("MULTIMATCH",     "Low MAPQ and this contig fragment maps well to multiple locations");
  h.addFilterField("LOWSPANDSCRD",   "Discordant-only cluster is too small given isize distribution to call confidently");
  h.addFilterField("SIMPLESEQUENCE", "Major portion of one contig mapping falls in a simple sequence, as given by -R flag. Assembly-only filter");

  // INFOs
  h.addInfoField("REPSEQ",    "1", "String",  "Repeat sequence near the event");
  h.addInfoField("NM",        "1", "Integer", "Number of mismatches of this alignment fragment to reference");
  h.addInfoField("MATENM",    "1", "Integer", "Number of mismatches of partner alignment fragment to reference");
  h.addInfoField("SVTYPE",    "1", "String",  "Type of structural variant");
  h.addInfoField("HOMSEQ",    "1", "String",  "Sequence of base pair identical micro-homology at event breakpoints. Plus strand sequence displayed.");
  h.addInfoField("IMPRECISE", "0", "Flag",    "Imprecise structural variation");
  h.addInfoField("SECONDARY", "0", "Flag",    "SV calls comes from a secondary alignment");
  h.addInfoField("HOMLEN",    "1", "Integer", "Length of base pair identical micro-homology at event breakpoints");
  h.addInfoField("MAPQ",      "1", "Integer", "Mapping quality (BWA-MEM) of this fragement of the contig (-1 if discordant only)");
  h.addInfoField("MATEMAPQ",  "1", "Integer", "Mapping quality of the partner fragment of the contig");
  h.addInfoField("MATEID",    "1", "String",  "ID of mate breakends");
  h.addInfoField("SUBN",      "1", "Integer", "Number of secondary alignments associated with this contig fragment");
  h.addInfoField("NUMPARTS",  "1", "Integer", "If detected with assembly, number of parts the contig maps to. Otherwise 0");
  h.addInfoField("EVDNC",     "1", "String",  "Evidence for variant. ASSMB assembly only, ASDIS assembly+discordant. DSCRD discordant only, TSI_L templated-sequence insertion (local, e.g. AB or BC of an ABC), TSI_G global (e.g. AC of ABC)");
  h.addInfoField("SCTG",      "1", "String",  "Identifier for the contig assembled by svaba to make the SV call");
  h.addInfoField("INSERTION", "1", "String",  "Sequence insertion at the breakpoint.");
  h.addInfoField("SPAN",      "1", "Integer", "Distance between the breakpoints. -1 for interchromosomal");
  h.addInfoField("DISC_MAPQ", "1", "Integer", "Mean mapping quality of discordant reads mapped here");

  // SvABA2.0 (v3): VCF 4.5-aligned SV INFO fields.
  h.addInfoField("SOMATIC",   "0", "Flag",    "Somatic call per svaba's somatic LOD model (see BreakPoint::score_somatic)");
  h.addInfoField("EVENT",     "1", "String",  "ID of the SV event this breakend belongs to. Paired BND records share the same EVENT; a symbolic SV is its own event. Populated from BreakPoint::id (col 52 of bps.txt) when available.");
  h.addInfoField("END",       "1", "Integer", "End position for symbolic alleles (<DEL>/<DUP>/<INV>). Not set for BND.");
  h.addInfoField("SVLEN",     "1", "Integer", "Length of the SV in bp. Negative for deletions. Only set on symbolic alleles.");
  h.addInfoField("SVCLAIM",   "1", "String",  "Evidence type supporting the SV call, per VCFv4.5: J=junction-level (split/assembly), D=depth, DJ=both.");
  h.addInfoField("MAXLOD",    "1", "Float",   "Maximum per-sample log-odds (variant vs error) across all samples (see BreakPoint::max_lod)");
  h.addInfoField("SOMLOD",    "1", "Float",   "Somatic log-odds (LLR of somatic vs non-somatic; see SvabaModels::SomaticLOD)");

  // FORMATs (SV)
  h.addFormatField("GT",   "1", "String",  "Most likely genotype");
  h.addFormatField("AD",   "1", "Integer", "Allele depth: Number of reads supporting the variant");
  h.addFormatField("DP",   "1", "Integer", "Depth of coverage: Number of reads covering site.");
  h.addFormatField("GQ",   "1", "String",  "Genotype quality (currently not supported. Always 0)");
  h.addFormatField("PL",   ".", "Float",   "Normalized likelihood of the current genotype");
  h.addFormatField("SR",   "1", "Integer", "Number of spanning reads for this variant");
  h.addFormatField("DR",   "1", "Integer", "Number of discordant-supported reads for this variant");
  h.addFormatField("LO",   "1", "Float",   "Log-odds that this variant is real vs artifact");
  h.addFormatField("LO_n", "1", "Float",   "Log-odds that this variant is AF=0 vs AF>=0.5");
}

void populate_indel_header(VCFHeader& h) {

  // FILTERs
  h.addFilterField("LOWMAPQ",        "Assembly contig has less than MAPQ 10");
  h.addFilterField("SHORTALIGNMENT", "Matched (M) contig frag to left or right of indel < 20 bp");
  h.addFilterField("LOWLOD",         "LOD score is less than the cutoff");
  h.addFilterField("MULTIMATCH",     "Low MAPQ and this contig fragment maps well to multiple locations");
  h.addFilterField("LOWAS",          "Less than 80% of contig length is covered by a supporting read");
  h.addFilterField("NONVAR",         "0/0 is the most likely genotype");
  h.addFilterField("PASS",           "LOD score pass");
  h.addFilterField("VLOWAF",         "allelic fraction < 0.05");
  h.addFilterField("REPVAR",         "Multiple conflicting variants at a highly repetitive region");

  // INFOs
  h.addInfoField("SCTG",     "1", "String",  "Identifier for the contig assembled by svaba to make the indel call");
  h.addInfoField("MAPQ",     "1", "Integer", "Mapping quality (BWA-MEM) of the assembled contig");
  h.addInfoField("SPAN",     "1", "Integer", "Size of the indel");
  h.addInfoField("REPSEQ",   "1", "String",  "Repeat sequence near the variant");
  h.addInfoField("NM",       "1", "Integer", "Number of mismatches of this alignment fragment to reference");
  h.addInfoField("READNAMES",".", "String",  "IDs of ALT reads");
  h.addInfoField("BX",       ".", "String",  "Table of BX tag counts for supporting reads");
  h.addInfoField("DBSNP",    "0", "Flag",    "Variant found in dbSNP");
  h.addInfoField("LOD",      "1", "Float",   "Log of the odds that variant is real vs artifact");

  // SvABA2.0 (v3): VCF 4.5-aligned indel INFO fields.
  h.addInfoField("SOMATIC",  "0", "Flag",    "Somatic indel per svaba's somatic LOD model");
  h.addInfoField("EVENT",    "1", "String",  "Stable indel identifier (populated from BreakPoint::id when available)");
  h.addInfoField("SVTYPE",   "1", "String",  "INS / DEL for indels emitted with explicit SVTYPE");
  h.addInfoField("SVLEN",    "1", "Integer", "Length of the indel (negative for deletions)");
  h.addInfoField("MAXLOD",   "1", "Float",   "Maximum per-sample log-odds (variant vs error) across samples");
  h.addInfoField("SOMLOD",   "1", "Float",   "Somatic log-odds (LLR of somatic vs non-somatic)");

  // FORMATs (indel)
  h.addFormatField("GT",   "1", "String",  "Most likely genotype");
  h.addFormatField("AD",   "1", "Integer", "Allele depth: Number of reads supporting the variant");
  h.addFormatField("DP",   "1", "Integer", "Depth of coverage: Number of reads covering site.");
  h.addFormatField("GQ",   "1", "String",  "Genotype quality (currently not supported. Always 0)");
  h.addFormatField("PL",   ".", "Float",   "Normalized likelihood of the current genotype");
  h.addFormatField("SR",   "1", "Integer", "Number of spanning reads for this variant");
  h.addFormatField("CR",   "1", "Integer", "Number of cigar-supported reads for this variant");
  h.addFormatField("LO",   "1", "Float",   "Log-odds that this variant is real vs artifact");
  h.addFormatField("LO_n", "1", "Float",   "Log-odds that this variant is AF=0 vs AF>=0.5");
}

} // namespace

VCFFile::VCFFile(std::string              file,
                 std::string              id,
                 const SvabaSharedConfig& sc,
                 const VCFHeader&         vheader,
                 bool                     nopass,
                 bool                     m_verbose,
                 bool                     skip_dedup_in_ctor) {

  filename        = file;
  analysis_id     = id;
  verbose         = m_verbose;
  include_nonpass = nopass;
  skip_dedup      = skip_dedup_in_ctor;  // captured before ctor's deduplicate() runs

  // Seed both headers from the caller-provided base (which carries contig
  // and sample field lines), then layer on the SV- and indel-specific
  // INFO/FILTER/FORMAT blocks.
  sv_header    = vheader;
  indel_header = vheader;
  populate_sv_header(sv_header);
  populate_indel_header(indel_header);

  // Open the gzipped bps file stream.
  igzstream infile(file.c_str(), std::ios::in);
  if (!infile) {
    std::cerr << "Can't read file " << file << " for parsing VCF\n";
    std::exit(EXIT_FAILURE);
  }

  std::cerr << "...vcf - reading breakpoints: " << file << '\n';

  // Per-contig cap: if more than VCF_SECONDARY_CAP pairs share the same
  // cname, drop additional ones to avoid runaway output at pileups.
  std::unordered_map<std::string, int> cname_count;

  std::string line;
  // Skip the first row — it's the column header emitted by BreakPoint::header().
  std::getline(infile, line, '\n');

  int line_count = 0;
  while (std::getline(infile, line, '\n')) {

    if (line.empty() || line[0] == '#')
      continue;

    // BreakPoint parser will fail if any chr is unknown. Skip defensively.
    if (line.find("Unknown") != std::string::npos)
      continue;

    // Parse one bps.txt.gz row into a heap BreakPoint. Must be shared_ptr
    // because BreakPoint is non-copyable/non-movable.
    std::shared_ptr<BreakPoint> bp;
    try {
      bp = std::make_shared<BreakPoint>(line, &sc);
    } catch (const std::exception& e) {
      std::cerr << "Warning: failed to parse bps line: " << e.what() << '\n';
      continue;
    }

    // Honor the unfiltered/filtered flag from the caller. Gate on
    // `confidence` (the FILTER column string) rather than bp->pass —
    // the parser doesn't repopulate the `pass` bool from the dump,
    // and the scorer only updates `confidence`, so `confidence` is
    // the single source of truth for "did this record PASS".
    if (!include_nonpass && bp->confidence != "PASS")
      continue;

    // Per-contig cap: cap the number of entries per contig name.
    if (++cname_count[bp->cname] >= VCF_SECONDARY_CAP)
      continue;

    auto vpair = std::make_shared<VCFEntryPair>(bp);

    if (bp->svtype == SVType::INDEL)
      indels.emplace(line_count, vpair);
    else
      entry_pairs.emplace(line_count, vpair);

    ++line_count;
  }

  std::cerr << "...vcf - read in "
            << SeqLib::AddCommas(indels.size())      << " indels and "
            << SeqLib::AddCommas(entry_pairs.size()) << " SVs\n";
  std::cerr << "...vcf - SV deduplicating "
            << SeqLib::AddCommas(entry_pairs.size()) << " events\n";

  deduplicate(sc.header);

  std::cerr << "...vcf - SV deduplicated down to "
            << SeqLib::AddCommas(entry_pairs.size() - dups.size())
            << " break pairs\n";
}

// ===========================================================================
// Dedupe — SVs via paired-breakend interval-tree matching, indels via hash.
// ===========================================================================

namespace {

// Extension of GenomicRegion carrying the entry_pairs key and a pass flag.
class GenomicRegionWithID : public SeqLib::GenomicRegion {
 public:
  GenomicRegionWithID(int32_t c, uint32_t p1, uint32_t p2, int i, int p)
      : SeqLib::GenomicRegion(c, p1, p2), id(i), pass(p) {}
  uint32_t id   : 30;
  uint32_t pass : 2;
};

} // namespace

void VCFFile::deduplicate(const SeqLib::BamHeader& header) {

  // SvABA2.0: `svaba tovcf` receives a bps.txt.gz that has already been
  // sorted + deduped by scripts/svaba_postprocess.sh (step 3). Rerunning
  // the internal paired-interval-tree dedup on that input is wasted work
  // and can introduce its own collisions via the high-overlap blacklist
  // logic below. Short-circuit when the caller opts in.
  if (skip_dedup) {
    if (verbose)
      std::cerr << "...vcf - skip_dedup set, bypassing internal dedupe\n";
    return;
  }


  // Build separate interval trees for left (b1) and right (b2) breakends.
  SeqLib::GenomicRegionCollection<GenomicRegionWithID> grv1;
  SeqLib::GenomicRegionCollection<GenomicRegionWithID> grv2;

  for (const auto& kv : entry_pairs) {
    const auto& bpp = kv.second->bp;
    // Use the confidence string (not the stale bp->pass bool) to
    // decide PASS status for the dedup grouping — see BreakPoint.cpp
    // parser note: pass never gets repopulated from the dumped
    // confidence column, so it's 0 across all parsed records.
    const int bp_pass_int = (bpp->confidence == "PASS") ? 1 : 0;
    grv1.add(GenomicRegionWithID(bpp->b1.gr.chr, bpp->b1.gr.pos1,
                                 bpp->b1.gr.pos2, kv.first,
                                 bp_pass_int));
    grv2.add(GenomicRegionWithID(bpp->b2.gr.chr, bpp->b2.gr.pos1,
                                 bpp->b2.gr.pos2, kv.first,
                                 bp_pass_int));
  }
  grv1.CreateTreeMap();
  grv2.CreateTreeMap();
  assert(grv1.size() == grv2.size());

  // Regions blacklisted mid-pass because they turned out to be SV pileups.
  SeqLib::GenomicRegionCollection<SeqLib::GenomicRegion> high_overlap_blacklist;

  size_t count = 0;

  for (const auto& kv : entry_pairs) {
    const int current_key = kv.first;
    const auto& pair = kv.second;
    const auto& bpp  = pair->bp;

    // Tighter pad for small intrachrom events, looser for interchrom / large.
    const int pad =
        (bpp->b1.gr.chr != bpp->b2.gr.chr ||
         std::abs(bpp->b1.gr.pos1 - bpp->b2.gr.pos1) > DEDUPEPAD * 2)
            ? DEDUPEPAD
            : 10;

    if (verbose && (count % 20000 == 0))
      std::cerr << "...deduping at " << SeqLib::AddCommas(count) << '\n';
    ++count;

    // Blacklist check for either breakend.
    SeqLib::GenomicIntervalVector bad1, bad2;
    auto fb1 = high_overlap_blacklist.GetTree()->find(bpp->b1.gr.chr);
    auto fb2 = high_overlap_blacklist.GetTree()->find(bpp->b2.gr.chr);
    if (fb1 != high_overlap_blacklist.GetTree()->end()) {
      fb1->second.findOverlapping(bpp->b1.gr.pos1 - pad,
                                  bpp->b1.gr.pos1 + pad, bad1);
      if (!bad1.empty()) { dups.insert(current_key); continue; }
    }
    if (fb2 != high_overlap_blacklist.GetTree()->end()) {
      fb2->second.findOverlapping(bpp->b2.gr.pos1 - pad,
                                  bpp->b2.gr.pos1 + pad, bad2);
      if (!bad2.empty()) { dups.insert(current_key); continue; }
    }

    // Find all entries whose b1/b2 fall within pad of this one's b1/b2.
    SeqLib::GenomicIntervalVector giv1, giv2;
    auto ff1 = grv1.GetTree()->find(bpp->b1.gr.chr);
    auto ff2 = grv2.GetTree()->find(bpp->b2.gr.chr);
    assert(ff1 != grv1.GetTree()->end());
    assert(ff2 != grv2.GetTree()->end());

    ff1->second.findContained(bpp->b1.gr.pos1 - pad,
                              bpp->b1.gr.pos1 + pad, giv1);
    ff2->second.findContained(bpp->b2.gr.pos1 - pad,
                              bpp->b2.gr.pos1 + pad, giv2);

    // If this breakend is in a region hit by more than HIGH_OVERLAP_LIMIT
    // other breakends, treat it as an SV pileup: blacklist the region and
    // drop the current entry.
    if (giv1.size() > HIGH_OVERLAP_LIMIT) {
      SeqLib::GenomicRegion gr_bad = bpp->b1.gr;
      std::cerr << "...added " << gr_bad << " with " << giv1.size()
                << " overlaps to SV-pileup blacklist at dedupe\n"
                << "NB: Any SV with one breakend occuring more than "
                << HIGH_OVERLAP_LIMIT
                << " times will be removed. To change this behavior, "
                << "adjust HIGH_OVERLAP_LIMIT in vcf.cpp and recompile/run\n";
      gr_bad.Pad(BAD_REGION_PAD);
      high_overlap_blacklist.add(gr_bad);
      high_overlap_blacklist.CreateTreeMap();
      dups.insert(current_key);
      continue;
    }
    if (giv2.size() > HIGH_OVERLAP_LIMIT) {
      SeqLib::GenomicRegion gr_bad = bpp->b2.gr;
      std::cerr << "...added " << gr_bad << " with " << giv2.size()
                << " overlaps to SV-pileup blacklist at dedupe\n";
      gr_bad.Pad(BAD_REGION_PAD);
      high_overlap_blacklist.add(gr_bad);
      high_overlap_blacklist.CreateTreeMap();
      dups.insert(current_key);
      continue;
    }

    assert(giv1.size() <= HIGH_OVERLAP_LIMIT);
    assert(giv2.size() <= HIGH_OVERLAP_LIMIT);

    // A hit is a dup if the same entry_pairs key appears on BOTH sides, with
    // matching strands and matching pass status. Source of truth for
    // PASS is confidence == "PASS" (see note at top of deduplicate()).
    const bool   is_pass  = (pair->e1.bp->confidence == "PASS");
    const char   s1_here  = bpp->b1.gr.strand;
    const char   s2_here  = bpp->b2.gr.strand;

    std::unordered_map<int, std::pair<size_t, size_t>> key_count;

    for (const auto& j : giv1) {
      const auto& hit = grv1.at(j.value);
      if (hit.id == current_key) continue;
      if ((is_pass == static_cast<bool>(hit.pass)) &&
          entry_pairs[hit.id]->bp->b1.gr.strand == s1_here)
        ++key_count[hit.id].first;
    }
    for (const auto& j : giv2) {
      const auto& hit = grv2.at(j.value);
      if (hit.id == current_key) continue;
      if ((is_pass == static_cast<bool>(hit.pass)) &&
          entry_pairs[hit.id]->bp->b2.gr.strand == s2_here)
        ++key_count[hit.id].second;
    }

    // If this entry has a "twin" with a worse (lower) score than ours, mark
    // the twin — not ours — as a dup. Ties break by `operator<` on BreakPoint.
    // Exception: TSI_LOCAL/TSI_GLOBAL annotate the same underlying event in
    // two ways and should NOT dup each other.
    for (const auto& kc : key_count) {
      if (kc.second.first == 0 || kc.second.second == 0)
        continue;

      const auto& other = entry_pairs[kc.first];

      if (*bpp < *other->bp) {
        const bool local_global_pair =
            (bpp->svtype          == SVType::TSI_LOCAL  &&
             other->bp->svtype    == SVType::TSI_GLOBAL) ||
            (bpp->svtype          == SVType::TSI_GLOBAL &&
             other->bp->svtype    == SVType::TSI_LOCAL);
        if (local_global_pair) {
          // don't dup: intentional double-annotation
        } else {
          dups.insert(kc.first);
        }
      }
    }
  }

  // --- Indel dedupe: cheap hash on (chr, pos, REF, ALT). ---
  if (verbose)
    std::cerr << "...vcf - hashing " << SeqLib::AddCommas(indels.size())
              << " indels for dedupe\n";

  std::unordered_set<std::string> seen;
  VCFEntryPairMap deduped_indels;
  deduped_indels.reserve(indels.size());

  for (const auto& kv : indels) {
    std::string key;
    try {
      const auto& b1 = kv.second->e1.bp->b1;
      key = std::to_string(b1.gr.chr) + ':' +
            std::to_string(b1.gr.pos1) + '_' +
            kv.second->e1.getRefString() + '_' +
            kv.second->e1.getAltString(header);
    } catch (...) {
      std::cerr << "indel dedupe: error building hash key for entry "
                << kv.first << '\n';
      continue;
    }
    if (seen.insert(key).second)
      deduped_indels.emplace(kv.first, kv.second);
  }
  indels = std::move(deduped_indels);
}

// ===========================================================================
// Writers
// ===========================================================================

void VCFFile::writeIndels(std::string              basename,
                          bool                     /*zip*/,
                          bool                     /*onefile*/,
                          const SeqLib::BamHeader& header) const {

  // The legacy writer emitted a single combined indel VCF regardless of
  // zip/onefile flags; we preserve that behavior to avoid surprising
  // downstream consumers.
  const std::string out_name = basename + "germline.indel.vcf";

  std::ofstream out(out_name);
  if (!out) {
    std::cerr << "vcf: failed to open " << out_name << " for writing\n";
    return;
  }

  out << indel_header << '\n';

  // Collect e1 of every indel pair into a vector, sort by genomic position.
  VCFEntryVec rows;
  rows.reserve(indels.size());
  for (const auto& kv : indels)
    rows.push_back(kv.second->e1);

  std::sort(rows.begin(), rows.end());

  for (const auto& r : rows) {
    // PASS/FILTER gate: use the confidence string as source of truth
    // (bp->pass is a stale bool; see parser note in BreakPoint.cpp).
    if (r.bp->confidence != "PASS" && !include_nonpass)
      continue;
    out << r.toFileString(header) << '\n';
  }
}

void VCFFile::writeSVs(std::string              basename,
                       bool                     /*zip*/,
                       bool                     /*onefile*/,
                       const SeqLib::BamHeader& header) const {

  const std::string out_name = basename + "germline.sv.vcf";

  std::ofstream out(out_name);
  if (!out) {
    std::cerr << "vcf: failed to open " << out_name << " for writing\n";
    return;
  }

  out << sv_header << '\n';

  // Pull both breakends of each surviving pair into a flat vector, sort,
  // then emit. Skip entries marked as dups during dedupe.
  VCFEntryVec rows;
  rows.reserve(entry_pairs.size() * 2);
  for (const auto& kv : entry_pairs) {
    if (dups.count(kv.first)) continue;
    rows.push_back(kv.second->e1);
    rows.push_back(kv.second->e2);
  }

  std::sort(rows.begin(), rows.end());

  for (const auto& r : rows) {
    // PASS/FILTER gate: use the confidence string as source of truth
    // (bp->pass is a stale bool; see parser note in BreakPoint.cpp).
    if (r.bp->confidence != "PASS" && !include_nonpass)
      continue;
    out << r.toFileString(header) << '\n';
  }
}

// ===========================================================================
// SvABA2.0: symbolic-allele classifier + single-file writers (tovcf path).
// ===========================================================================

namespace {

// Classify a BreakPoint pair to decide whether it can be emitted as a
// single symbolic-allele record (<DEL>/<DUP>/<INV>) instead of two paired
// BND records. Conservative:
//   - both breakends must live on the same chromosome
//   - INDEL records keep their existing REF/ALT shape (not this path)
//   - strand pattern determines the symbolic kind:
//       +/+  or  -/-  -> <DEL>
//       -/+            -> <DUP>
//       +/-            -> <INV>
//   - anything ambiguous / inter-chrom falls through to BND.
// Returns empty string when the pair should stay as paired BND.
std::string classify_symbolic_kind(const BreakPoint& bp) {
  if (bp.svtype == SVType::INDEL)        return "";
  if (bp.b1.gr.chr != bp.b2.gr.chr)      return "";
  const char s1 = bp.b1.gr.strand;
  const char s2 = bp.b2.gr.strand;
  if (s1 == '+' && s2 == '+')            return "DEL";
  if (s1 == '-' && s2 == '-')            return "DEL";
  if (s1 == '-' && s2 == '+')            return "DUP";
  if (s1 == '+' && s2 == '-')            return "INV";
  return "";
}

// Sort key used by the single-file writers: for a BND entry use its
// represented breakend; for a symbolic entry use the lower of the two
// endpoints (matching VCF convention for the POS of a symbolic SV).
bool entry_lt_for_output(const VCFEntry& a, const VCFEntry& b) {
  const BreakEnd* ba = a.symbolic_rep
      ? (a.bp->b1.gr.pos1 < a.bp->b2.gr.pos1 ? &a.bp->b1 : &a.bp->b2)
      : (a.id_num == 1 ? &a.bp->b1 : &a.bp->b2);
  const BreakEnd* bb = b.symbolic_rep
      ? (b.bp->b1.gr.pos1 < b.bp->b2.gr.pos1 ? &b.bp->b1 : &b.bp->b2)
      : (b.id_num == 1 ? &b.bp->b1 : &b.bp->b2);
  return ba->gr < bb->gr;
}

// Open path for writing (gzipped or plain). Returns empty unique_ptr on
// failure, with the error already logged to stderr.
std::unique_ptr<std::ostream> open_vcf_out(const std::string& path, bool gzip) {
  if (gzip) {
    auto g = std::make_unique<ogzstream>(path.c_str(), std::ios::out);
    if (!g->good()) {
      std::cerr << "vcf: failed to open " << path << " for writing\n";
      return {};
    }
    return g;
  }
  auto f = std::make_unique<std::ofstream>(path);
  if (!f->good()) {
    std::cerr << "vcf: failed to open " << path << " for writing\n";
    return {};
  }
  return f;
}

} // namespace

void VCFFile::writeSvsSingleFile(const std::string&       path,
                                 bool                     gzip,
                                 const SeqLib::BamHeader& header) const {

  auto out_owner = open_vcf_out(path, gzip);
  if (!out_owner) return;
  std::ostream& out = *out_owner;

  out << sv_header << '\n';

  // Build the emit list. For each pair:
  //   - skip dups (if dedup ran)
  //   - classify symbolic eligibility when sv_format == SYMBOLIC_WHEN_OBVIOUS
  //   - emit one e1 (symbolic) or two (BND) entries into rows
  VCFEntryVec rows;
  rows.reserve(entry_pairs.size() * 2);
  for (const auto& kv : entry_pairs) {
    if (dups.count(kv.first)) continue;

    VCFEntryPair pair_copy = *kv.second;  // shallow copy to stamp symbolic_rep
    if (sv_format == SvFormat::SYMBOLIC_WHEN_OBVIOUS) {
      const std::string kind = classify_symbolic_kind(*pair_copy.bp);
      if (!kind.empty()) {
        pair_copy.e1.symbolic_rep  = true;
        pair_copy.e1.symbolic_kind = kind;
        rows.push_back(pair_copy.e1);
        continue;
      }
    }
    rows.push_back(pair_copy.e1);
    rows.push_back(pair_copy.e2);
  }

  std::sort(rows.begin(), rows.end(), entry_lt_for_output);

  for (auto& r : rows) {
    // PASS/FILTER gate: use the confidence string as source of truth
    // (bp->pass is a stale bool; see parser note in BreakPoint.cpp).
    if (r.bp->confidence != "PASS" && !include_nonpass)
      continue;
    // Stamp the per-entry somlod threshold so fillInfoFields applies
    // the right cutoff to the SOMATIC flag.
    r.somatic_threshold = somatic_threshold;
    out << r.toFileString(header, qual_mode) << '\n';
  }
}

void VCFFile::writeIndelsSingleFile(const std::string&       path,
                                    bool                     gzip,
                                    const SeqLib::BamHeader& header) const {

  auto out_owner = open_vcf_out(path, gzip);
  if (!out_owner) return;
  std::ostream& out = *out_owner;

  out << indel_header << '\n';

  // Indels are always a single REF/ALT record per event (no BND pair),
  // so this is simpler than the SV path: collect e1 of every indel,
  // sort, emit.
  VCFEntryVec rows;
  rows.reserve(indels.size());
  for (const auto& kv : indels)
    rows.push_back(kv.second->e1);

  std::sort(rows.begin(), rows.end());

  for (auto& r : rows) {
    // PASS/FILTER gate: use the confidence string as source of truth
    // (bp->pass is a stale bool; see parser note in BreakPoint.cpp).
    if (r.bp->confidence != "PASS" && !include_nonpass)
      continue;
    r.somatic_threshold = somatic_threshold;
    out << r.toFileString(header, qual_mode) << '\n';
  }
}
