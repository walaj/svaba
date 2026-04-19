#include "AlignedContig.h"
#include "SvabaUtils.h"
#include "ContigAlignmentScore.h"

#include <cctype>
#include <iomanip>
#include <sstream>
#include <unordered_set>

// SvABA2.0: the old ASCII emitter (printToAlignmentsFile → alignments.txt.gz)
// and its support header PlottedRead.h have been removed. All alignment-dump
// data now goes through the structured r2c TSV (printToR2CTsv /
// r2cTsvHeader), which viewers parse and re-plot on demand.

// bav is all of the alignments of the contig to the reference
AlignedContig::AlignedContig(BamRecordPtrVector& bav,
			     const GenomicRegion& region, 
			     const SvabaSharedConfig* sc_) : sc(sc_) {
  
  if (!bav.size())
    return;
  
  // find the longest sequence, taking the first one.
  // make sure sequence dir is set to same as first alignment
  for (const auto& i : bav) {
    if (i->Sequence().length() > m_seq.length()) {
      if (i->ReverseFlag() == bav.front()->ReverseFlag()) {
	m_seq = i->Sequence();
      } else {
	m_seq = i->Sequence();
	SeqLib::rcomplement(m_seq);
      }
    }
  }
  
  // set the sequence. Convention is store as it came off
  // the assembler for first alignment
  if (bav.front()->ReverseFlag()) 
    SeqLib::rcomplement(m_seq);

  // instantiate coverage of reads on the contig
  aligned_coverage = std::vector<int>(m_seq.length(), 0);
  
  // find the number of primary alignments
  size_t num_align = 0;
  for (const auto& i : bav)
    if (!i->SecondaryFlag())
      ++num_align;
  
  // make the individual alignments and add
  for (auto& i : bav) {
    bool flip = (m_seq != i->Sequence()); // if the seq was flipped, need to flip the AlignmentFragment
    if (!i->SecondaryFlag()) {
      m_frag_v.emplace_back(AlignmentFragment(i, flip, region, sc));
      m_frag_v.back().num_align = num_align;
    } else if (m_frag_v.size()) {
      m_frag_v.back().secondaries.emplace_back(AlignmentFragment(i, flip, region, sc));
    } else {
      throw std::runtime_error("Secondary alignment encountered with no primary.");
    }

    // set the sequence for the aligned contig in the flipped (or not) orientation
    //m_frag_v.back().m_seq = m_seq;
    
    // set the aligned coverage
    SeqLib::Cigar cig = i->GetCigar();
    size_t p = 0;
    if (!i->SecondaryFlag())
      for (auto& c : cig) {
	for (size_t j = 0; j < c.Length(); ++j) {
	  if (c.Type() == 'M' || c.Type() == 'I')  // consumes contig and not clip
	    ++aligned_coverage[p];
	  if (c.ConsumesQuery() || c.Type() == 'H') // consumes contig, move iterator
	    ++p;
	}
      }
  }
  
  if (!m_frag_v.size())
    throw std::runtime_error("AlignedContig - no AlignmentFragments made");

  // extract the indels on primary alignments
  for (auto& f : m_frag_v) 
    f.SetIndels();

  // sort fragments by order on fwd-strand contig
  // underneath this calls AlignmentFragment::operator<
  // which is defined based on the "flip-convention"
  // alignment position on the contig, NOT on the genome
  std::sort(m_frag_v.begin(), m_frag_v.end());
  
  // get rearrangement breaks out of it
  setMultiMapBreakPairs();
  
  // filter indels that land too close to a multi-map break
  //filterIndelsAtMultiMapSites(5);
}

// void AlignedContig::filterIndelsAtMultiMapSites(size_t buff) {
  
//   if (m_frag_v.size() < 2)
//     return;
  
//   // make the ranges ON CONIG for the multimaps
//   SeqLib::GRC grc;
//   if (m_global_bp->b1.cpos > m_global_bp->b2.cpos) {
//     grc.add(SeqLib::GenomicRegion(0, m_global_bp->b2.cpos-buff, m_global_bp->b1.cpos+buff)); // homology
//   }
//     else {
//       grc.add(SeqLib::GenomicRegion(0, m_global_bp->b1.cpos-buff, m_global_bp->b2.cpos+buff)); // insertion	
//     }    
    
//     for (auto& i : m_local_breaks) {
//       if (i.b1.cpos > i.b2.cpos) {
// 	grc.add(SeqLib::GenomicRegion(0, i.b2.cpos-buff, i.b1.cpos+buff)); // homology
//       }
//       else {
// 	grc.add(SeqLib::GenomicRegion(0, i.b1.cpos-buff, i.b2.cpos+buff)); // insertion	
//       }
//     }
//     grc.CreateTreeMap();

//     // check if 
//     for (auto& i : m_frag_v) {
//       BPVec new_indel_vec;
//       for (auto& b : i.m_indel_breaks) {
// 	if (!grc.CountOverlaps(SeqLib::GenomicRegion(0, b.b1.cpos, b.b2.cpos)))
// 	  new_indel_vec.push_back(b);
//       }
//       i.m_indel_breaks = new_indel_vec;
//     }
    
//   }
  
SeqLib::GenomicRegionVector AlignedContig::getAsGenomicRegionVector() const {
  SeqLib::GenomicRegionVector g;
  for (auto& i : m_frag_v)
    g.push_back(i.m_align->AsGenomicRegion());
  return g;
}

// void AlignedContig::printContigFasta(std::ofstream& os) const {
//   os << ">" << getContigName() << std::endl;
//   os << getSequence() << std::endl;
// }

void AlignedContig::blacklist(SeqLib::GRC &grv) {
  
  // loop through the indel breaks and blacklist
  for (auto& i : m_frag_v) 
    for (auto& j : i.m_indel_breaks) 
      j->checkBlacklist(grv);
}

void AlignedContig::splitCoverage() { 
  
  // for (auto& i : m_local_breaks_secondaries) 
  //   i.splitCoverage(m_bamreads);
  
  // for (auto& i : m_global_bp_secondaries) 
  //   i->splitCoverage(m_bamreads);
  
  for (auto& i : m_frag_v) 
    for (auto& j : i.m_indel_breaks) 
      j->splitCoverage(m_bamreads);
  
  for (auto& i : m_local_breaks) 
    i->splitCoverage(m_bamreads);
  
  if (m_global_bp) 
    m_global_bp->splitCoverage(m_bamreads);
  
}

// ---------------------------------------------------------------------------
// SvABA2.0: structured re-plot format.
//
// printToR2CTsv emits this contig's info in a tab-separated form that
// downstream tools (e.g. viewer/bps_explorer.html) parse and re-plot on
// demand. It's the successor to the old printToAlignmentsFile() /
// alignments.txt.gz pre-rendered ASCII output, which has been removed.
//
// The schema is record-type discriminated: one "contig" row per contig,
// then one "read" row per r2c-aligned read, all sharing a `contig_name`
// key. Contig-only fields (sequence, fragments, BPs, read count) are
// populated on the contig row and "NA" on read rows; read-only fields
// are populated on read rows and "NA" on the contig row.
//
// Anything that was encoded only implicitly in the old ASCII art (e.g.
// fragment orientation via '>' vs '<', or leading/trailing soft-clip
// bases via lowercase letters) is here available as an explicit field
// so viewers don't have to reverse-engineer it.
// ---------------------------------------------------------------------------

namespace {
// TSV-escape: replace tab/CR/LF with spaces so a single column value can't
// break row parsing. Used everywhere we emit free-form text into a cell.
std::string r2cEscape(const std::string& s) {
  std::string out;
  out.reserve(s.size());
  for (char c : s) {
    if (c == '\t' || c == '\n' || c == '\r') out.push_back(' ');
    else out.push_back(c);
  }
  return out;
}

// Serialize a SeqLib::Cigar to the standard string form (e.g. "30M1D15M")
// without relying on the overload of operator<< (which is fine but writes
// through a stream; we want a plain std::string for cell composition).
std::string r2cCigarString(const SeqLib::Cigar& c) {
  std::ostringstream oss;
  for (auto it = c.begin(); it != c.end(); ++it)
    oss << it->Length() << static_cast<char>(it->Type());
  return oss.str();
}
}  // namespace

// Header line for the r2c TSV. Callers should write this once before any
// rows. Not a row itself — do not prefix with a comment char.
std::string AlignedContig::r2cTsvHeader() {
  static const char* const kCols[] = {
    "record_type",        // "contig" or "read"
    "contig_name",        // shared key: group reads back to their contig
    "contig_len",         // set on both row types for convenience
    "contig_seq",         // assembled sequence     [contig only; NA on read]
    "frags",              // BWA hits of contig vs. ref [contig only]
    "bps",                // breakpoints on this contig  [contig only]
    "n_reads",            // # read rows that follow this contig [contig only]
    // read-only:
    "read_id",            // svabaRead UniqueName (sample-prefixed qname)
    "read_chrom",         // chromosome the read originally mapped to (BAM)
    "read_pos",           // BAM alignment position
    "read_flag",          // 0 if r2c forward, 16 if reverse-complemented
    "r2c_cigar",          // read-to-contig CIGAR (e.g. "24M1D15M7I1M7D103M")
    "r2c_start",          // start_on_contig (0-based)
    "r2c_rc",             // 0/1: was r2c alignment reverse-complemented?
    "r2c_nm",             // NM (mismatch count) from r2c alignment
    "support",            // "split" | "disc" | "both" | "none" (categorical)
    // SvABA2.0: unambiguous per-BP support lists. The categorical
    // `support` above is the OR over these — use these columns when
    // you need to know *which* BP(s) a read supports on this contig.
    // Both are comma-joined lists of bp_id; "NA" if no match.
    "split_bps",          // bp_ids where this read is a split-supporter
    "disc_bps",           // bp_ids where this read (or its mate) is in a
                          // discordant cluster associated with the BP
    "r2c_score",          // alignment score under r2c CIGAR
    "native_score",       // alignment score under BAM CIGAR (or corrected)
    "read_seq"            // raw read sequence (pre-rc, pre-gap-expand)
  };
  std::ostringstream oss;
  for (size_t i = 0; i < sizeof(kCols) / sizeof(kCols[0]); ++i) {
    if (i) oss << '\t';
    oss << kCols[i];
  }
  return oss.str();
}

std::string AlignedContig::printToR2CTsv(const SeqLib::BamHeader& h) const {
  std::ostringstream out;
  const std::string cname = this->getContigName();
  const std::string NA    = "NA";

  // --- per-contig fields ------------------------------------------------

  // frags: "|"-separated; within a frag, ":"-separated:
  //   chr:pos:strand:cigar:mapq:cpos_break1:cpos_break2:gpos_break1:gpos_break2:flipped
  std::string frags;
  for (size_t i = 0; i < m_frag_v.size(); ++i) {
    const auto& frag = m_frag_v[i];
    const auto& a = frag.m_align;
    if (!a) continue;
    if (!frags.empty()) frags.push_back('|');
    frags += r2cEscape(a->ChrName(sc->header));
    frags += ':'; frags += std::to_string(a->Position());
    frags += ':'; frags += (a->ReverseFlag() ? '-' : '+');
    frags += ':'; frags += r2cEscape(a->CigarString());
    frags += ':'; frags += std::to_string(a->MapQuality());
    frags += ':'; frags += std::to_string(frag.break1);
    frags += ':'; frags += std::to_string(frag.break2);
    frags += ':'; frags += std::to_string(frag.gbreak1);
    frags += ':'; frags += std::to_string(frag.gbreak2);
    frags += ':'; frags += (frag.flipped ? '1' : '0');
  }
  if (frags.empty()) frags = NA;

  // bps: "|"-separated; within a bp, ":"-separated (10 fields):
  //   kind:bp_id:chr1:pos1:strand1:chr2:pos2:strand2:span:insertion
  // kind ∈ {global, multi, indel}. bp_id is the unique BP identifier
  // assigned by svabaThreadUnit::next_bp_id() ("." if unset — should
  // not happen in normal runs). insertion is "." if none. All fields
  // fixed-position so the field count stays at 10 regardless of row.
  auto bpCell = [&](const char* kind, const BreakPoint& bp) -> std::string {
    std::ostringstream os;
    os << kind << ':'
       << (bp.id.empty() ? std::string(".") : r2cEscape(bp.id)) << ':'
       << r2cEscape(bp.b1.gr.ChrName(h)) << ':' << bp.b1.gr.pos1 << ':'
       << (bp.b1.gr.strand ? bp.b1.gr.strand : '.') << ':'
       << r2cEscape(bp.b2.gr.ChrName(h)) << ':' << bp.b2.gr.pos1 << ':'
       << (bp.b2.gr.strand ? bp.b2.gr.strand : '.') << ':'
       << bp.getSpan() << ':'
       << (bp.insertion.empty() ? std::string(".") : r2cEscape(bp.insertion));
    return os.str();
  };

  std::string bps;
  auto appendBp = [&](const std::string& cell) {
    if (!bps.empty()) bps.push_back('|');
    bps += cell;
  };
  if (m_global_bp) appendBp(bpCell("global", *m_global_bp));
  for (const auto& b : m_local_breaks) if (b) appendBp(bpCell("multi", *b));
  for (const auto& frag : m_frag_v)
    for (const auto& b : frag.m_indel_breaks)
      if (b) appendBp(bpCell("indel", *b));
  if (bps.empty()) bps = NA;

  const size_t n_reads = m_bamreads.size();

  // --- contig row -------------------------------------------------------
  //
  // Fields in order (21 total, mirroring r2cTsvHeader):
  //   record_type contig_name contig_len contig_seq frags bps n_reads
  //   read_id read_chrom read_pos read_flag r2c_cigar r2c_start r2c_rc
  //   r2c_nm support split_bps disc_bps r2c_score native_score read_seq
  // On the contig row, all read-specific fields are NA.
  out << "contig\t"
      << r2cEscape(cname)   << '\t'
      << m_seq.length()     << '\t'
      << r2cEscape(m_seq)   << '\t'
      << frags              << '\t'
      << bps                << '\t'
      << n_reads            << '\t'
      // read-only columns, 14 of them:
      << NA << '\t' << NA << '\t' << NA << '\t' << NA << '\t' << NA << '\t'
      << NA << '\t' << NA << '\t' << NA << '\t' << NA << '\t' << NA << '\t'
      << NA << '\t' << NA << '\t' << NA << '\t' << NA << '\n';

  // --- per-BP support tracking -----------------------------------------
  //
  // Replaces the older flat split_supporters / disc_supporters sets.
  // For each BP on this contig we collect:
  //   - bp_id: the unique identifier (may be "." if unset — shouldn't
  //     happen in a normal svaba run, but we don't block emission on it)
  //   - split_reads: UniqueNames credited as split-supporters for this
  //     specific BP by BreakPoint::splitCoverage
  //   - disc_reads:  UniqueNames (of reads or mates) appearing in this
  //     BP's own DiscordantCluster (populated by
  //     CombineWithDiscordantClusterMap). Note the scope is per-BP, not
  //     per-contig — a contig can have multiple BPs, and a discordant
  //     cluster may be attributed to only some of them.
  //
  // Per-read output then walks this vector and concatenates the bp_ids
  // this read shows up under into `split_bps` / `disc_bps`. A read
  // supporting two different BPs on the same contig (e.g. a multi-map
  // contig with a global + a local BP) will show both ids, comma-
  // joined, instead of the previous unambiguous "split" tag that hid
  // which BP it actually supported.
  struct BpSupportEntry {
    std::string bp_id;
    std::unordered_set<std::string> split_reads;
    std::unordered_set<std::string> disc_reads;
  };
  std::vector<BpSupportEntry> per_bp_support;
  per_bp_support.reserve(8);
  for (const auto& bp : getAllBreakPoints()) {
    if (!bp) continue;
    BpSupportEntry e;
    e.bp_id = bp->id.empty() ? std::string(".") : bp->id;
    for (const auto& un : bp->getAllSupportingReads()) e.split_reads.insert(un);
    for (const auto& kv : bp->dc.reads) e.disc_reads.insert(kv.first);
    for (const auto& kv : bp->dc.mates) e.disc_reads.insert(kv.first);
    per_bp_support.push_back(std::move(e));
  }

  auto fmt_score = [](double s) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(1) << s;
    return os.str();
  };

  // --- read rows --------------------------------------------------------
  for (const auto& i : m_bamreads) {
    const std::string sr       = i->UniqueName();
    const r2c         this_r2c = i->GetR2C(getContigName());

    // r2c CIGAR: empty if no r2c alignment (shouldn't happen for reads
    // listed here, but guard anyway).
    std::string r2c_cig_str = r2cCigarString(this_r2c.cig);
    if (r2c_cig_str.empty()) r2c_cig_str = NA;

    // Scores: same source-of-truth precedence BreakPoint::splitCoverage
    // uses for its "r2c better than native" gate — we reproduce it here
    // so the TSV lets you diagnose split-coverage decisions by field
    // inspection alone.
    std::string r2c_score_str = NA;
    if (this_r2c.cig.size() > 0) {
      r2c_score_str = fmt_score(
          svaba::readAlignmentScore(this_r2c.cig, this_r2c.nm));
    }
    std::string native_score_str = NA;
    {
      SeqLib::Cigar native_cig;
      int32_t       native_nm = -1;
      if (i->corrected_native_cig.size() > 0) {
        native_cig = i->corrected_native_cig;
        native_nm  = i->corrected_native_nm;
      } else {
        native_cig = i->GetCigar();
        i->GetIntTag("NM", native_nm);
      }
      if (native_cig.size() > 0) {
        native_score_str = fmt_score(
            svaba::readAlignmentScore(native_cig, native_nm));
      }
    }

    // Build the per-BP support lists for this specific read. Comma-
    // joined bp_ids, "NA" if this read doesn't support any BP in that
    // category on this contig. The categorical `support` is derived
    // from whether either list is non-empty, so the new columns are
    // a strict superset of the old info — downstream code keying on
    // `support == "split"` still works.
    std::string split_bps_str;
    std::string disc_bps_str;
    bool any_split = false, any_disc = false;
    for (const auto& e : per_bp_support) {
      if (e.split_reads.count(sr)) {
        if (!split_bps_str.empty()) split_bps_str.push_back(',');
        split_bps_str += e.bp_id;
        any_split = true;
      }
      if (e.disc_reads.count(sr)) {
        if (!disc_bps_str.empty()) disc_bps_str.push_back(',');
        disc_bps_str += e.bp_id;
        any_disc = true;
      }
    }
    if (split_bps_str.empty()) split_bps_str = NA;
    if (disc_bps_str.empty())  disc_bps_str  = NA;

    const char* support_kind = any_split && any_disc ? "both"
                             : any_split             ? "split"
                             : any_disc              ? "disc"
                             :                         "none";

    // read_chrom: BAM-native chromosome; may be empty for unmapped mates.
    std::string read_chrom =
        (i->ChrID() >= 0) ? i->ChrName(h) : std::string("*");

    out << "read\t"
        << r2cEscape(cname)     << '\t'
        << m_seq.length()       << '\t'
        << NA                   << '\t'   // contig_seq
        << NA                   << '\t'   // frags
        << NA                   << '\t'   // bps
        << NA                   << '\t'   // n_reads
        // read fields:
        << r2cEscape(sr)            << '\t'
        << r2cEscape(read_chrom)    << '\t'
        << i->Position()            << '\t'
        << (this_r2c.rc ? 16 : 0)   << '\t'
        << r2c_cig_str              << '\t'
        << this_r2c.start_on_contig << '\t'
        << (this_r2c.rc ? 1 : 0)    << '\t'
        << this_r2c.nm              << '\t'
        << support_kind             << '\t'
        << split_bps_str            << '\t'
        << disc_bps_str             << '\t'
        << r2c_score_str            << '\t'
        << native_score_str         << '\t'
        << r2cEscape(i->Sequence()) << '\n';
  }

  return out.str();
}

void AlignedContig::setMultiMapBreakPairs() {
  
  // if single mapped contig, nothing to do here
  if (m_frag_v.size() == 1)
    return;
   
  // walk along the ordered contig list and make the breakpoint pairs
  // these are the so called "local" breakpoints, meaning they are
  // AB rearrangemnts without an intervening sequence
  // so this loop will do AB, BC, CD, etc
  for (AlignmentFragmentVector::const_iterator it = m_frag_v.begin();
       it != m_frag_v.end() - 1;
       it++) {
    
    AlignmentFragmentVector bwa_hits_1, bwa_hits_2;
    
    // first add the primary alignments
    bwa_hits_1.push_back(*it);
    bwa_hits_2.push_back(*(it+1));
    
    // then the secondaries 
    // bwa_hits_1.insert(bwa_hits_1.end(),
    // 		      it->secondaries.begin(),
    // 		      it->secondaries.end());
    
    // bwa_hits_2.insert(bwa_hits_2.end(),
    // 		      (it+1)->secondaries.begin(),
    // 		      (it+1)->secondaries.end());
    
    // make all of the local breakpoints
    // OK so if we turn off secondary alignments, this will
    // just be one iteration through each loop (A and B)
    // but that's OK
    for (auto& a : bwa_hits_1) {
      for (auto& b : bwa_hits_2) {
	
	// initialize the breakpoint
	BreakPointPtr bp = std::make_shared<BreakPoint>(sc, &a, &b, this);

	// add the the vector of breakpoints
	m_local_breaks.push_back(bp);
      }
    }
  } // end frag iterator loop
  
  // if this is a double mapping, we are done
  if (m_frag_v.size() == 2) {
    return;
  }
  
  //////////////
  // 3+ mappings
  //////////////
  
  // go through alignments and find start and end that reach mapq
  // bstart and bend are indices of local alignments
  // that should contribute to a global
  // e.g. bstart = 0 = A
  int bstart = -1; //MAX_CONTIG_SIZE; // dummy value to be overwritten
  int bend = -1; 
  
  // loop and find contigs with strong support
  for (size_t i = 0; i < m_frag_v.size(); i++) {
    if (m_frag_v[i].m_align->MapQuality() >= 50) {
      // first good alignment fragment, that's the start
      if (bstart < 0)
	bstart = i;
      bend = i;
    }
  }
  
  // if nothing good found, then just pick first and last (A-C or whatever)
  if (bstart == bend || bstart < 0) { 
    bstart = 0;
    bend = m_frag_v.size() -1;
  }
  
  assert(bend <= m_frag_v.size() && bend > 0);
  assert(bstart <= m_frag_v.size());
  assert(bstart >= 0); 
  
  // if we had ABC but only AB or BC is good connection
  // then don't treat as complex
  if (bend - bstart == 1) {
    return; 
  }

  /////
  // there are 3+ quality mappings. Set a global break
  /////
  // set some values in the breakpoint for this contig
  m_global_bp = std::make_shared<BreakPoint>(sc,
					     &m_frag_v[bstart],
					     &m_frag_v[bend  ],
					     this);
  
  m_global_bp->svtype = SVType::TSI_GLOBAL; // complex=true;
  
  // now, each local breakpoint is a TSI local not a simple rearrangement
  for (auto& i : m_local_breaks) {
    i->svtype = SVType::TSI_LOCAL;
  }    
  
}

std::string AlignedContig::printDiscordantClusters(const SeqLib::BamHeader& h) const {
  
  std::stringstream out;
  if (m_dc.size() == 0)
    return "none";

  for (const auto& d : m_dc) {
    out << d.print(h) << " ";
  }
  return out.str();
  
}

// bool AlignedContig::checkLocal(const SeqLib::GenomicRegion& window)
// {
//   bool has_loc = false;
//   for (auto& i : m_frag_v) 
//     if (i.checkLocal(window))
//       has_loc = true;
  
//   // check for all of the breakpoints of non-indels (already handling indels)
//   for (auto& i : m_local_breaks) 
//     i.checkLocal(window);
//   // for (auto& i : m_global_bp_secondaries) 
//   //   i.checkLocal(window);
//   m_global_bp->checkLocal(window);
  
//   return has_loc;
  
// }


BreakPointPtrVector AlignedContig::getAllBreakPoints() const {
  
  BreakPointPtrVector out;
  
  // get the indel BreakPoints
  for (auto& i : m_frag_v) {
    out.insert(out.end(), i.m_indel_breaks.begin(), i.m_indel_breaks.end());
  }
  
  // get the global breakpoint
  if (m_global_bp) //->isEmpty())
    out.push_back(m_global_bp);
  
  // get the local breakpoints
  out.insert(out.end(), m_local_breaks.begin(), m_local_breaks.end());
  
  return out;
}

// std::vector<BreakPoint> AlignedContig::getAllBreakPointsSecondary() const {

//   std::vector<BreakPoint> out;

//   for (auto& i : m_global_bp_secondaries)
//     if (!i.isEmpty())
// 	out.push_back(i);

//   return out;
// }

  
bool AlignedContig::hasVariant() const { 
  
  if (m_global_bp) //!m_global_bp->isEmpty())
    return true;
  
  if (m_local_breaks.size())
    return true; 
  
  // if (m_global_bp_secondaries.size())
  //   return true; 
  
  // if (m_local_breaks_secondaries.size())
  //     return true; 

  for (const auto& i : m_frag_v)
    for (const auto& j : i.m_indel_breaks)
      if (j->hasMinimal()) 
	return true;
  
  return false;
  
}

// void AlignedContig::addDiscordantCluster(DiscordantClusterMap& dmap)
// {

//   // loop through the breaks and compare with the map
//   for (auto& i : m_local_breaks)
//     i.__combine_with_discordant_cluster(dmap);

//   if (!m_global_bp->isEmpty())
//     m_global_bp->__combine_with_discordant_cluster(dmap);
  
//   if (m_global_bp->hasDiscordant())
//     m_dc.push_back(m_global_bp->dc);
  
//   for (auto& i : m_global_bp_secondaries)
//     i.__combine_with_discordant_cluster(dmap);
// }

  // void AlignedContig::assessRepeats() {

  //   for (auto& a : m_frag_v)
  //     for (auto& i : a.m_indel_breaks)
  // 	i.repeatFilter();
    
  //   for (auto& b : m_local_breaks)
  //     b.repeatFilter();
  //   for (auto& b : m_local_breaks_secondaries)
  //     b.repeatFilter();
  //   for (auto& b : m_global_bp_secondaries)
  //     b.repeatFilter();

  // }


// void AlignedContig::refilterComplex() {
  
//   if (m_global_bp && m_global_bp->num_align <= 2)
//     return;
  
//   // initialize split counts vec
//   std::vector<int> scounts(m_local_breaks.size(), 0);
  
//   for (size_t i = 0; i < scounts.size() ; ++i) {
//     for (auto& j : m_local_breaks[i]->allele)
//       scounts[i] += j.second.split;
//   }

//   // every local BP needs to have 4+ split reads
//   // to be a quality TSI event
//   // if not, then forget the global / TSI altogether
//   for (auto& i : scounts) {
//     if (i < 4) {
//       m_global_bp = nullptr;
//       for (auto& i : m_local_breaks)
// 	i->svtype = SVType::ASSMB; 
//       return;
//     }
//   }
  
// }

std::string AlignedContig::getContigName() const {
  if (!m_frag_v.size()) 
    return "";  
    return m_frag_v[0].m_align->Qname(); 
  }

int AlignedContig::getMaxMapq() const { 
  int m = -1;
  for (auto& i : m_frag_v)
    if (i.m_align->MapQuality() > m)
      m = i.m_align->MapQuality();
  return m;
  
}

int AlignedContig::getMinMapq() const {
  int m = 1000;
  for (auto& i : m_frag_v)
    if (i.m_align->MapQuality() < m)
      m = i.m_align->MapQuality();
  return m;
}

// bool AlignedContig::hasLocal() const { 
//   for (auto& i : m_frag_v) 
//     if (i.local) 
//       return true; 
//   return false; 
// }

// void AlignedContig::writeAlignedReadsToBAM(SeqLib::BamWriter& bw) { 
//   for (const auto& i : m_bamreads)
//     bw.WriteRecord(*i);
// } 

// void AlignedContig::writeToBAM(SeqLib::BamWriter& bw) const { 
//   for (auto& i : m_frag_v) {
//     i.writeToBAM(bw);
//   }
// } 

// std::string AlignedContig::getSequence() const { 
//   assert(m_seq.length()); 
//   return m_seq; 
// }

void AlignedContig::AddAlignedRead(svabaReadPtr& br) {
  m_bamreads.push_back(br);
}
