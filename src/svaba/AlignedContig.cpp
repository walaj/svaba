#include "AlignedContig.h"
#include "PlottedRead.h"
#include "SvabaUtils.h"
#include "ContigAlignmentScore.h"

#include <cctype>
#include <iomanip>
#include <sstream>
#include <unordered_set>

// SvABA2.0: rendering toggle for alignments.txt.gz read display.
//
//   true  -> one read per output line (easier to grep, easier to inspect a
//            specific read; what we want when debugging variant-support
//            decisions). Each read still gets the same per-read info token
//            with r2c/native scores etc, so the file stays diff-able.
//   false -> the original space-packed mode where multiple non-overlapping
//            reads share a line (compact view, fewer lines for big contigs).
//
// Flip and rebuild to switch between modes. The packing code path (when
// false) is preserved verbatim so we can switch back without code changes.
namespace {
constexpr bool kAlignmentsOneReadPerLine = true;
}  // namespace

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

std::string AlignedContig::printToAlignmentsFile(const SeqLib::BamHeader& h) const {
  
  std::stringstream out;
  std::string cname = this->getContigName();
  
  // print the global complex breakpoint
  if (m_global_bp)
    out << "Global BP: " << m_global_bp->printSimple(h) << 
      //      " ins_against_contig " << insertion_against_contig_read_count << 
      //" del_against_contig " << deletion_against_contig_read_count << "  " << 
      cname << std::endl;       
  
  // print the global breakpoint for secondaries
  // if (m_global_bp_secondaries.size())
  //   out << "SECONDARY Global BP: " << m_global_bp->print(h) << 
  //     " ins_against_contig " << insertion_against_contig_read_count << 
  //     " del_against_contig " << deletion_against_contig_read_count << "  " << 
  //     cname << std::endl;       
  
  // print the rearrangement breakpoints
  for (auto& i : m_local_breaks) {
    if (i) {
      out << "Multi-map BP: " << i->printSimple(h) <<
	" -- " << cname << std::endl;
    }
  }
  
  // print the indel breakpoints
  for (auto& i : m_frag_v) {
    for (auto& j : i.m_indel_breaks) {
      if (j)
	out << "Indel: " << j->printSimple(h) << " -- " << cname << "\n";
	  // " ins_a_contig " << insertion_against_contig_read_count << 
	  // " del_a_contig " << deletion_against_contig_read_count <<
    }
  }
    
  // print the AlignmentFragments alignments
  for (auto& i : m_frag_v) 
    out << i.printToAlignmentsFile() << " Disc: " <<
      this->printDiscordantClusters(h) << " -- " << cname <<
      "\n";
  
  // print the secondary alignments of the contigs
  // bool draw_divider = true;
  // for (auto& i : m_frag_v) {
  //   for (auto& j : i.secondaries) {
  //     if (draw_divider) {
  // 	out << std::string(m_seq.length(), 'S') << std::endl;
  // 	draw_divider = false;
  //     }
  //     out << j.print() << " Disc: " << printDiscordantClusters(h) << " -- " << getContigName() << std::endl;
  //   }
  // }
  
  // print the break locations for indel deletions
  for (auto& i : m_frag_v) {
    for (auto& j : i.m_indel_breaks) {
      if (j->isIndel() && j->insertion.empty()) // deletion
	out << j->printDeletionMarksForAlignmentsFile() << "\n";
    }
  }
  
  out << m_seq << "    " << cname << std::endl;
  PlottedReadVector plot_vec;

  // SvABA2.0: build the set of split-supporting and discordant-supporting
  // read UniqueNames so we can tag each rendered read line with its
  // support kind (split / disc / both / none). Split-supporting reads are
  // those that landed in any BreakPoint::SampleInfo::supporting_reads on
  // any BP for this contig (populated by BreakPoint::splitCoverage).
  // Discordant-supporting reads are those keyed in any associated
  // DiscordantCluster::reads or ::mates (keyed by UniqueName, see
  // DiscordantCluster.cpp:134).
  std::unordered_set<std::string> split_supporters;
  for (const auto& bp : getAllBreakPoints()) {
    if (!bp) continue;
    for (const auto& un : bp->getAllSupportingReads())
      split_supporters.insert(un);
  }
  std::unordered_set<std::string> disc_supporters;
  for (const auto& dc : m_dc) {
    for (const auto& kv : dc.reads)  disc_supporters.insert(kv.first);
    for (const auto& kv : dc.mates)  disc_supporters.insert(kv.first);
  }

  // print out the individual reads
  for (const auto& i : m_bamreads) {

    std::string seq = i->Sequence();
    std::string sr = i->UniqueName();

    // get the read to contig alignment information
    r2c this_r2c = i->GetR2C(getContigName());

    int pos = this_r2c.start_on_contig;

    if (this_r2c.rc)
      SeqLib::rcomplement(seq);

    // SvABA2.0: build the rendered string with three regions:
    //   leading-S bases  -> lowercase, rendered immediately to the LEFT of
    //                       the matched portion (so they appear at columns
    //                       [start_on_contig - leading_S, start_on_contig)).
    //   M / D bases      -> uppercase / '-' as before.
    //   trailing-S bases -> lowercase, rendered to the RIGHT of the match.
    // I (insertion wrt contig) is dropped since it has no contig column.
    // If the leading clip would render before column 0 we truncate it so
    // PlottedReadLine's `pos >= last_loc` invariant is preserved.
    size_t lead_clip_len  = 0;
    size_t trail_clip_len = 0;
    if (this_r2c.cig.size() > 0) {
      const auto& fop = this_r2c.cig.front();
      const auto& lop = this_r2c.cig.back();
      if (fop.Type() == 'S') lead_clip_len  = fop.Length();
      // protect against a CIGAR that is a single S op (front == back)
      if (lop.Type() == 'S' && this_r2c.cig.size() > 1) trail_clip_len = lop.Length();
    }

    std::string lead_clip;
    std::string trail_clip;
    if (lead_clip_len > 0 && lead_clip_len <= seq.length()) {
      lead_clip = seq.substr(0, lead_clip_len);
      for (auto& ch : lead_clip) ch = static_cast<char>(std::tolower(ch));
    }
    if (trail_clip_len > 0 && trail_clip_len <= seq.length()) {
      trail_clip = seq.substr(seq.length() - trail_clip_len);
      for (auto& ch : trail_clip) ch = static_cast<char>(std::tolower(ch));
    }

    // edit the string to reflect gapped and clipped alignments
    // that is, only show the match portion, and put "-" on reads for deletions wrt contig
    size_t p = 0; // move along on sequence, starting at first non-clipped base
    std::string gapped_seq;
    for (SeqLib::Cigar::const_iterator c = this_r2c.cig.begin(); c != this_r2c.cig.end(); ++c) {
      if (c->Type() == 'M') { //
	assert(p + c->Length() <= seq.length());
	gapped_seq += seq.substr(p, c->Length());
      } else if (c->Type() == 'D') {
	gapped_seq += std::string(c->Length(), '-');
      }

      if (c->ConsumesQuery()) //c->Type() == 'I' || c->Type() == 'M' || c->Type() == 'S') // consumes query
	p += c->Length();
    }
    seq = lead_clip + gapped_seq + trail_clip;

    pos = abs(pos);

    // shift the render anchor left by the leading-clip length so the
    // lowercase clip bases appear in the correct contig columns.
    pos -= static_cast<int>(lead_clip_len);
    if (pos < 0) {
      // leading clip would extend past the contig left edge; chop the
      // overhanging lowercase bases so we don't violate the renderer's
      // monotonic-position invariant.
      const size_t to_drop = static_cast<size_t>(-pos);
      if (to_drop >= seq.size()) seq.clear();
      else seq = seq.substr(to_drop);
      pos = 0;
    }

    int padlen = m_seq.size() - pos - seq.size() + 5;
    padlen = std::max(5, padlen);

    // SvABA2.0: variant-support kind for this read on this contig.
    // The viewer parses this token with /\bsupport: (split|disc|both)\b/
    // and uses it to drive the split-only / disc-only / both toggles.
    const bool is_split = split_supporters.count(sr) > 0;
    const bool is_disc  = disc_supporters.count(sr)  > 0;
    const char* support_kind = is_split && is_disc ? "both"
                             : is_split            ? "split"
                             : is_disc             ? "disc"
                             :                       "none";

    // SvABA2.0: per-read SV-focused alignment scores. These are the exact
    // scores BreakPoint::splitCoverage uses for its "r2c better than
    // native" gate (see svaba::readAlignmentScore in
    // ContigAlignmentScore.h). Surfacing them here makes it possible to
    // diagnose, by inspection of alignments.txt.gz alone, why a particular
    // read was or wasn't credited as a variant supporter.
    //
    //   r2cScore    = score(r2c CIGAR, r2c NM)               always present
    //                                                        when r2c.cig
    //                                                        is non-empty
    //   nativeScore = score(BAM CIGAR, BAM NM tag)           "NA" if either
    //                                                        is missing
    //
    // The gate is: variant-supporter iff r2cScore > nativeScore.
    auto fmt_score = [](double s) {
      std::ostringstream os;
      os << std::fixed << std::setprecision(1) << s;
      return os.str();
    };

    std::string r2c_score_str = "NA";
    if (this_r2c.cig.size() > 0) {
      r2c_score_str = fmt_score(
          svaba::readAlignmentScore(this_r2c.cig, this_r2c.nm));
    }

    // Match BreakPoint::splitCoverage's source-of-truth precedence: prefer
    // the corrected-read realignment (apples-to-apples with r2c), fall
    // back to the original BAM record only when corrected_native is empty.
    // Without this preference the dump would mislead — the gate decides
    // on one number and the .txt would print a different one.
    std::string native_score_str = "NA";
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

    std::stringstream rstream;
    rstream << sr << "--" << (i->ChrID()+1) << ":" << i->Position()
	    << " r2c POS: " << this_r2c.start_on_contig
	    << " FLAG: " << (this_r2c.rc ? 16 : 0)
	    << " NM: " << this_r2c.nm
	    << " support: " << support_kind
	    << " r2cScore: "    << r2c_score_str
	    << " nativeScore: " << native_score_str
	    << " CIGAR: " << this_r2c.cig;

    plot_vec.push_back({pos, seq, rstream.str()});
  }
  
  std::sort(plot_vec.begin(), plot_vec.end());

  PlottedReadLineVector line_vec;

  // SvABA2.0: when kAlignmentsOneReadPerLine is true, skip the "fits in an
  // existing line?" test and always start a new line per read. The packed
  // mode (else branch) is the original behavior, kept verbatim so flipping
  // the toggle restores it byte-for-byte.
  if (kAlignmentsOneReadPerLine) {
    for (auto& i : plot_vec) {
      PlottedReadLine prl;
      prl.contig_len = m_seq.length();
      prl.addRead(&i);
      line_vec.push_back(prl);
    }
  } else {
    // plot the reads from the ReadPlot vector
    for (auto& i : plot_vec) {
      bool found = false;
      for (auto& j : line_vec) {
        if (j.readFits(i)) { // it fits here
          j.addRead(&i);
          found = true;
          break;
        }
      }
      if (!found) { // didn't fit anywhere, so make a new line
        PlottedReadLine prl;
        prl.contig_len = m_seq.length();
        prl.addRead(&i);
        line_vec.push_back(prl);
      }
    }
  }

  // plot the lines. Add contig identifier to each
  for (auto& i : line_vec)
    out << i << " " << cname << std::endl;

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
