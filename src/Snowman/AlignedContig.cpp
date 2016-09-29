#include "AlignedContig.h"
#include "PlottedRead.h"

AlignedContig::AlignedContig(const SeqLib::BamRecordVector& bav, const std::set<std::string>& pref) {
    
    if (!bav.size())
      return;

    // find the longest sequence, taking the first one. 
    // make sure sequence dir is set to same as first alignment
    for (auto& i : bav) {
      if (i.Sequence().length() > m_seq.length()) {
	if (i.ReverseFlag() == bav.begin()->ReverseFlag()) {
	  m_seq = i.Sequence();
	} else {
	  m_seq = i.Sequence();
	  SeqLib::rcomplement(m_seq);
	}
      }
    }

    // set the sequence. Convention is store as it came off assembler for first alignment
    if (bav.begin()->ReverseFlag()) 
      SeqLib::rcomplement(m_seq);

    prefixes = pref;

    for (size_t i = 0; i < m_seq.length(); ++i)
      aligned_coverage.push_back(0);
    
    // find the number of primary alignments
    size_t num_align = 0;
    for (auto& i : bav)
      if (!i.SecondaryFlag())
	++num_align;

    // make the individual alignments and add
    for (auto& i : bav) {
      if (!i.SecondaryFlag()) {
	bool flip = (m_seq != i.Sequence()); // if the seq was flipped, need to flip the AlignmentFragment
	m_frag_v.push_back(AlignmentFragment(i, flip));
	m_frag_v.back().num_align = num_align;
      } else {
	bool flip = (m_seq != i.Sequence()); // if the seq was flipped, need to flip the AlignmentFragment
	if (m_frag_v.size())
	  m_frag_v.back().secondaries.push_back(AlignmentFragment(i, flip));
      }      

      // set the aligned coverage
      SeqLib::Cigar cig = i.GetCigar();
      size_t p = 0;
      if (!i.SecondaryFlag())
      for (auto& c : cig) {
	for (size_t j = 0; j < c.Length(); ++j) {
	  if (c.Type() == 'M' || c.Type() == 'I')  // consumes contig and not clip
	    ++aligned_coverage[p];
	  if (c.ConsumesQuery() || c.Type() == 'H') // consumes contig, move iterator
	    ++p;
	}
      }
    }

    // find the total aligned coverage
    for (const auto& i : aligned_coverage)
      if (i) // 0 for no cov at this base, 1 is single cov, 2 double, etc
	++aligned_covered;
 
    // extract the indels on primary alignments
    for (auto& f : m_frag_v) 
      f.SetIndels(this);
   
    // sort fragments by order on fwd-strand contig
    if (m_frag_v.size() > 1) 
      std::sort(m_frag_v.begin(), m_frag_v.end());

    if (!m_frag_v.size())
      return;

    // get breaks out of it
    setMultiMapBreakPairs();

    // filter indels that land too close to a multi-map break
    filterIndelsAtMultiMapSites(5);
  }

  void AlignedContig::filterIndelsAtMultiMapSites(size_t buff) {
    
    if (m_frag_v.size() < 2)
      return;

    // make the ranges ON CONIG for the multimaps
    SeqLib::GRC grc;
    if (m_global_bp.b1.cpos > m_global_bp.b2.cpos) {
      grc.add(SeqLib::GenomicRegion(0, m_global_bp.b2.cpos-buff, m_global_bp.b1.cpos+buff)); // homology
    }
    else {
      grc.add(SeqLib::GenomicRegion(0, m_global_bp.b1.cpos-buff, m_global_bp.b2.cpos+buff)); // insertion	
    }    
    
    for (auto& i : m_local_breaks) {
      if (i.b1.cpos > i.b2.cpos) {
	grc.add(SeqLib::GenomicRegion(0, i.b2.cpos-buff, i.b1.cpos+buff)); // homology
      }
      else {
	grc.add(SeqLib::GenomicRegion(0, i.b1.cpos-buff, i.b2.cpos+buff)); // insertion	
      }
    }
    grc.CreateTreeMap();

    // check if 
    for (auto& i : m_frag_v) {
      BPVec new_indel_vec;
      for (auto& b : i.m_indel_breaks) {
	if (!grc.CountOverlaps(SeqLib::GenomicRegion(0, b.b1.cpos, b.b2.cpos)))
	  new_indel_vec.push_back(b);
      }
      i.m_indel_breaks = new_indel_vec;
    }
    
  }
  
  SeqLib::GenomicRegionVector AlignedContig::getAsGenomicRegionVector() const {
    SeqLib::GenomicRegionVector g;
    for (auto& i : m_frag_v)
      g.push_back(i.m_align.AsGenomicRegion());
    return g;
  }

  void AlignedContig::printContigFasta(std::ofstream& os) const {
    os << ">" << getContigName() << std::endl;
    os << getSequence() << std::endl;
  }
  
void AlignedContig::blacklist(SeqLib::GRC &grv) {
    
    // loop through the indel breaks and blacklist
    for (auto& i : m_frag_v) 
      for (auto& j : i.m_indel_breaks) 
	j.checkBlacklist(grv);
    
  }
  
  void AlignedContig::splitCoverage() { 
    
    for (auto& i : m_local_breaks_secondaries) 
      i.splitCoverage(m_bamreads);

    for (auto& i : m_global_bp_secondaries) 
      i.splitCoverage(m_bamreads);

    for (auto& i : m_frag_v) 
      for (auto& j : i.m_indel_breaks) 
      j.splitCoverage(m_bamreads);
    
    for (auto& i : m_local_breaks) 
      i.splitCoverage(m_bamreads);
    
    if (!m_global_bp.isEmpty()) 
      m_global_bp.splitCoverage(m_bamreads);
    
  }
  
  std::ostream& operator<<(std::ostream& out, const AlignedContig &ac) {
    
    // print the global breakpoint
    if (!ac.m_global_bp.isEmpty())
      out << "Global BP: " << ac.m_global_bp << 
	" ins_aginst_contig " << ac.insertion_against_contig_read_count << 
	" del_against_contig " << ac.deletion_against_contig_read_count << "  " << 
	ac.getContigName() << std::endl;       
    
    // print the global breakpoint for secondaries
    if (ac.m_global_bp_secondaries.size())
      out << "SECONDARY Global BP: " << ac.m_global_bp << 
	" ins_aginst_contig " << ac.insertion_against_contig_read_count << 
	" del_against_contig " << ac.deletion_against_contig_read_count << "  " << 
	ac.getContigName() << std::endl;       
    
    // print the multi-map breakpoints
    for (auto& i : ac.m_local_breaks)
      if (!i.isEmpty())
	out << "Multi-map BP: " << i << " -- " << ac.getContigName() << std::endl;       
    // print the multi-map breakpoints for secondary
    for (auto& i : ac.m_local_breaks_secondaries)
      if (!i.isEmpty())
	out << "SECONDARY Multi-map BP: " << i << " -- " << ac.getContigName() << std::endl;       
    
    // print the indel breakpoints
    for (auto& i : ac.m_frag_v)
      for (auto& j : i.getIndelBreaks()) 
	if (!j.isEmpty())
	  out << "Indel: " << j << " -- " << ac.getContigName() << " ins_a_contig " << ac.insertion_against_contig_read_count << 
	    " del_a_contig " << ac.deletion_against_contig_read_count << std::endl;       
    
    // print the AlignmentFragments alignments
    for (auto& i : ac.m_frag_v) 
      out << i << " Disc: " << ac.printDiscordantClusters() << " -- " << ac.getContigName() << std::endl;
    bool draw_divider = true;
    for (auto& i : ac.m_frag_v) {
      for (auto& j : i.secondaries) {
	if (draw_divider) {
	  out << std::string(ac.m_seq.length(), 'S') << std::endl;
	  draw_divider = false;
	}
	out << j << " Disc: " << ac.printDiscordantClusters() << " -- " << ac.getContigName() << std::endl;
      }
    }

    // print the break locations for indel deletions
    for (auto& i : ac.m_frag_v) {
      for (auto& j : i.getIndelBreaks()) {
	if (j.num_align == 1 && j.insertion == "") // deletion
	  //std::cerr << j.b1.cpos << " " << j.b2.cpos << " name " << ac.getContigName() << " ins " << j.insertion << std::endl;
	  out << std::string(j.b1.cpos, ' ') << "|" << std::string(j.b2.cpos-j.b1.cpos-1, ' ') << '|' << "   " << ac.getContigName() << std::endl;	
      }
    }
    
    out << ac.getSequence() << "    " << ac.getContigName() << std::endl; 
    PlottedReadVector plot_vec;
    
    // print out the individual reads
    for (auto& i : ac.m_bamreads) {
      
      int pos = -1;
      int aln = -1;
      int rc = 0;
      std::string this_cig;
      std::string seq = i.QualitySequence();
      std::string sr = i.GetZTag("SR");
      
      
      // get the more complex tags (since there can be multiple annotations per tag)
      std::vector<int> posvec = i.GetSmartIntTag("SL"); // start positions ON CONTIG
      std::vector<int> alnvec = i.GetSmartIntTag("TS"); // start positions ON READ
      std::vector<int> rcvec = i.GetSmartIntTag("RC"); // read reverse complmented relative to contig
      std::vector<std::string> cigvec = i.GetSmartStringTag("SC"); // read against contig CIGAR
      std::vector<std::string> cnvec = i.GetSmartStringTag("CN");
      
      if (posvec.size() != alnvec.size() ||
	  posvec.size() != rcvec.size() ||
	  cigvec.size() != posvec.size() ||
	  cnvec.size() != posvec.size())
	continue;
      
      assert(cnvec.size() == posvec.size());
      size_t kk = 0;
      for (; kk < cnvec.size(); kk++) 
	if (cnvec[kk] == ac.getContigName()) {
	  pos = posvec[kk];
	  aln = alnvec[kk];
	  rc = rcvec[kk];
	  this_cig = cigvec[kk];
	  break;
	}
      
      // reverse complement if need be
      if (rc)
	SeqLib::rcomplement(seq);      
      
      if (aln > 0)
	try {
	  seq = seq.substr(aln, seq.length() - aln);
	} catch (...) {
	  std::cerr << "AlignedContig::operator<< error: substring out of bounds. seqlen " << 
	    seq.length() << " start " << aln << " length " << (seq.length() - aln) << std::endl;
	}
      
      if ( (pos + seq.length() ) > ac.getSequence().length()) 
	try { 
	  seq = seq.substr(0, ac.getSequence().length() - pos);
	} catch (...) {
	  std::cerr << "AlignedContig::operator<< (2) error: substring out of bounds. seqlen " << 
	    seq.length() << " start " << 0 << " pos " << pos << " ac.getSequence().length() " << 
	    ac.getSequence().length() << std::endl;

	}
      
      
      assert(kk != cnvec.size()); // assure that we found something
      pos = abs(pos);
      int padlen = ac.getSequence().size() - pos - seq.size() + 5;
      padlen = std::max(5, padlen);
      
      std::stringstream rstream;
      assert(pos < MAX_CONTIG_SIZE && padlen < MAX_CONTIG_SIZE); // bug, need to check
      rstream << sr << "--" << (i.ChrID()+1) << ":" << i.Position() << " r2c CIGAR: " << this_cig;
      
      plot_vec.push_back({pos, seq, rstream.str()});
    }
    
    std::sort(plot_vec.begin(), plot_vec.end());
    
    PlottedReadLineVector line_vec;
    
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
	prl.contig_len = ac.getSequence().length();
	prl.addRead(&i);
	line_vec.push_back(prl);
      }
    }
    
    // plot the lines. Add contig identifier to each
    for (auto& i : line_vec) 
      out << i << " " << ac.getContigName() << std::endl;
    
    return out;
  }
  
  void AlignedContig::setMultiMapBreakPairs() {
     
    // if single mapped contig, nothing to do here
    if (m_frag_v.size() == 1)
      return;
    
    // if single mapped contig, nothing to do here  
    // initialize the breakpoint, fill with basic info
    BreakPoint bp;
    for (auto& i : prefixes) 
      bp.allele[i].indel = false;

    bp.seq = getSequence();
    bp.aligned_covered = aligned_covered;
    bp.num_align = m_frag_v.size();
    assert(bp.num_align > 0);
    
    bp.cname = getContigName(); 
    assert(bp.cname.length());
    
    bp.has_local_alignment = m_frag_v[0].m_align.GetIntTag("LA"); // has local alignment?

    // walk along the ordered contig list and make the breakpoint pairs  
    for (AlignmentFragmentVector::const_iterator it = m_frag_v.begin(); it != m_frag_v.end() - 1; it++) {
      
      AlignmentFragmentVector bwa_hits_1, bwa_hits_2;
      bwa_hits_1.push_back(*it);
      bwa_hits_2.push_back(*(it+1));
      bwa_hits_1.insert(bwa_hits_1.end(), it->secondaries.begin(), it->secondaries.end());
      bwa_hits_2.insert(bwa_hits_2.end(), (it+1)->secondaries.begin(), (it+1)->secondaries.end());
      
      // make all of the local breakpoints
      for (auto& a : bwa_hits_1) {
	for (auto& b : bwa_hits_2) {
	  
	  bp.b1 = a.makeBreakEnd(true);
	  bp.b2 = b.makeBreakEnd(false); 

	  // set the insertion / homology
	  bp.__set_homologies_insertions();
	  
	  // order the breakpoint
	  bp.order();

	  bp.secondary = a.m_align.SecondaryFlag() || b.m_align.SecondaryFlag();

	  assert(bp.valid());

	  // add the the vector of breakpoints
	  if (!bp.secondary) {
	    m_local_breaks.push_back(bp);
	  } else {
	    m_local_breaks_secondaries.push_back(bp);	  
	  }
	}
      }
    } // end frag iterator loop
    
    // if this is a double mapping, we are done
    if (m_frag_v.size() == 2) {
      m_global_bp = m_local_breaks[0];
      m_local_breaks.clear();
      for (auto& i : m_local_breaks_secondaries) {
	m_global_bp_secondaries.push_back(i);
      }
      m_local_breaks_secondaries.clear();
      return;
    }
    
    // 3+ mappings. If all good, then don't make the "global"
    // Actually, do make the global
    //bool make_locals = true;
    //for (size_t i = 1; i < m_frag_v.size() - 1; ++i)
    //  if (m_frag_v[i].m_align.MapQuality() < 50)
    //	make_locals = false;
    
    //if (make_locals) { // intermediates are good, so just leave locals as-is
    //  return;
    //}
    
    // TODO support 3+ mappings that contain secondary
    
    // go through alignments and find start and end that reach mapq 
    size_t bstart = MAX_CONTIG_SIZE; //1000 is a dummy
    size_t bend = m_frag_v.size() - 1;

    for (size_t i = 0; i < m_frag_v.size(); i++)
      if (m_frag_v[i].m_align.MapQuality() >= 60) {
	bend = i;
	if (bstart == MAX_CONTIG_SIZE)
	  bstart = i;
      }
    if (bstart == bend || bstart==MAX_CONTIG_SIZE) {
      bstart = 0;
      bend = m_frag_v.size() -1 ;
    }

    assert(bend <= m_frag_v.size());
    assert(bstart <= m_frag_v.size());
    assert(bstart != MAX_CONTIG_SIZE);
    
    // there are 3+ mappings, and middle is not great. Set a global break
    m_global_bp = bp;
    m_global_bp.b1 = m_frag_v[bstart].makeBreakEnd(true);
    m_global_bp.b2 = m_frag_v[bend].makeBreakEnd(false);
    m_global_bp.complex=true;
    for (auto& i : m_local_breaks) {
      i.complex = true;
      i.complex_local = true; // this is a sub-piece
    }
    
    // set the homologies
    m_global_bp.__set_homologies_insertions();
    
    // order the breakpoint
    m_global_bp.order();
    
    assert(m_global_bp.valid());
  }
  
  std::string AlignedContig::printDiscordantClusters() const {
    
    std::stringstream out;
    if (m_dc.size() == 0)
      return "none";
    
    for (std::vector<DiscordantCluster>::const_iterator it = m_dc.begin(); it != m_dc.end(); it++)
      out << *it << " ";
    return out.str();
    
  }
  
bool AlignedContig::checkLocal(const SeqLib::GenomicRegion& window)
  {
    bool has_loc = false;
    for (auto& i : m_frag_v) 
      if (i.checkLocal(window))
	has_loc = true;

    // check for all of the breakpoints of non-indels (already handling indels)
    for (auto& i : m_local_breaks) 
      i.checkLocal(window);
    for (auto& i : m_global_bp_secondaries) 
      i.checkLocal(window);
    m_global_bp.checkLocal(window);

    return has_loc;
    
  }
  
  
void AlignedContig::checkAgainstCigarMatches(const std::unordered_map<std::string, SeqLib::CigarMap>& cmap) {

    for (auto& i : m_frag_v)
      i.indelCigarMatches(cmap);
    
  }

  
  
  std::vector<BreakPoint> AlignedContig::getAllBreakPoints(bool local_restrict) const {
    
    std::vector<BreakPoint> out;
    for (auto& i : m_frag_v) {
      if (i.local || !local_restrict) // only output if alignment is local for indels
	for (auto& k : i.m_indel_breaks)
	  out.push_back(k);
    }
 
    if (!m_global_bp.isEmpty())
      out.push_back(m_global_bp);
    
    out.insert(out.end(), m_local_breaks.begin(), m_local_breaks.end());
    out.insert(out.end(), m_global_bp_secondaries.begin(), m_global_bp_secondaries.end());
    
    return out;
  }
  
  std::vector<BreakPoint> AlignedContig::getAllBreakPointsSecondary() const {
    
    std::vector<BreakPoint> out;
    
    for (auto& i : m_global_bp_secondaries)
      if (!i.isEmpty())
	out.push_back(i);
    
    return out;
  }
  
  
  bool AlignedContig::hasVariant() const { 
    
    if (!m_global_bp.isEmpty())
      return true;

    if (m_local_breaks.size())
      return true; 

    if (m_global_bp_secondaries.size())
      return true; 

    if (m_local_breaks_secondaries.size())
      return true; 

    for (auto& i : m_frag_v)
      if (i.local && i.m_indel_breaks.size())
	return true;
    
    return false;
    
  }
  
  void AlignedContig::addDiscordantCluster(DiscordantClusterMap& dmap)
  {
    
    // loop through the breaks and compare with the map
    for (auto& i : m_local_breaks)
      i.__combine_with_discordant_cluster(dmap);

    if (!m_global_bp.isEmpty())
      m_global_bp.__combine_with_discordant_cluster(dmap);
    
    if (m_global_bp.hasDiscordant())
      m_dc.push_back(m_global_bp.dc);
    
    for (auto& i : m_global_bp_secondaries)
      i.__combine_with_discordant_cluster(dmap);
  }
  
  void AlignedContig::assessRepeats() {

    for (auto& a : m_frag_v)
      for (auto& i : a.m_indel_breaks)
	i.repeatFilter();
    
    for (auto& b : m_local_breaks)
      b.repeatFilter();
    for (auto& b : m_local_breaks_secondaries)
      b.repeatFilter();
    for (auto& b : m_global_bp_secondaries)
      b.repeatFilter();

  }


  void AlignedContig::refilterComplex() {

    if (m_global_bp.num_align <= 2)
      return;

    // initialize split counts vec
    std::vector<int> scounts(m_local_breaks.size(), 0);

    for (size_t i = 0; i < scounts.size() ; ++i) {
      for (auto& j : m_local_breaks[i].allele)
	scounts[i] += j.second.split;
    }

    bool bad = false;
    for (auto& i : scounts)
      if (i < 4)
	bad = true;

    if (bad) 
      m_global_bp = BreakPoint();

  }
  

std::string AlignedContig::getContigName() const { 
    if (!m_frag_v.size()) 
      return "";  
    return m_frag_v[0].m_align.Qname(); 
  }

int AlignedContig::getMaxMapq() const { 
  int m = -1;
  for (auto& i : m_frag_v)
    if (i.m_align.MapQuality() > m)
      m = i.m_align.MapQuality();
  return m;
  
}

int AlignedContig::getMinMapq() const {
  int m = 1000;
  for (auto& i : m_frag_v)
    if (i.m_align.MapQuality() < m)
      m = i.m_align.MapQuality();
  return m;
}

bool AlignedContig::hasLocal() const { 
  for (auto& i : m_frag_v) 
    if (i.local) 
      return true; 
  return false; 
}

void AlignedContig::writeAlignedReadsToBAM(SeqLib::BamWriter& bw) { 
  for (auto& i : m_bamreads)
    bw.WriteRecord(i);
} 


void AlignedContig::writeToBAM(SeqLib::BamWriter& bw) const { 
  for (auto& i : m_frag_v) {
    i.writeToBAM(bw);
  }
} 

std::string AlignedContig::getSequence() const { 
  assert(m_seq.length()); 
  return m_seq; 
}

void AlignedContig::AddAlignedRead(const SeqLib::BamRecord& br) {
  m_bamreads.push_back(br);
}
