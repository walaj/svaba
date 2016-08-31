#include "AlignedContig.h"

#include <unordered_map>

#define MAX_CONTIG_SIZE 5000000
#define MIN_INDEL_MATCH_BRACKET 1
#define MAX_INDELS 10000
#define MAX_INDEL_PER_CONTIG 6 

using namespace SeqLib;

  static std::vector<std::string> repr = {"AAAAA", "TTTTT", "CCCCC", "GGGG", 
					  "TATATATA", "ATATATAT", 
					  "GCGCGCGC", "CGCGCGCG", 
					  "TGTGTGTG", "GTGTGTGT", 
					  "TCTCTCTC", "CTCTCTCT", 
					  "CACACACA", "ACACACAC", 
					  "GAGAGAGA", "AGAGAGAG"};
  

  AlignedContig::AlignedContig(const BamRecordVector& bav, const std::set<std::string>& pref) {
    
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
    if (bav.begin()->ReverseFlag()) {
      SeqLib::rcomplement(m_seq);
    }

    prefixes = pref;

    // zero the coverage
    for (auto& i : cov)
      i.second = std::vector<int>(m_seq.length(), 0);
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
	m_frag_v.push_back(AlignmentFragment(i, flip, pref));
	m_frag_v.back().num_align = num_align;
      } else {
	bool flip = (m_seq != i.Sequence()); // if the seq was flipped, need to flip the AlignmentFragment
	if (m_frag_v.size())
	  m_frag_v.back().secondaries.push_back(AlignmentFragment(i, flip, pref));
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
    for (auto& i : aligned_coverage)
      if (i)
	++aligned_covered;

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
    GRC grc;
    if (m_global_bp.b1.cpos > m_global_bp.b2.cpos) {
      grc.add(GenomicRegion(0, m_global_bp.b2.cpos-buff, m_global_bp.b1.cpos+buff)); // homology
    }
    else {
      grc.add(GenomicRegion(0, m_global_bp.b1.cpos-buff, m_global_bp.b2.cpos+buff)); // insertion	
    }    
    
    for (auto& i : m_local_breaks) {
      if (i.b1.cpos > i.b2.cpos) {
	grc.add(GenomicRegion(0, i.b2.cpos-buff, i.b1.cpos+buff)); // homology
      }
      else {
	grc.add(GenomicRegion(0, i.b1.cpos-buff, i.b2.cpos+buff)); // insertion	
      }
    }
    grc.CreateTreeMap();

    // check if 
    for (auto& i : m_frag_v) {
      BPVec new_indel_vec;
      for (auto& b : i.m_indel_breaks) {
	if (!grc.CountOverlaps(GenomicRegion(0, b.b1.cpos, b.b2.cpos)))
	  new_indel_vec.push_back(b);
      }
      i.m_indel_breaks = new_indel_vec;
    }
    
  }
  
  SeqLib::GenomicRegionVector AlignedContig::getAsGenomicRegionVector() const {
    SeqLib::GenomicRegionVector g;
    for (auto& i : m_frag_v)
      g.push_back(i.m_align.asGenomicRegion());
    return g;
  }

  void AlignedContig::printContigFasta(std::ofstream& os) const {
    os << ">" << getContigName() << std::endl;
    os << getSequence() << std::endl;
  }
  
  void AlignedContig::blacklist(GRC &grv) {
    
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
  
  bool AlignedContig::checkLocal(const GenomicRegion& window)
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
  
  bool AlignmentFragment::checkLocal(const GenomicRegion& window)
  {
    // make a region for this frag
    GenomicRegion gfrag(m_align.ChrID(), m_align.Position(), m_align.PositionEnd());
    
    if (window.GetOverlap(gfrag)) {
      local = true;
      return true;
    }
    
    return false;
  }

  AlignmentFragment::AlignmentFragment(const BamRecord &talign, bool flip, const std::set<std::string>& prefixes) {

    m_align = talign;

    sub_n = talign.GetIntTag("SQ");

    // orient cigar so it is on the contig orientation. 
    // need to do this to get the right ordering of the contig fragments below
    // We only flip if we flipped the sequence, and that was determined
    // by the convention set in AlignedContig::AlignedContig, so we need
    // to pass that information explicitly
    if (flip/*m_align.ReverseFlag()*/) {
      m_cigar = m_align.GetReverseCigar();
    } else { 
      m_cigar = m_align.GetCigar();
    }

    // find the start position of alignment ON CONTIG
    start = 0; 
    for (auto& i : /*m_align.GetCigar()*/ m_cigar) {
      if (i.Type() != 'M')
	start += i.Length();
      else
	break;
    }

    // set the left-right breaks
    unsigned currlen  = 0; 
    
    // CPOS is zero based

    // cigar is oriented to as is from aligner
    for (auto& i : m_cigar/*m_align.GetCigar()*/ /*align.CigarData*/) { //CigarOpVec::const_iterator j = align.cigar.begin(); j != align.cigar.end(); j++) {
      
      // SET THE CONTIG BREAK (treats deletions and leading S differently)
      // the first M gets the break1, pos on the left
      if (i.Type() == 'M' && break1 == -1)
	break1 = currlen;
      if (i.Type() != 'D') // m_skip deletions but not leading S, but otherwise update
	currlen += i.Length();
      if (i.Type() == 'M') // keeps triggering every M, with pos at the right
	break2 = currlen;
    }

    // assign the genomic coordinates of the break
    if (m_align.ReverseFlag()) {
      gbreak2 = m_align.Position() + 1;
      gbreak1 = m_align.PositionEnd();
    } else {
      gbreak1 = m_align.Position() + 1;
      gbreak2 = m_align.PositionEnd();
    }

    if (break1 >= MAX_CONTIG_SIZE || break2 >= MAX_CONTIG_SIZE || break1 < 0 || break2 < 0) 
      std::cerr << " break1 " << break1 << " break2 " << break2 << " " << (*this) << std::endl;
    assert(break1 < MAX_CONTIG_SIZE);
    assert(break2 < MAX_CONTIG_SIZE);

    if (break1 < 0 || break2 < 0) {
      std::cout << (*this) << std::endl;
    }
    assert(break1 >= 0);
    assert(break2 >= 0);

    // parse right away to see if there are indels on this alignment
    BreakPoint bp;
    size_t fail_safe_count = 0;
    while (parseIndelBreak(bp) && fail_safe_count++ < 100 && !m_align.SecondaryFlag()) {
      //std::cerr << bp << " " << (*this) << std::endl;
      assert(bp.valid());
      for (auto& i : prefixes)
	bp.allele[i].indel = true;
      m_indel_breaks.push_back(bp);
      assert(bp.num_align == 1);
    }

    assert(fail_safe_count != MAX_INDELS);

}

  std::ostream& operator<<(std::ostream &out, const AlignmentFragment &c) {
    
    // sets the direction to print
    char jsign = '>'; 
    if (c.m_align.ReverseFlag())
      jsign = '<';
    
    // print the cigar value per base
    for (auto& j : /*c.m_align.GetCigar()*/ c.m_cigar) { //c.align.CigarData) { // print releative to forward strand
      if (j.Type() == 'M')
	out << std::string(j.Length(), jsign);
      else if (j.Type() == 'I') 
	out << std::string(j.Length(), 'I');
      else if (j.Type() == 'S' || j.Type() == 'H')
	out << std::string(j.Length(), '.');
    }
    
    // print contig and genome breaks
    out << "\tC[" << c.break1 << "," << c.break2 << "] G[" << c.gbreak1 << "," << c.gbreak2 << "]";
    
    // add local info
    std::string chr_name = c.m_align.GetZTag("MC");
    if (!chr_name.length())
      chr_name = std::to_string(c.m_align.ChrID()+1);
    out << "\tLocal: " << c.local << "\tAligned to: " << chr_name << ":" << c.m_align.Position() << "(" << (c.m_align.ReverseFlag() ? "-" : "+") << ") CIG: " << c.m_align.CigarString() << " MAPQ: " << c.m_align.MapQuality() << " SUBN " << c.sub_n;
  
    return out;
  }
  
  void AlignedContig::checkAgainstCigarMatches(const std::unordered_map<std::string, CigarMap>& cmap) {

    for (auto& i : m_frag_v)
      i.indelCigarMatches(cmap);
    
  }

  
  void AlignmentFragment::indelCigarMatches(const std::unordered_map<std::string, CigarMap>& cmap) {

    // loop through the indel breakpoints
    for (auto& i : m_indel_breaks) {
      
      assert(i.getSpan() > 0);

      // get the hash string in same formate as cigar map (eg. pos_3D)
      std::string st = i.getHashString();

      for (auto& c : cmap) {
	CigarMap::const_iterator ff = c.second.find(st);
	// if it is, add it
	if (ff != c.second.end()) {
	  i.allele[c.first].cigar = ff->second;	  
	}
      }
    }      
    
  }

  // convention is that cpos and gpos for deletions refer to flanking REF sequence.
  // eg a deletion of 1 bp of base 66 will have gpos1 = 65 and gpos2 = 67
  bool AlignmentFragment::parseIndelBreak(BreakPoint &bp) {
    
     // make sure we have a non-zero cigar
    if (m_cigar.size() == 0) {
      std::cerr << "CIGAR of length 0 on " << *this << std::endl;
      return false;
    }

    // reject if too many mismatches
    if (di_count == 0) {
      for (auto& i : m_cigar)
	if (i.Type() == 'D' || i.Type() == 'I')
	  ++di_count;
    }
    if (di_count > MAX_INDEL_PER_CONTIG && m_align.Qname().substr(0,2) == "c_") // only trim for snowman assembled contigs
      return false;
    
     // reject if first alignment is I or D or start with too few M
    if (m_cigar.begin()->Type() == 'I' || m_cigar.begin()->Type() == 'D' || m_cigar.back().Type() == 'D' || m_cigar.back().Type() == 'I') {
      return false;
    }

    // use next available D / I by looping until can increment idx
    size_t loc = 0; // keep track of which cigar field
    for (auto& i : m_cigar) {
      ++loc;
      if ( (i.Type() == 'D' || i.Type() == 'I')) {
	assert (loc != 1 && loc != m_cigar.size()); // shuldn't start with I or D
	bool prev_match = (m_cigar[loc-2].Type() == 'M' && m_cigar[loc-2].Length() >= MIN_INDEL_MATCH_BRACKET);
	bool post_match = (m_cigar[loc].Type() == 'M' && m_cigar[loc].Length() >= MIN_INDEL_MATCH_BRACKET);
	if (loc > idx && prev_match && post_match) { // require 15M+ folowoing I/D
	  idx = loc;
	  break;
	}
      }
    }

    // we made it to the end, no more indels to report
    if (loc == m_cigar.size())
      return false;

    // clear out the old bp just in case
    bp = BreakPoint();
    bp.num_align = 1;
    
    int curr = 0;
    int gcurrlen = 0; 
    
    // make break ends with dummy positions
    bp.b1 = BreakEnd(m_align);
    bp.b2 = BreakEnd(m_align);

    // assign the contig-wide properties
    bp.cname = m_align.Qname();
    bp.seq = m_align.Sequence();
    assert(bp.cname.length());
    
    size_t count = 0; // count to make sure we are reporting the right indel
    for (auto& i : m_cigar) { // want breaks in CONTIG coordnates, so use oriented cigar
      ++count;
      
      // set the contig breakpoint
      if (i.Type() == 'M' || i.Type() == 'I' || i.Type() == 'S') 
	curr += i.Length();

      // update the left match side
      if (i.Type() == 'M' && count < idx)
	bp.left_match += i.Length();

      // update the right match side
      if (i.Type() == 'M' && count > idx)
	bp.right_match += i.Length();

      // if deletion, set the indel
      if (i.Type() == 'D' && bp.b1.cpos == -1 && count == idx) {
	
	bp.b1.cpos = curr-1;
	bp.b2.cpos = curr;
	
      } 
      
      // if insertion, set the indel
      if (i.Type() == 'I' && bp.b1.cpos == -1 && count == idx) {
	bp.b1.cpos = curr - i.Length() - 1; // -1 because cpos is last MATCH
	bp.b2.cpos = curr/* - 1*/; 
	bp.insertion = m_align.Sequence().substr(bp.b1.cpos+1, i.Length()); // +1 because cpos is last MATCH.
      }
      
      // set the genome breakpoint
      if (bp.b1.cpos > 0 && count == idx) {
	if (i.Type() == 'D') {
	  if (!m_align.ReverseFlag()) {
	    bp.b1.gr.pos1 =  m_align.Position() + gcurrlen; // dont count this one//bp.b1.cpos + align.Position; //gcurrlen + align.Position;
	    bp.b2.gr.pos1 = bp.b1.gr.pos1 + i.Length() + 1;
	  } else {
	    bp.b2.gr.pos1 =  (m_align.PositionEnd()) - gcurrlen + 1; // snowman81 removed a -1, set to +1  //bp.b1.cpos + align.Position; //gcurrlen + align.Position;
	    bp.b1.gr.pos1 =  bp.b2.gr.pos1 - i.Length()- 1; 

	  }
	} else if (i.Type() == 'I') {
	  if (!m_align.ReverseFlag()) {
	    bp.b1.gr.pos1 = m_align.Position() + gcurrlen; //gcurrlen + align.Position;
	    bp.b2.gr.pos1 = bp.b1.gr.pos1 + 1;	
	  } else {
	    // GetEndPosition is 1 too high
	    bp.b2.gr.pos1 = (m_align.PositionEnd()-1) - gcurrlen; //gcurrlen + align.Position;
	    bp.b1.gr.pos1 = bp.b2.gr.pos1 - 1;	
	  }
	}
	//break; // already got it, so quit cigar loop
      }
      
      // update the position on the genome
      if (i.Type() == 'M' || i.Type() == 'D') {
	gcurrlen += i.Length();
      } 
      
      
    } // end cigar loop
    
    // set the dummy other end
    bp.b1.gr.pos2 = bp.b1.gr.pos1; 
    bp.b2.gr.pos2 = bp.b2.gr.pos1;
   
    // should have been explicitly ordered in the creation above
    if (!(bp.b1.gr < bp.b2.gr)) {
      return false;
    }

    bp.b1.gr.strand = '+';
    bp.b2.gr.strand = '-';

    assert(bp.valid());
    return true;
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

  BreakEnd AlignmentFragment::makeBreakEnd(bool left) {
    
    BreakEnd b;

    b.sub_n = sub_n;
    b.chr_name = m_align.GetZTag("MC");
    assert(b.chr_name.length());
    b.mapq = m_align.MapQuality();
    b.matchlen = m_align.NumMatchBases();
    b.local = local;
    b.nm = std::max(m_align.GetIntTag("NM") - m_align.MaxDeletionBases(), (uint32_t)0);
    b.simple = m_align.GetIntTag("SZ");

    b.as_frac = (double)m_align.GetIntTag("AS") / (double) m_align.NumMatchBases();


    // if alignment is too short, zero the mapq
    // TODO get rid of this by integrating matchlen later
    if (b.matchlen < 30)
      b.mapq = 0;

    if (left) {
      b.gr = SeqLib::GenomicRegion(m_align.ChrID(), gbreak2, gbreak2);
      b.gr.strand = m_align.ReverseFlag() ? '-' : '+'; 
      b.cpos = break2; // take the right-most breakpoint as the first
    } else {
      b.gr = SeqLib::GenomicRegion(m_align.ChrID(), gbreak1, gbreak1);
      b.gr.strand = m_align.ReverseFlag() ? '+' : '-';
      b.cpos = break1;  // take the left-most of the next one
    }

    assert(b.cpos < MAX_CONTIG_SIZE);
    
    return b;
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
  
