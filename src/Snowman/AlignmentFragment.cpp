#include "AlignmentFragment.h"
#include "AlignedContig.h"

#define MIN_INDEL_MATCH_BRACKET 6
#define MAX_INDELS 10000
#define MAX_INDEL_PER_CONTIG 6 

// write the alignment record to a BAM file
void AlignmentFragment::writeToBAM(SeqLib::BamWriter& bw) const { 
  bw.WriteRecord(m_align); 
} 


bool AlignmentFragment::checkLocal(const SeqLib::GenomicRegion& window) {

    // make a region for this frag
    SeqLib::GenomicRegion gfrag(m_align.ChrID(), m_align.Position(), m_align.PositionEnd());
    
    if (window.GetOverlap(gfrag)) {
      local = true;
      return true;
    }
    
    return false;
  }

AlignmentFragment::AlignmentFragment(const SeqLib::BamRecord &talign, bool flip) {

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

}

void AlignmentFragment::SetIndels(const AlignedContig * c) {

    BreakPoint bp;
    size_t fail_safe_count = 0;
    while (parseIndelBreak(bp) && fail_safe_count++ < MAX_INDELS && !m_align.SecondaryFlag()) {
      bp.aligned_covered = c->aligned_covered;
      assert(bp.aligned_covered);
      assert(bp.valid());
      for (auto& i : c->prefixes)
	bp.allele[i].indel = true;
      m_indel_breaks.push_back(bp);
      assert(bp.num_align == 1);
    }
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


void AlignmentFragment::indelCigarMatches(const std::unordered_map<std::string, SeqLib::CigarMap>& cmap) {

    // loop through the indel breakpoints
    for (auto& i : m_indel_breaks) {
      
      assert(i.getSpan() > 0);

      // get the hash string in same formate as cigar map (eg. pos_3D)
      std::string st = i.getHashString();

      for (auto& c : cmap) {
	SeqLib::CigarMap::const_iterator ff = c.second.find(st);
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

