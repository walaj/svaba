#include "AlignmentFragment.h"
#include "AlignedContig.h"

#define MIN_INDEL_MATCH_BRACKET 0
#define MAX_INDELS 10000
#define MAX_INDEL_PER_CONTIG 6 

// write the alignment record to a BAM file
// void AlignmentFragment::writeToBAM(SeqLib::BamWriter& bw) const { 
//   bw.WriteRecord(*m_align); 
// } 

using SeqLib::CigarField;

bool AlignmentFragment::operator<(const AlignmentFragment& str) const {
  return (start < str.start);
}


// bool AlignmentFragment::checkLocal() const {

//   assert(region_.pos1 >= 0);
  
//   // make a region for this frag
//   SeqLib::GenomicRegion gfrag(m_align->ChrID(),
// 			      m_align->Position(),
// 			      m_align->PositionEnd());
  
//   if (region_.GetOverlap(gfrag)) {
//     return true;
//   }
  
//   return false;
// }

AlignmentFragment::AlignmentFragment(BamRecordPtr &talign,
				     bool flip,
				     const GenomicRegion& local_region,
				     const SvabaSharedConfig* sc_) :
  region_(local_region), m_align(talign), sc(sc_) {
  
  // orient cigar so it is on the contig orientation. 
  // need to do this to get the right ordering of the contig fragments below
  // We only flip if we flipped the sequence, and that was determined
  // by the convention set in AlignedContig::AlignedContig, so we need
  // to pass that information explicitly
  // if (flip/*m_align->ReverseFlag()*/) {
  //   m_cigar = m_align->GetReverseCigar();
  // } else { 
  //   m_cigar = m_align->GetCigar();
  // }

  flipped = flip;
  
  // set the left-right breaks
  unsigned currlen  = 0; 
  
  // NB: CPOS is zero based
  // NB: for break1 and break2, we use these later in the
  //     for either printing to alignments file, or for
  //     getting the *contig* coordinate breakpoints
  //     So this to account for "flip" convention across
  //     AlignmentFragments
  const auto& cig = flipped ? m_align->GetReverseCigar() : m_align->GetCigar();
  for (auto& i : cig) { 
    
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
  if (m_align->ReverseFlag()) {
    gbreak2 = m_align->Position() + 1;
    gbreak1 = m_align->PositionEnd();
  } else {
    gbreak1 = m_align->Position() + 1;
    gbreak2 = m_align->PositionEnd();
  }
  
  // shouldn't hit either of the next two conditions
  // if (break1 >= MAX_CONTIG_SIZE || break2 >= MAX_CONTIG_SIZE || break1 < 0 || break2 < 0) {
  //   throw std::runtime_error(
  // 			     "Invalid breakpoints: break1 = " + std::to_string(break1) +
  // 			     ", break2 = " + std::to_string(break2) +
  // 			     ", context = " + this->print()
  // 			     );
  // }
  
  if (break1 < 0 || break2 < 0) {
    throw std::runtime_error("Negative breakpoint detected in AlignmentFragment "); //  + printToAlignmentsFile());
  }
  
  // find the start position of alignment ON CONTIG
  start = 0; 
  for (auto& i : cig) { //m_align->GetCigar()) { //m_cigar) {
    if (i.Type() != 'M')
      start += i.Length();
    else
      break;
  }
  
}

void AlignmentFragment::SetIndels() {
  
  if (m_align->SecondaryFlag())
    return; // ignore these

  const auto& cig = m_align->GetCigar();
  // loops the alignment fragment to extract indels
  // start at 1 and end at len-1 because we don't consider indels if they
  //   start the cigar string e.g. 1D5M, not reliable
  for (size_t i = 1; i < cig.size() - 1; ++i) { 
    const CigarField& c = cig[i];
    const CigarField& pre = cig[i-1];
    const CigarField& post = cig[i+1];
    if ( (c.Type() == 'D' || c.Type() == 'I')) {
      bool prev_match = pre.Type() == 'M'; /* && cig[loc-2].Length() >= MIN_INDEL_MATCH_BRACKET*/
      bool post_match = post.Type() == 'M';/* && cig[loc].Length() >= MIN_INDEL_MATCH_BRACKET*/
      
      // only consider indels if 
      if (prev_match && post_match) {
	// convention is that cpos and gpos for deletions refer to flanking REF sequence.
	// eg a deletion of 1 bp of base 66 will have gpos1 = 65 and gpos2 = 67
	BreakPointPtr bp = std::make_shared<BreakPoint>(this, i, sc);	
	assert(bp);
	// add the indel    
	m_indel_breaks.push_back(bp);
      }
    }
  }
}

std::string AlignmentFragment::printToAlignmentsFile() const {
  
  std::stringstream out;
  // sets the direction to print
  char jsign = '>'; 
  if (m_align->ReverseFlag())
    jsign = '<';

  // NB: for printing, this should actually finally
  //   consider if this fragment is "flipped" relative
  //   to the BWA alignment, since we want to sync
  //   orientations across all of the fragments.
  const auto& cig = flipped ? m_align->GetReverseCigar() : m_align->GetCigar(); 
  
  // print the cigar value per base
  for (auto& j : cig) {
    if (j.Type() == 'M')
      out << std::string(j.Length(), jsign);
    else if (j.Type() == 'I') 
      out << std::string(j.Length(), 'I');
    else if (j.Type() == 'S' || j.Type() == 'H')
      out << std::string(j.Length(), '.');
  }
  
  // print contig and genome breaks
  out << "\tC[" << break1 << "," << break2 << "] G[" << gbreak1 << "," << gbreak2 << "]";

  // print it
  std::string chr_name = m_align->ChrName(sc->header);
  assert(chr_name.length());
  out << "\tAligned to: " <<
    chr_name << ":" << m_align->Position() << "(" <<
    (m_align->ReverseFlag() ? "-" : "+") << ") CIG: " <<
    m_align->CigarString() << " MAPQ: " <<
    m_align->MapQuality(); // << " SUBN " << sub_n;
  
  return out.str();
}

// void AlignmentFragment::fillRearrangementBreakEnd(bool left,
// 					 BreakEnd& b) {
    
//     b.transferContigAlignmentData(m_align);

//     if (left) {
//       b.gr = SeqLib::GenomicRegion(m_align->ChrID(), gbreak2, gbreak2);
//       b.gr.strand = m_align->ReverseFlag() ? '-' : '+'; 
//       b.cpos = break2; // take the right-most breakpoint of the left  as the first
//     } else {
//       b.gr = SeqLib::GenomicRegion(m_align->ChrID(), gbreak1, gbreak1);
//       b.gr.strand = m_align->ReverseFlag() ? '+' : '-';
//       b.cpos = break1;  // take the left-most of the next one
//     }

//     assert(b.cpos < MAX_CONTIG_SIZE);

//   }
//       assert(i.getSpan() > 0);

//       // get the hash string in same formate as cigar map (eg. pos_3D)
//       std::string st = i.getHashString();

//       for (auto& c : cmap) {
// 	SeqLib::CigarMap::const_iterator ff = c.second.find(st);
// 	// if it is, add it
// 	if (ff != c.second.end()) {
// 	  i.allele[c.first].cigar = ff->second;	  
// 	}
//       }
//     }      
//   }
