#include "BreakPoint.h"

#include <getopt.h>
#include <iomanip>
#include <cassert>
#include <tuple>
#include <optional>

#include "gzstream.h"
#include "svabaUtils.h"
#include "STCoverage.h"
#include "svabaOptions.h"

#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/BamHeader.h"
#include "AlignmentFragment.h"
#include "AlignedContig.h"
#include "svabaModels.h"

// n is the max integer given the int size (e.g. 255). x is string with int
#define INTNSTOI(x,n) std::min((int)n, std::stoi(x));

static inline std::string to_string(SVType t) {
  switch(t) {
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

static inline std::string to_string(SomaticState s) {
  switch(s) {
  case SomaticState::NOTSET:       return "NOTSET";
  case SomaticState::NORMAL_LOD:   return "0";
  case SomaticState::SOMATIC_LOD:  return "1";
  case SomaticState::FAILED:       return "0";
  }
  return "UNKNOWN_SOMATICSTATE";
}

std::ostream& operator<<(std::ostream& out, const BreakPoint::SampleInfo& a) {
  out << " split: " << a.split << " cigar " << a.cigar << " alt " << a.alt << " cov " << a.cov << " disc " << a.disc;
  return out;
}

GenomicRegion BreakPoint::BreakEndAsGenomicRegionLeft() const {
  return b1.gr;
}

GenomicRegion BreakPoint::BreakEndAsGenomicRegionRight() const {
  return b2.gr;
}

bool BreakPoint::isIndel() const {
  if (svtype == SVType::NOTSET)
    throw std::runtime_error("BreakPoint without an svtype");

  return svtype == SVType::INDEL;
}

std::string BreakPoint::printDeletionMarksForAlignmentsFile() const {
  if (!isIndel() || !insertion.empty())
    return std::string();
  
  std::stringstream out;
  out << std::string(b1.cpos, ' ') <<
    "|" << std::string(b2.cpos - b1.cpos-1, ' ') <<
    "|   " << cname;

  return out.str();
  
}


HashVector BreakPoint::getBreakEndHashes() {
  HashVector vec;

  assert(!confidence.empty());
  assert(svtype != SVType::NOTSET);
  
  if (svtype != SVType::INDEL)
    return vec;

  if (b1.gr.pos1 > 0) {
    {
      vec.push_back(b1.hash(0));
      vec.push_back(b1.hash(1));
      vec.push_back(b1.hash(-1));
    }
  }

  return vec;
}

BreakPoint::SampleInfo operator+(const BreakPoint::SampleInfo& a1, const BreakPoint::SampleInfo& a2) {

  BreakPoint::SampleInfo a;
  
  a.disc = a1.disc + a2.disc;
  a.split = a1.split + a2.split;
  a.cigar = a1.cigar + a2.cigar;
  a.cigar_near = a1.cigar_near + a2.cigar_near;
  a.cov = a1.cov + a2.cov;
  
  // add the reads
  for (auto& i : a1.supporting_reads)
    a.supporting_reads.insert(i);
  for (auto& i : a2.supporting_reads)
    a.supporting_reads.insert(i);
  
  //a.alt = std::max((int)a.supporting_reads.size(), a.cigar);
  
  if (a.supporting_reads.size()) // we have the read names, so do that (non-refilter run)
    a.UpdateAltCounts();
  else // no read names stored, so just get directly
    a.alt = a1.alt + a2.alt;
  
  return a;
}


static constexpr double scale_factor = 10.0;
static std::unordered_map<int, double> ERROR_RATES =
  {{0, scale_factor * 1e-4}, {1, scale_factor * 1e-4},
   {2,  scale_factor * 1e-4}, {3,  scale_factor * 1e-4},
   {4,  scale_factor * 1e-4}, {5, scale_factor * 2e-4},
   {6, scale_factor * 5e-4}, {7, scale_factor * 1e-3},
   {8, scale_factor * 2e-3}, {9, scale_factor * 3e-3},
   {10, scale_factor * 1e-2}, {11, scale_factor * 2e-2},
   {12, scale_factor * 3e-2}};

double __myround(double x) { return std:: floor(x * 10) / 10; }

BreakPoint::BreakPoint(const SvabaSharedConfig* _sc,
		       const AlignmentFragment* left,
		       const AlignmentFragment* right,
		       const AlignedContig* alc
		       ) : sc(_sc) {

  // this is at least of type ASSMB (can upgrade later to ASDIS or TSI)
  svtype = SVType::ASSMB;

  // instantiate the SampleInfo objects for each BAM
  for (auto const& p : sc->prefixes) {
    // if p isn't already present, constructs SampleInfo(this) in-place
    allele.try_emplace(p);
  }

  assert(alc);

  // transfer AlignedContig level information
  seq = alc->m_seq;
  num_align = alc->m_frag_v.size();
  cname = alc->getContigName();

  // set the local alignment
  // int la = 0;
  // alc->m_frag_v[0].m_align->GetIntTag("LA", la);
  // if (la > 0)
  //   local = LocalAlignment::NONVAR_LOCAL_REALIGNMENT;
  
  // set the break coordinates 
  bool isleft = true;
  b1.transferContigAlignmentData(left,    isleft);
  b2.transferContigAlignmentData(right, !isleft);

  // set the insertion / homology 
  set_homologies_insertions();

  // set seconary flag
  secondary = left->m_align->SecondaryFlag() ||
    right->m_align->SecondaryFlag();
  
}

// make the file string
std::string BreakPoint::toFileString(const BamHeader& header) const {
  
  // make sure we already ran scoring
  assert(confidence.length());
  assert(svtype != SVType::NOTSET);
  
  std::string sep = "\t";
  std::stringstream ss;
  
  double max_lod = 0;
  for (const auto&  [_,al] : allele) 
    max_lod = std::max(max_lod, al.LO);
  
  std::string evidence = to_string(svtype);

  std::string somatic_string = "NA";
  std::string somatic_lod_string = "NA";
  for (const auto& [pref,_] : allele) {
    if (pref.at(0) == 'n') { // need to have one normal to call somatic
      somatic_string = to_string(somatic);
      somatic_lod_string = std::to_string(std::min((double)99,LO_s));
    }
  }
  

  //NB: we are adding +1 to the positions here to be consistent
  // with SAM/VCF which is 1-indexed, while we are using 0-indexed
  // positions internally to be consistent with htslib
  ss << b1.gr.ChrName(header) << sep
     << (b1.gr.pos1 +1) << sep
     << b1.gr.strand << sep //1-3
     << b2.gr.ChrName(header) << sep
     << (b2.gr.pos1 +1) << sep
     << b2.gr.strand << sep //4-6
     << ref << sep << alt << sep //7-8
     << getSpan() << sep //9 
     << a.split << sep
     << a.alt << sep << a.cov << sep
     << a.cigar << sep << a.cigar_near << sep    
     << dc.mapq1 << sep << dc.mapq2 << sep
     << dc.ncount << sep << dc.tcount << sep
     << b1.mapq << sep << b2.mapq << sep 
     << b1.nm << sep << b2.nm << sep
     << b1.as << sep << b2.as << sep
     << b1.sub << sep << b2.sub << sep      
     << (homology.length() ? homology : "x") << sep 
     << (insertion.length() ? insertion : "x") << sep
     << (repeat_seq.length() ? repeat_seq : "x") << sep
     << cname << sep
     << num_align << sep 
     << confidence << sep
     << evidence << sep
     << quality << sep
     << secondary << sep
     << somatic_string << sep
     << somatic_lod_string << sep 
     << max_lod << sep 
    //<< pon << sep
    // << (repeat_seq.length() ? repeat_seq : "x") << sep 
     << (rs.length() ? rs : "x");
  
  for (const auto& [_,al] : allele)
    ss << sep << al.toFileString(svtype);
  
  return ss.str();
  
}

// make the file string
std::string BreakPoint::printSimple(const SeqLib::BamHeader& h) const {
  
  std::stringstream out;

  // make sure we set everything forthis
  assert(svtype != SVType::NOTSET);
  assert(!cname.empty());
  assert(b1.gr.chr >= 0);
  assert(svtype == SVType::DSCRD || b1.nm >= 0);
  assert(svtype == SVType::DSCRD || b2.nm >= 0);  
  assert(svtype == SVType::DSCRD || b1.as >= 0);
  assert(svtype == SVType::DSCRD || b2.as >= 0);  
    
  if (svtype == SVType::INDEL) {
    out << ">" << (insertion.size() ? "INS: " : "DEL: ") << getSpan() << " " << 
      b1.gr.ToString(h) << " " << cname;
    for (const auto& [pref, al] : allele)
      out << " " << pref << ":" << al.split;  
  } else {
    out << ": " << b1.gr.PointString(h) << " to " << b2.gr.PointString(h) << " SPAN " << getSpan() << " " << cname;
    for (const auto& [pref, al] : allele)
      out << " " << pref << ":" << al.split;  
  }
  
  return out.str();
  
}
  

/*std::ostream& operator<<(std::ostream& out, const BreakPoint& b) {
  
  if (b.isindel) {
	out << ">" << (b.insertion.size() ? "INS: " : "DEL: ") << b.getSpan() << " " << 
	  b.b1.gr << " " << b.cname << " " << b.evidence;
	//<< " T/N split: " << b.t.split << "/" << b.n.split << " T/N cigar: " 
          //  << b.t.cigar << "/" << b.n.cigar << " T/N Cov " << b.t.cov << "/" << b.n.cov << " DBSNP: " << rs_t;
	for (auto& i : b.allele)
	  out << " " << i.first << ":" << i.second.split;  
      } else {
	out << ": " << b.b1.gr.PointString() << " to " << b.b2.gr.PointString() << " SPAN " << b.getSpan() << " " << b.cname << " " << b.evidence;
	//<< " T/N split: " << b.t.split << "/" << b.n.split << " T/N disc: " 
	  //  << b.dc.tcount << "/" << b.dc.ncount << " " << b.evidence;
	for (auto& i : b.allele)
	  out << " " << i.first << ":" << i.second.split;  
      }
      
      return out;
      


      }*/

void BreakPoint::splitCoverage(svabaReadPtrVector& bav) {

  // dummy left-most and right-most positions that have an
  // overlapping read2contig alignment
  split_cov_bounds = {std::numeric_limits<int>::max(), -1};
  
  // track if first and second mate covers same split. fishy and remove them both
  std::unordered_map<std::string, bool> qname_and_num;
  
  // keep track of which reads already added
  std::unordered_set<std::string> qnames;
  
  // keep track of reads to reject
  std::set<std::string> reject_qnames;
  
  // keep track of which SR tags are valid splits
  std::unordered_set<std::string> valid_reads;
  
  // get the homology length. useful bc if read alignment ends in homologous region, it is not split
  int homlen = b1.cpos - b2.cpos;
  if (homlen < 0)
    homlen = 0;
  
  // loop all of the read to contig alignments for this contig
  for (const auto& j : bav) {
    
    r2c this_r2c = j->GetR2C(cname);
    
    bool read_should_be_skipped = false;
    if (num_align == 1) {
      
      std::vector<int> del_breaks;
      std::vector<int> ins_breaks;
      int pos = 0;
      
      // if this is a nasty repeat, don't trust non-perfect alignmentx on r2c alignment
      // if (checkHomopolymer(j->Sequence())) 
      // 	read_should_be_skipped = true;
      
      // loop through r2c cigar and see positions
      for(auto& i : this_r2c.cig) {
	
	if (i.Type() == 'D') 
	  del_breaks.push_back(pos);
	else if (i.Type() == 'I')
	  ins_breaks.push_back(pos);
	
	// update position on contig
	if (i.ConsumesReference())
	  pos += i.Length(); // update position on contig
	
      }
      
      size_t buff = std::max((size_t)3, repeat_seq.length() + 3);
      for (auto& i : del_breaks)
	if (i > b1.cpos - buff || i < b1.cpos + buff) // if start of insertion is at start of a del of r2c
	  read_should_be_skipped = true;
      for (auto& i : ins_breaks)
	if (i > b1.cpos - buff || i < b1.cpos + buff) // if start of insertion is at start of a del of r2c
	  read_should_be_skipped = true;
    } 
    
    if (read_should_be_skipped)  // default is r2c does not support var, so don't amend this_r2c
      continue;
    
    // get read ID
    std::string sample_id = j->Prefix(); //substr(0,4); // maybe just make this prefix
    std::string sr = j->UniqueName(); 
    
    // need read to cover past variant by some buffer. If there is a repeat,
    // then this needs to be even longer to avoid ambiguity
    int this_tbuff = T_SPLIT_BUFF + repeat_seq.length();
    int this_nbuff = N_SPLIT_BUFF + repeat_seq.length();      
    
    int rightbreak1 = b1.cpos + (j->Tumor() ? this_tbuff : this_nbuff); // read must extend this far right of break1
    int leftbreak1  = b1.cpos - (j->Tumor() ? this_tbuff : this_nbuff); // read must extend this far left of break1
    int rightbreak2 = b2.cpos + (j->Tumor() ? this_tbuff : this_nbuff);
    int leftbreak2  = b2.cpos - (j->Tumor() ? this_tbuff : this_nbuff);
    
    std::string contig_qname; // for sanity checking
    // get the alignment position on contig
    int pos = this_r2c.start_on_contig;
    int te  = this_r2c.end_on_contig;
    
    int rightend = te; 
    int leftend  = pos;
    bool issplit1 = (leftend <= leftbreak1) && (rightend >= rightbreak1);
    bool issplit2 = (leftend <= leftbreak2) && (rightend >= rightbreak2);
    
    bool both_split = issplit1 && issplit2;
    bool one_split = issplit1 || issplit2;
    
    
    // be more permissive for NORMAL, so keep out FPs
    bool valid  = (both_split && (j->Tumor() || homlen > 0)) ||
      (one_split && !j->Tumor() && homlen == 0) || (one_split && insertion.length() >= INSERT_SIZE_TOO_BIG_SPAN_READS);
    
    // requiring both break ends to be split for homlen > 0 is for situation beow
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>A..........................
    // ............................B>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // if we set break2 at pos B, then reads that end at A will count as splits on B,
    // when they should count as non-splits on A. For any read, if it is non-split on one 
    // segment, it should be non-split on all, for overlapping alignments like above. Not true of 
    // insertions at junctions, where one can split at one and not the other because of the intervening sequence buffer
    
    // check that deletion (in read to contig coords) doesn't cover break point
    size_t p = pos; // move along on contig, starting at first non-clipped base
    for (SeqLib::Cigar::const_iterator c = this_r2c.cig.begin(); c != this_r2c.cig.end(); ++c) {
      if (c->Type() == 'D') { 
	if ( (p >= leftbreak1 && p <= rightbreak1) || (p >= leftbreak2 && p <= rightbreak2) )
	  read_should_be_skipped = true;
      }
      if (c->ConsumesReference()) // if it moves it along the contig
	p += c->Length();
    }
    
    // add the split reads for each end of the break
    // a read is split if it is spans both break ends for tumor, one break end for normal (to
    // be more sensitive to germline) and if it spans both ends for deletion (should be next to 
    // each other), or one end for insertions larger than 10, or this is a complex breakpoint
    
    if (valid) { 
      std::string qn = j->Qname();
      
      
      // if read seen and other read was other mate, then reject
      if (num_align > 1 && qname_and_num.count(qn) && qname_and_num[qn] != j->FirstFlag()) {
	
	// need to reject all reads of this qname
	// because we saw both first and second in pair hit same split
	
	// 2023 - turning this off -- why? Beacuse for very short
	// insert size distributions, this can still be possible
	//reject_qnames.insert(qn);
	
      } else {
	
	// if haven't seen this read, add here. If have, then dont because want to avoid re-setting first/second convention
	// e.g. if we see a split read with first designation, we want to reject all reads with same qname if at any time
	// we see one with second mate designation. If we don't have this conditional, we can get fooled if the order is 
	// 1, 2, 1
	if (!qname_and_num.count(qn)) 
	  qname_and_num[qn] = j->FirstFlag();
	
	// this is a valid read
	this_r2c.supports_var = true;
	valid_reads.insert(sr);
	
	// how much of the contig do these span
	// for a given read QNAME, get the coverage that 
	split_cov_bounds.first = std::min(split_cov_bounds.first, pos);
	split_cov_bounds.second = std::max(split_cov_bounds.second, te);
      }
    }
    
    // update the counters for each break end
    if (issplit1 && valid)
      ++b1.split[sample_id];
    if (issplit2 && valid)
      ++b2.split[sample_id];	
    
    // add r2c back, but as amended
    //j.AddR2C(cname, this_r2c);
    
  } // end read loop
  
    // process valid reads
  for (auto& i : bav) {
    
    r2c this_r2c = i->GetR2C(cname);
    if (valid_reads.count(i->UniqueName())) {
      
      std::string qn = i->Qname();
      if (qnames.count(qn))
	continue; // don't count support if already added and not a short event
      // check that it's not a bad 1, 2 split
      if (reject_qnames.count(qn)) {
	this_r2c.supports_var = false;
	i->AddR2C(cname, this_r2c); // update that this actually does not support
	continue; 
      }
      
      reads.push_back(i);
      
      // keep track of qnames of split reads
      qnames.insert(qn);
      auto prefix = i->Prefix();
      if (auto it = allele.find(prefix); it != allele.end()) {
	auto& sample = it->second;
	sample.supporting_reads.insert(i->UniqueName());
	++sample.split;
      } else {
	throw std::runtime_error("SampleInfo not instantiated with this bam prefix: " + prefix);
      }      
    }
  }
  
  // adjust the alt count
  for (auto& [_,al] : allele) {
    al.UpdateAltCounts();
  }
  
}

void BreakPoint::checkBlacklist(GRC &grv) {
  if (grv.CountOverlaps(b1.gr) || grv.CountOverlaps(b2.gr)) 
    confidence = "BLACKLIST";
}

void BreakPoint::set_homologies_insertions() {
  try { 
    if (b1.cpos > b2.cpos)
      homology = seq.substr(b2.cpos, b1.cpos-b2.cpos);
    else if (b2.cpos > b1.cpos)
      insertion = seq.substr(b1.cpos, b2.cpos-b1.cpos);
    //if (insertion.length() == 0)
    //	;//insertion = "x";
    //if (homology.length() == 0)
    //;//homology = "x";
  } catch (...) {
    std::cerr << "cname: " << cname << " b1.cpos " << b1.cpos << " b2.cpos " << b2.cpos << " seq.length " << seq.length() << std::endl;
    std::cerr << "Caught error with contig on global-getBreakPairs: " << cname << std::endl;
    std::cerr << b1.cpos << " " << b2.cpos << " seq.length() " << seq.length() << " num_align " << num_align << std::endl;
  }
}

// isleft - this is the left-sided AlignmentFragment, so take the right break
//          as the breakpoint (and vice versa)
void BreakEnd::transferContigAlignmentData(const AlignmentFragment* f,
					   bool isleft) {

  const BamRecordPtr &r = f->m_align;
  
  // get number of suboptimal alignments
  r->GetIntTag("XS", sub);  
  // sub could be 0
  
  // get number of !SNV! mismatches
  // this is a svaba internal trick, since we are interested in
  // indels, so just track number of point mutation differences
  r->GetIntTag("NM", nm);
  assert(nm >= 0);
  nm = std::max(nm
              - static_cast<int>(r->MaxInsertionBases())
              - static_cast<int>(r->MaxDeletionBases()),
              0);

  // get alignment score
  r->GetIntTag("AS", as);
  assert(as >= 0);

  // transfer mapq and match length
  mapq     = r->MapQuality();
  matchlen = r->NumMatchBases();

  // set the coordinates as right-most (for an isleft=true fragment)
  // or left-most (for an isleft = false fragment)
  if (isleft) {
    gr = GenomicRegion(r->ChrID(), f->gbreak2, f->gbreak2);
    gr.strand = r->ReverseFlag() ? '-' : '+';     
    cpos = f->break2; 
  } else {
    gr = GenomicRegion(r->ChrID(), f->gbreak1, f->gbreak1);
    gr.strand = r->ReverseFlag() ? '+' : '-';
    cpos = f->break1;
  }
}
 
// construct an indel BreakPoint
BreakPoint::BreakPoint(const AlignmentFragment* f,
		       const int idx, 
		       const SvabaSharedConfig* _sc) : sc(_sc)
{

  svtype = SVType::INDEL;
  
  // instantiate the SampleInfo objects for each BAM
  for (auto const& p : sc->prefixes) {
    allele.try_emplace(p);
  }
  
  b1.transferContigAlignmentData(f, true);
  b2.transferContigAlignmentData(f, true);

  // re-zero the positions, not set yet
  b1.gr.pos1 = b1.gr.pos2 = -1;
  b2.gr.pos1 = b2.gr.pos2 = -1;
  
  // this alignment fragemnt underlying alignment
  const BamRecordPtr &r = f->m_align;

  // assign the contig-wide properties
  cname = r->Qname();
  seq = r->Sequence();
  assert(cname.length());
  assert(seq.length());
  
  num_align = 1;
  
  b1.gr.strand = '+'; // always the case for indels
  b2.gr.strand = '-';
  
  // set seconary flag
  secondary = r->SecondaryFlag();

  int curr = 0;     // current position on the contig (starting at zero)
  int gcurrlen = 0; // current position on the genome (starting at zero)

  //debugprint
  // if (cname == "c_1_49001_74001_14C")
  //   std::cerr << "...constructing... " << idx <<
  //     " m_align " << *r << " flipped " << f->flipped << std::endl;

  //NB: curr is the counter for number of based traversed on the contig
  // this is actually 1-based in how we use it, since e.g. if we have
  // a cigar like 1M... then right away curr is 1 and we are looking at the
  // first based of the contig. A bit confusing since C is 0-based, but
  // this is effectively how this is used below. 1-based is good since
  // this is what we want the output to be for VCF / SAM spec

  // NB: contig position and genome position internally are 0-based
  // similar to HTSlib
  size_t count = 0; // count to make sure we are reporting the right indel
  // NB: don't worry about "flip"" stuff from AlignedContig,
  // that's only important for comparing across different alignments/
  // Since we are parsing indels from a single alignment, which is always
  // referenced to the forward strand, just go with the raw BamRecord data
  for (auto& i : f->m_align->GetCigar()) { 

    // update the left match side
    if (i.Type() == 'M' && count < idx)
      left_match += i.Length();
    
    // update the right match side
    if (i.Type() == 'M' && count > idx)
      right_match += i.Length();

    // std::cerr << "LM " << left_match << " RM " << right_match << " count " <<
    //   count << " idx " << idx << std::endl; //debug
    
    //////
    // set the breakpoint positions on the contig
    //////
    
    // for either insertion or deletion, first base is the one before the
    // actual variant
    // NB: recall that curr is the number of bases thus far
    //     that have consumed the contig/query e.g. for 140M5D, curr will be 140 since
    //     the 5D does not consume the query/contig. So in a 0-based system like
    //     we use here internally to be consistent with htslib and C-rules,
    //     we add a -1 so that we count the last base of the contig before the deletion
    //     (e.g. the 140th base, or position 139 in 0-based system) as the first
    //     base on the contig for the deletion (or insertion).
    if (count == idx && (i.Type() == 'D' || i.Type() == 'I')) {
      b1.cpos = curr - 1; // 0-based position of first base before the deletion/insertion
    }
    
    // if deletion
    if (i.Type() == 'D' && count == idx) {
      b2.cpos = curr;     // 0-based position of first base after the deletion
    } 
    
    // if insertion
    if (i.Type() == 'I' && count == idx) {
      b2.cpos = b1.cpos + i.Length();
      // here, we add the +1 to start so that the actual insertion sequence is
      // at the first base of the insertion. For VCF reporting in the ALT
      // column, I'll add the leading reference base later
      // e.g. for a 2 bp insertion, insertion becomes e.g. AG while
      // in the ALT column of the VCF/bps it is GAG (if G is last reference match)
      insertion = f->m_align->Sequence().substr(b1.cpos + 1, i.Length()); 
    }
    
    // update the contig position traversal -- consumes query
    // M,I,S
    if (i.ConsumesQuery()) 
      curr += i.Length();
    
    //////
    // set the breakpoint positions on the genome
    //////
    
    // set the genome breakpoint

    if (b1.cpos >= 0 && count == idx) {

      if (i.Type() != 'I' && i.Type() != 'D')
	throw std::runtime_error("BreakPoint indel parsing - unexpected to get non del/ins CIGAR field");
      
      // here again the -1 is to get htslib position 0-based to SAM format (1-based), since e.g. for
      //        140M5D, we want the genome position of the start based to be the first
      //        matching base before the deletion, so again have to -1
      b1.gr.pos1 = r->Position() + gcurrlen - 1;

      // for deletions, it's the opposite to contig, so deletion consumes the genome
      if (i.Type() == 'D')
	b2.gr.pos1 = b1.gr.pos1 + i.Length();
      
      if (i.Type() == 'I') 
	b2.gr.pos1 = b1.gr.pos1 + 1;
    }
    
    // update the position on the genome
    // D,M
    if (i.ConsumesReference()) { 
      gcurrlen += i.Length();
    }
    
    // if (cname == "c_1_122501_141530_11C") {
    //   std::cerr << " count " << count << "    curr " << curr << " gcurr " << gcurrlen <<
    // 	" C " << i << " insertion " << insertion << " cname " << cname << std::endl;
    //   std::cerr << f->m_seq << std::endl;
    //   std::cerr << "FIPPED " << f->flipped << std::endl;
    //   std::cerr << " b1.gr " << b1.gr << " b2.gr " << b2.gr << std::endl;
    // }

    ++count;
    
  } // end cigar loop
  
    // set the dummy other end
  b1.gr.pos2 = b1.gr.pos1; 
  b2.gr.pos2 = b2.gr.pos1;
  
  // should have been explicitly ordered in the creation above
  if (!(b1.gr < b2.gr)) {
    std::cerr << "B1 " << b1.gr << " - " << b2.gr << " " << cname << " b1.cpos " << b1.cpos << " b2.cpos " << b2.cpos << std::endl;
    std::cerr << f->printToAlignmentsFile() << std::endl;
    throw std::runtime_error("invalid BreakPoint creation in indels in AlignmentFragment");
  }


  //// add the repeat sequence
  const std::string& seq = f->m_align->Sequence();
  const int REPBUFF = 3;
  for (const auto& [start, end, rep] : svabaUtils::find_long_homopolymers(seq)) {
    bool crosses_breakpoint = start < (b1.cpos + REPBUFF) && (b2.cpos - REPBUFF) < end;
    if (rep.length() > repeat_seq.length() && crosses_breakpoint) 
      repeat_seq = rep;
  }
  
  for (const auto& [start, end, rep] : svabaUtils::find_long_dinuc_repeats(seq)) {
    bool crosses_breakpoint = start < (b1.cpos + REPBUFF) && (b2.cpos - REPBUFF) < end;
    if (rep.length() > repeat_seq.length() && crosses_breakpoint)
      repeat_seq = rep;
  }
  ///////////
}

void BreakPoint::CombineWithDiscordantClusterMap(DiscordantClusterMap& dmap)
{

  // don't operate on indels
  if (svtype == SVType::INDEL)
    return;

  // since discordant clusters are ordered with the
  // m_reg1 being more left on the genome, we want to
  // match that convention for the breakpoint. However,
  // we can't actually change the bp1 and bp2 ordering
  // since this is set by the alignment orderings on the contig
  // not on the genome. So can make a dummy swap just for purposes of
  // finding the matching discorantcluster
  // on the GENOME
  const int PAD = 100;
  GenomicRegion bp1 = (b1.gr < b2.gr) ? b1.gr : b2.gr;
  GenomicRegion bp2 = (b1.gr < b2.gr) ? b2.gr : b1.gr;
  bp1.Pad(PAD);
  bp2.Pad(PAD);
  assert(! (bp1 > bp2) );
    
  for (auto& [_,d] : dmap) {

    // checks basic things about discordant cluster
    if (!d.valid())
      continue;


    assert(d.m_reg1 < d.m_reg2);
      
    // do either breakends overlap
    bool bp1reg1 = bp1.GetOverlap(d.m_reg1) > 0;
    bool bp2reg2 = bp2.GetOverlap(d.m_reg2) > 0;
    
    bool s1 = bp1.strand == d.m_reg1.strand;
    bool s2 = bp2.strand == d.m_reg2.strand;
    
    // get the edge of the cluster
    // if its a forward alignment discordant read, then get "right" most
    // if its a reverse alignment discordant read, then get "left" most
    int pos1 = d.m_reg1.strand == '+' ? d.m_reg1.pos2 : d.m_reg1.pos1;
    int pos2 = d.m_reg2.strand == '+' ? d.m_reg2.pos2 : d.m_reg2.pos1;

    // both sides overlap and both orientions agree
    bool pass = bp1reg1 && bp2reg2 && s1 && s2;

    // check that the ends are not way off
    // if (std::abs(pos1 - b1.gr.pos1) > PAD ||
    // 	std::abs(pos2 - b2.gr.pos1) > PAD)
    //   pass = false;
    
    //std::cerr << " cname " << cname << " pad " << PAD << " pass " << pass << " DC pos1 " << pos1 << " DC pos2 " << pos2 << 
    //  " s1 " << s1 << " s2 " << s2 << " bp1reg1 " << bp1reg1 << " bp2reg2 " << bp2reg2 << " diff1 " << std::abs(pos1 - b1.gr.pos1) <<
    //  " diff2 " << std::abs(pos2 - b2.gr.pos1) << " bp pos1 " << b1.gr.pos1 << " bp pos2 " << b2.gr.pos1 << std::endl;

    if (!pass)
      return;
    
    // check that we haven't already added a cluster to this breakpoint
    // if so, chose the one with more tumor support
    if (dc.isEmpty() || dc.tcount < d.tcount) {
      
      dc = d;
      d.m_contig = cname;
      
      // add the read counts from DiscordantCluster to this BreakPoint
      // the DiscordantCluster "counts" are created at DiscordantCluster constructor
      assert(!d.counts.empty());      
      for (auto& c : d.counts) {
	allele.at(c.first).disc = c.second;
      }

      // add the discordant reads names to supporting reads for each sampleinfo
      // reads
      for (auto& [_, read] : d.reads) {
	allele.at(read->Prefix())
          .supporting_reads
          .insert(read->UniqueName());
      }
      
      // mates (and fix the typo in your assert)
      for (auto& [_, mate] : d.mates) {
	allele.at(mate->Prefix())
          .supporting_reads
          .insert(mate->UniqueName());
      }      

      // update the variant support read counts per bam
      // (i.e. per SampleInfo object in our BreakPoint "alleles" structure)
      for (auto& [_,al] : allele)
       	al.UpdateAltCounts();
    } 
  }

  svtype = SVType::ASDIS;
}

// void BreakPoint::set_evidence() {
  
//   // if we are in refilter, then this is already set
//   if (!evidence.empty())
//     return;
  
//   bool isdisc = (dc.tcount + dc.ncount) != 0;
  
//   if (num_align == 1)
//     evidence = "INDEL";
//   else if ( isdisc && svtype == SVType::REARRANGEMENT) //!complex && num_align > 0)
//     evidence = "ASDIS";
//   else if ( isdisc && num_align < 3)
//     evidence = "DSCRD";
//   else if (svtype == SVType::REARRANGEMENT)
//     evidence = "ASSMB";
//   else if (svtype == SVType::TSI_GLOBAL) ///complex && !complex_local) // is A-C of an ABC
//     evidence = "TSI_G";
//   else if (svtype == SVType::TSI_LOCAL)
//     //else if (complex && complex_local) // is AB or BC of an ABC 
//     evidence = "TSI_L";
  
//   assert(evidence.length());
  
// }


void BreakPoint::score_assembly_only() {

  assert(local != LocalAlignment::NOTSET);
  assert(secondary != -1);
  
  int span = getSpan();
  int num_split = t.split + n.split;
  int cov_span = split_cov_bounds.second - split_cov_bounds.first ;

  // check for high repeats
  // bool hi_rep = false;
  // for (auto& rr : hirepr)
  //   if (seq.find(rr) != std::string::npos)
  //     hi_rep = true;

  float as_frac1  = static_cast<float>(b1.as)  / static_cast<float>(b1.matchlen);
  float sub_frac1 = static_cast<float>(b1.sub) / static_cast<float>(b1.matchlen);  
  float as_frac2  = static_cast<float>(b2.as)  / static_cast<float>(b2.matchlen);
  float sub_frac2 = static_cast<float>(b2.sub) / static_cast<float>(b2.matchlen);  
    
  if (local == LocalAlignment::FROM_DISTANT_REGION && svtype != SVType::TSI_LOCAL)  // added this back in v71
    // issue is that if a read is secondary aligned, it could be 
    // aligned to way off region. Saw cases where this happend in tumor
    // and not normal, so false-called germline event as somatic.
    confidence = "NOLOCAL";
  else if (local == LocalAlignment::NONVAR_LOCAL_REALIGNMENT)
    confidence = "LOCALMATCH";
  else if ( num_split > 1 && ( (cov_span <= (sc->readlen + 5 ) && cov_span > 0) || cov_span < 0) )
    confidence = "DUPREADS"; // the same sequences keep covering the split
  else if (homology.length() >= 20 && (span > 1500 || span == -1) && std::max(b1.mapq, b2.mapq) < 60)
    confidence = "NODISC";
  else if ((int)seq.length() < sc->readlen + 30)
    confidence = "TOOSHORT";
  else if (a.split < 7 && (span > 1500 || span == -1))  // large and inter chrom need 7+
    confidence = "NODISC";
  else if (std::max(b1.mapq, b2.mapq) <= 40 || std::min(b1.mapq, b2.mapq) <= 10) 
    confidence = "LOWMAPQ";
  else if ( std::min(b1.mapq, b2.mapq) <= 30 && a.split <= 8 ) 
    confidence = "LOWMAPQ";
  else if (std::max(b1.nm, b2.nm) >= 10 || std::min(as_frac1, as_frac2) < 0.8) 
    confidence = "LOWAS";
  else if ( (std::max(b1.nm, b2.nm) >= 3 || std::min(as_frac1, as_frac2) < 0.85) && getSpan() < 0 )
    confidence = "LOWAS";      
  //else if ((double)aligned_covered / (double)seq.length() < 0.80) // less than 80% of read is covered by some alignment
  //  confidence = "LOWAS";        
  else if ( (b1.matchlen < 50 && b1.mapq < 60) || (b2.matchlen < 50 && b2.mapq < 60) )
    confidence = "LOWMAPQ";
  else if ( std::min(b1.nm, b2.nm) >= 10)
    confidence = "LOWMAPQ";
  else if (a.split <= 3 && span <= 1500 && span != -1) // small with little split
    confidence = "LOWSPLITSMALL";
  else if (b1.gr.chr != b2.gr.chr && std::min(b1.matchlen, b2.matchlen) < 60) // inter-chr, but no disc reads, weird alignment
    confidence = "LOWICSUPPORT";
  else if (b1.gr.chr != b2.gr.chr && std::max(b1.nm, b2.nm) >= 3 && std::min(b1.matchlen, b2.matchlen) < 150) // inter-chr, but no disc reads, and too many nm
    confidence = "LOWICSUPPORT";
  else if (std::min(b1.matchlen, b2.matchlen) < 0.6 * sc->readlen)
    confidence = "LOWICSUPPORT";      
  else if (std::min(b1.mapq, b2.mapq) < 50 && b1.gr.chr != b2.gr.chr) // interchr need good mapq for assembly only
    confidence = "LOWMAPQ";
  else if (std::min(b1.matchlen, b2.matchlen) < 40 ||
	   (svtype == SVType::TSI_LOCAL && std::min(b1.matchlen, b2.matchlen) < 100)) // not enough evidence
    confidence = "LOWMATCHLEN";    
  else if (std::min(b1.matchlen - homology.length(), b2.matchlen - homology.length()) < 40)
    confidence = "LOWMATCHLEN";          
  else if ((/*b1.sub_n && */b1.mapq < 30) || (/*b2.sub_n && */b2.mapq < 30)) 
    confidence = "LOWMAPQ";
  else if (secondary && std::min(b1.mapq, b2.mapq) < 30)
    confidence = "SECONDARY";
  // else if ((repeat_seq.length() >= 10 && std::max(t.split, n.split) < 7) || hi_rep)
  //   confidence = "WEAKSUPPORTHIREP";
  else if (num_split < 6 && getSpan() < 300 && b1.gr.strand==b2.gr.strand) 
    confidence = "LOWQINVERSION";
  // else if ( (b1.matchlen - b1.simple < 15 || b2.matchlen - b2.simple < 15) )
  //   confidence = "SIMPLESEQUENCE";
  else if ((int)homology.length() * HOMOLOGY_FACTOR > sc->readlen) // if homology is too high, tough to tell from mis-assemly
    confidence = "HIGHHOMOLOGY";
  else
    confidence = "PASS";
  
  assert(confidence.length());
  
}

void BreakPoint::score_somatic(double error_fwd) {
  
  // this is LOD of normal being REF vs AF = 0.5+
  // We want this to be high for a somatic call
  // NB: this is almost the flip of LOD alt vs ref.
  //     where high log-odds is more likely there is
  //     a variant. It makes sense, you always think of
  //     log-odds as "higher-is-better", but are asking
  //     two different things:
  //        LO - is this variant versus ref
  //        LO_n - is this ref *in the normal* = somatic

  // this is a fully Bayesian approach -- if all normals strongly support ref, the sum is huge.
  // Cons: one bad normal (low or negative LO_n) can be overwhelmed by many good normals.
  // Commented out since favor conservative approach
  // double somatic_lod = 0.0;
  // for (const auto& [pref, al] : allele) {
  //   if (pref[0] == 'n') {
  //     somatic_lod += al.LO_n;
  //   }
  //}

  // this is a more conservative approach - if any normal BAM shows evidence for
  // variant, then get the somatic score from that
  int thiscov_n = n.cov;
  if (n.alt >= n.cov)  
    thiscov_n = n.alt;
  int thiscov_t = t.cov;
  if (t.alt >= t.cov)  
    thiscov_t = t.alt;

  // now to be conservative, also count normal near-cigar matches
  // as alt reads
  int normal_cigar = n.cigar + n.cigar_near;
  
  // adjust the alt count
  if (n.alt < normal_cigar)
    n.alt = normal_cigar;
  if (t.alt < t.split)
    t.alt = t.split;
  if (t.alt < t.cigar)
    t.alt = t.cigar;
  if (t.alt < t.split)
    t.alt = t.split;
  
  double a_cov_n = (double)thiscov_n * (double)(sc->readlen - 2 * T_SPLIT_BUFF)/sc->readlen;
  double a_cov_t = (double)thiscov_t * (double)(sc->readlen - 2 * T_SPLIT_BUFF)/sc->readlen;  
  double scaled_alt_n = std::min((double)n.alt, a_cov_n);  
  double scaled_ref_n = a_cov_n- scaled_alt_n;
  double scaled_alt_t = std::min((double)t.alt, a_cov_t);  
  double scaled_ref_t = a_cov_t- scaled_alt_t;
  double error_rev = 1e-6;
  
  LO_s =
    SvabaModels::SomaticLOD(scaled_alt_n, a_cov_n,
			    scaled_alt_t, a_cov_t,
			    error_fwd, error_rev);

  if (false)
  std::cerr << std::fixed << std::setprecision(6)
	    << "[DEBUG] somatic_lod=" << LO_s
	    << " | n.cov=" << n.cov
	    << " | n.alt=" << n.alt
	    << " | n.cigar=" << n.cigar
	    << " | t.cov=" << t.cov
	    << " | t.alt=" << t.alt
	    << " | t.cigar=" << t.cigar
	    << " | t.split=" << t.split
	    << " | thiscov_n=" << thiscov_n
	    << " | thiscov_t=" << thiscov_t
	    << " | readlen=" << sc->readlen
	    << " | T_SPLIT_BUFF=" << T_SPLIT_BUFF
	    << " | a_cov_n=" << a_cov_n
	    << " | a_cov_t=" << a_cov_t
	    << " | scaled_alt_n=" << scaled_alt_n
	    << " | scaled_ref_n=" << scaled_ref_n
	    << " | scaled_alt_t=" << scaled_alt_t
	    << " | scaled_ref_t=" << scaled_ref_t
	    << " | error_fwd=" << error_fwd
	    << " | error_rev=" << error_rev
	    << " | somatic_LOD=" << LO_s
	    << std::endl;
  
  /*for (const auto& [pref, al] : allele) {
    if (pref[0] == 'n') {
      LO_s = std::min(somatic_lod, al.LO_n);
    }
    }*/
  
  // find the somatic to normal ratio
  double ratio = n.alt > 0 ?
    (double)t.alt / (double)n.alt : 100;    
  
  if (svtype == SVType::INDEL) {
    
    // somatic score is just true or false for now
    // use the specified cutoff for indels, taking into account whether at dbsnp site
    double cutoff = (rs.empty() || rs=="x") ? sc->opts.lodSomatic : sc->opts.lodSomaticDb;
    somatic = LO_s > cutoff ? SomaticState::SOMATIC_LOD : SomaticState::NORMAL_LOD;
    
    // can't call somatic with 5+ normal reads or <5x more tum than norm ALT
    //if ((ratio <= 12 && n.cov > 10) || n.alt > 5)
    //if (n.alt > 5)
    // for SVs, use LOD and then also a hard cutoff
    // for gauging somatic vs germline
  } else {
    
    // passes
    somatic = (LO_s > sc->opts.lodSomatic) ? SomaticState::SOMATIC_LOD : SomaticState::NORMAL_LOD;
    
    // require no reads in normal or MAX 1 read and tons of tumor reads
    if (ratio < MIN_SOMATIC_RATIO ||
	n.split > 2 && dc.ncount > 1)
      somatic = SomaticState::FAILED;
  }
  
  // set germline if single normal read in discordant clsuter
  if (svtype == SVType::DSCRD && n.alt > 0 && somatic == SomaticState::SOMATIC_LOD)
    somatic = SomaticState::FAILED;
  
}

void BreakPoint::score_assembly_dscrd() {

  int this_mapq1 = b1.mapq;
  int this_mapq2 = b2.mapq;
  int span = getSpan();
  bool germ = dc.ncount > 0 || n.split > 0;
  
  int max_a_mapq = std::max(this_mapq1, dc.mapq1);
  int max_b_mapq = std::max(this_mapq2, dc.mapq2);
  
  // how much of contig is covered by split reads
  int cov_span = split_cov_bounds.second - split_cov_bounds.first ;
  
  // set the total number of supporting reads for tumor / normal
  // these alt counts should already be one qname per alt (ie no dupes)
  int t_reads = 0;
  int n_reads = 0;
  for (auto& [pref,al] : allele) {
    if (pref.at(0) == 't')
      t_reads += al.alt;
    else
      n_reads += al.alt;
  }
  
  int total_count = t_reads + n_reads; //n.split + t.split + dc.ncount + dc.tcount;
  int disc_count = dc.tcount + dc.ncount;
  //int hq = dc.tcount_hq + dc.ncount_hq;
  
  if ( (max_a_mapq < 30 && b1.local == LocalAlignment::FROM_DISTANT_REGION) ||
       (max_b_mapq < 30 && b2.local == LocalAlignment::FROM_DISTANT_REGION) ||
       (/*b1.sub_n > 7 && */b1.mapq < 10 && b1.local == LocalAlignment::FROM_DISTANT_REGION) ||
       (/*b2.sub_n > 7 && */b2.mapq < 10 && b2.local == LocalAlignment::FROM_DISTANT_REGION) )
    confidence = "LOWMAPQ";
  else if ( std::min(b1.nm, b2.nm) >= 10)
    confidence = "LOWMAPQ";
  else if ( std::min(b1.mapq, b2.mapq) < 10/* && std::min(dc.mapq1, dc.mapq2) < 10 */)
    confidence = "LOWMAPQ";      
  else if ( total_count < 4 || (std::max(t.split, n.split) <= 5 && cov_span < (sc->readlen + 5) && disc_count < 7) )
    confidence = "LOWSUPPORT";
  else if ( total_count < 15 && germ && span == -1) // be super strict about germline interchrom
    confidence = "LOWSUPPORT";
  else if ( std::min(b1.matchlen, b2.matchlen) < 50 && b1.gr.chr != b2.gr.chr ) 
    confidence = "LOWICSUPPORT";
  else if (secondary && getSpan() < 1000) // local alignments are more likely to be false for alignemnts with secondary mappings
    confidence = "SECONDARY";	
  /*else if (dc.tcount_hq + dc.ncount_hq < 3) { // multimathces are bad if we don't have good disc support too
    if ( ((b1.sub_n && dc.mapq1 < 1) || (b2.sub_n && dc.mapq2 < 1))  )
    confidence = "MULTIMATCH";
    else if ( ( (secondary || b1.sub_n > 1) && !b1.local) && ( std::min(max_a_mapq, max_b_mapq) < 30 || std::max(dc.tcount, dc.ncount) < 10)) 
    confidence = "SECONDARY";
    else 
    confidence = "PASS";
    }*/
  else
    confidence = "PASS";
}
 
void BreakPoint::score_dscrd() { 
  
  t.alt = dc.tcount;
  n.alt = dc.ncount;

  int nm_count = std::max(dc.nm1, dc.nm2);
  int disc_count = dc.ncount + dc.tcount;
  //int hq_disc_count = dc.ncount_hq + dc.tcount_hq;
  int disc_cutoff = 8;
  //int hq_disc_cutoff = disc_count >= 10 ? 3 : 5; // reads with both pair-mates have high MAPQ

  const int min_dscrd_size = 2000;
  if (getSpan() > 0 && (getSpan() < min_dscrd_size && b1.gr.strand == '+' && b2.gr.strand == '-')) // restrict span for del (FR) type 
    confidence = "LOWSPANDSCRD";
  else if ((disc_count < disc_cutoff || std::min(dc.mapq1, dc.mapq2) < 15))
    confidence = "LOWMAPQDISC";
  else if (nm_count >= 4)
    confidence = "HIGHNM"; 
  else if (!dc.m_id_competing.empty())
    confidence = "COMPETEDISC";
  else if ( disc_count < disc_cutoff)
    confidence = "WEAKDISC";
  else 
    confidence = "PASS";
  
  assert(confidence.length());
}

void BreakPoint::score_indel() {
  
  assert(b1.mapq == b2.mapq);
  
  bool is_refilter = !confidence.empty(); // act differently if this is refilter run
  
  // for refilter, only consider ones that were low lod or PASS
  // ie ones that with a different lod threshold may be changed
  // if confidence is empty, this is original run so keep going
  if (confidence != "LOWLOD" &&
      confidence != "PASS" &&
      is_refilter)
    return;
  
  double max_lod = 0;
  for (const auto& [_,al] : allele) 
    max_lod = std::max(max_lod, al.LO);

  // check if homozygous reference is most likely GT
  // bool homozygous_ref = true;
  // for (auto& [_,al] : allele) {
  //   if (al.genotype_likelihoods[0] > al.genotype_likelihoods[1] ||
  // 	al.genotype_likelihoods[0] > al.genotype_likelihoods[2])
  //     homozygous_ref = false;
  // }
  
  // get the allelic ractions, just for VLOWAF filter
  double af_t = t.cov > 0 ? (double)t.alt / (double)t.cov : 0;
  double af_n = n.cov > 0 ? (double)n.alt / (double)n.cov : 0;
  double af = std::max(af_t, af_n);

  assert(b1.local != LocalAlignment::NOTSET);  
  if (b1.local != LocalAlignment::FROM_LOCAL_REGION)
    confidence="NOLOCAL";
  if (b1.mapq < 10) 
    confidence="LOWMAPQ";
  //else if (!is_refilter && (double)aligned_covered / (double)seq.length() < 0.80) // less than 80% of read is covered by some alignment
  //  confidence = "LOWAS";  
  else if ((/*b1.sub_n && */b1.mapq < 30) || (/*b2.sub_n && */b2.mapq < 30)) 
    confidence = "LOWMAPQ";      
  else if (max_lod < sc->opts.lod && rs.empty())        // non db snp site
    confidence = "LOWLOD";
  else if (max_lod < sc->opts.lodDb && !rs.empty()) // be more permissive for dbsnp site
    confidence = "LOWLOD";
  else if (af < 0.05) // if really low AF, get rid of 
    confidence = "VLOWAF";
  else if (!is_refilter && std::min(left_match, right_match) < 20) 
    confidence = "SHORTALIGNMENT"; // no conf in indel if match on either side is too small
  else
    confidence="PASS";

  //debugprint
  // std::cerr <<" MAPQ " << b1.mapq << " b1.sub_n " << b1.sub_n <<
  //   " b2.sub_n " << b2.sub_n << " max_lod " << max_lod << " af " << af << " is_refilter " << is_refilter <<
  //   std::endl;
  // for (const auto& [pref,al] : allele) {
  //   std::cerr << " LOD SCORES " << pref << " - " << al.LO << std::endl;
  //   std::cerr << " PL likelihood 0/0 " << al.phred_likelihoods[0] <<
  //     " PL likelihood 0/1 " << al.phred_likelihoods[1] <<
  //     " PL likelihood 1/1 " << al.phred_likelihoods[2] << std::endl;
  // }
  // std::cerr << "AF T " << af_t << " AF N " << af_n << " AF " << af << std::endl;
  // std::cerr << " confidene " << confidence << std::endl;
  
}

// LOD cutoff is just whether this is ref / alt
// LOD cutoff dbsnp is same, but at DBSNP site (should be lower threshold since we have prior)
// LOD SOM cutoff is cutoff for whether NORMAL BAM(s) are AF = 0
// LOD SOM DBSNP cutoff same as above, but at DBSNP site (should be HIGHER threshold,
//   since we have prior that it's NOT somatic
void BreakPoint::scoreBreakpoint() {

  // some sanity checking
  assert(svtype != SVType::NOTSET);
  if (!confidence.empty() && confidence != "BLACKLIST")
    throw std::runtime_error("BreakPoint confidence already set (unexpected)");

  // set the error rate
  // depends on teh repeat context
  double error_rate = (repeat_seq.length() > 10) ?
    MAX_ERROR : ERROR_RATES[repeat_seq.length()];

  // first calculate log-odds and genotype likelihoods for each individual bam
  std::vector<double> lods;
  for (auto& [pref,al] : allele) {
    al.modelSelection(error_rate, sc->readlen);
    lods.push_back(al.LO);
  }

  /////
  // Calculate a site-level quality score
  ///// 
  // Sum the LODs across all your BAMs (assuming independence) to get a combined log-odds
  double lod_site = std::accumulate(lods.begin(), lods.end(), 0.0);
  // Convert that back to an error probability
  double p_err = 1.0 / (1.0 + std::pow(10.0, lod_site));
  // Phred-scale (-log10) to get a QUAL score
  int qual_calc  = static_cast<int>(std::lround(-10.0 * std::log10(p_err)));
  // Cap at 99 
  quality = std::min({qual_calc, 99});
  
  // combine the all tumor bam support and normal bam support
  for (auto& [pref, al] : allele) {
    if (pref.at(0) == 't') {
      t = t + al;
    } else {
      n = n + al;
    }
  }

  // combine all read support regardless of tumor / normal
  a = t + n;

  // std::cerr << " T MODEL " << error_rate << std::endl;
  // t.modelSelection(error_rate, sc->readlen);
  // std::cerr << " N MODEL " << error_rate << std::endl;  
  // n.modelSelection(error_rate, sc->readlen);
  
  // kludge. make sure we have included the DC counts (should have done this arleady...)
  if (svtype == SVType::DSCRD || svtype == SVType::ASDIS) {
    t.disc = dc.tcount;
    n.disc = dc.ncount;
  }
  
  // provide a scaled LOD that accounts for MAPQ. Heuristic, not really used
  //int mapqr1 = b1.local ? std::max(30, b1.mapq) : b1.mapq; // if local, don't drop below 30
  //int mapqr2 = b2.local ? std::max(30, b2.mapq) : b2.mapq; // if local (aligns to within window), don't drop below 30
  //double scale = (double)( std::min(mapqr1, mapqr2) - 2 * b1.nm) / (double)60;
  //for (auto& i : allele) 
  //  i.second.SLO = i.second.LO * scale;
  //t.SLO = t.LO * scale;
  //n.SLO = n.LO * scale;
  //a.SLO = a.LO * scale;
  
  // // sanity check
  // int split =0;
  // for (auto& i : allele) 
  //   split += i.second.split;
  // assert( (split == 0 && t.split == 0 && n.split==0) ||
  // 	  (split > 0 && (t.split + n.split > 0)));

  // if already overlapping blacklist, then we're done
  if (confidence == "BLACKLIST")
    return;
  
  // do the scoring
  bool iscomplex =
    svtype == SVType::TSI_LOCAL ||
    svtype == SVType::TSI_GLOBAL; 

  // score rearrangements of various kinds
  if (svtype == SVType::ASSMB || (iscomplex  && (dc.ncount + dc.tcount)==0))
    score_assembly_only();
  
  else if (svtype == SVType::ASDIS || (iscomplex && (dc.ncount + dc.tcount))) 
    score_assembly_dscrd();
  
  else if (svtype == SVType::DSCRD) 
    score_dscrd();
  
  // it failed assembly filters, but might pass discordant filters
  // if (evidence == "ASDIS" && confidence != "PASS") { 
  //   evidence = "DSCRD";
  //   score_dscrd(min_dscrd_size);
  // }
  else if (svtype == SVType::INDEL) 
    score_indel(); 
  else
    throw std::runtime_error("no scoring method called in BreakPoint::score_breakpoints");
  
  // set the somatic_score field to true or false
  score_somatic(error_rate); 
  
  // quality score is odds that read is
  // non-homozygous reference (max 99)
  // quality = 0;
  // for (auto& a : allele)
  //   quality = std::max(a.second.NH_GQ, (double)quality); 
  
}

void BreakPoint::setRefAlt(const RefGenome* main_rg,
			   const BamHeader& header) {
  
  assert(!main_rg->IsEmpty());
  assert(ref.empty());
  assert(alt.empty());
  assert(svtype != SVType::NOTSET);
  
  if (svtype != SVType::INDEL) {
    
    try {
      // get the reference for BP1
      ref = main_rg->QueryRegion(b1.gr.ChrName(header), b1.gr.pos1, b1.gr.pos1);
    } catch (const std::invalid_argument& ia) {}
    
    try {
      alt = main_rg->QueryRegion(b2.gr.ChrName(header), b2.gr.pos1, b2.gr.pos1);
    } catch (const std::invalid_argument& ia) {}
    
    
  } else {

    // insertion
    if (!insertion.empty()) {
      
      try {
	// we store b1.gr internally as zero-indexed,
	// to be consistent with htslib, so since this passed
	// to htslib query, keep as zero-indexed
	ref = main_rg->QueryRegion(b1.gr.ChrName(header),
				   b1.gr.pos1,
				   b1.gr.pos1);
      } catch (const std::invalid_argument& ia) {
	ref = "N";
	std::cerr << "Caught exception in BreakPoint:setRefAlt for indel ref: " << ia.what() << std::endl;
      }

      // alt 
      alt = ref + insertion;
      
    // deletion
    } else {	
      
      // reference
      assert(b2.gr.pos1 - b1.gr.pos1 - 1 >= 0);
      try {

	// again, 0-based, see above comment
	ref = main_rg->QueryRegion(b1.gr.ChrName(header),
				   b1.gr.pos1,
				   b2.gr.pos2);

	// if (cname == "c_1_49001_74001_14C")
	//   std::cerr << " REF : " << ref << " b1 " << b1.gr <<
	//     " b2 " << b2.gr << " cname " << cname << std::endl;
	
      } catch (const std::invalid_argument& ia) {
	ref = std::string(std::abs(b1.gr.pos1 - b2.gr.pos1), 'N');
	std::cerr << "Caught exception in BreakPoint:setRefAlt for indel ref: " << ia.what() << std::endl;	  
      }
      alt = ref.substr(0,1);
    }
  }
  
  if (ref.empty()) {
    ref = "N";
  }
  if (alt.empty()) {
    alt = "N";
  }
  
}

int BreakPoint::getSpan() const {

  assert(svtype != SVType::NOTSET);
  if (svtype == SVType::INDEL && insertion.empty())
    return b2.gr.pos1 - b1.gr.pos1;

  if (svtype == SVType::INDEL && !insertion.empty())
    return insertion.length();
  
  if (b1.gr.chr == b2.gr.chr)
    return abs(b2.gr.pos1-b1.gr.pos1);

  // interchromosomal
  else
    return -1;
}

// ReducedBreakPoint::ReducedBreakPoint(const std::string &line, const SeqLib::BamHeader& h) {

//   if (h.isEmpty()) {
//     std::cerr << "ReducedBreakPoint::ReducedBreakPoint - Must supply non-empty header" << std::endl;
//     exit(EXIT_FAILURE);
//   }

//   std::istringstream iss(line);
//   std::string val;
//   size_t count = 0;
  
//   ref = nullptr;
//   alt = nullptr;
//   cname = nullptr;
//   evidence = nullptr;
//   confidence = nullptr;
//   insertion = nullptr;
//   homology = nullptr;
  
//   //float afn, aft;
//   std::string ref_s, alt_s, cname_s, insertion_s, homology_s, evidence_s, confidence_s, read_names_s, bxtable_s;
  
//   std::string chr1, pos1, chr2, pos2, repeat_s; 
//   char strand1 = '*', strand2 = '*';
//   while (std::getline(iss, val, '\t')) {
//     try{
//       switch(++count) {
//       case 1: chr1 = val; break;
//       case 2: pos1 = val; break; 
//       case 3: assert(val.length()); strand1 = val.at(0); break;
//       case 4: chr2 = val; break;
//       case 5: pos2 = val; break; 
//       case 6: assert(val.length()); strand2 = val.at(0); break;
//       case 7: 
// 	ref_s = val;
// 	break; 
//       case 8: 
// 	alt_s = val;
// 	break;
//       case 9: break; //span = stoi(val); break; // automatically calculated
//       case 10: //mapq1
// 	b1 = ReducedBreakEnd(GenomicRegion(chr1, pos1, pos1, h), std::stoi(val)); b1.gr.strand = strand1; break;
//       case 11: //mapq2
// 	b2 = ReducedBreakEnd(GenomicRegion(chr2, pos2, pos2, h), std::stoi(val)); b2.gr.strand = strand2; break;
//       case 12: b1.nm = INTNSTOI(val,255); break;
//       case 13: b2.nm = INTNSTOI(val,255); break;
//       case 14: dc.mapq1 = INTNSTOI(val, 255); break;
//       case 15: dc.mapq2 = INTNSTOI(val, 255); break;
//       case 16: break; //split (not needed, since in genotype) 
//       case 17: break; //cigar (not needed, since in genotype) 
//       case 18: break; //alt (not needed, since in genotype)
//       case 19: cov = INTNSTOI(val,65535); break;
//       case 20: b1.sub_n = INTNSTOI(val,255); break;
//       case 21: b2.sub_n = INTNSTOI(val,255); break;
//       case 22: 
// 	homology_s = val;
// 	break; 
//       case 23: 
// 	insertion_s = val;
// 	break; 
//       case 24: cname_s = val; break;
//       case 25: num_align = std::min((int)31, std::stoi(val)); break;
//       case 26: 
// 	pass = val == "PASS";
// 	confidence_s = val;
// 	break;
//       case 27: 
// 	evidence_s = val;
// 	indel = val == "INDEL"; 
// 	imprecise = val == "DSCRD"; 
// 	break; 
//       case 28: quality = std::stod(val); break; //std::min((int)255,std::stoi(val)); break;
//       case 29: secondary = val == "1" ? 1 : 0;
//       case 30: somatic_score = std::stod(val); break;
//       case 31: LO_s = std::stod(val); break;
//       case 32: true_lod = std::stod(val); break;
//       case 33: pon = std::min(255,std::stoi(val)); break;
//       case 34: repeat_s = val; break; // repeat_seq
//       case 35: blacklist = (val=="1" ? 1 : 0); break;
//       case 36: dbsnp = val != "x"; break;
//       case 37: read_names_s = val; break; //reads
//       case 38: bxtable_s = val; break; //bx tags
//       default:
// 	format_s.push_back(val);
//       }
      
//     } catch(...) {
//       std::cerr << "caught stoi/stod/stof error on: " << val << " for count " << count << std::endl;
//       std::cerr << line << std::endl;
//       exit(1);
//     }
//   }
  
//   confidence = __string_alloc2char(confidence_s, confidence);
//   evidence   = __string_alloc2char(evidence_s, evidence);
//   insertion  = __string_alloc2char(insertion_s, insertion);
//   homology   = __string_alloc2char(homology_s, homology);
//   cname      = __string_alloc2char(cname_s, cname);
//   ref        = __string_alloc2char(ref_s, ref);
//   alt        = __string_alloc2char(alt_s, alt);
//   repeat     = repeat_s.empty() ? nullptr : __string_alloc2char(repeat_s, repeat);
//   if (somatic_score && confidence_s == "PASS")
//     read_names     = read_names_s; // == "x" ? nullptr : __string_alloc2char(read_names_s, read_names);
//   bxtable = bxtable_s;
  
// }


void BreakPoint::SampleInfo::modelSelection(double er, int readlen) {
  
  // can't have more alt reads than total reads
  // well you can for SVs...

  
  // adjust the coverage to be more in line with restrictions on ALT.
  // namely that ALT reads must overlap the variant site by more than T_SPLIT_BUFF
  // bases, but the raw cov calc does not take this into account. Therefore, adjust here
  //  if (readlen) {
  double a_cov = (double)cov * (double)(readlen - 2 * T_SPLIT_BUFF)/readlen;
  a_cov = a_cov < 0 ? 0 : a_cov;
  //  } else {
  //    a_cov = thiscov;
  //  }
  af = a_cov > 0 ? (double)alt / (double)a_cov : 1;
  af = af > 1 ? 1 : af;
  
  // mutect log liklihood against error
  // how likely to see these ALT counts if true AF is af
  // vs how likely to see these ALT counts if true AF is 0
  // and all the ALTs are just errors. 
  // The higher the error rate, the more negative ll_alt will go,
  // which will drive LO lower and decrease our confidence.
  // A high er will also drive ll_err up (or less negative), since
  // it will be more likely to see ALT reads generated by errors
  // As ALT and COV go higher, we should see LO go higher because
  // the indiviual calcs will have more confidence. Ultimately,
  // LO will represent the log likeihood that the variant is AF = af
  // vs AF = 0
  double scaled_alt = std::min((double)alt, a_cov);  
  double scaled_ref = a_cov - scaled_alt;

  // assume indels don't back-mutate as a sequencing error back to ref
  // but it's good to have this non-zero, so that with the scaling
  // tricks above, if we do get a non-zero ref count but var af = 1
  // (which technically isn't mathematically allowed), then having some
  // non-zero chance that the ref is a back-mutate makes the LO not
  // become uncomputable (which would revert to -1e12 which is not really what we want)
  double er_back = 0.000001; 
  double ll_alt = SvabaModels::LogLikelihood(scaled_ref, scaled_alt, af, er, er_back);
  double ll_err = SvabaModels::LogLikelihood(scaled_ref, scaled_alt, 0 , er, er_back); // likelihood that variant is actually true REF

  //std::cerr << "scaled_ref " << scaled_ref << " scaled_alt " << scaled_alt << " af " << af << " er " << er << " er back " << er_back << " ll_alt " << ll_alt << " ll_err " << ll_err << std::endl;
  // std::cerr << "SCALD REF " << scaled_ref << " scaled alt " <<
  //   scaled_alt << " af " << af << " ll_ALT " << ll_alt <<
  //   " ll_err " << ll_err << " scaled af " << scaled_af << std::endl;
  
  LO = ll_alt - ll_err; 
  
  // mutetct log likelihood normal
  // er = 0.0005; // make this low, so that ALT in REF is rare and NORM in TUM gives low somatic prob
  // actually, dont' worry about it too much. 3+ alt in ref excludes somatic anyways.
  // so a high LO_n for the normal means that we are very confident that site is REF only in 
  // normal sample. This is why LO_n from the normal can be used as a discriminant for 
  // germline vs somatic. eg if somatic_lod = normal.LO_n is above a threshold X
  // or above a larger threshold XX if at DBSNP site, then we accept as somatic
  // LO_n should not be used for setting the confidence that something is real, just the 
  // confidence that it is somatic
  // NB: the likelihood that the variant is actually a ref, and any alt reads are explained by errors, is ll_err.
  //     so ll_err = ll_ref_norm
  // NB: so it's a bit confusing here, but even though we are calculating the LO_n for
  //     any sample (tumor, normal, combined) that is passed, we are only interested in (and will use)
  //     normal.LO_n. This becomes the somatic score, as higher LO_n means it is more likely that the
  //     variant, specifically in the normal, is explained by sequencing errors with a ground-truth of not-present
  //     than by a germline indel. For the tumor, LO_n is not really meaningful, so we look at LO in general for
  //     whether the mutation is around AT ALL.
  double ll_alt_norm = SvabaModels::LogLikelihood(scaled_ref, scaled_alt, std::max(af, 0.5), er, er_back); // likelihood that variant is at AF >= 0.5 (as expected for germline)
  LO_n = ll_err - ll_alt_norm; // higher number means more likely to be AF = 0 (ref) than AF = 0.5 (alt). 

  // genotype calculation as provided in 
  // http://bioinformatics.oxfordjournals.org/content/early/2011/09/08/bioinformatics.btr509.full.pdf+html
  //int scaled_cov = std::floor((double)cov * 0.90);
  //int this_alt = std::min(alt, scaled_cov);
  genotype_likelihoods[0] = SvabaModels::GenotypeLikelihoods(2, er, scaled_alt, a_cov); // 0/0
  genotype_likelihoods[1] = SvabaModels::GenotypeLikelihoods(1, er, scaled_alt, a_cov); // 0/1
  genotype_likelihoods[2] = SvabaModels::GenotypeLikelihoods(0, er, scaled_alt, a_cov); // 1/1
  
  //debugprint
  //std::cerr << " ALT " << alt << " scaled alt " << scaled_alt << " ER " << er << " A_COV " << a_cov << 
  //  " 0/0 " << genotype_likelihoods[0] << " 0/1 " << genotype_likelihoods[1] << 
  //  " 1/1 " << genotype_likelihoods[2] << " LOD " << LO << " LO_n " << LO_n << 
  //  " af " << af << std::endl;

  double max_likelihood = *std::max_element(genotype_likelihoods.begin(), genotype_likelihoods.end());
  if (max_likelihood == genotype_likelihoods[0])
    genotype = "0/0";
  else if (max_likelihood == genotype_likelihoods[1])
    genotype = "0/1";
  else 
    genotype = "1/1";
  
  // compute PLs (phred-scaled likelihoods) 
  for(int i = 0; i < 3; ++i) {
    // PL = 10 * (GL_i  maxGL), rounded to nearest integer
    phred_likelihoods[i] = int(std::round( -10.0 * (genotype_likelihoods[i] - max_likelihood) ));
  }
  
  //debugprint
  // std::cerr << " thiscov " << thiscov << " alt " << alt << " cigar " << cigar << " split " << split <<
  //   " a_cov " << a_cov << " af " << af << " scaled_alt " << scaled_alt << " ll_alt " << ll_alt <<
  //   " ll_err " << ll_err << " LO " << LO << " ll_alt_norm " <<
  //   ll_alt_norm << " LO_n " << LO_n << " GT 0/0 " <<
  //   genotype_likelihoods[0] << " GT 0/1 " << genotype_likelihoods[1] <<
  //   " GT 1/1 " << genotype_likelihoods[2] << std::endl;
  
  // get the genotype quality
  GQ = SvabaModels::GenotypeQuality(phred_likelihoods);
  
  // get the genotype quality that it is not hom ref
  NH_GQ   = std::min(phred_likelihoods[0], 99);  
  
}

void BreakPoint::addCovs(const std::unordered_map<std::string, STCoverage*>& covs) {
  
  for (auto& [pref,i] : covs)  {
    int c = 0;
    for (int j = -COVERAGE_AVG_BUFF; j <= COVERAGE_AVG_BUFF; ++j) {
      c +=  i->getCoverageAtPosition(b1.gr.chr, b1.gr.pos1 + j);
      c +=  i->getCoverageAtPosition(b2.gr.chr, b2.gr.pos1 + j);
    }
    allele.at(pref).cov = c / 2 / (COVERAGE_AVG_BUFF*2 + 1);
  }
  
}

// std::string BreakEnd::print(const SeqLib::BamHeader& h) const {
//   return gr.ToString(h) + " - " + id + " mapq " + std::to_string(mapq) + " subn " + std::to_string(sub_n);
// }

/*  std::ostream& operator<<(std::ostream& out, const BreakEnd& b) {
    out << b.gr << " - " << b.id << " mapq " << b.mapq << " subn " << b.sub_n << std::endl;
    return out;
    }*/
  

std::string BreakPoint::SampleInfo::toFileString(SVType svtype) const {

  std::stringstream pl;
  pl << phred_likelihoods[0] << ',' 
     << phred_likelihoods[1] << ',' 
     << phred_likelihoods[2];
  
  // std::stringstream gt;
  // gt << genotype_likelihoods[0] << ',' 
  //    << genotype_likelihoods[1] << ',' 
  //    << genotype_likelihoods[2];
  
  std::stringstream ss;

  //GT:AD:DP:SR:DR:GQ:PL:LO:LO_n
  if (svtype == SVType::INDEL)
    ss << std::setprecision(4)
       << genotype << ":"
       << std::max(alt, cigar) << ":"
       << cov << ":"
       << split << ":"
       << cigar << ":"
       << GQ << ":"
       << pl.str() << ":"
       << LO << ":"
       << LO_n;
  else
    ss << std::setprecision(4)
       << genotype << ":"
       << alt << ":"
       << cov << ":"
       << split << ":"
       << disc << ":"
       << GQ << ":"
       << pl.str() << ":"
       << LO << ":"
       << LO_n;
  
  return ss.str();
  
}

void BreakPoint::SampleInfo::FillFromString(const std::string& s,
					    SVType svtype) { 
  
  std::string val;
  int count = 0;
  std::istringstream input(s);

  while (std::getline(input, val, ':')) {
    
    // tmp fix
    if (val.find("nan") != std::string::npos)
      val = "0";
    
    switch(++count) {

    case 1: genotype = val; break;      
    case 2: alt = std::stoi(val); break;
    case 3: cov = std::stoi(val); break;
    case 4: split = std::stoi(val); break;
    case 5: 
      if (svtype == SVType::INDEL)
	cigar = std::stoi(val);
      else
	disc = std::stoi(val);
      break;
    case 6: GQ = std::stod(val);  break;
    case 7: phred_likelihoods = svabaUtils::parsePLString(val);  break;
    case 8: LO =   std::stod(val); break;
    case 9: LO_n = std::stod(val); break;
    }
  }

}

// void BreakPoint::__rep(int rep_num, std::string& rseq, bool fwd) {
  
//   // move left and right from breakend 1
//   //int replen = 0;
//   //int curr_replen = 1;
//   rseq = "";
  
//   // return if we are too close to the edge for some reason
//   if (b1.cpos > (int)seq.length() - 5)
//     return;
  
//   const int REP_BUFF = 50;
  
//   std::string sss = seq;
//   if (!fwd)
//     std::reverse(sss.begin(), sss.end());
//   int cpos = b1.cpos + 1;
//   if (!fwd)
//     cpos = seq.length() - b1.cpos;
  
//   int i = cpos; //b1.cpos + 1;
//   int stop = std::min((int)seq.length() - 1, cpos + REP_BUFF); //fwd ? std::min((int)seq.length() - 1, b1.cpos + REP_BUFF) : std::max(0, b1.cpos - REP_BUFF);
//   int end = -1;
  
//   bool no_rep = true;
  
//   assert(stop - i < 80);
  
//   std::string c = sss.substr(cpos, std::min(rep_num, (int)seq.length() -  cpos));
  
//   i += rep_num; //*fr;
  
//   for (; i < stop; i+= rep_num) {
    
//     if (i >= (int)seq.length() - 1 || i + rep_num > (int)seq.length()) // shouldn't happen, but ensures no substr errors
//       break;
    
//     // terminating
//     if (c != sss.substr(i,rep_num)) { // || i >= (stop - rep_num)) {
      
//       end = i;
//       break;
      
//     } else {
      
//       no_rep = false;
//       //curr_replen += rep_num;
//     }
//   }
  
//   if (end == -1)
//     end = stop;
  
//   if (!no_rep)
//     rseq = sss.substr(cpos, std::min(end - cpos, (int)seq.length() - cpos));
// }

std::string BreakEnd::hash(int offset) const {
  
  return (std::to_string(gr.chr) + ":" + std::to_string(gr.pos1 + offset));
  
}

bool BreakPoint::hasMinimal() const {
  int count;
  count = dc.tcount + dc.ncount + t.split + n.split;
  
  if (local == LocalAlignment::FROM_DISTANT_REGION)
    return false;
  
  if (count >= 2)
    return true;
  
  return false;
}

std::string BreakPoint::getHashString() const {

  if (!isIndel())
    return "";

  bool isdel = insertion.length() == 0;
  std::string st = std::to_string(b1.gr.chr) +
    "_" + std::to_string(b1.gr.pos1) +  
    "_" + std::to_string(this->getSpan()) +
    (isdel ? "D" : "I");
  
  return st;
  
}

void BreakPoint::setLocal(const GenomicRegion& window) {
  
  b1.setLocal(window);
  b2.setLocal(window);

  assert(b1.local != LocalAlignment::NOTSET);
  assert(b2.local != LocalAlignment::NOTSET);

  // if this was already set as a non-variant local re-alignment
  // during alignedContig construction, then that has precedent
  if (local == LocalAlignment::NONVAR_LOCAL_REALIGNMENT)
    return;
  
  // set the BreakPoint level local flag
  if (b1.local == LocalAlignment::FROM_LOCAL_REGION ||
      b2.local == LocalAlignment::FROM_LOCAL_REGION)
    local = LocalAlignment::FROM_LOCAL_REGION;
  else
    local = LocalAlignment::FROM_DISTANT_REGION;
}

void BreakEnd::setLocal(const GenomicRegion& window) {

  assert(local == LocalAlignment::NOTSET);
  if (gr.GetOverlap(window))
    local = LocalAlignment::FROM_LOCAL_REGION;
  else
    local = LocalAlignment::FROM_DISTANT_REGION;
}

void BreakPoint::SampleInfo::UpdateAltCounts() {
  
  // get unique qnames
  std::unordered_set<std::string> qn;
  for (auto& r : supporting_reads) {
    size_t posr = r.find("_", 5) + 1; // extract the qname from tXXX_XXX_QNAME
    qn.insert(r.substr(posr, r.length() - posr));
  }
  
  // alt count is max of cigar or unique qnames (includes split and discordant)
  alt = std::max((int)qn.size(), cigar);
  
}


// bool ReducedBreakPoint::operator<(const ReducedBreakPoint& bp) const { 

//   //ASDIS > ASSMB > COMP > DSCRD
//   if (std::strcmp(evidence,bp.evidence) < 0) // <
//     return true;
//   else if (std::strcmp(evidence, bp.evidence) > 0) // >
//     return false;
  
//   if (cov > bp.cov)
//     return true;
//   else if (cov < bp.cov)
//     return false;

//   //if (nsplit > bp.nsplit) 
//   //  return true;
//   //else if (nsplit < bp.nsplit)
//   //  return false;
  
//   //if (tsplit > bp.tsplit)
//   //  return true;
//   //else if (tsplit < bp.tsplit)
//   //  return false;
  
//   //if (dc.ncount > bp.dc.ncount)
//   //  return true;
//   //else if (dc.ncount < bp.dc.ncount)
//   //  return false;
  
//   //if (dc.tcount > bp.dc.tcount)
//   // return true;
//   //else if (dc.tcount < bp.dc.tcount)
//   // return false;

//   // break the tie somehow
//   if (cname > bp.cname)
//     return true;
//   else if (cname < bp.cname)
//     return false;
  
//   return false;
// }

// std::ostream& operator<<(std::ostream& os, const ReducedBreakEnd& rbe) {
//     // Example format: [Chr Name: genomic_region (MapQ, SubN, NM)]
//     // Adjust formatting as needed
//      os << rbe.gr << " ("
//        << "MapQ: " << static_cast<int>(rbe.mapq) << ", "
//        << "SubN: " << static_cast<int>(rbe.sub_n) << ", "
//        << "NM: " << static_cast<int>(rbe.nm) << ")]";
//     return os;
// }

// make a breakpoint from a discordant cluster 
BreakPoint::BreakPoint(DiscordantCluster& tdc,
		       DiscordantClusterMap& dmap, 
		       const GenomicRegion& region,
		       const SvabaSharedConfig* sc_) : sc(sc_) {

  svtype = SVType::DSCRD;

  // instantiate the SampleInfo objects for each BAM
  for (auto const& p : sc->prefixes) {
    allele.try_emplace(p);
  }  
  
  num_align = 0;
  dc = tdc;
  
  std::string chr_name1, chr_name2;

  // not based on aligmnent of contig, so set as NA (-1)
  secondary = -1;
  
  try {
    // if this throw error, it means that there are more chr in 
    // reads than in reference
    chr_name1 = sc->header.IDtoName(dc.m_reg1.chr);
    chr_name2 = sc->header.IDtoName(dc.m_reg2.chr);      
    //	chr_name1 = bwa->ChrIDToName(dc.m_reg1.chr); //bwa->ChrIDToName(tdc.reads.begin()->second.ChrID());
    //      chr_name2 = bwa->ChrIDToName(dc.m_reg2.chr); //bwa->ChrIDToName(tdc.reads.begin()->second.ChrID());
  } catch (...) {
    std::cerr << "Warning: Found mismatch between reference genome and BAM genome for discordant cluster " <<
      dc.print(sc->header) << std::endl;
    chr_name1 = "Unknown";
    chr_name2 = "Unknown";
  }
  
  assert(chr_name1.length());
  assert(chr_name2.length());
  
  int pos1 = (dc.m_reg1.strand == '+') ? dc.m_reg1.pos2 : dc.m_reg1.pos1;
  int pos2 = (dc.m_reg2.strand == '+') ? dc.m_reg2.pos2 : dc.m_reg2.pos1;

  b1.gr = GenomicRegion(dc.m_reg1.chr, pos1, pos1);
  b2.gr = GenomicRegion(dc.m_reg2.chr, pos2, pos2);
  b1.mapq = dc.mapq1;
  b2.mapq = dc.mapq2;
  b1.gr.strand = dc.m_reg1.strand;
  b2.gr.strand = dc.m_reg2.strand;
  
  // set the alt counts, counting only unique qnames
  std::unordered_map<std::string, std::unordered_set<std::string>> alt_counts;
  for (auto const& p : sc->prefixes) {
    alt_counts.try_emplace(p);
}
  
  // add the supporting read info to alleles
  auto process = [this, &alt_counts](auto& container){
    for (auto& [_, rec] : container) {
      const auto& prefix = rec->Prefix();
      
      auto ait = allele.find(prefix);
      ait->second.supporting_reads.insert(rec->UniqueName());
      
      auto act = alt_counts.find(prefix);
      act->second.insert(rec->Qname());
    }
  };
  process(dc.reads);
  process(dc.mates);  
  
  // add the alt counts
  for (auto& [pref,al] : allele) {
    al.disc = alt_counts[pref].size(); //allele[i.first].supporting_reads.size();
    al.UpdateAltCounts();
    //i.second.alt = i.second.disc;
  }
  
  // give a unique id (OK to give an empty header here, since used internally)
  cname = dc.toRegionString(sc->header) + "__" + std::to_string(region.chr+1) + "_" + std::to_string(region.pos1) + 
    "_" + std::to_string(region.pos2) + "D";
  
  // check if another cluster overlaps, but different strands
  if (getSpan() > 800 || getSpan() == -1) { // only check for large events.
    for (auto& [_,d] : dmap) {
      
      // don't overlap if on different chr, or same event
      if (dc.m_reg1.chr != d.m_reg1.chr || dc.m_reg2.chr != d.m_reg2.chr || d.ID() == dc.ID())
	continue;
      
      // isolate and pad
      GenomicRegion gr1 = d.m_reg1;
      GenomicRegion gr2 = d.m_reg2;
      gr1.Pad(100);
      gr2.Pad(100);
      
      if (dc.m_reg1.GetOverlap(gr1) && dc.m_reg2.GetOverlap(gr2))
	if (dc.m_reg1.strand != d.m_reg1.strand || dc.m_reg2.strand != d.m_reg2.strand) {
	  dc.m_id_competing = d.ID();
	  tdc.m_id_competing = d.ID();
	  d.m_id_competing = dc.ID();
	}
    }
  }
}


void BreakPoint::indelCigarCheck(const CigarMapMap& cmap) {

  if (!isIndel())
    return;

  // get the hash count for this
  std::string hash = getHashString();

  // tally the number of reads that have cigar that support this
  // NB: the totalling for tumor, normal occur late in scoreBreakpoints
  for (auto& [sample,si] : allele) {
    size_t val = cmap.at(sample).count(hash) ? cmap.at(sample).at(hash) : 0;
    //std::cerr << " OG HASH VAL " << val << " hash " << hash << " sample " << sample << std::endl;
    si.cigar += val;
  }

  int span = getSpan();
  // for somatic, be extra careful and score "nearby" cigars
  // so e.g. if the normal has
  // this is really weird and inefficient way to do this, but should be so fast it's inconsequential
  char in_del[2] = {'D', 'I'};
  const int CIGAR_NEAR_BUFFER = 4;
  std::string thischr = std::to_string(b1.gr.chr) + "_";
  for (int i = -CIGAR_NEAR_BUFFER; i <= CIGAR_NEAR_BUFFER; i++) {
    for (int len = 1; len < (span + 5); len++) {
      std::string lenst = std::to_string(len);
      for (int ii = 0; ii < 2; ii++) { // loop I or D
	char c = in_del[ii];
	std::string st = thischr + std::to_string(b1.gr.pos1 + i) + "_" + lenst + c;
	//std::cerr << " cname " << cname << " st " << st << " hash " << hash << " i " << i << " len " << len << " ii " << ii << " span " << span << std::endl;

	if (hash == st)
	  continue; // don't redo this
	for (auto& [sample, si] : allele) {
	  size_t val = cmap.at(sample).count(st) ? cmap.at(sample).at(st) : 0;
	  si.cigar_near += val;
	  //	  std::cerr << " sample " << sample << " cname " << cname << " st " << st << " hash " << hash << " val " << val << std::endl;	  
	}
      }
    }
  }
  
}
