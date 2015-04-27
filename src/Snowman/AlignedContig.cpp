#include "AlignedContig.h"
#include <regex>
#include <unordered_map>
#include "SeqanTools.h"

using namespace std;
using namespace BamTools;

/* HENG LI CODE FROM SAMTOOLS -> BAM_IMPORT.C  */
#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

/* HENG LI CODE FROM SAMTOOLS -> BAM_IMPORT.C  */
static inline uint8_t *alloc_data(bam1_t *b, int size)
{
  if (b->m_data < size) {
    b->m_data = size;
    kroundup32(b->m_data);
    b->data = (uint8_t*)realloc(b->data, b->m_data);
  }
  return b->data;
}

#ifdef HAVE_SEQAN_BASIC_H
//#define SWALIGN 1
#endif

void AlignedContig::addAlignment(const BamTools::BamAlignment &align, const SnowTools::GenomicRegion &window,
				 const CigarMap &nmap, const CigarMap &tmap) { 

  AlignmentFragment tal(align, window, nmap, tmap); 
  m_align.push_back(tal);
  
  if (m_align.size() > 1)
    sort(m_align.begin(), m_align.end());

}

void AlignedContig::printContigFasta(std::ofstream& os) const {
  os << ">" << getContigName() << endl;
  os << getSequence() << endl;
}

void AlignedContig::blacklist(SnowTools::GenomicRegionCollection<SnowTools::GenomicRegion> &grv) {

 if (m_skip)
    return;

 // loop through the indel breaks and blacklist
 for (auto& i : m_align) 
   for (auto& j : i.indel_breaks) 
     j.checkBlacklist(grv);

}

void AlignedContig::splitCoverage() { 
  
  if (m_skip)
    return;

  // require that alignReadsToContig already ran
  assert(m_tried_align_reads);
  
  for (auto& i : m_align) {
    for (auto& j : i.indel_breaks) {
      //cout << "calling indel split coverage " << endl;
      j.splitCoverage(m_bamreads);
    }
  }
  
  for (auto& i : m_local_breaks) {
    //cout << "calling local split coverage " << endl;
    i.splitCoverage(m_bamreads);
  }

  if (!m_global_bp.isEmpty()) {
    //cout << "calling global split coverage " << endl;
    m_global_bp.splitCoverage(m_bamreads);
  }

}

ostream& operator<<(ostream& out, const AlignedContig &ac) {

  // print the global breakpoint
  if (!ac.m_global_bp.isEmpty())
    out << "Global BP: " << ac.m_global_bp << " contig " << ac.getContigName() << endl;       

  // print the multi-map breakpoints
  for (auto& i : ac.m_local_breaks)
    if (!i.isEmpty())
      out << "Multi-map BP: " << i << " contig " << ac.getContigName() << endl;       

  // print the indel breakpoints
  for (auto& i : ac.m_align)
    for (auto& j : i.indel_breaks) 
      if (!j.isEmpty())
	out << "Indel: " << j << " contig " << ac.getContigName() << endl;       

  // print the AlignmentFragments alignments
  for (auto& i : ac.m_align) 
    out << i << " Discordant: " << ac.printDiscordantClusters() << " contig " << ac.getContigName() << endl;

  // print the break locations for indel deletions
  for (auto& i : ac.m_align) {
    for (auto& j : i.indel_breaks) {
      if (j.isindel && j.insertion == "") // deletion
	out << string(j.cpos1, ' ') << "|" << string(j.cpos2-j.cpos1-1, ' ') << '|' << endl;	
    }
  }
    
   // print the contig base-pairs
  out << ac.getSequence() << "    " << ac.getContigName() << endl; 
   
   PlottedReadVector plot_vec;

   // print out the individual reads
   for (auto& i : ac.m_bamreads) {

     string seq ="", sr ="";
     int pos =-1, sw = -1;

     r_get_trimmed_seq(i, seq);
     r_get_Z_tag(i, "SR", sr);

     // reverse complement if need be
     int32_t rc;
     r_get_int32_tag(i, "RC", rc);
     if (rc)
       SnowTools::rcomplement(seq);

     // get the more complex tags (since there can be multiple annotations per tag)
     vector<int> posvec = SnowTools::GetIntTag(i, "AL");
     vector<int> swvec = SnowTools::GetIntTag(i, "SW");
     vector<string> cnvec = SnowTools::GetStringTag(i, "CN");
     assert(posvec.size() == swvec.size());
     assert(cnvec.size() == posvec.size());
     size_t kk = 0;
     for (; kk < cnvec.size(); kk++) 
       if (cnvec[kk] == ac.getContigName()) {
	 pos = posvec[kk];
	 sw = swvec[kk];
	 break;
       }

     assert(kk != cnvec.size()); // assure that we found something
     pos = abs(pos);
     int padlen = ac.getSequence().size() - pos - seq.size() + 5;
     padlen = max(5, padlen);

     stringstream rstream;
     assert(pos < 1e4 && padlen < 1e4); // bug, need to check
     rstream << sr << "--" << (r_id(i)+1) << ":" << r_pos(i)  << " SW: " << sw;

     plot_vec.push_back({pos, seq, rstream.str()});

   }

   sort(plot_vec.begin(), plot_vec.end());
  
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

   // plot the lines
   for (auto& i : line_vec) 
     out << i << " contig " << ac.getContigName() << endl;
   
   return out;
}

void AlignedContig::setMultiMapBreakPairs() {
   
  if (m_skip)
    return;

  // if single mapped contig, nothing to do here
  if (m_align.size() == 1)
    return;

  // initialize the breakpoint, fill with basic info
  BreakPoint bp;
  bp.seq = getSequence();
  bp.num_align = m_align.size();
  bp.cname = getContigName(); 
  assert(bp.cname.length());
  
  // set the discovar support if it's there
  parseDiscovarName(bp.disco_norm, bp.disco_tum);

  // walk along the ordered contig list and make the breakpoint pairs  
  for (AlignmentFragmentVector::iterator it = m_align.begin(); it != m_align.end() - 1; it++) {
    
    bp.gr1 = SnowTools::GenomicRegion(it->align.RefID, it->gbreak2, it->gbreak2);
    bp.gr2 = SnowTools::GenomicRegion((it+1)->align.RefID, (it+1)->gbreak1, (it+1)->gbreak1);
    //bp.gr1.strand = it->align.IsReverseStrand() ? '-' : '+';
    //bp.gr2.strand = (it+1)->align.IsReverseStrand() ? '+' : '-';
    bp.gr1.strand = !it->align.IsReverseStrand();
    bp.gr2.strand = (it+1)->align.IsReverseStrand();

    
    bp.cpos1 = it->break2; // take the right-most breakpoint as the first
    bp.cpos2 = (it+1)->break1;  // take the left-most of the next one
    
    assert(bp.cpos1 < 10000);
    assert(bp.cpos2 < 10000);

    // set the mapq
    bp.mapq1 = it->align.MapQuality;
    bp.mapq2 = (it+1)->align.MapQuality;
    
    assert(bp.mapq1 < 1000 && bp.mapq2 < 1000);
    
    // set the local
    bp.local1 = it->local;
    bp.local2 = (it+1)->local;
    
    // set the match length
    for (CigarOpVec::const_iterator cc = it->align.CigarData.begin(); cc != it->align.CigarData.end(); cc++)
      if (cc->Type == 'M')
	bp.matchlen1 += cc->Length;
    for (CigarOpVec::const_iterator cc = (it+1)->align.CigarData.begin(); cc != (it+1)->align.CigarData.end(); cc++)
      if (cc->Type == 'M')
	bp.matchlen2 += cc->Length;
    
    // set the NM
    int nmtag;
    if (!it->align.GetTag("NM", nmtag))
      nmtag = 0;
    bp.nm1 = nmtag;
    if (!(it+1)->align.GetTag("NM", nmtag))
      nmtag = 0;
    bp.nm2 = nmtag;
    
    // set the insertion / homology
    try {
      if (bp.cpos1 >= bp.cpos2)
	bp.homology = m_seq.substr(bp.cpos2, bp.cpos1-bp.cpos2);
      else if (bp.cpos2 >= bp.cpos1)
	bp.insertion = m_seq.substr(bp.cpos1, bp.cpos2-bp.cpos1);
      if (bp.insertion.length() == 0)
	bp.insertion = "";
      if (bp.homology.length() == 0)
	bp.homology = "";
    } catch (...) {
      unsigned hom = abs(bp.cpos1 - bp.cpos2);
      cout << "Caught error with contig on fine-getBreakPairs: " << it->align.Name << endl; 
      cout << "m_seq length: " << m_seq.length() << " bp.cpos1: " << bp.cpos1 << " bp.cpos2: " << bp.cpos2 << " bp.cpos1-bp.cpos2: " << hom << " m_seq: " << m_seq << endl;
    }

    m_local_breaks.push_back(bp);
    
  }
  
  // if this is a double mapping, we are done
  if (m_align.size() == 2) {
    m_global_bp = bp;
    m_local_breaks.clear();
    return;
  }
  
  // go through alignments and find start and end that reach mapq 
  size_t bstart = 1000; //1000 is a dummy
  size_t bend = m_align.size() - 1;
  for (size_t i = 0; i < m_align.size(); i++)
    if (m_align[i].align.MapQuality >= 60) {
      bend = i;
      if (bstart == 1000)
	bstart = i;
    }
  if (bstart == bend || bstart==1000) {
    bstart = 0;
    bend = m_align.size() -1 ;
  }
  assert(bend <= m_align.size());
  assert(bstart <= m_align.size());
  assert(bstart != 1000);
  
  // there are 3+ mappings
  m_global_bp = bp;
  m_global_bp.cpos1 = m_align[bstart].break2; // first mapping
  m_global_bp.gr1.pos1 = m_align[bstart].gbreak2;
  m_global_bp.gr1.pos2 = m_global_bp.gr1.pos1;
  m_global_bp.gr2.pos1 = m_align[bend].gbreak1;
  m_global_bp.gr2.pos2 = m_global_bp.gr2.pos1;
  m_global_bp.gr1.chr = m_align[bstart].align.RefID;
  m_global_bp.gr2.chr = m_align[bend].align.RefID;
  //m_global_bp.pos1  = m_align[bstart].gbreak2;
  m_global_bp.cpos2 = m_align[bend].break1; // last mapping
   //m_global_bp.pos2  = m_align[bend].gbreak1;
   //m_global_bp.refID1 = m_align[bstart].align.RefID; 
   //m_global_bp.refID2 = m_align[bend].align.RefID; 

   // set the strands
   //m_global_bp.gr1.strand = m_align[bstart].align.IsReverseStrand() ? '-' : '+';
   //m_global_bp.gr2.strand = m_align[bend].align.IsReverseStrand()   ? '+' : '-';
   m_global_bp.gr1.strand = !m_align[bstart].align.IsReverseStrand();
   m_global_bp.gr2.strand = m_align[bend].align.IsReverseStrand();


   // set the splits
   //m_global_bp.nsplit1 = m_align[bstart].nsplit2;
   //m_global_bp.tsplit1 = m_align[bstart].tsplit2;
   //m_global_bp.nsplit2 = m_align[bend].nsplit1;
   //m_global_bp.tsplit2 = m_align[bend].tsplit1;

   //m_global_bp.nsplit = min(m_global_bp.nsplit1, m_global_bp.nsplit2);
   //m_global_bp.tsplit = min(m_global_bp.tsplit1, m_global_bp.tsplit2);

   // set the mapping quality
   m_global_bp.mapq1 = m_align[bstart].align.MapQuality;
   m_global_bp.mapq2 = m_align[bend].align.MapQuality;

   if (m_global_bp.mapq1 > 60 || m_global_bp.mapq2 > 60) {
     cerr << "bad mapq GLOBAL" << endl;
     cerr << m_global_bp.toString() << endl;
     cerr << " m_align size " << m_align.size() << " bstart " << bstart << " bend " << bend << endl;
     exit(EXIT_FAILURE); 
   }

   // set the homologies
   try {
     if (m_global_bp.cpos1 >= m_global_bp.cpos2)
       m_global_bp.homology = m_seq.substr(m_global_bp.cpos2, m_global_bp.cpos1-m_global_bp.cpos2);
     else if (m_global_bp.cpos2 >= m_global_bp.cpos1)
       m_global_bp.insertion = m_seq.substr(m_global_bp.cpos1, m_global_bp.cpos2-m_global_bp.cpos1);
     if (m_global_bp.insertion.length() == 0)
       m_global_bp.insertion = "N";
     if (m_global_bp.homology.length() == 0)
       m_global_bp.homology = "N";
     m_global_bp.order();
   } catch (...) {
     cout << "Caught error with contig on global-getBreakPairs: " << getContigName() << endl;
   }

}
 
string AlignedContig::printDiscordantClusters() const {

  stringstream out;
  if (m_dc.size() == 0)
    return "No Discordant Clusters";

  for (vector<DiscordantCluster>::const_iterator it = m_dc.begin(); it != m_dc.end(); it++)
    out << *it << " ";
  return out.str();

}

// make an aligned contig from a sam record
AlignedContig::AlignedContig(const string &sam, const BamReader * reader, const SnowTools::GenomicRegion &twindow,
			     const CigarMap &nmap, const CigarMap &tmap) {
  
  window = twindow;
   
  samrecord = sam;

  std::istringstream iss(sam);
  std::string val, line;

  int i = 0;
  string cigar;

  std::regex reg_xp("^XA:Z:(.*)");
  std::regex reg_nm("^NM:[A-Za-z]:(.*)");

  //tryit
  //bam1_t *b = bam_init1();
  //bam1_core_t *c = &b->core;
  //size_t doff = 0;

  while (getline(iss, line, '\n')) {
    std::istringstream issv(line);

    BamTools::BamAlignment a;
    while(getline(issv, val, '\t')) {

      switch(i) {
      case 0 : a.Name = val; 
	//c->l_qname = val.length() + 1; 
	//memcpy(alloc_data(b, doff + c->l_qname) + doff, /*str->s*/ val.c_str(), c->l_qname);
	break;
      case 1 : a.AlignmentFlag = std::stoi(val); break;
      case 2 : try {a.RefID = (val == "*") ? -1 : reader->GetReferenceID(val);} catch (...) { a.RefID = -1; } break;
      case 3 : a.Position = std::stoi(val); break;
      case 4 : a.MapQuality = std::stoi(val); break;
      case 5 : cigar = val; break;
	//case 7 -8 are for paired end reads
      case 9 : a.QueryBases = val; break;
      case 10 : (val == "*") ? a.Qualities = string(a.QueryBases.length(), 'I') : a.Qualities = val; break;
      }

      // set the sequence, always on positive strand
      if (m_seq.length() == 0 && i == 9) 
	m_seq = val;

      // deal with the cigar
      if (i == 5 && val != "*") 
	a.CigarData = BamToolsUtils::stringToCigar(val);

      // process the tags
      //if (i > 10) 
      //  parseTags(val, a);

      i++;
	
    }

    assert(a.MapQuality <= 60);
    if (a.CigarData.size() > 0) // ensure that it is at least mapped
      addAlignment(a, window, nmap, tmap);
    
    i = 0;
  }

  // reject if not at last partially in window
  m_skip = true;
  for (auto& ca : m_align)
    if (ca.local)
      m_skip = false;

  // set break pairs for multi mapped reads
  setMultiMapBreakPairs();

  // set the min alignment end length
  //if (m_align.size() > 1) {
  //  m_min_align_len = min(m_align[0].align.Length, m_align.back().align.Length);
  //}

}

void AlignedContig::alignReadsToContigs(ReadVec &bav) {

  // BOWTIE ATTEMPT
  /*typedef String<Dna> TString;
  typedef StringSet<TString> TStringSet;
  typedef Index<StringSet<TString>, FMIndex<> > TIndex;

  TStringSet stringSet, readSet;
  TString str0 = "TATAGTACGTGCTATATCGGCGATATCCGATCGATTACTGCGGACTACTATCGAGCGACGATCTACGGCGATCATCGATCTACTAGC";
  appendValue(stringSet, str0);
  TString read0 = "TCGGCGATATCCGATCGATTACTGCGGACTACTATCGAGCGACGATCT";
  appendValue(readSet, read0);

  seqan::CharString haystack = "Simon, send more money!";
  seqan::CharString needle = "more";
  */
  //FragmentStore<> frag;
  
  
  ////////////////////
  
  if (m_skip)
    return;

  m_tried_align_reads = true; // set the try flag

  // MATCHING BY FIND
  int buff = 8;
  //int pad = 10;

#ifdef SWALIGN
  TSequence contig = m_seq;
#endif

  // 
  for (auto& j : bav) {

    string QB;
    r_get_trimmed_seq(j, QB);
    
    string RQB;
    int seqlen = QB.length();
    int32_t score = seqlen * 4;;
    
    // 
    if (seqlen > 35 && r_cig_size(j) > 1 /* don't align 101M, dont believe that they could be split */) {

      string read_name;
      r_get_Z_tag(j, "SR", read_name);
      string short_name = read_name.substr(0,2);
      
      size_t pos = m_seq.find(QB);
      int32_t aligned_pos;
      
      bool addread = false;
      bool isrev = false;
      
      // make some substrings to match. Dont proceed if none of these match anywhere
      string sub1 = QB.substr(10, buff); // 10-(10+buff)
      string sub2 = QB.substr(seqlen-20, buff); // (end-20) - (end-(20-buff))
      
      // matched completely, add
      if (pos != string::npos) {
	aligned_pos = pos;
	addread = true;
      } 
      if (!addread) {
	RQB = QB;
	SnowTools::rcomplement(RQB);
	pos = m_seq.find(RQB);
	if (pos != string::npos) {
	  aligned_pos = pos;
	  addread = true;
	  isrev = true;
	}
      }
#ifdef SWALIGN
      int cutoff = score - 25;
      // didn't match completely, SW align
      if (!addread) {
	if ((m_seq.find(sub1) != string::npos || m_seq.find(sub2) != string::npos || true) ) {
	  if (SeqanTools::SWalign(contig, aligned_pos, QB, score, cutoff, false, indel)) {
	    addread = true;
	  }
	}
      }
      
      // reverse complement the attempts
      if (!addread) {
	SnowTools::rcomplement(sub1);
	SnowTools::rcomplement(sub2);
      }
      
      // forwards SW didn't make it, try reverse
      if (!addread)
	if ((m_seq.find(sub1) != string::npos || m_seq.find(sub2) != string::npos || true)) {
	  if (SeqanTools::SWalign(contig, aligned_pos, RQB, score, cutoff, false, indel)) {
	    isrev = true;
	    addread = true;
	    //aligned_pos = QB.length() - 1 - aligned_pos; // because it is rev comp
	    //assert(aligned_pos >= 0);
	    //assert(aligned_pos < QB.length());
	  }
	}
#endif
      
      // add some tags. remove others
      if (addread) {
	
	SnowTools::SmartAddTag(j, "AL", to_string(aligned_pos));
	SnowTools::SmartAddTag(j, "CN", getContigName());
	SnowTools::SmartAddTag(j, "SW", to_string(score));
	
	int tt = isrev ? 1 : 0;
	r_add_int32_tag(j, "RC", tt); // flag for rev comp
	//j->EditTag("TS", "Z", QB); // stores reverse comp if need be, or shortened read
	
	m_bamreads.push_back(j); // make a copy of the data
      }
    } // end if > 20
    
  } // end read loop
  
}

AlignmentFragment::AlignmentFragment(const BamTools::BamAlignment &talign, const SnowTools::GenomicRegion &window, 
				     const CigarMap &nmap, const CigarMap &tmap) {
  
  
  m_seq = talign.QueryBases;
  m_name = talign.Name;

  align = talign;

  assert(align.CigarData.size());

  cigar = align.CigarData;

  SnowTools::GenomicRegion tmpw = window;
  tmpw.pad(1000);

  if (tmpw.getOverlap(SnowTools::GenomicRegion(align.RefID, align.Position, align.Position)))
    local = true;

  // orient cigar so it is on the contig orientation. 
  // need to do this to get the right ordering of the contig fragments below
  if (talign.IsReverseStrand())
    BamToolsUtils::flipCigar(cigar);

  // find the start position of alignment ON CONTIG
  start = 0; 
  for (auto& i : cigar) {
    if (i.Type != 'M')
      start += i.Length;
    else
      break;
  }

  // set the left-right breaks
  unsigned currlen  = 0; 

  // cigar is oriented to as is from aligner
  for (auto& i : cigar /*align.CigarData*/) { //CigarOpVec::const_iterator j = align.cigar.begin(); j != align.cigar.end(); j++) {
    
    // SET THE CONTIG BREAK (treats deletions and leading S differently)
    // the first M gets the break1, pos on the left
    if (i.Type == 'M' && break1 == -1)
      break1 = currlen;
    if (i.Type != 'D') // m_skip deletions but not leading S, but otherwise update
      currlen += i.Length;
    if (i.Type == 'M') // keeps triggering every M, with pos at the right
      break2 = currlen;
    
  }
  
  gbreak1 = align.Position;
  gbreak2 = align.GetEndPosition() - 1;
    
  assert(break1 < 10000);
  assert(break2 < 10000);
  if (break1 < 0)
    cout << "break " << break1 << " cigar "  << BamToolsUtils::cigarToString(cigar) << " frag " << *this << endl;
  assert(break1 >= 0);
  assert(break2 >= 0);
  
  // parse right away to see if there are indels on this alignment
  BreakPoint bp;
  size_t fail_safe_count = 0;
  while (parseIndelBreak(bp) && fail_safe_count++ < 100) 
    indel_breaks.push_back(bp);

  assert(fail_safe_count != 100);

  // set the cigar matchesx
  indelCigarMatches(nmap, tmap);
  
}

ostream& operator<<(ostream &out, const AlignmentFragment &c) {
  
  // sets the direction to print
  char jsign = '>'; 
  if (c.align.IsReverseStrand())
    jsign = '<';
  
  // print the cigar value per base
  for (auto& j : c.cigar) { //c.align.CigarData) { // print releative to forward strand
    if (j.Type == 'M')
      out << string(j.Length, jsign);
    else if (j.Type == 'I') 
      out << string(j.Length, 'I');
    else if (j.Type == 'S' || j.Type == 'H')
      out << string(j.Length, '.');
  }
  
  // print the info
  out << "    " << c.align.RefID + 1 << ":" << c.align.Position 
      << " MAPQ: " << c.align.MapQuality << " OrientedCigar: " << BamToolsUtils::cigarToString(c.align.CigarData)
	  << " Start: " << c.start
	  << " Breaks: " << c.break1 << " " << c.break2 
      << " GenomeBreaks: " << c.gbreak1 << " " << c.gbreak2 << " cname " << c.m_name;
    //	  << " Tsplits: " << c.tsplit1 << " " << c.tsplit2 
    //    << " Nsplits: " << c.nsplit1 << " " << c.nsplit2;
  
  return out;
}


void AlignmentFragment::indelCigarMatches(const CigarMap &nmap, const CigarMap &tmap) { 

  for (auto& i : indel_breaks) {

    assert(i.isindel);
    if (i.getSpan() <= 0) {
      cerr << "weird span detected " << i.getSpan();
      cerr << i << endl;
      cerr << *this << endl;
    }
    //assert(i.getSpan() > 0);
    
    string st = i.getHashString();

    CigarMap::const_iterator ffn = nmap.find(st);
    CigarMap::const_iterator fft = tmap.find(st);

    if (ffn != nmap.end())
      i.ncigar = ffn->second;
    if (fft != tmap.end())
      i.tcigar = fft->second;

  }
}

bool AlignedContig::isWorse(const AlignedContig &ac) const {

  // TODO.  right now is TRUE if any more than 1 bp per contig
  vector<const BreakPoint*> bpv1 = getAllBreakPointPointers();
  vector<const BreakPoint*> bpv2 = ac.getAllBreakPointPointers();
  if (bpv1.size() != 1 || bpv2.size() != 1)
    return false;
  
  bool same = (*(bpv2.back())) == (*(bpv1.back()));
  if (!same)
    return false;
    
  int mapq_this = max(bpv1.back()->mapq1, bpv1.back()->mapq2);
  int mapq_that = max(bpv2.back()->mapq1, bpv2.back()->mapq2);

  if (mapq_this < mapq_that) {
    return true;
  } 

  if ( (m_seq.length() < ac.m_seq.length()) && (mapq_this == mapq_that)) {
    return true;
  }

  return false;
  

}

bool AlignmentFragment::parseIndelBreak(BreakPoint &bp) {
  
  // reject out of hand if not in interval
  if (!local)
    return false;

  //cout << "parsing indel break" << endl;

  assert(cigar.size());

  // reject if it has small matches, could get confused. Fix later
  for (auto& i : cigar) 
    if (i.Type == 'M' && i.Length < 5)
      return false;

  // reject if first alignment is I or D
  for (auto& i : cigar) {
    if (i.Type == 'D' || i.Type == 'I') {
      cerr << "rejcting cigar for starting in I or D" << endl;
      return false;
    }
    if (i.Type == 'M')
      break;
  }

  // reject if last alignment is I or D
  CigarOpVec tmpcig = cigar;
  BamToolsUtils::flipCigar(tmpcig);
  for (auto& i : tmpcig) {
    if (i.Type == 'D' || i.Type == 'I') {
      cerr << "rejcting cigar for ending in I or D" << endl;
      return false;
    }
    if (i.Type == 'M')
      break;
  }

  // use next available largest D / I
  size_t loc = 0; // keep track of which cigar field
  for (auto& i : cigar) {
    loc++;
    if (i.Type == 'D' || i.Type == 'I') {
      
      // if it starts / stop with I, D, reject
      if (loc == 1 || loc == cigar.size()) 
	return false;
      if (loc > idx) {
	idx = loc;
	break;
      }
	
    }
  }

  // we made it to the end, no more indels to report
  //if (idx == align.CigarData.size())
  //  return false;

  // we made it to the end, no more indels to report
  if (loc == cigar.size())
    return false;

  // clear out the old bp just in case
  bp.insertion = "";
  bp.homology = "";

  bp.seq = m_seq;
  bp.isindel = true;

  int curr = 0;
  int gcurrlen = -1;

  bp.gr1.pos1 = -1;
  bp.gr1.pos2 = -1;
  bp.gr2.pos1 = -1;
  bp.gr2.pos2 = -1;
  bp.gr1.chr = align.RefID;
  bp.gr2.chr = align.RefID;

  bp.cpos1 = -1;
  bp.cpos2 = -1;
  bp.num_align = 1;
  bp.cname = m_name;
  bp.mapq1 = align.MapQuality;
  bp.mapq2 = align.MapQuality;
  bp.gr1.strand = true;
  bp.gr2.strand = false;

  assert(bp.cname.length());

  size_t count = 0; // count to make sure we are reporting the right indel

  for (auto& i : cigar) { // want breaks in CONTIG coordaintes, so use oriented cigar
    count++;

    // set the contig breakpoint
    if (i.Type == 'M' || i.Type == 'I') 
      curr += i.Length;
    if (i.Type == 'D' && bp.cpos1 == -1 && count == idx) {

      if (m_name == "c_1_8456000_8462000_45")
	cout << "Deletion found of " << i.Length << i.Type << " and insertion " << bp.insertion << endl;
	
      bp.cpos1 = curr-1;
      bp.cpos2 = curr;
    } 
    if (i.Type == 'I' && bp.cpos1 == -1 && count == idx) {
      bp.cpos1 = curr - i.Length - 1;
      bp.cpos2 = curr - 1;
      bp.insertion = align.QueryBases.substr(curr, i.Length);

      if (m_name == "c_1_8456000_8462000_45")
	cout << "Insertion found of " << i.Length << i.Type << endl;

    }

    // set the genome breakpoint
    if (bp.cpos1 > 0) {
      if (i.Type == 'D') {
	if (!align.IsReverseStrand()) {
	  bp.gr1.pos1 =  align.Position + gcurrlen; // dont count this one//bp.cpos1 + align.Position; //gcurrlen + align.Position;
	  bp.gr2.pos1 =  bp.gr1.pos1 + i.Length + 1;
	} else {
	  bp.gr2.pos1 =  (align.GetEndPosition()-1) - gcurrlen; //bp.cpos1 + align.Position; //gcurrlen + align.Position;
	  bp.gr1.pos1 =  bp.gr2.pos1 - i.Length - 1;
	}
      } else if (i.Type == 'I') {
	if (!align.IsReverseStrand()) {
	  bp.gr1.pos1 = align.Position + gcurrlen; //gcurrlen + align.Position;
	  bp.gr2.pos1 = bp.gr1.pos1 + 1;	
	} else {
	  // GetEndPosition is 1 too high
	  bp.gr2.pos1 = (align.GetEndPosition()-1) - gcurrlen; //gcurrlen + align.Position;
	  bp.gr1.pos1 = bp.gr2.pos1 - 1;	
	}
      }
      break; // already got it, so quit cigar loop
    }
    
    // update the position on the genome
    if (i.Type == 'M' || i.Type == 'D') {
      gcurrlen += i.Length;
    } 


  } // end cigar loop

  // set the dummy other end
  bp.gr1.pos2 = bp.gr1.pos1;
  bp.gr2.pos2 = bp.gr2.pos1;

  // error check the length
  int seq_len = m_seq.length();
  if (bp.cpos1 > seq_len || bp.cpos2 > seq_len) {
    cerr << "bp " << bp << endl;
    cerr << "CA " << *this << endl;
    cerr << "seq len " << seq_len << endl;
  }
  assert(bp.cpos1 <= seq_len);
  assert(bp.cpos2 <= seq_len);

  bp.order();

  if (m_name == "c_1_8456000_8462000_45")
    cout << "FINAL found of " <<  bp.insertion << endl;
  
  return true;
}

bool AlignedContig::parseDiscovarName(size_t &tumor, size_t &normal) {

  string s;
  bool valid = m_align[0].align.GetTag("TN", s);
  if (!valid)
    return false;
  
  // set the tumor support
  std::regex reg("^t([0-9]+)_*");
  std::smatch match;
  //if (std::regex_search(s.begin(), s.end(), match, reg))
  if (std::regex_search(s, match, reg))
    tumor = std::stoi(match[1].str());
  
  // set the normal support
  std::regex reg2("^t[0-9]+n([0-9]+)");
  std::smatch match2;
  if (std::regex_search(s, match2, reg2))
    normal = std::stoi(match2[1].str());
  
  return false;
}

vector<const BreakPoint*> AlignedContig::getAllBreakPointPointers() const  {

  vector<const BreakPoint*> out;
  for (auto& i : m_align) {
    for (auto& k : i.indel_breaks)
      out.push_back(&k);
  }
  
  if (!m_global_bp.isEmpty())
    out.push_back(&m_global_bp);
  
  return out;
}

vector<BreakPoint> AlignedContig::getAllBreakPoints() const {

  vector<BreakPoint> out;
  for (auto& i : m_align) {
    for (auto& k : i.indel_breaks)
      out.push_back(k);
  }
  
  if (!m_global_bp.isEmpty())
    out.push_back(m_global_bp);
  
  return out;
}


bool AlignedContig::hasVariant() const { 
  
  if (!m_global_bp.isEmpty())
    return true;

  for (auto& i : m_align)
    if (i.indel_breaks.size())
      return true;

  return false;
  
}
