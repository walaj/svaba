#include "AlignedContig.h"
#include <regex>
#include "GenomicRegion.h"
#include <unordered_map>
#include "SVBamReader.h"
#include "api/algorithms/Sort.h"
#include "SnowUtils.h"
#include "SeqanTools.h"

using namespace std;
using namespace BamTools;

#define SPLIT_BUFF 12


// constructor taking in a BAM record from the contig BAM
AlignedContig::AlignedContig(const BamTools::BamAlignment align) { 

  m_window = GenomicRegion(align.Name);
  addAlignment(align);

}

// if we are keeping this AlignedContig, then do more heavy processing
void AlignedContig::settleContigs() {

  m_window = GenomicRegion(m_tmpalign[0].Name);

  for (BAVec::const_iterator it = m_tmpalign.begin(); it != m_tmpalign.end(); it++)
    addAlignment(*it);

  m_tmpalign.clear();
}

// add a new contig alignment
void AlignedContig::addAlignment(const BamTools::BamAlignment align) { 

  CigarOpVec cig = orientCigar(align);

  string cigstring = cigarToString(cig);
    
  CAlignment tal(align, cig, cigstring);

  // add the alignment to the seq-field
  if (align.QueryBases.length() > m_seq.length()) {
    m_seq = align.QueryBases;
    if (align.IsReverseStrand()) // make sure on right strand
      SnowUtils::rcomplement(m_seq);
  }
    
  // find the start position
  tal.start = 0; 
  for (CigarOpVec::const_iterator it = cig.begin(); it != cig.end(); it++) 
    if (it->Type != 'M')
      tal.start += it->Length;
    else
      break;
  
  // determine if the current alignment is local
  int buff = 8000;
  GenomicRegion buff_region(align.RefID, align.Position - buff, align.Position + buff);
  if (buff_region.getOverlap(m_window) != 0) {
    m_local = true; // overall AlignedContig has a local
    tal.ca_local = true; // this particular alignment is local
  }
  //OLD: if (align.RefID == m_window.chr && align.Position > m_window.pos1 - buff && align.Position < m_window.pos2 + buff)

  // set the breakpoints
  setBreaks(tal);

  // set the max mapq
  m_maxmapq = max(m_maxmapq, (int)align.MapQuality);
  m_minmapq = min(m_minmapq, (int)align.MapQuality);

  // add to the overall struct
  //if (align.IsPrimaryAlignment())
  m_align.push_back(tal); 
  //else
  //  m_align_second.push_back(tal);       

  // ensure that they are sorted
  if (m_align.size() > 1)
    sort(m_align.begin(), m_align.end());

}

// make work function for getting PER-ALIGNMENT breaks from BWA-MEM
void AlignedContig::setBreaks(CAlignment &align) { 

    unsigned currlen  = 0; 
    unsigned gcurrlen = 0;
    BAVec::const_iterator i;
    CigarOpVec::const_iterator j, jend;

    align.break1  = -1; // dummy initialize
    align.break2  = -1; 
    align.gbreak1 = -1; // dummy initialize
    align.gbreak2 = -1; 

    bool trig = false;
    for (CigarOpVec::const_iterator j = align.cigar.begin(); j != align.cigar.end(); j++) {

      // SET THE CONTIG BREAK (treats deletions and leading S differently)
      // the first M gets the break1, pos on the left
      if (j->Type == 'M' && align.break1 == -1)
        align.break1 = currlen;
      if (j->Type != 'D') // skip deletions but not leading S, but otherwise update
	  currlen += j->Length;
      if (j->Type == 'M') // keeps triggering every M, with pos at the right
        align.break2 = currlen;

      // SET THE GENOME BREAK (treats insertions differently)
      //if (j->Type == 'M' && align.gbreak1 == -1
      //if (j->Type == 'M' && !trig)
      //	trig = true;
      trig = j->Type == 'M' || trig;
      //        align.gbreak1 = gcurrlen;
      //if (j->Type != 'I' && ( (j->Type != 'S' || j->Type != 'H') || trig)) // skip insertions and leading S, but otherwise update
      if (j->Type != 'I' && j->Type != 'S' && j->Type != 'H' && trig) 
	gcurrlen += j->Length;
      //if (j->Type == 'M') // keeps triggering every M, with pos at the right
      //  align.gbreak2 = gcurrlen;
    }

    // convert genome break in contig coords to reference coords
    if (align.align.IsReverseStrand()) {
      //align.gbreak1 = gcurrlen + align.align.Position - align.gbreak1 + 1;
      align.gbreak1 = align.align.Position + gcurrlen;
      align.gbreak2 = align.align.Position + 1; //gcurrlen + align.align.Position - align.gbreak2;
    } else {
      align.gbreak1  = align.align.Position + 1;
      align.gbreak2  = align.align.Position + gcurrlen;
      //align.gbreak1 += align.align.Position + 1;
      //align.gbreak2 += align.align.Position + 1 + gcurrlen;         
    }

    
    /*if (align.align.Name.compare("contig_8:94302502-94310502_1")==0)
      cout << "GBREAK1: " << align.gbreak1 << " GBREAK2: " << align.gbreak2 << 
	" gcurrlen: " << gcurrlen << " IsReverseStrand: " << align.align.IsReverseStrand() << 
	" cigar: " << align.cigstring << " Position: " << align.align.Position << endl;*/
}

// print to a fasta
void AlignedContig::printContigFasta(ofstream &ostream) const {
  ostream << ">" << m_align[0].align.Name << endl;
  ostream << m_seq << endl;
}

void AlignedContig::sortReads() { 
  std::sort( m_bamreads.begin(), m_bamreads.end(), BamTools::Algorithms::Sort::ByTag<int32_t>("AL", BamTools::Algorithms::Sort::AscendingOrder) );
}

void AlignedContig::splitCoverage() { 
  
   for (AlignVec::iterator it = m_align.begin(); it != m_align.end(); it++) {

     it->nsplit1 = 0;
     it->tsplit1 = 0;
     it->nsplit2 = 0;
     it->tsplit2 = 0;

     for (vector<BamAlignment>::const_iterator j = m_bamreads.begin(); j != m_bamreads.end(); j++) {

       string jw;
       assert(j->GetTag("JW", jw));
       bool tumor_read = (jw.at(0) == 't');
       string seq;
       assert(j->GetTag("TS", seq));
       int pos;
       assert(j->GetTag("AL", pos));

       int rightend = pos + seq.length();
       int leftend  = pos;
       int rightbreak1= it->break1 + SPLIT_BUFF;
       int leftbreak1 = it->break1 - SPLIT_BUFF;
       int rightbreak2= it->break2 + SPLIT_BUFF;
       int leftbreak2 = it->break2 - SPLIT_BUFF;
       bool issplit1 = (leftend <= leftbreak1) && (rightend >= rightbreak1);
       bool issplit2 = (leftend <= leftbreak2) && (rightend >= rightbreak2);

       //bool issplit1 = (j->pos <= (it->break1 - SPLIT_BUFF)) && ((j->pos + j->seq.length() >= (it->break1 + SPLIT_BUFF))
       //bool issplit1 = (j->pos >= (it->break1 - j->seq.length() + SPLIT_BUFF) ) && (j->pos <= it->break1 + j->seq.length() - SPLIT_BUFF);
       //bool issplit2 = (j->pos >= (it->break2 - j->seq.length() + SPLIT_BUFF) ) && (j->pos <= it->break2 + j->seq.length() - SPLIT_BUFF);

       if (issplit1 && tumor_read)
	 it->tsplit1++;
       if (issplit1 && !tumor_read)
	 it->nsplit1++;
       if (issplit2 && tumor_read)
	 it->tsplit2++;
       if (issplit2 && !tumor_read)
	 it->nsplit2++;
       
     }
   }

}

// determine whether this AlignedContig intersects with another AlignmentContig
bool AlignedContig::intersect(const AlignedContig *al) const {
  
  AlignVec av = al->getAlignments();

  for (AlignVec::const_iterator it = m_align.begin(); it != m_align.end(); it++)
    for (AlignVec::const_iterator jt = av.begin(); jt != av.end(); jt++)
      if (jt->align.RefID == it->align.RefID)
	if (abs(it->align.Position - jt->align.Position) < 2000)
	  return true;
  return false;

}

void AlignedContig::printAlignments(ofstream &ostream) const {

  if (m_align.size() == 0) {
    cerr << "in printAlignments: no alignments" << endl;
    return;
  }
   // print the contig
   //ostream << m_align[0].align.Name << endl;

  // print the local breakpoints only if there are 3+ alignments
  //if (m_breaks.size() > 2) {
  //  for (BPVec::const_iterator it = m_breaks.begin(); it != m_breaks.end(); it++)
  //    ostream << "Local:  " << it->toString() << endl;
  //}

  ostream << "Global: " << m_farbreak.toString() << endl; 

   // print the BWA alignments
   for (AlignVec::const_iterator i = m_align.begin(); i != m_align.end(); i++) {

     // sets the direction to print
     char jsign = '>'; 
     if (i->align.IsReverseStrand())
       jsign = '<';

     // print the cigar value per base
     for (CigarOpVec::const_iterator j = i->cigar.begin(); j != i->cigar.end(); j++) {
       if (j->Type == 'M')
	 ostream << string(j->Length, jsign);
       else if (j->Type == 'I') 
	 ostream << string(j->Length, 'I');
       else if (j->Type == 'S' || j->Type == 'H')
	 ostream << string(j->Length, '.');
     }

     // print the info
     ostream << "    " << i->align.RefID + 1 << ":" << i->align.Position 
	       << " MAPQ: " << i->align.MapQuality << " OrientedCigar: " << i->cigstring
               << " Start: " << i->start << " has_local: " << m_local
               << " Breaks: " << i->break1 << " " << i->break2 
               << " GenomeBreaks: " << i->gbreak1 << " " << i->gbreak2 
               << " Tsplits: " << i->tsplit1 << " " << i->tsplit2 
               << " Nsplits: " << i->nsplit1 << " " << i->nsplit2 
	     << " DiscordantCluster: " << printDiscordantClusters() << endl;
   }

   // print the masked contig base-pairs
   ostream << m_masked_seq << "    RepeatMasked contig" << endl;
   for (RepeatMaskerVec::const_iterator it = m_rep_vec.begin(); it != m_rep_vec.end(); it++) {
     int len = max(it->query_end - it->query_begin + 1, 0);
     int pre = max(it->query_begin - 1, 0);
     int post= min((int)m_seq.length() - len - pre, 1000);
     ostream << string(pre, '.') << string(len, '*') << string(post, '.') << 
       "    Repeat Element: " << it->repeat_name << " Repeat Class: " << it->repeat_class << " Percent Divergence: " << it->perc_div << " SW Score: " << it->sw << endl;
   }

   // print the contig base-pairs
   ostream << m_seq << "    " << m_align[0].align.Name << endl; 

   // print the individual reads
   /*
   for (R2CVec::const_iterator i = m_reads.begin(); i != m_reads.end(); i++) {
     int padlen = m_seq.size() - i->pos - i->seq.size() + 5;
     padlen = max(5, padlen);

     // print it out if it is reasonable. If not, probably a bug. Protect against massive output
     if (i->pos > 1e4 || padlen > 1e4) { // big problem, don't print
       cout << "Overflow for contig: " << m_align[0].align.Name << " Padlen: " << padlen << " pos: " << i->pos << endl;
       return;
     } else {
       ostream << string(i->pos, ' ') << i->seq << string(padlen, ' ') << i->rname << " " << i->a.RefID+1 << ":" << i->a.Position << " SW: " << i->sw_score << endl;
     }

   }
   */

   for (vector<BamAlignment>::const_iterator i = m_bamreads.begin(); i != m_bamreads.end();  i++) {

     string seq, jw;
     int pos, sw;
     i->GetTag("TS", seq);
     i->GetTag("AL", pos);
     i->GetTag("SW", sw);
     i->GetTag("JW", jw);
     
     int padlen = m_seq.size() - pos - seq.size() + 5;
     padlen = max(5, padlen);

     // print it out if it is reasonable. If not, probably a bug. Protect against massive output
     if (pos > 1e4 || padlen > 1e4) { // big problem, don't print
       cout << "Overflow for contig: " << m_align[0].align.Name << " Padlen: " << padlen << " pos: " << pos << endl;
       return;
     } else {
       ostream << string(pos, ' ') << seq << string(padlen, ' ') << jw << "--" << i->Name << " " << i->RefID+1 << ":" << i->Position  << " SW: " << sw << endl;
     }

   }

}

// flips the cigar if the contig is aligned to the opposite strand
CigarOpVec AlignedContig::orientCigar(const BamTools::BamAlignment align) {

    CigarOpVec cig; 
    
    // reverse the cigar if necesssary
    if (align.IsReverseStrand()) 
      for (CigarOpVec::const_iterator it = align.CigarData.end() - 1; it != align.CigarData.begin() - 1; it--) 
        cig.push_back(*it);
    else 
      cig = align.CigarData;
 
    return cig;
}

// converts the BamTools cigar format to string
string AlignedContig::cigarToString(CigarOpVec cig) {
    stringstream cigstring;
    for (CigarOpVec::const_iterator it = cig.begin(); it != cig.end(); it++) 
      cigstring << it->Length << it->Type;
    return cigstring.str();
}

// parses the contig file name to determine where the anchor window was
/*void AlignedContig::setWindow(const string s) {
  GenomicRegion gr(s);
  
  
   // set the chromosome
   int chrom = -1;
   try {
     regex reg("^c_([0-9XYM]+):*");
     smatch match;
     if (regex_search(s.begin(), s.end(), match, reg))
       chrom = stoi(match[1].str()) - 1;
   } catch (...) {
     cerr << "Caught error with contig: " << s << endl;
     chrom = 0;
   }

   // set pos1
  int pos1 = -1; 
   try { 
     regex creg("^c_[0-9XYM]+:([0-9]+)-*");
     smatch cmatch;
     if (regex_search(s.begin(), s.end(), cmatch, creg))
       pos1 = stoi(cmatch[1].str());
   } catch (...) {
     cerr << "Caught error with contig: " << s << endl;
     pos1 = 1;
   }

   // set pos2
   int pos2 = -1;
   try {
     regex creg2("^c_[0-9XYM]+:[0-9]+-([0-9]+).*");
     smatch cmatch2;
     if (regex_search(s.begin(), s.end(), cmatch2, creg2))
       pos2 = stoi(cmatch2[1].str());
   } catch (...) {
     cerr << "Caught error with contig: " << s << endl;
     pos2 = 1;
   }
     
   Window w(chrom, pos1, pos2);
   m_window = w;


   }*/

// find the breakpoint pairs by looping through ALL the alignments
void AlignedContig::getBreakPairs() {
   
   if (m_align.size() < 2)
     return;

   //int tall = getNumReads('t');
   //int nall = getNumReads('n');

   // walk along the ordered contig list and make the breakpoint pairs
   BreakPoint bp;
   //bp.tall = tall;
   //bp.nall = nall;
   bp.seq = m_seq;
   bp.num_align = m_align.size();
   bp.cname = m_align[0].align.Name;

   if (bp.cname=="") {
     cerr << "Contig name is empty" << endl;
     exit(EXIT_FAILURE); 
   }

   bp.window = m_window;
   bp.part_of_local = m_local; // set whether breakpoint is part of contig with local

   // set the discovar support if it's there
   string s;
   bool valid = m_align[0].align.GetTag("TN", s);
   if (valid) {
     
     bp.discovar = true;

     // set the tumor support
     std::regex reg("^t([0-9]+)_*");
     std::smatch match;
     //if (std::regex_search(s.begin(), s.end(), match, reg))
     if (std::regex_search(s, match, reg))
       bp.disco_tum = std::stoi(match[1].str());

     // set the normal support
     std::regex reg2("^t[0-9]+n([0-9]+)");
     std::smatch match2;
     if (std::regex_search(s, match2, reg2))
       bp.disco_norm = std::stoi(match2[1].str());

   }

   for (AlignVec::iterator it = m_align.begin(); it != m_align.end() - 1; it++) {

     // set the strand
     //bp.strand1 = it->align.IsReverseStrand() ? '-' : '+';
     //bp.strand2 = it->align.IsReverseStrand() ? '+' : '-';

     bp.cpos1 = it->break2; // take the right-most breakpoint as the first
     bp.pos1  = it->gbreak2; 
     bp.refID1 = it->align.RefID;
     bp.cpos2 = (it+1)->break1;  // take the left-most of the next one
     bp.pos2  = (it+1)->gbreak1; 
     bp.refID2 = (it+1)->align.RefID;

     if (bp.refID1 == bp.refID2)
       bp.span = abs((int)bp.pos1-(int)bp.pos2);
     else
       bp.span = -1;
    
     // set the splits
     bp.nsplit1 = it->nsplit2;
     bp.tsplit1 = it->tsplit2;
     bp.nsplit2 = (it+1)->nsplit1;
     bp.tsplit2 = (it+1)->tsplit1;

     bp.nsplit = min(bp.nsplit1, bp.nsplit2);
     bp.tsplit = min(bp.tsplit1, bp.tsplit2);

     // set the mapq
     bp.mapq1 = it->align.MapQuality;
     bp.mapq2 = (it+1)->align.MapQuality;

    if (bp.mapq1 > 60 || bp.mapq2 > 60) {
       cerr << "bad mapq" << endl;
       cerr << bp.toString() << endl;
       exit(EXIT_FAILURE); 
    }

     // set the local
     bp.local1 = it->ca_local;
     bp.local2 = (it+1)->ca_local;

     // set the match length
     for (CigarOpVec::const_iterator cc = it->align.CigarData.begin(); cc != it->align.CigarData.end(); cc++)
       if (cc->Type == 'M')
	 bp.matchlen1 += cc->Length;
     for (CigarOpVec::const_iterator cc = (it+1)->align.CigarData.begin(); cc != (it+1)->align.CigarData.end(); cc++)
       if (cc->Type == 'M')
	 bp.matchlen2 += cc->Length;
     
     // set the NM
     int nmtag;
     it->align.GetTag("NM", nmtag);
     bp.nm1 = nmtag;
     (it+1)->align.GetTag("NM", nmtag);
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

     // send the strandedness 
     bp.strand1 =     it->align.IsReverseStrand() ? '-' : '+';
     bp.strand2 = (it+1)->align.IsReverseStrand() ? '+' : '-';
     //bp.strand2 = it->align.IsReverseStrand() ? '+' : '-';

     // make sure it is in order
     bp.order();

     m_breaks.push_back(bp);   

   }

   // skip the middle contigs, to make a master breakpoint pair (for 3+ mappings)
   if (m_align.size() < 3) {
     m_farbreak = bp;
     m_farbreak_filt = bp;
     return;
   }

   // go through alignments and find start and end that reach mapq 
   size_t bstart = 1000; //1000 is a dummy
   size_t bend = m_align.size() - 1;
   for (size_t i = 0; i < m_align.size(); i++)
     if (m_align[i].align.MapQuality >= mapq_threshold) {
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
   m_farbreak = bp;
   m_farbreak.cpos1 = m_align[bstart].break2; // first mapping
   m_farbreak.pos1  = m_align[bstart].gbreak2;
   m_farbreak.cpos2 = m_align[bend].break1; // last mapping
   m_farbreak.pos2  = m_align[bend].gbreak1;
   m_farbreak.refID1 = m_align[bstart].align.RefID; 
   m_farbreak.refID2 = m_align[bend].align.RefID; 
   if (m_farbreak.refID1 == m_farbreak.refID2)
     m_farbreak.span = abs((int)m_farbreak.pos1-(int)m_farbreak.pos2);
   else
     m_farbreak.span = -1;

   // set the strands
   m_farbreak.strand1 = m_align[bstart].align.IsReverseStrand() ? '-' : '+';
   m_farbreak.strand2 = m_align[bend].align.IsReverseStrand()   ? '+' : '-';

   // set the splits
   m_farbreak.nsplit1 = m_align[bstart].nsplit2;
   m_farbreak.tsplit1 = m_align[bstart].tsplit2;
   m_farbreak.nsplit2 = m_align[bend].nsplit1;
   m_farbreak.tsplit2 = m_align[bend].tsplit1;

   m_farbreak.nsplit = min(m_farbreak.nsplit1, m_farbreak.nsplit2);
   m_farbreak.tsplit = min(m_farbreak.tsplit1, m_farbreak.tsplit2);

   // set the mapping quality
   m_farbreak.mapq1 = m_align[bstart].align.MapQuality;
   m_farbreak.mapq2 = m_align[bend].align.MapQuality;

   if (m_farbreak.mapq1 > 60 || m_farbreak.mapq2 > 60) {
     cerr << "bad mapq GLOBAL" << endl;
     cerr << m_farbreak.toString() << endl;
     cerr << " m_align size " << m_align.size() << " bstart " << bstart << " bend " << bend << endl;
     exit(EXIT_FAILURE); 
   }


   // set the homologies
   try {
     if (m_farbreak.cpos1 >= m_farbreak.cpos2)
       m_farbreak.homology = m_seq.substr(m_farbreak.cpos2, m_farbreak.cpos1-m_farbreak.cpos2);
     else if (m_farbreak.cpos2 >= m_farbreak.cpos1)
       m_farbreak.insertion = m_seq.substr(m_farbreak.cpos1, m_farbreak.cpos2-m_farbreak.cpos1);
     if (m_farbreak.insertion.length() == 0)
       m_farbreak.insertion = "N";
     if (m_farbreak.homology.length() == 0)
       m_farbreak.homology = "N";
   m_farbreak.order();
   } catch (...) {
       cout << "Caught error with contig on global-getBreakPairs: " << m_align[0].align.Name << endl;
   }
 
 }
 
/*void fillExtraReads() {

  vector<BamAlignment> tbav, nbav, pbav, bav, tbav_d, nbav_d, pbav_d;
  int isize = 1000;
  int mapq = 0;
  int qualthresh = 8;
  int minOverlap = 32;

  SVBamReader t_reader(tumor_bam, "te", isize, mapq, qualthresh, minOverlap, true, false, 2);
  t_reader.findBamIndex();
  SVBamReader n_reader(tumor_bam, "ne", isize, mapq, qualthresh, minOverlap, true, false, 2);
  n_reader.findBamIndex();
  SVBamReader p_reader(tumor_bam, "pe", isize, mapq, qualthresh, minOverlap, true, false, 2);
  p_reader.findBamIndex();

  
  for (AlignVec::const_iterator it = m_align.begin(); it != m_align.end(); it ++) {

    int refID = it->align.RefID;
    int pos1 = it->align.Position - 700;
    int pos2 = it->align.Position + 700

    if (t_reader.setBamRegion(refID, refID, pos1, pos2)) 
      t_reader.bamToBAVec(tbav);
    if (n_reader.setBamRegion(refID, refID, pos1, pos2)) 
      n_reader.bamToBAVec(nbav);
    if (p_reader.setBamRegion(refID, refID, pos1, pos2)) 
      p_reader.bamToBAVec(pbav);

  }
  
  SVBamReader::deduplicateReads(tbav, tbav_d);
  SVBamReader::deduplicateReads(nbav, nbav_d);
  SVBamReader::deduplicateReads(pbav, pbav_d);
  tbav.clear();
  nbav.clear();
  pbav.clear();
  SVBamReader::deduplicateReadsPos(tbav_d, tbav);
  SVBamReader::deduplicateReadsPos(nbav_d, nbav);
  SVBamReader::deduplicateReadsPos(pbav_d, pbav);

  m_reads.clear();
  for (vector<BamAlignment>::const_iterator it = tbav.begin(); it != tbav.end(); it++) {
    string tag;
    it->GetTag("JW", tag);
    string ts;
    it->GetTag("TS", ts);
    Read2Contig rc(m_align[0].align.Name, tag, 0, *it, ts);
  }
  for (vector<BamAlignment>::const_iterator it = nbav.begin(); it != nbav.end(); it++) {
    string tag;
    it->GetTag("JW", tag);
    string ts;
    it->GetTag("TS", ts);
    Read2Contig rc(m_align[0].align.Name, tag, 0, *it, ts);
  }
  for (vector<BamAlignment>::const_iterator it = pbav.begin(); it != pbav.end(); it++) {
    string tag;
    it->GetTag("JW", tag);
    string ts;
    it->GetTag("TS", ts);
    Read2Contig rc(m_align[0].align.Name, tag, 0, *it, ts);
  }

  }*/

// doing a more stringent alignment of reads to the contig
// this is to remove normals that ruin somatic calls
void AlignedContig::realignReads() { 

   vector<BamAlignment> rr;
   // declare the readers
   int isize = 1000;
   int mapq = 0;
   int qualthresh = 4;
   int minOverlap = 35;
   bool skip_supp = false;
   bool skip_r2 = false;
   int min_clip = 8;
   int verbose = 0;

   SVBamReader treader = SVBamReader(tbam,  "tb", isize, mapq, qualthresh, minOverlap, skip_supp, skip_r2, min_clip, verbose);
   SVBamReader nreader = SVBamReader(nbam,  "nb", isize, mapq, qualthresh, minOverlap, skip_supp, skip_r2, min_clip, verbose);
   SVBamReader preader = SVBamReader(pbam,  "pb", isize, mapq, qualthresh, minOverlap, skip_supp, skip_r2, min_clip, verbose);
   SVBamReader rreader = SVBamReader(rbam,  "pb", 1, 0, qualthresh, minOverlap, skip_supp, skip_r2, min_clip, verbose);
   treader.clipOnly = true;
   nreader.clipOnly = true;
   preader.clipOnly = true;
   rreader.clipOnly = true;

   // read in the reads
   bool tgood = treader.getBamFilename().length() > 1; // && false;
   bool ngood = nreader.getBamFilename().length() > 1; // && false;
   bool pgood = preader.getBamFilename().length() > 1; // && false;
   bool rgood = rreader.getBamFilename().length() > 1;;

   size_t tvecsize = 0;

   // set the limit on number of reads to read in
   treader.setReadLimit(5000);
   nreader.setReadLimit(5000);
   preader.setReadLimit(5000);
   rreader.setReadLimit(5000);

   if (tgood) 
     if (!treader.findBamIndex())
       cerr << "Failed to open BAM index in Tumor" << endl;
   if (ngood) 
     if (!nreader.findBamIndex())
       cerr << "Failed to open BAM index in Normal" << endl;
   if (pgood) 
     if (!preader.findBamIndex())
       cerr << "Failed to open BAM index in Panel" << endl;
   if (rgood) 
     if (!rreader.findBamIndex())
       cerr << "Failed to open BAM index in R2C BAM" << endl;

   for (AlignVec::const_iterator it = m_align.begin(); it != m_align.end(); it++) {

     BamAlignmentVector this_vec;
     
     GenomicRegion rg(it->align.RefID, it->align.Position, it->align.Position + it->align.Length);

     // set the BAM region
     if (tgood) {
       if (!treader.setBamRegion(rg))
	 cerr << "Failed to set BAM position in Tumor" << endl;
       if (!treader.bamToBAVec(this_vec)) {
	 //cerr << "Failed to get BAM reads for tumor" << endl;
       }
     }

     // set the BAM region
     if (ngood) {
       if (!nreader.setBamRegion(rg))
	 cerr << "Failed to set BAM position in Normal" << endl;
       if (!nreader.bamToBAVec(this_vec)) {
	 //cerr << "Failed to get BAM reads for Normal" << endl;
       }
     }

     // set the BAM region
     if (pgood) {
       if (!preader.setBamRegion(rg))
	 cerr << "Failed to set BAM position in Panel" << endl;
       if (!preader.bamToBAVec(this_vec)) {
	 //cerr << "Failed to get BAM reads for Panel" << endl;
       }
     }

     // set the BAM region
     if (rgood) {
       if (!rreader.setBamRegion(rg))
	 cerr << "Failed to set BAM position in R2C BAM" << endl;
       if (!rreader.bamToBAVec(this_vec)) {
	 cerr << "Failed to get BAM reads for R2C BAM" << endl;
       }
     }
     
     tvecsize += this_vec.size();

     // MATCHING BY FIND
     int buff = 12;
     int pad = 10;

    for (BamAlignmentVector::iterator j = this_vec.begin(); j != this_vec.end(); j++) {
      
      string QB;
      
      if (!j->GetTag("TS", QB))
	QB = j->QueryBases;
      int seqlen = QB.length();
      size_t posa = m_seq.find(QB.substr(pad, buff)); // try first part of read
      size_t posb = m_seq.find(QB.substr(max(seqlen-buff-pad,0),buff)); // try second part of read
      bool hit1 = posa != string::npos;
      bool hit2 = posb != string::npos;
      string read_name;
      //j->GetTag("JW", read_name);
      
      // PROCEED IF ALIGNS TO FORWARD
      if (hit1 || hit2) {
	rr.push_back(*j);
	//int tpos = posa - pad;
	//i->addRead(*j, std::max(tpos, 0), true);
	//name_map.insert(pair<string, unsigned>(read_name, 0));
      }
      
      //OTHERWISE TRY REVERSE
      
      else {
	string rstring = QB;
	SnowUtils::rcomplement(rstring); 
	posa = m_seq.find(rstring.substr(pad,buff)); 
	posb = m_seq.find(rstring.substr(max(seqlen-buff-pad,0),buff)); 
	hit1 = posa != string::npos;
	hit2 = posb != string::npos;
	
	if (hit1 || hit2) {
	  j->EditTag("TS", "Z", rstring); // edit the tag to be reverseComplemented
	  rr.push_back(*j);
	}
      }
      ////////////////////////////
    } // end read loop
   
    } // end alignment loop

    if (rr.size() > 10000) {
      cout << "Greater than 10k reads at contig: " << m_align[0].align.Name << " with # reads: " << rr.size() << endl;
     return;
   }
   
   vector<BamAlignment> rr_name_dd, rr_pos_dd;
   SVBamReader::deduplicateReads(rr, rr_name_dd);
   SVBamReader::deduplicateReadsPos(rr_name_dd, rr_pos_dd);
   rr = rr_pos_dd;

   double mseqlen = 0;

   // set the contig strings    
   TSequence contig = m_seq;

   // set the parametersa
   int lencutoff = floor( mseqlen * 0.5 );

   // loop through the reads and do alignment
   vector<BamAlignment>::iterator it = rr.begin();
   for (; it != rr.end(); it++) {
     string read_seq;
     it->GetTag("TS", read_seq);
     double seqlen = static_cast<double>(read_seq.length());
     double cutoff = 4 * seqlen - 25; 

     if (seqlen > lencutoff) {

       // check for a 100% match, to save compute
       int32_t pos = m_seq.find(read_seq);
       int32_t score = seqlen * 4;

       // matched completely, add
       if (pos != string::npos) {
	 it->AddTag("AL", "i", pos);
	 it->AddTag("CN", "Z", m_align[0].align.Name);
	 it->AddTag("SW", "i", score);
	 m_bamreads.push_back(*it);
       // didn't match completely, SW align
       } else if (SeqanTools::SWalign(contig, pos, read_seq, score, false) > cutoff) {
	 it->AddTag("AL", "i", pos);
         it->AddTag("CN", "Z", m_align[0].align.Name);
	 it->AddTag("SW", "i", score);
	 m_bamreads.push_back(*it);
       }
	 // forwards SW didn't make it, try reverse
	 /*else if (SWalign(contig, true, pos, read_seq, score) > cutoff) {
	 it->AddTag("AL", "i", pos);
	 it->AddTag("CN", "Z", m_align[0].align.Name);
	 it->EditTag("TS", "Z", read_seq);
	 it->AddTag("SW", "i", score);
	 m_bamreads.push_back(*it);

	 }*/

     }
   
   }

}

// print out the alignments to R
string AlignedContig::printForR() const {

  stringstream ss;
  //for (R2CVec::const_iterator it = m_reads.begin(); it != m_reads.end(); it++) 
  //  ss << it->toString() << endl;
  return ss.str();

}

// read in an R2C bam with TS, CN, AL and SW tags already supplied
void AlignedContig::readR2Creads() {

  if (rbam == "")
    return;

   vector<BamAlignment> rr;

   SVBamReader rreader = SVBamReader(rbam);

   if (!rreader.findBamIndex())
     cerr << "Failed to open BAM index in R2C BAM" << endl;

   //debug this
   for (AlignVec::const_iterator it = m_align.begin(); it != m_align.end(); it++) {
     if (!rreader.setBamRegion(it->align.RefID, it->align.Position - 500, it->align.Position + it->align.Length + 500))
       cerr << "Failed to set BAM position in R2C BAM" << endl;
     if (!rreader.R2CbamToBAVec(rr)) 
       cerr << "Failed to get BAM reads for R2C BAM" << endl;
     
   }
  
   vector<BamAlignment> rr_name_dd, rr_pos_dd;
   SVBamReader::deduplicateReads(rr, rr_name_dd);
   SVBamReader::deduplicateReadsPos(rr_name_dd, rr_pos_dd);
   m_bamreads = rr_pos_dd;
   
}

// print the discordant clusters that align to this contig
string AlignedContig::printDiscordantClusters() const {

  stringstream out;
  if (m_dc.size() == 0)
    return "No Discorant Clusters";

  for (vector<DiscordantCluster>::const_iterator it = m_dc.begin(); it != m_dc.end(); it++)
    out << *it << " ";
  return out.str();

}

/*void AlignedContig::discordantNormalCheck() {

  SVBamReader nreader = SVBamReader(nbam); 
  if (!nreader.findBamIndex())
    cerr << "Failed to open BAM index in Normal" << endl;

  for (AlignVec::const_iterator it = m_align.begin(); it != m_align.end(); it++) { 
    if (ndisc > 0)
      break;
    GenomicRegion gr(it->align.RefID, it->align.Position, it->align.Position);
    gr.pad(600);
    int span = (it->align.RefID != it->align.)
    nreader.setBamRegion(gr);
    ndisc += nreader.discordantCount();
  }

  }*/
