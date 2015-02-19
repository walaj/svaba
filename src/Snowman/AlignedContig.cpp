#include "AlignedContig.h"
#include <regex>
#include "GenomicRegion.h"
#include <unordered_map>
//#include "SVBamReader.h"
#include "api/algorithms/Sort.h"
#include "SnowUtils.h"
#include "SeqanTools.h"
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace BamTools;

#define SPLIT_BUFF 8

// add a new contig alignment
void AlignedContig::addAlignment(const BamTools::BamAlignment align) { 

  CAlignment tal(align, orientCigar(align));
  m_align.push_back(tal);
  
  if (m_align.size() > 1)
    sort(m_align.begin(), m_align.end());

}

void AlignedContig::setIndelBreaks(CAlignment &ca) {

  assert(m_align.size());

  // reject some
  bool reject = true;
  
  // reject if multiple D / I
  for (auto& i : ca.align.CigarData) {
    if (!reject && i.Type=='D' || i.Type == 'I') { // TODO support mulitple indels
      reject = true;
      break;
    }
    if (i.Type == 'D' || i.Type == 'I') {
      reject = false;
    }
  }

  // reject if not in window
  GenomicRegion tmpw = window;
  tmpw.pad(1000);
  if (tmpw.getOverlap(GenomicRegion(ca.align.RefID, ca.align.Position, ca.align.Position)) == 0)
    reject = true;

  if (reject)
    return;

  hasvariant = true;
  int curr = 0;
  int gcurrlen = 0;
  char strand1 = '*', strand2 = '*';
  string insertion = "";
  for (auto& i : ca.align.CigarData) {

    // set the contig breakpoint
    if (i.Type == 'M' || i.Type == 'S' || i.Type == 'I') {
      curr += i.Length;
    } 
    if (i.Type == 'D' && ca.break1 == -1) {
      ca.break1 = curr;
      ca.break2 = curr + 1; //i.Length;
    } 
    if (i.Type == 'I' && ca.break1 == -1) {
      ca.break1 = curr - i.Length - 1;
      ca.break2 = curr - 1;
      insertion = ca.align.QueryBases.substr(curr, i.Length);
    }
   
    // set the genome length diferently, skip insertions
    if (i.Type != 'I' && i.Type != 'S' && i.Type != 'H') 
      gcurrlen += i.Length;
    if (ca.gbreak1 == -1 && ca.break1 > 0) {
      /* if (ca.align.IsReverseStrand() && i.Type == 'D') {
	ca.gbreak1 = ca.align.Position + ca.align.Length - ca.break1;
	ca.gbreak2 = ca.align.Position + ca.align.Length - ca.break1 - i.Length;
	strand1 = '-';
	strand2 = '+';
	} else */ if (/*!ca.align.IsReverseStrand() && */i.Type == 'D') {
	ca.gbreak1 =  ca.break1 + ca.align.Position;
	ca.gbreak2 = ca.gbreak1 + i.Length;
	strand1 = '+';
	strand2 = '-';
      } /*else if (ca.align.IsReverseStrand() && i.Type == 'I') {
	ca.gbreak1 = ca.align.Position + ca.align.Length - ca.break1;
	ca.gbreak2 = ca.gbreak1 - 1; 
	strand1 = '-';
	strand2 = '+';	
	} */ else if (/*!ca.align.IsReverseStrand() && */i.Type == 'I') {
	ca.gbreak1 = ca.break1 + ca.align.Position;
	ca.gbreak2 = ca.gbreak1 + 1;	
	strand1 = '+';
	strand2 = '-';
      }
    }


    /*     cout << i.Length << i.Type << " gcurrlen " << gcurrlen << " curr " << curr << 
      " ca.break1 " << ca.break1 << " ca.break2 " << ca.break2 << " ca.gbreak1 " << 
       ca.gbreak1 << " ca.gbreak2 " << ca.gbreak2 <<" cigar " << ca.cigstring << endl;
    */
      }

//   cout << "D_GBREAK1: " << ca.dgbreak1 << " D_GBREAK2: " << ca.dgbreak2 
//        << "D_BREAK1: " << ca.dbreak1 << " D_BREAK2: " << ca.dbreak2 << 
//     " gcurrlen: " << gcurrlen << " IsReverseStrand: " << ca.ca.IsReverseStrand() << 
//     " cigar: " << ca.cigstring << " Position: " << ca.align.Position << endl;
  
 
  // immediately set the global breakpoint
  m_farbreak.gr1 = GenomicRegion(ca.align.RefID, ca.gbreak1, ca.gbreak1, strand1, ca.align.MapQuality);
  m_farbreak.gr2 = GenomicRegion(ca.align.RefID, ca.gbreak2, ca.gbreak2, strand2, ca.align.MapQuality);
  m_farbreak.cpos1 = ca.break1;
  m_farbreak.cpos2 = ca.break2;
  m_farbreak.span = abs(ca.gbreak2 - ca.gbreak1);
  m_farbreak.insertion = insertion;
  m_farbreak.seq = m_seq;
  m_farbreak.num_align = m_align.size();
  m_farbreak.cname = m_align[0].align.Name;

  m_farbreak.order();

}

// make work function for getting PER-ALIGNMENT breaks from BWA-MEM
void AlignedContig::setBreaks(CAlignment &align) { 

  hasvariant = true;
  unsigned currlen  = 0; 
  unsigned gcurrlen = 0;

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
    trig = j->Type == 'M' || trig;
    if (j->Type != 'I' && j->Type != 'S' && j->Type != 'H' && trig) 
      gcurrlen += j->Length;
  }
  
  // convert genome break in contig coords to reference coords
  if (align.align.IsReverseStrand()) {
    align.gbreak1 = align.align.Position + gcurrlen;
    align.gbreak2 = align.align.Position + 1; //gcurrlen + align.align.Position - align.gbreak2;
  } else {
    align.gbreak1  = align.align.Position + 1;
    align.gbreak2  = align.align.Position + gcurrlen;
  }
    
  /*If (align.align.Name.compare("contig_8:94302502-94310502_1")==0)
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
  //  std::sort( m_bamreads.begin(), m_bamreads.end(), ByALTag());
}

void AlignedContig::splitCoverage() { 
  
   for (AlignVec::iterator it = m_align.begin(); it != m_align.end(); it++) {

     it->nsplit1 = 0;
     it->tsplit1 = 0;
     it->nsplit2 = 0;
     it->tsplit2 = 0;

     int rightbreak1= it->break1 + SPLIT_BUFF;
     int leftbreak1 = it->break1 - SPLIT_BUFF;
     int rightbreak2= it->break2 + SPLIT_BUFF;
     int leftbreak2 = it->break2 - SPLIT_BUFF;

     for (auto& j : m_bamreads) {

       string jw;
       assert(j->GetTag("SR", jw));
       bool tumor_read = (jw.at(0) == 't');
       string seq;
       if (!j->GetTag("TS",seq))
	 seq = j->QueryBases;
       assert(seq.length() > 0);
       //assert(j->GetTag("TS", seq));

       int pos = SnowUtils::GetIntTag(j, "AL").back();

       int rightend = pos + seq.length();
       int leftend  = pos;
       bool issplit1 = (leftend <= leftbreak1) && (rightend >= rightbreak1);
       bool issplit2 = (leftend <= leftbreak2) && (rightend >= rightbreak2);

       // add the reads
       if (issplit1)
	 it->reads_b1.push_back(j);
       if (issplit2)
	 it->reads_b2.push_back(j);

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

ostream& operator<<(ostream& out, const AlignedContig &ac) {
  out << ac.printAlignments();
  return out;
}

string AlignedContig::printAlignments() const {

  if (m_align.size() == 0) {
    cerr << "in printAlignments: no alignments" << endl;
    return "";
  }

  stringstream ostream;
  
  if (!m_farbreak.isEmpty())
    ostream << "Global: " << m_farbreak << endl; 

   // print the BWA alignments
  for (auto& i : m_align) 
    ostream << i << endl;

   // print the masked contig base-pairs
   if (m_masked_seq.length() > 0)
     ostream << m_masked_seq << "    RepeatMasked contig" << endl;
   for (RepeatMaskerVec::const_iterator it = m_rep_vec.begin(); it != m_rep_vec.end(); it++) {
     int len = max(it->query_end - it->query_begin + 1, 0);
     int pre = max(it->query_begin - 1, 0);
     int post= min((int)m_seq.length() - len - pre, 1000);
     ostream << string(pre, '.') << string(len, '*') << string(post, '.') << 
       "    Repeat Element: " << it->repeat_name << " Repeat Class: " << it->repeat_class << " Percent Divergence: " << it->perc_div << " SW Score: " << it->sw << endl;
   }

   // print the break locations
   if ( (m_align[0].break1 >= 0 && m_align[0].break1 < 3000) && (m_align[0].break2 >= 0 && m_align[0].break2 < 3000) && m_align.size() == 1) {
     ostream << string(m_align[0].break1, ' ') << "|" << string(m_align[0].break2-m_align[0].break1, ' ') << '|' << endl;
   }

   // print the contig base-pairs
   ostream << m_seq << "    " << m_align[0].align.Name << endl; 

   // print out the individual reads
   for (auto& i : m_bamreads) {

     string seq, sr;
     int pos, sw;
     if (!i->GetTag("TS", seq))
       seq = i->QueryBases;
     assert(i->GetTag("SR",sr));

     // get the more complex tags (since there can be multiple annotations per tag)
     vector<int> posvec = SnowUtils::GetIntTag(i, "AL");
     vector<int> swvec = SnowUtils::GetIntTag(i, "SW");
     vector<string> cnvec = SnowUtils::GetStringTag(i, "CN");
     assert(posvec.size() == swvec.size());
     assert(cnvec.size() == posvec.size());
     size_t kk = 0;
     for (; kk < cnvec.size(); kk++) 
       if (cnvec[kk] == m_align[0].align.Name) {
	 pos = posvec[kk];
	 sw = swvec[kk];
	 break;
       }
     assert(kk != cnvec.size()); // assure that we found something
     pos = abs(pos);
     int padlen = m_seq.size() - pos - seq.size() + 5;
     padlen = max(5, padlen);

     // print it out if it is reasonable. If not, probably a bug. Protect against massive output
     if (pos > 1e4 || padlen > 1e4) { // big problem, don't print
       cout << "Overflow for contig: " << m_align[0].align.Name << " Padlen: " << padlen << " pos: " << pos << endl;
       return ostream.str();
     } else {
       ostream << string(pos, ' ') << seq << string(padlen, ' ') << sr << "--" << i->RefID+1 << ":" << i->Position  << " SW: " << sw << endl;
     }

   }

   return ostream.str();

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
string CAlignment::cigarToString(CigarOpVec cig) {
    stringstream cigstring;
    for (CigarOpVec::const_iterator it = cig.begin(); it != cig.end(); it++) 
      cigstring << it->Length << it->Type;
    return cigstring.str();
}

// find the breakpoint pairs by looping through ALL the alignments
void AlignedContig::getBreakPairs() {
   
  if (m_align.size() < 2) {
    if (m_align.size() == 1) {
      m_farbreak.reads.insert(m_farbreak.reads.begin(), m_align.back().reads_b1.begin(), m_align.back().reads_b1.end());
      m_farbreak.reads.insert(m_farbreak.reads.begin(), m_align.back().reads_b2.begin(), m_align.back().reads_b2.end());
    }
    return; 
  }

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

   //bp.window = m_window;
   //bp.part_of_local = m_local; // set whether breakpoint is part of contig with local

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

     bp.gr1 = GenomicRegion(it->align.RefID, it->gbreak2, it->gbreak2);
     bp.gr2 = GenomicRegion((it+1)->align.RefID, (it+1)->gbreak1, (it+1)->gbreak1);
     bp.gr1.strand = it->align.IsReverseStrand() ? '-' : '+';
     bp.gr2.strand = (it+1)->align.IsReverseStrand() ? '+' : '-';

     bp.cpos1 = it->break2; // take the right-most breakpoint as the first
     //bp.pos1  = it->gbreak2; 
     //bp.refID1 = it->align.RefID;
     bp.cpos2 = (it+1)->break1;  // take the left-most of the next one
     //bp.pos2  = (it+1)->gbreak1; 
     //bp.refID2 = (it+1)->align.RefID;
     bp.reads.insert(bp.reads.end(), it->reads_b2.begin(), it->reads_b2.end());
     bp.reads.insert(bp.reads.end(), (it+1)->reads_b1.begin(), (it+1)->reads_b1.end());

     if (bp.gr1.chr == bp.gr2.chr)
       bp.span = abs(bp.gr1.pos1-bp.gr2.pos1);
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
     bp.gr1.mapq = it->align.MapQuality;
     bp.gr2.mapq = (it+1)->align.MapQuality;

    if (bp.gr1.mapq > 60 || bp.gr2.mapq > 60) {
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
     bp.gr1.nm = nmtag;
     (it+1)->align.GetTag("NM", nmtag);
     bp.gr2.nm = nmtag;

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
     //bp.strand1 =     it->align.IsReverseStrand() ? '-' : '+';
     //bp.strand2 = (it+1)->align.IsReverseStrand() ? '+' : '-';
     //bp.strand2 = it->align.IsReverseStrand() ? '+' : '-';

     // make sure it is in order
     //bp.order();

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
   //m_farbreak.pos1  = m_align[bstart].gbreak2;
   m_farbreak.cpos2 = m_align[bend].break1; // last mapping
   //m_farbreak.pos2  = m_align[bend].gbreak1;
   //m_farbreak.refID1 = m_align[bstart].align.RefID; 
   //m_farbreak.refID2 = m_align[bend].align.RefID; 
   if (m_farbreak.gr1.chr == m_farbreak.gr2.chr)
     m_farbreak.span = abs(m_farbreak.gr1.pos1-m_farbreak.gr2.pos1);
   else
     m_farbreak.span = -1;

   // set the strands
   m_farbreak.gr1.strand = m_align[bstart].align.IsReverseStrand() ? '-' : '+';
   m_farbreak.gr2.strand = m_align[bend].align.IsReverseStrand()   ? '+' : '-';

   // set the splits
   m_farbreak.nsplit1 = m_align[bstart].nsplit2;
   m_farbreak.tsplit1 = m_align[bstart].tsplit2;
   m_farbreak.nsplit2 = m_align[bend].nsplit1;
   m_farbreak.tsplit2 = m_align[bend].tsplit1;

   m_farbreak.nsplit = min(m_farbreak.nsplit1, m_farbreak.nsplit2);
   m_farbreak.tsplit = min(m_farbreak.tsplit1, m_farbreak.tsplit2);

   // set the mapping quality
   m_farbreak.gr1.mapq = m_align[bstart].align.MapQuality;
   m_farbreak.gr2.mapq = m_align[bend].align.MapQuality;

   if (m_farbreak.gr1.mapq > 60 || m_farbreak.gr2.mapq > 60) {
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
   //m_farbreak.order();
   } catch (...) {
       cout << "Caught error with contig on global-getBreakPairs: " << m_align[0].align.Name << endl;
   }
 
 }
 

// doing a more stringent alignment of reads to the contig
// this is to remove normals that ruin somatic calls
/*
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
	 else if (SWalign(contig, true, pos, read_seq, score) > cutoff) {
	 it->AddTag("AL", "i", pos);
	 it->AddTag("CN", "Z", m_align[0].align.Name);
	 it->EditTag("TS", "Z", read_seq);
	 it->AddTag("SW", "i", score);
	 m_bamreads.push_back(*it);

	 }

     }
   
   }

}

*/

// print out the alignments to R
string AlignedContig::printForR() const {

  stringstream ss;
  //for (R2CVec::const_iterator it = m_reads.begin(); it != m_reads.end(); it++) 
  //  ss << it->toString() << endl;
  return ss.str();

}

// print the discordant clusters that align to this contig
string AlignedContig::printDiscordantClusters() const {

  stringstream out;
  if (m_dc.size() == 0)
    return "No Discordant Clusters";

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

// make an aligned contig from a sam record
AlignedContig::AlignedContig(string sam, const BamReader * reader, GenomicRegion &twindow) {
  
  window = twindow;
  
  m_farbreak.gr1 = GenomicRegion(-1,-1,-1); // set a dummy

  samrecord = sam;

  std::istringstream iss(sam);
  std::string val, line;

  int i = 0;
  string cigar;

  std::regex reg_xp("^XA:Z:(.*)");
  std::regex reg_nm("^NM:[A-Za-z]:(.*)");

  while (getline(iss, line, '\n')) {
    std::istringstream issv(line);

    BamTools::BamAlignment a;
    while(getline(issv, val, '\t')) {

      switch(i) {
      case 0 : a.Name = val; break;
      case 1 : a.AlignmentFlag = std::stoi(val); break;
      case 2 : try {a.RefID = (val == "*") ? -1 : reader->GetReferenceID(val);} catch (...) { a.RefID = -1; } break;
      case 3 : a.Position = std::stoi(val); break;
      case 4 : a.MapQuality = std::stoi(val); break;
      case 5 : cigar = val; break;
	//case 7 : 
      case 9 : a.QueryBases = val; break;
      case 10 : (val == "*") ? a.Qualities = string(a.QueryBases.length(), 'I') : a.Qualities = val; break;
      }

      // set the sequence, always on positive strand
      if (m_seq.length() == 0 && i == 9) {
	m_seq = val;
	//if (a.IsReverseStrand())
	//  SnowUtils::rcomplement(m_seq);
      }

      // deal with the cigar
      if (i == 5 && val != "*") 
	a.CigarData = stringToCigar(val);

      // process the tags
      if (i > 10) 
	parseTags(val, a);

      i++;
	
    }

    assert(a.MapQuality <= 60);
    addAlignment(a);
    aln.push_back(a);
    i = 0;
    //break; // only do the first record of multi-part
  }

  assert(aln.size());

  // set the breaks
  if (m_align.size() > 1) {
    for (auto& i : m_align)
      setBreaks(i);
  } else {
    setIndelBreaks(m_align[0]);
  }

}

// align sequencing reads to the contig using SW alignment
void AlignedContig::alignReadsToContigs(BamAlignmentUPVector &bav) {

  hasvariant = true;
  // MATCHING BY FIND
  int buff = 1;
  //int pad = 10;
  TSequence contig = m_seq;

  // 
  for (auto& j : bav) {
      
      string QB;
      
      if (!j->GetTag("TS", QB))
	QB = j->QueryBases;

      int seqlen = QB.length();
      int32_t score = seqlen * 4;;
      int cutoff = score - 25;

      if (seqlen > 20) {

	string read_name;
	j->GetTag("SR",read_name);
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
	  // didn't match completely, SW align
	} else if (m_seq.find(sub1) != string::npos || m_seq.find(sub2) != string::npos) {
	  if (SeqanTools::SWalign(contig, aligned_pos, QB, score, cutoff, false))
	    addread = true;
	  else {
	    SnowUtils::rcomplement(sub1);
	    SnowUtils::rcomplement(sub2);
	  }
	}
	// forwards SW didn't make it, try reverse
	else if (!addread && (m_seq.find(sub1) != string::npos || m_seq.find(sub2) != string::npos)) {
	  if (SeqanTools::SWalign(contig, aligned_pos, QB, score, cutoff, true)) {
	    addread = true;
	    isrev = true;
	  }
	}
	
	// add some tags. remove others
	if (addread) {
	  
	  SnowUtils::SmartAddTag(j, "AL", to_string(aligned_pos));
	  SnowUtils::SmartAddTag(j, "CN", m_align[0].align.Name);
	  SnowUtils::SmartAddTag(j, "SW", to_string(score));
	  
	  if (isrev)
	    j->EditTag("TS", "Z", QB); // stores reverse comp if need be

	  
	  m_bamreads.push_back(j); // make a copy of the data
	}
      } // end if > 20
	
  } // end read loop

  // 
  sortReads();

  // if its an indel call, add reads directly
  // otherwise, get added in getBreakPairs
  if (m_align.size() == 1) {
    for (auto& i : m_bamreads) {
      m_farbreak.reads.push_back(i);
    }
  }
  
}

int AlignedContig::maxSplit() const {
  
  int out = max(max(m_farbreak.tsplit1, m_farbreak.tsplit2), max(m_farbreak.nsplit1, m_farbreak.nsplit2));
  return out;
  
}

/** 
 * Parse a cigar string into a vector<CigarOp>
 *
 * @param val CIGAR string to be parsed
 * @return The parsed CIGAR string
 */
CigarOpVec AlignedContig::stringToCigar(const string& val) {
  
  string v = val;
  vector<string> str_vec; // #2: Search for tokens
  boost::split(str_vec, v, boost::is_any_of("0123456789"), boost::token_compress_on ); // SplitVec == { "hello abc","ABC","aBc goodbye" }
  
  vector<string> len_vec; // #2: Search for tokens
  boost::split(len_vec, v, boost::is_any_of("MIDSHPN"), boost::token_compress_on ); // SplitVec == { "hello abc","ABC","aBc goodbye" }
  
  assert(len_vec.size() == str_vec.size());
  assert(len_vec.size());

  CigarOpVec cigop;
  for (size_t kk = 0; kk < (str_vec.size() - 1); kk++) // first strvec is empty and last len_vec is empty (due to token orderingin cigar)
    cigop.push_back(CigarOp(str_vec[kk+1].at(0), stoi(len_vec[kk])));

  assert(cigop.size());

  return cigop;
  
}

/**
 * Parse tag data from a SAM record, add to BamAlignment
 *
 * @param val SAM field to parse
 * @param a BamAlignment to add the tag to
 */
void AlignedContig::parseTags(const string& val, BamAlignment &a) {

  regex reg_xp("^XA:Z:(.*)");
  regex reg_nm("^NM:[A-Za-z]:(.*)");

  smatch match;
  if (regex_search(val, match, reg_xp)) {
    string xp = match[1].str();
    a.AddTag("XP","Z",xp);
  } else if (regex_search(val, match, reg_nm)) {
    int nm = stoi(match[1].str());
    a.AddTag("NM","i",nm);
  }

}



// CONSTRUCTOR
CAlignment::CAlignment(BamTools::BamAlignment talign, CigarOpVec tcigar) : align(talign), cigar(tcigar) {
  cigstring = cigarToString(tcigar);
  
  // find the start position of alignment ON CONTIG
  start = 0; 
  for (auto& i : cigar) {
    if (i.Type != 'M')
      start += i.Length;
    else
      break;
  }
}

ostream& operator<<(ostream &out, const CAlignment &c) {
  
  // sets the direction to print
  char jsign = '>'; 
  if (c.align.IsReverseStrand())
    jsign = '<';
  
  // print the cigar value per base
  for (auto& j : c.align.CigarData) { // print releative to forward strand
    if (j.Type == 'M')
      out << string(j.Length, jsign);
    else if (j.Type == 'I') 
      out << string(j.Length, 'I');
    else if (j.Type == 'S' || j.Type == 'H')
      out << string(j.Length, '.');
  }
  
  // print the info
  out << "    " << c.align.RefID + 1 << ":" << c.align.Position 
	  << " MAPQ: " << c.align.MapQuality << " OrientedCigar: " << c.cigstring
	  << " Start: " << c.start
	  << " Breaks: " << c.break1 << " " << c.break2 
	  << " GenomeBreaks: " << c.gbreak1 << " " << c.gbreak2 
	  << " Tsplits: " << c.tsplit1 << " " << c.tsplit2 
	  << " Nsplits: " << c.nsplit1 << " " << c.nsplit2 
      << " DiscordantCluster: "; // << printDiscordantClusters();
  
  return out;
}


void AlignedContig::splitIndelCoverage() { 
  
  int rightbreak1= m_farbreak.cpos1 + SPLIT_BUFF;
  int leftbreak1 = m_farbreak.cpos1 - SPLIT_BUFF;
  int rightbreak2= m_farbreak.cpos2 + SPLIT_BUFF;
  int leftbreak2 = m_farbreak.cpos2 - SPLIT_BUFF;
  
  for (auto& j : m_bamreads) {
    
    string sr;
    assert(j->GetTag("SR", sr));
    assert(sr.at(0) == 't' || sr.at(0) == 'n');
    bool tumor_read = (sr.at(0) == 't');
    string seq;
    if (!j->GetTag("TS",seq))
      seq = j->QueryBases;
    assert(seq.length() > 0);
    //assert(j->GetTag("TS", seq));
    
    int pos = SnowUtils::GetIntTag(j, "AL").back();
    
    int rightend = pos + seq.length();
    int leftend  = pos;
    bool issplit1 = (leftend <= leftbreak1) && (rightend >= rightbreak1);
    bool issplit2 = (leftend <= leftbreak2) && (rightend >= rightbreak2);
    
    // add the reads
    if (issplit1 || issplit2)
      m_farbreak.reads.push_back(j);
    if ( (issplit1 || issplit2) && tumor_read)
      m_farbreak.tsplit++;
    if ( (issplit1 || issplit2) && !tumor_read)
      m_farbreak.nsplit++;

    //    cout << "issplit1 " << issplit1 << " issplit2 "  << issplit2 << " leftedn " << leftend << " rightend " << rightend << " leftbrea1k " << leftbreak1 << " leftbreak2 " << leftbreak2 << " rightbreak1 " << rightbreak1 << " righrbreak2 " << rightbreak2 << endl;
    // cout << "checkng coverage tslpit " << m_farbreak.tsplit << " (isspli1 || issplit2) " << (issplit1 || issplit2) << " tumor_read " << tumor_read <<  endl;
    //cout << "checkng coverage nslpit " << m_farbreak.nsplit << endl;
  }

}
