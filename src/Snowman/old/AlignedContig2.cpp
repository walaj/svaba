#include "AlignedContig2.h"
#include <boost/algorithm/string.hpp>
#include <regex>
#include "SeqanTools.h"

using namespace std;

/** 
 * Constructor from a set of SAM records
 */
AlignedContig2::AlignedContig2(const string& sam, const BamTools::BamReader *treader) {

  istringstream iss(sam);
  string val, line;

  int i = 0;
  string cigar;

  // loop through each of the sam alignments
  while (getline(iss, line, '\n')) {

    istringstream issv(line);

    BamTools::BamAlignment a;
    
    // loop through each of the fields of the SAM line
    while(getline(issv, val, '\t')) {

      switch(i) {
        case 0 : a.Name = val; break;
        case 1 : a.AlignmentFlag = stoi(val); break;
        case 2 : try {a.RefID = (val == "*") ? -1 : treader->GetReferenceID(val);} catch (...) { a.RefID = -1; } break;
        case 3 : a.Position = stoi(val); break;
        case 4 : a.MapQuality = stoi(val); break;
        case 5 : cigar = val; break; // 6-8 are Mate info related
        case 9 : a.QueryBases = val; break;
        case 10 : (val == "*") ? a.Qualities = string(a.QueryBases.length(), 'I') : a.Qualities = val; break;
      }

      // set the sequence, always on positive strand
      if (seq.length() == 0 && i == 9) {
	seq = val;
	if (a.IsReverseStrand())
	  SnowUtils::rcomplement(seq);
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
    aln.push_back(a);
    i = 0;

  }

  setBreaks();
  
}

/** 
 * Parse a cigar string into a vector<CigarOp>
 *
 * @param val CIGAR string to be parsed
 * @return The parsed CIGAR string
 */
CigarOpVec AlignedContig2::stringToCigar(const string& val) {
  
  vector<string> str_vec; // #2: Search for tokens
  boost::split(str_vec, val, boost::is_any_of("0123456789"), boost::token_compress_on ); // SplitVec == { "hello abc","ABC","aBc goodbye" }
  
  vector<string> len_vec; // #2: Search for tokens
  boost::split(len_vec, val, boost::is_any_of("MIDSHPN"), boost::token_compress_on ); // SplitVec == { "hello abc","ABC","aBc goodbye" }
  
  assert(len_vec.size() == str_vec.size());
  
  CigarOpVec cigop;
  for (size_t kk = 0; kk < (str_vec.size() - 1); kk++) // first strvec is empty and last len_vec is empty (due to token orderingin cigar)
    cigop.push_back(CigarOp(str_vec[kk+1].at(0), stoi(len_vec[kk])));
  
  return cigop;
  
}

/**
 * Parse tag data from a SAM record, add to BamAlignment
 *
 * @param val SAM field to parse
 * @param a BamAlignment to add the tag to
 */
void AlignedContig2::parseTags(const string& val, BamAlignment &a) {

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


/**
 * Aligned sequencing reads to a contig
 *
 * This function performs two types of alignment: Smith-Waterman
 * and exact matching. Exact matching is performed first (for efficiency). 
 * If this fails, checks that at least 8 bases in front or back of read
 * are on the contig. If so, does full SW, forward and if need be, reverse.
 * This function will add AL, CN, SW tags to reads which align.
 *
 * @param bav A vector of pointers to sequencing reads that are to be aligned. 
 */
void AlignedContig2::alignReadsToContigs(BamAlignmentUPVector &bav) {

  // MATCHING BY FIND
  //int buff = 12;
  //int pad = 10;
  TSequence contig = seq;

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
	
	size_t pos = seq.find(QB);
	int32_t aligned_pos;
	
	bool addread = false;
	bool isrev = false;
	
	// make some substrings to match. Dont proceed if none of these match anywhere
	string sub1 = QB.substr(10, 6); // 10-16
	string sub2 = QB.substr(seqlen-20, 6); // (end-20) - (end-14)

	// matched completely, add
	if (pos != string::npos) {
	  aligned_pos = pos;
	  addread = true;
	  // didn't match completely, SW align
	} else if (seq.find(sub1) != string::npos || seq.find(sub2) != string::npos) {
	  if (SeqanTools::SWalign(contig, aligned_pos, QB, score, cutoff, false))
	    addread = true;
	  else {
	    SnowUtils::rcomplement(sub1);
	    SnowUtils::rcomplement(sub2);
	  }
	}
	// forwards SW didn't make it, try reverse
	else if (!addread && (seq.find(sub1) != string::npos || seq.find(sub2) != string::npos)) {
	  if (SeqanTools::SWalign(contig, aligned_pos, QB, score, cutoff, true)) {
	    addread = true;
	    isrev = true;
	  }
	}
	
	// add some tags. remove others
	if (addread) {
	  
	  SnowUtils::SmartAddTag(j, "AL", to_string(aligned_pos));
	  SnowUtils::SmartAddTag(j, "CN", name);
	  SnowUtils::SmartAddTag(j, "SW", to_string(score));
	  
	  if (isrev)
	    j->EditTag("TS", "Z", QB); // stores reverse comp if need be

	  reads.push_back(j); // make a copy of the data

	  // add the coverage
	  vector<size_t> tmp;
	  if (cov.count(short_name)) 
	    tmp = cov[short_name];
	  else 
	    tmp = vector<size_t> (seq.length(), 0);
	  for (size_t kk = aligned_pos; kk <= min(aligned_pos + QB.length(), seq.length() - 1); kk++) 
	    tmp[kk]++;
	  cov[short_name] = tmp;


	}
      } // end if > 20
	
  } // end read loop
  
}

void AlignedContig2::setBreaks() {

  assert(aln.size());

  // deal with multi-mappings differently
  if (aln.size() > 1) {

    vector<ContigFragment> cfv;

    for (auto& i : aln) {
      //
      CigarOpVec cig = orientCigar(i);
      unsigned currlen  = 0; 
      unsigned gcurrlen = 0;

      int b1 = -1;
      int b2 = -1;
      int g1 = -1;
      int g2 = -1;

      for (auto& c : cig) { 

	// SET THE CONTIG BREAK (treats deletions and leading S differently)
	// the first M gets the break1, pos on the left
	if (c.Type == 'M' && b1 ==  -1)
	  b1 = currlen;
	if (c.Type != 'D') // skip deletions but not leading S, but otherwise update
	  currlen += c.Length;
	if (c.Type == 'M') // keeps triggering every M, with pos at the right
	  b2 = currlen;
	
	// SET THE GENOME BREAK (treats insertions differently)
	if (c.Type != 'I' && c.Type != 'S' && c.Type != 'H') 
	  gcurrlen += c.Length;
      }
      
      // convert genome break in contig coords to reference coords
      if (i.IsReverseStrand()) {
        g1 = i.Position + gcurrlen;
	g2 = i.Position + 1; //gcurrlen + align.align.Position - align.gbreak2;
      } else {
	g1  = i.Position + 1;
	g2  = i.Position + gcurrlen;
      }

      ContigFragment cf(b1, b2, g1, g2);
      cfv.push_back(cf);
    }
   
    // sort the contig fragment vector so left most fragment is left most on genome
    sort(cfv.begin(), cfv.end());
 
    for (auto c = cfv.begin(); c != cfv.end(); c++) {
      Break br();
    }
    //splitCoverage(br, )
    
    
  }

}

CigarOpVec AlignedContig2::orientCigar(const BamTools::BamAlignment& align) {

    CigarOpVec cig; 
    
    // reverse the cigar if necesssary
    if (align.IsReverseStrand()) 
      for (CigarOpVec::const_iterator it = align.CigarData.end() - 1; it != align.CigarData.begin() - 1; it--) 
        cig.push_back(*it);
    else 
      cig = align.CigarData;
 
    return cig;
}
