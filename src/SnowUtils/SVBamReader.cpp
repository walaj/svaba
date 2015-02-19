#include "SVBamReader.h"
#include "unistd.h"
#include "api/algorithms/Sort.h"
#include <time.h> // for now
#include "SnowUtils.h"
#include <set>

#include "ClusterReads.h"

//#define DEBUG_SVREADS 2

#define MIN_READ_LENGTH

// Phred score transformations
inline int char2phred(char b) {
  uint8_t v = b;
  assert(v >= 33);
  return v - 33;
}

inline std::string string_numf(int n, int len)
{
  std::string result(len--, '0');
  for (int val=(n<0)?-n:n; len>=0&&val!=0; --len,val/=10)
    result[len]='0'+val%10;
  if (len>=0&&n<0) result[0]='-';
  return result;
}

// Perform a soft-clipping of the sequence by removing low quality bases from the
// 3' end using Heng Li's algorithm from bwa
void SVBamReader::softClip(int qualTrim, std::string &seq, std::string const &qual, unsigned clipnum, 
			   bool &tooshort, bool &not_real_clip, int &startpoint, int &endpoint)
{
    assert(seq.size() == qual.size());

    endpoint = seq.length()-1;
    startpoint = 0;
    int i = 0; 
 
    // get the start point (loop forward)
    while(i < (int)seq.length()) {
        int ps = char2phred(qual[i]);
        if (ps >= qualTrim) {
          startpoint = i;
          break;
	}
	i++;
    }

    // get the end point (loop backwards)
    i = seq.length() - 1;
    while(i >= 0) {
        int ps = char2phred(qual[i]);
        if (ps >= qualTrim) {
          endpoint = i;
          break;
	}
	i--;
    }

    // check that they aren't all bad
    if (startpoint == 0 && endpoint == (int)seq.length() - 1 && char2phred(qual[0]) < qualTrim) {
      seq = "";
      tooshort = true;
      startpoint = 0;
      endpoint = -1;
      //qual = "";
      return;
    }

    int readlen = seq.length();

    // Clip the read
    seq = seq.substr(startpoint, endpoint);
    //qual = qual.substr(startpoint, endpoint);

    int new_clipnum = clipnum - (readlen - seq.length() - 1);
    
    // assign the values
    tooshort = seq.length() < std::floor(m_mol * 1.5);
    int int_minclip = m_minclip; // get rid of int compare warning
    not_real_clip = new_clipnum < int_minclip && int_minclip != 0;
    
    
    //bool lowq = seq.length() < m_mol || (new_clipnum < 10); 
    //std::cout << "NEW CLIP: " << new_clipnum << " SEQLENGTH: " << seq.length() << std::endl;
    //if (lowq) { // num clipped that aren't removed must be bigger than 10
    //  low_qual++;
    //  return false; // don't include read
    //} else {
    //  return true;
    //}
}

// Perform a soft-clipping of the sequence by removing low quality bases from the
// 3' end using Heng Li's algorithm from bwa
void SVBamReader::qualityTrimRead(int qualTrim, std::string &seq, std::string &qual) {

    assert(seq.size() == qual.size());

    int endpoint = -1; //seq.length();
    int startpoint = 0;
    int i = 0; 
 
    // get the start point (loop forward)
    while(i < (int)seq.length()) {
        int ps = char2phred(qual[i]);
        if (ps >= qualTrim) {
          startpoint = i;
          break;
	}
	i++;
    }

    // get the end point (loop backwards)
    i = seq.length() - 1;
    while(i >= 0) {
        int ps = char2phred(qual[i]);
        if (ps >= qualTrim) {
          endpoint = i + 1; // endpoint is one past edge
          break;
	}
	i--;
    }
    // check that they aren't all bad
    if (startpoint == 0 && endpoint == -1) {
      seq = "";
      qual = "";
      //tooshort = true;
      //startpoint = 0;
      //endpoint = -1;
      //qual = "";
      return;
    }

    int readlen = seq.length();

    // Clip the read
    seq =   seq.substr(startpoint, endpoint);
    qual = qual.substr(startpoint, endpoint);
    return;
    
    //int new_clipnum = clipnum - (readlen - seq.length() - 1);
    
    // assign the values
    //tooshort = seq.length() < std::floor(m_mol * 1.5);
    //int int_minclip = m_minclip; // get rid of int compare warning
    //not_real_clip = new_clipnum < int_minclip && int_minclip != 0;
    
    
    //bool lowq = seq.length() < m_mol || (new_clipnum < 10); 
    //std::cout << "NEW CLIP: " << new_clipnum << " SEQLENGTH: " << seq.length() << std::endl;
    //if (lowq) { // num clipped that aren't removed must be bigger than 10
    //  low_qual++;
    //  return false; // don't include read
    //} else {
    //  return true;
    //}
}


bool SVBamReader::bamToBAVec(BamAlignmentVector &bav) { 

  // setup a map of reads to make sure we aren't duplicating
  if (!m_reader.Open(m_bam)) {
    std::cerr << "FAILED TO OPEN BAM: " << m_bam << std::endl;
    std::cerr << m_bam << std::endl;
    return false;
  }

  if (m_bai.length() == 0) {
    std::cerr << "BAM index not set. Run reader.findBamIndex()" << endl;
    return false;
  }

  if (!m_reader.OpenIndex(m_bai)) {
    if (!m_reader.OpenIndex(m_bam_bai)) {
      std::cerr << "FAILED TO OPEN INDEX: " << m_bam_bai << endl;
      return false;
    }
  }

  if (!m_reader.SetRegion(m_region.chr, m_region.pos1, m_region.chr, m_region.pos2)) {
    std::cerr << "FAILED TO SET REGION\n"; 
    return false;
  }

  unsigned counter = 0;

  BamTools::BamAlignment a;

  // loop through the reads
  while (m_reader.GetNextAlignmentCore(a)) {
 
    total_reads++;

    // get the number of soft-clipped bases BEFORE TRIMMING
    unsigned clipnum = getClipCount(a);
    
    // get some other info about the alignment before filling string fields
    //bool single_align = a.MateRefID == -1;
    bool unmap = !a.IsMapped() && !disc_cluster_only;
    bool discr = (m_isize > 0) && (std::abs(a.InsertSize) > m_isize || a.RefID != a.MateRefID); 
    discr = discr && !clipOnly;
    //discr = discr || single_align;
    bool mapqr = (a.MapQuality >= m_mapq  || unmap) && !a.IsFailedQC() && !a.IsDuplicate();
    bool supp  = a.IsPrimaryAlignment() || !m_skip_supp; // if m_skip_supp is true, only keep primary alignments
    mapqr = mapqr && supp;

    if ( (discr || unmap || (clipnum >= m_minclip || m_minclip == 0)) && mapqr ) {

      a.BuildCharData(); //populate the rest of the string data

      // clip bases with low optical quality
      std::string trimmed_bases = a.QueryBases;
      std::string trimmed_quals = a.Qualities;
      bool tooshort = false;
      bool not_real_clip = false;
      int startpoint;
      int endpoint;
      
      // only need to clip if SE and SA tags not there
      int32_t se, ss;
      if (!a.GetTag("SS", ss) || !a.GetTag("SE", se)) {
	softClip(m_qualthresh, trimmed_bases, trimmed_quals, clipnum, tooshort, not_real_clip, startpoint, endpoint);
	a.AddTag("SS","i",ss);
	a.AddTag("SE","i",se);
	//a.Qualities = "";
      } else {
	startpoint = ss;
	endpoint = se;
	trimmed_bases = a.QueryBases.substr(ss, se - ss + 1);
	tooshort = (se - ss) < 60; 
      }
	
      // check that there are no N bases
      bool npass = passNtest(a.QueryBases);
      
      // get the number of matches
      uint32_t nm;
      if (a.GetTag("NM", nm)) {} else { nm = 0; }
      
      // update the counters and determine if clipped/disc/mapq/unmap pass
      // check what properties read satisfies 
      bool nmpass = nm <= m_nmlim;
      
      bool clipp = clipnum >= m_minclip && !not_real_clip; // make sure raw clip is sufficient, plus processed
      clipp = clipp || m_minclip == 0;
      clipp = clipp && !disc_cluster_only;
      // update the counters
      if (discr && !clipp && mapqr)
	disc_num++;
      if (unmap)
	unmap_num++;
      if (clipp && !discr && mapqr)
	clip_num++;
      if (clipp && discr && mapqr)
      disc_clip_num++;
      if (mapqr)
	mapq_num++;
      
      // pass CLIPPED, DISC or UNMAP reads, provided they satisfy mapqr and are not too short
      bool cont = (discr || unmap || clipp) && mapqr && !tooshort && npass && nmpass;
      
      if ( cont ) {

	// add the trimmed sequences to the read
	a.AddTag("TS", "Z", trimmed_bases);
	
	// want to index by position of anchor read
	std::string hp = string_numf(a.RefID, 2) + "_" + string_numf(a.Position, 9);
	a.AddTag("HP", "Z", hp);
	
	//set the prefix for this read name 
	std::ostringstream s_first;
	//s << m_prefix << "_" << a.AlignmentFlag << "_" << a.Name;
	s_first << m_prefix << a.AlignmentFlag << "_" << a.Name;
	a.AddTag("JW", "Z", s_first.str());
	
	std::string r2, q2, trim_seq;
	bool r2exists = a.GetTag("R2", r2) && !m_skip_r2;
	bool has_good_r2 = false;
	
	if (r2exists && !disc_cluster_only) {
	  
	  //std::string trimmed_bases = a.QueryBases;
	  //std::string trimmed_quals = a.Qualities;
	  trim_seq = r2; // assign the trimmed sequence the original
	  has_good_r2 = GetR2Read(a, r2, q2, trim_seq); //trimmed_bases, trimmed_quals);
	  
	  if (has_good_r2) {

	    BamTools::BamAlignment ba_r2;
	    ba_r2.Name = a.Name;
	    ba_r2.Length = r2.length();   
	    ba_r2.QueryBases = r2;
	    ba_r2.Qualities = q2;
	    ba_r2.InsertSize = -a.InsertSize;
	    BamTools::CigarOp cig_op('M', r2.length());
	    std::vector<BamTools::CigarOp> tmp_vec;
	    tmp_vec.push_back(cig_op);
	    ba_r2.CigarData = tmp_vec;
	    
	    ba_r2.RefID = a.MateRefID;
	    ba_r2.MateRefID = a.RefID;
	    
	    // new way with HP tag
	    ba_r2.Position = a.MatePosition;
	    ba_r2.MatePosition = a.Position;          
	    
	    // want to index by position of anchor read
	    ba_r2.AddTag("HP", "Z", hp);
	    
	    // add the trimmed sequences to the read
	    ba_r2.AddTag("TS", "Z", trim_seq);
	    
	    // add a flag to hold the anchor position
	    ba_r2.AddTag<int32_t>("IR", "i", 1);      // set a tag to let know is R2
	    
	    // set the flag
	    ba_r2.SetIsFirstMate(a.IsSecondMate());
	    ba_r2.SetIsSecondMate(a.IsFirstMate());
	    ba_r2.SetIsMateMapped(true);
	    ba_r2.SetIsMateReverseStrand(a.IsReverseStrand());
	    ba_r2.SetIsPaired(true);
	    ba_r2.SetIsPrimaryAlignment(true);
	    ba_r2.SetIsReverseStrand(a.IsMateReverseStrand());
	    ba_r2.SetIsFailedQC(false);     
	    ba_r2.SetIsMapped(true);
	    ba_r2.SetIsDuplicate(false);
	    ba_r2.SetIsProperPair(false);
	    
	    std::ostringstream s_second;
	    s_second << m_prefix <<  ba_r2.AlignmentFlag << "_R2_" << a.Name;
	    ba_r2.AddTag("JW", "Z", s_second.str());
	    
	    // add the pairmate name tag
	    ba_r2.AddTag("J2", "Z", s_first.str());
	    a.AddTag("J2", "Z", s_second.str());
	    
	    bav.push_back(ba_r2);
	    used_num++;
	  } // end has good r2
	  
	} // end if R2 tag
	
	// add the anchor read if R2 not there, or R2 is good, or set to skip R2
	if (!r2exists || (r2exists && has_good_r2) || disc_cluster_only) {
	  used_num++;

	  // remove tags for memory saving
	  a.RemoveTag("Q2"); 
	  a.RemoveTag("R2"); 
	  a.RemoveTag("XM"); 
	  a.RemoveTag("SA"); 
	  a.RemoveTag("XS"); 
	  a.RemoveTag("RG"); 
	  a.RemoveTag("AS"); 
	  a.RemoveTag("XS");
	  a.RemoveTag("OQ");

	  //a.Qualities=""; // extreme
	  //a.QueryBases = ""; // extreme	
	  bav.push_back(a);
	  counter++;

	  if (counter > limit) { // hit the read limit
	    bav.clear();
	    //cout << "Limit of " << limit << " hit on region: " << 
	    //  m_region.toStringOffset()  << 
	    //  " in " << m_bam << endl;
	    return true;
	  }

	}
	
      } // end the if (cont)
      
    } // end buildCharData
    
  } // end the while loop

  bool dothesort = true;
  if (m_prefix.size() > 1)
    if (m_prefix.at(1) == 'd')
      dothesort = false;

  // ensure that it is sorted by HP
  if (dothesort) // don't need to sort if is discordant lookup
    std::sort( bav.begin(), bav.end(), BamTools::Algorithms::Sort::ByTag<std::string>("HP", BamTools::Algorithms::Sort::AscendingOrder) );
    
  if (m_verbose > 1) 
    std::cout << m_prefix << " Total Reads: " << total_reads << " Used reads: " << bav.size() << 
      "\tUnmapped: " << unmap_num << "\tDiscOnly: " << disc_num << "\tClipNum: " << clip_num << 
      "\tDisc-Clip Num: " << disc_clip_num << "\tMapQ pass: " << mapq_num << 
      "\tLow-Qual Remove: " << low_qual << "\tNon-standard Chr Remove: " << non_num << std::endl;
  
  return true;
  
}

bool SVBamReader::findBamIndex() {

  m_bai = m_bam;
  m_bai = m_bai.substr(0, m_bam.size()-3).append("bai");
  m_bam_bai = m_bam;
  m_bam_bai = m_bam_bai.append(".bai");
  return true;

}

bool SVBamReader::setBamRegion(int refid, int pos1, int pos2) {

  m_region = GenomicRegion(refid, pos1, pos2);

  return true;

}

bool SVBamReader::setBamRegion(GenomicRegion gp) {

  m_region = gp;
  
  return true;

}

// count the number of soft clips
// check that the cigar is good
/*unsigned SVBamReader::getClipCount(BamTools::BamAlignment a) const {

      std::vector<int> clipSize;
      std::vector<int> readPos;
      std::vector<int> genPos;
      a.GetSoftClips(clipSize, readPos, genPos, false);
  
      // get the clip number
      unsigned clipnum = 0;
      for(std::vector<int>::iterator j=clipSize.begin();j!=clipSize.end();++j)
	clipnum += *j;

      return clipnum;
}
*/

bool SVBamReader::updateCounters(BamTools::BamAlignment a, unsigned clipnum) {
 
    // check what properties read satisfies 
    bool discr = m_isize > 0 && (std::abs(a.InsertSize) > m_isize || a.RefID != a.MateRefID); 
    bool unmap = !a.IsMapped();
    bool clipp = clipnum >= m_minclip || m_minclip == 0;
    bool mapqr = a.MapQuality >= m_mapq && !a.IsFailedQC() && !a.IsDuplicate();
    
    //debug
    //if (a.Name.compare("C2MEGACXX131207:2:1316:2411:64457")==0) 
    //  std::cerr << "DISCR: " << discr << " UNMAP: " << unmap << " CLIPP: " << clipnum << " mapqr: " << mapqr << " FLAG: " << a.AlignmentFlag << " ISIZE: " << a.InsertSize << 
    //	" m_isize: " << m_isize << std::endl;

    // update the counters
    if (discr && !clipp && mapqr)
      disc_num++;
    if (unmap)
      unmap_num++;
    if (clipp && !discr && mapqr)
      clip_num++;
    if (clipp && discr && mapqr)
      disc_clip_num++;
    if (mapqr)
      mapq_num++;

    bool output = ( (discr || unmap || clipp)) && mapqr;
    return output;
 
}

bool SVBamReader::passNtest(std::string str) const {

  for (size_t i = 0; i < str.size(); i++)
    if (str[i] == 'N')
      return false;
  return true;

}

bool SVBamReader::GetR2Read(BamTools::BamAlignment const &a, std::string const &r2, std::string &q2, std::string &trim_seq) {

  //if (!a.GetTag("R2", r2)) 
  //  return false;

   total_reads++;

   // get the cont
   if (!a.GetTag("Q2", q2))
     q2 = std::string(r2.length(), 'I'); // fill with dummy, being optimistic for sensitivity

   // clip bases with low optical quality
   bool tooshort = false;
   bool not_real_clip = false; // doesn't matter here really, because it's already discordant, so doesn't matter if clipped
   int startpoint;
   int endpoint;
   softClip(m_qualthresh, trim_seq, q2, trim_seq.length(), tooshort, not_real_clip, startpoint, endpoint); // pretend all bases clipped, so accepts it

   // check if it is a non_standard chromosom
   bool non_stan = a.MateRefID > 24;
   if (non_stan)
     non_num++;

   // discard if there are Ns in the sequence
   bool npassr2 = passNtest(r2);

   return npassr2 && !tooshort && !non_stan;

}

bool SVBamReader::contigBamToBAVec(BamAlignmentVector &bav) {
  
    // setup a map of reads to make sure we aren't duplicating
  if (!m_reader.Open(m_bam)) {
    std::cerr << "FAILED TO OPEN BAM: " << m_bam << std::endl;
    std::cerr << m_bam << std::endl;
    return false;
  }

  if (m_bai.length() == 0) {
    std::cerr << "BAM index not set. Run reader.findBamIndex()" << endl;
    return false;
  }

  if (!m_reader.OpenIndex(m_bai)) {
    if (!m_reader.OpenIndex(m_bam_bai)) {
      std::cerr << "FAILED TO OPEN INDEX: " << m_bam_bai << endl;
      return false;
    }
  }

  if (!m_reader.SetRegion(m_region.chr, m_region.pos1, m_region.chr, m_region.pos2)) {
    std::cerr << "FAILED TO SET REGION\n"; 
    return false;
  }

  BamTools::BamAlignment a;

  // loop through the reads
  while (m_reader.GetNextAlignment(a)) {
    total_reads++;
    std::string hp = string_numf(a.RefID, 2) + "_" + string_numf(a.Position, 9);
    a.AddTag("HP", "Z", hp);

    // set the prefix
    std::ostringstream s_first;
    s_first << m_prefix << a.AlignmentFlag << "_" << a.Name;
    a.AddTag("JW", "Z", s_first.str());

    // add it
    bav.push_back(a);
  }

  // ensure that it is sorted by HP
  std::sort( bav.begin(), bav.end(), BamTools::Algorithms::Sort::ByTag<std::string>("HP", BamTools::Algorithms::Sort::AscendingOrder) );

  return true;

}


bool SVBamReader::R2CbamToBAVec(BamAlignmentVector &bav) { 

  // setup a map of reads to make sure we aren't duplicating
  if (!m_reader.Open(m_bam)) {
    std::cerr << "FAILED TO OPEN BAM: " << m_bam << std::endl;
    std::cerr << m_bam << std::endl;
    return false;
  }

  if (m_bai.length() == 0) {
    std::cerr << "BAM index not set. Run reader.findBamIndex()" << endl;
    return false;
  }

  if (!m_reader.OpenIndex(m_bai)) {
    if (!m_reader.OpenIndex(m_bam_bai)) {
      std::cerr << "FAILED TO OPEN INDEX: " << m_bam_bai << endl;
      return false;
    }
  }

  if (!m_reader.SetRegion(m_region.chr, m_region.pos1, m_region.chr, m_region.pos2)) {
    std::cerr << "FAILED TO SET REGION\n"; 
    return false;
  }

  BamTools::BamAlignment a;

  // loop through the reads
  while (m_reader.GetNextAlignment(a)) 
    bav.push_back(a);

  return true;
    
}

// check whether a read is a tumor read (must have a JW tag)
bool SVBamReader::IsTumorRead(const BamAlignment &a) {
  
  string tmp_name;
  if (!a.GetTag("JW", tmp_name)) {
    cerr << "IsTumorRead: No JW tag" << endl;
    return false;
  }
  
  if (tmp_name.at(0) == 't')
    return true;
  return false;
}

// remove duplicate reads by name and alignment flag
void SVBamReader::deduplicateReads(const BamAlignmentVector &inbav, BamAlignmentVector &outbav) {
  
  unordered_map<string, bool> name_map;
  unordered_map<string, bool> seq_map;


  for (BamAlignmentVector::const_iterator it = inbav.begin(); it != inbav.end(); it++) {

    // deduplicate by query-bases / position
    string sname = to_string(it->RefID) + "_" + to_string(it->Position) + "_" + to_string(it->MateRefID) + "_" + to_string(it->MatePosition) + it->QueryBases;    
    // deduplicate by Name
    string uname = it->Name + "_" + to_string(it->IsFirstMate());

    // its not already add, insert
    if (name_map.count(uname) == 1 && seq_map.count(sname) == 1) {  
      name_map.insert(pair<string, int>(uname, true));
      seq_map.insert(pair<string, int>(sname, true));
      outbav.push_back(*it);
    } 

  }
}

// clean out the tags
void SVBamReader::clearFinalTags(BamAlignment * align) {
    align->RemoveTag("J2");
    //align->RemoveTag("RP");
    align->RemoveTag("HP");
}

// deduplicate the reads by position (deals with reads which are not correctly marked as duplicates)
void SVBamReader::deduplicateReadsPos(const BamAlignmentVector &inbav, BamAlignmentVector &outbav) {
    
  unordered_map<string, int> name_map;
  for (BamAlignmentVector::const_iterator it = inbav.begin(); it != inbav.end(); it++) {
    
    string uname = to_string(it->RefID) + "_" + to_string(it->Position) + "_" + to_string(it->MateRefID) + "_" + to_string(it->MatePosition) + it->QueryBases;
    
    if (name_map.find(uname) == name_map.end()) {  // its not already added
      name_map.insert(pair<string, int>(uname, 0));
      outbav.push_back(*it);
    } // else it is a dupe, don't addd
  }
}

unsigned SVBamReader::getClipCount(BamTools::BamAlignment a) {
  
  std::vector<int> clipSize;
  std::vector<int> readPos;
  std::vector<int> genPos;
  a.GetSoftClips(clipSize, readPos, genPos, false);
  
  // get the clip number
  unsigned clipnum = 0;
  for(std::vector<int>::iterator j=clipSize.begin();j!=clipSize.end();++j)
    clipnum += *j;
  return clipnum;
}

string SVBamReader::findBamIndex(string bam) {
  
  // try with bai extension
  string bai = bam;
  bai = bai.substr(0, bam.size()-3).append("bai");
  if (SnowUtils::existTest(bai))
    return bai;
  
  // try with .bam.bai extension
  bai = bam;
  bai = bam.append(".bai");
  if (SnowUtils::existTest(bai))
    return bai;
  
  return "";
  
}

// get the reference vector from a BAM header
void SVBamReader::getRefVector(string bamfile, RefVector &ref) {

  BamReader read;
  if (!read.Open(bamfile)) {
    cerr << "Failed to open contig BAM to get reference on bam " << bamfile << endl;
    cerr << "   Setting default of 1,2,..., Y" << endl;
    for (int i = 0; i < 25; i++) {
      RefData rf(CHR_NAME[i], CHR_LEN[i]);
      ref.push_back(rf);      
    }
    read.Close();
    return;
  }

  ref = read.GetReferenceData();
  read.Close();
}

// get the SamHeader from a bamfile
string SVBamReader::getSamHeader(string bamfile, SamHeader &sam) {

  BamReader read;
  if (!read.Open(bamfile)) {
    cerr << "Failed to open contig BAM to get header on bam " << bamfile << endl;
    cerr << "   Setting default of 'none'" << endl;
    sam = SamHeader ("none");
    read.Close();
    return "none";
  }


  sam = read.GetHeader();

  string out = read.GetHeaderText();
  read.Close();
  return out;

}

// count the number of discorant reads matching the region pair
size_t SVBamReader::discordantCount(GenomicRegion const &gr1, GenomicRegion const &gr2, int span) {
 
  // setup a map of reads to make sure we aren't duplicating
  if (!m_reader.Open(m_bam)) {
    std::cerr << "FAILED TO OPEN BAM: " << m_bam << std::endl;
    std::cerr << m_bam << std::endl;
    return false;
  }

  if (m_bai.length() == 0) {
    std::cerr << "BAM index not set. Run reader.findBamIndex()" << endl;
    return false;
  }

  if (!m_reader.OpenIndex(m_bai)) {
    if (!m_reader.OpenIndex(m_bam_bai)) {
      std::cerr << "FAILED TO OPEN INDEX: " << m_bam_bai << endl;
      return false;
    }
  }

  if (!m_reader.SetRegion(m_region.chr, m_region.pos1, m_region.chr, m_region.pos2)) {
    std::cerr << "FAILED TO SET REGION\n"; 
    return false;
  }
 
  BamAlignment a;
  bool intrachrom = gr1.chr == gr2.chr;
  
  size_t hitcount = 0;

  while (m_reader.GetNextAlignmentCore(a)) {
    if ((a.InsertSize > span - 1000 && a.InsertSize < span + 1000 && intrachrom) || (a.RefID != a.MateRefID && !intrachrom)) {
      GenomicRegion gr_read(a.RefID, a.Position, a.Position);
      if (gr1.getOverlap(gr_read) > 0 && gr2.getOverlap(gr_read) > 0)
	hitcount++;
    }
  }

  return hitcount;

}

bool SVBamReader::preprocessBam(BamTools::BamWriter &writer, BamQC &qc, bool qc_only) {

  // setup a map of reads to make sure we aren't duplicating
  if (!m_reader.Open(m_bam)) {
    std::cerr << "FAILED TO OPEN BAM: " << m_bam << std::endl;
    std::cerr << m_bam << std::endl;
    return false;
  }
  
  if (m_bai.length() == 0) {
    std::cerr << "BAM index not set. Run reader.findBamIndex()" << endl;
    return false;
  }
  
  if (!m_reader.OpenIndex(m_bai)) {
    if (!m_reader.OpenIndex(m_bam_bai)) {
      std::cerr << "FAILED TO OPEN INDEX: " << m_bam_bai << endl;
      return false;
    }
  }
  
  if (!m_reader.SetRegion(m_region.chr, m_region.pos1, m_region.chr, m_region.pos2)) {
    std::cerr << "FAILED TO SET REGION\n"; 
    return false;
  }
 
  int keep_counter = 0;
  int total = 0;
  int keep_counter_MAIN = 0;
  int total_MAIN = 0;

  int mapq0_keep_counter = 0;
  int discordant_keep_counter = 0;
  int n_keep_counter = 0;
  int clipped_keep_counter = 0;

  int mapq0_counter = 0;
  int discordant_counter = 0;
  int n_counter = 0;
  int clipped_counter = 0;

  int pileup = 0;
  
  BamTools::BamAlignment a;

  BamAlignmentVector bam_buffer;
  vector<int> mapq_buffer;

  /*
  GenomicIntervalTree tree;
  bool havetree = true;
  if (grm.count(m_region.chr) == 1)
    tree = grm[m_region.chr];
  else
    havetree = false;
  */

  // start the full keep check
  /*
  GenomicRegionVector::iterator gt;
  if (grv_fullkeep.size() > 0) {
    gt = grv_fullkeep.begin();
    while (gt->getOverlap(m_region) == 0)
      gt++;
  }
  */

  // loop through the reads
  //bool in_full_region = false;
  //bool grv_blank = grv_fullkeep.size() == 0;

  while (m_reader.GetNextAlignment(a)) {

    //bool old_full_region = in_full_region;
    //GenomicRegion region_old;
    //if (grv_fullkeep.size() > 0)
    //  region_old = (*gt);

    // check if we need to move the reigon pointer, and if in region
    /*
    if (gt != grv_fullkeep.end() && !grv_blank) { // it's not at the end
      if (a.Position > gt->pos2) { // we are beyond the region
	in_full_region = false;
	gt++;

	if (gt->chr != m_region.chr)
	  gt = grv_fullkeep.end();
      } else if (a.Position >= gt->pos1) {
	in_full_region = true;
      }
    }
    */
    /*
    // check that we just left a full region. if so, processes buffer 
    if (old_full_region && !in_full_region && grv_fullkeep.size() > 0) {
      // transfer mini-buffer to region-buffer
      bav_full_region.insert(bav_full_region.begin(), bam_buffer.begin(), bam_buffer.end());

      // loop and find pairs that are weird or doubled and set to keep
      set<string> name_set;
      //unordered_map<string, size_t> name_map;
      for (auto read : bav_full_region) {
	int32_t wp;
	if (read.GetTag("WP",wp))
	  name_set.insert(a.Name);
      }

      // loop through and actually write the keepers
      size_t keep_region_count = 0;
      int32_t wr;
      for (auto read : bav_full_region) {
	if (name_set.count(read.Name) == 1 || read.GetTag("WR", wr)) {
	  writer.SaveAlignment(read);
	  keep_region_count++;
	}
      }

      //
      if (m_verbose > 1) 
	cout << "Processed keep-all " << region_old << " with read count: " << region_old.tcount << " keeping " << keep_region_count << " weird / region-paired reads" << endl;


      // clear the buffers
      bav_full_region.clear();
      bam_buffer.clear();
      pileup = 0;
    }

    */

    //if (in_full_region)
    //  gt->tcount++;

    total++;
    total_MAIN++;
        
    if (m_verbose > 0 && total % 500000 == 0) {

      char buffer[100];
      int perc  = SnowUtils::percentCalc<int>(keep_counter_MAIN, total_MAIN); 
      string posstring = SnowUtils::AddCommas<int>(a.Position);
      sprintf (buffer, "Reading read %11s at position %2s:%-11s. Kept %11s (%2d%%) [running count across whole BAM]",  SnowUtils::AddCommas<int>(total_MAIN).c_str(), GenomicRegion::chrToString(a.RefID).c_str(), posstring.c_str(),  SnowUtils::AddCommas<int>(keep_counter_MAIN).c_str(), perc);
      printf ("%s\n",buffer);
      char buffer2[100];
      sprintf(buffer2, "   Filter (%% of kept)  -- Reads with N (%2d%%), Mapq0 (%2d%%), Discordant (%2d%%), Clipped (%2d%%)", 
	      SnowUtils::percentCalc<int>(n_keep_counter, keep_counter), 
	      SnowUtils::percentCalc<int>(mapq0_keep_counter, keep_counter), 
	      SnowUtils::percentCalc<int>(discordant_keep_counter, keep_counter), 
	      SnowUtils::percentCalc<int>(clipped_keep_counter, keep_counter));
      char buffer3[100];
      sprintf(buffer3, "   Filter (%% of total) -- Reads with N (%2d%%), Mapq0 (%2d%%), Discordant (%2d%%), Clipped (%2d%%)", 
	      SnowUtils::percentCalc<int>(n_counter, total), 
	      SnowUtils::percentCalc<int>(mapq0_counter, total), 
	      SnowUtils::percentCalc<int>(discordant_counter, total), 
	      SnowUtils::percentCalc<int>(clipped_counter, total));
      if (m_verbose > 1) {
	printf("%s\n", buffer2);
	printf("%s\n", buffer3);
      }

      // zero the counters
      mapq0_keep_counter = 0;
      discordant_keep_counter = 0;
      n_keep_counter = 0;
      clipped_keep_counter = 0;
      
      mapq0_counter = 0;
      discordant_counter = 0;
      n_counter = 0;
      clipped_counter = 0;
      keep_counter = 0;
      total = 0;

      // kill if seen 50m reads, and it's looking bad
      if (perc >= perclimit && total > 25000000) { 
	cerr << "This is a a really bad BAM after checking out 25m+ reads. Killing job. Percent weird reads: " << perc << " is above limit of " << perclimit << endl;
	cerr << "Reading in region" << m_region << endl;
	exit(EXIT_FAILURE);
      }
    }
    string rule_pass = m_mr->isValid(a);
    //    bool in_full_region = m_mr->isValid(a);
    
    /*
    // get the number of soft-clipped bases BEFORE TRIMMING
    unsigned clipnum = getClipCount(a);
  
    // get some other info about the alignment before filling string fields
    //bool single_align = a.MateRefID == -1;
    m_skip_supp = true;
    bool unmap = !a.IsMapped() || !a.IsMateMapped();
    bool discr = (m_isize > 0) && (std::abs(a.InsertSize) > m_isize || a.RefID != a.MateRefID); 
    bool mapqr = (a.MapQuality >= m_mapq || unmap) && !a.IsFailedQC() && !a.IsDuplicate();
    bool supp  = a.IsPrimaryAlignment() || !m_skip_supp; // if m_skip_supp is true, only keep primary alignments
    bool clip  = clipnum >= m_minclip || m_minclip == 0;
    mapqr = mapqr && supp;

    // FR read
    bool FR_f = !a.IsReverseStrand() && (a.Position < a.MatePosition) && (a.RefID == a.MateRefID) &&  a.IsMateReverseStrand();
    bool FR_r =  a.IsReverseStrand() && (a.Position > a.MatePosition) && (a.RefID == a.MateRefID) && !a.IsMateReverseStrand();
    bool FR = FR_f || FR_r;
    discr = discr || !FR;

    //    if ( (discr || unmap || clip) && (mapqr || discr)  && a.Length >= 60) {
    
    //a.BuildCharData(); //populate the rest of the string data
    // clip bases with low optical quality

    std::string trimmed_bases = a.QueryBases;
    std::string trimmed_quals = a.Qualities;
    bool tooshort = false;
    bool not_real_clip = false;
    int startpoint;
    int endpoint;
    
    softClip(m_qualthresh, trimmed_bases, trimmed_quals, clipnum, tooshort, not_real_clip, startpoint, endpoint);
    //cout << startpoint << " " << endpoint << " too short " << tooshort << endl;
    // check that there are no N bases
    bool npass = passNtest(a.QueryBases);
    
    // get the number of matches
    uint32_t nm;
    if (a.GetTag("NM", nm)) {} else { nm = 0; }
    bool nmpass = nm <= m_nmlim;

    // get the AS tag
    uint32_t as;
    if (a.GetTag("AS", as)) {} else { as = 0; }
    // get the XP tag
    uint32_t xp;
    if (a.GetTag("XP", as)) {} else { xp = 0; }

    // make sure raw clip is sufficient, plus processed
    bool clipp = (clip && !not_real_clip) || m_minclip == 0;

    // make sure the read is long enough
    bool lenpass = (endpoint - startpoint + 1) >= min_length;

    bool hardclip = false;
    // remove hard clips
    for (auto cig : a.CigarData)
      if (cig.Type == 'H')
	hardclip = true;
    
    // get the mean phred quality
    size_t i = 0;
    int phred = 0;
    while(i < a.Qualities.length()) {
        phred += char2phred(a.Qualities[i]);
	i++;
    }
    if (a.Qualities.length() > 0)
      phred = static_cast<int>(floor(static_cast<float>(phred) / a.Qualities.length()));

    */

    // build the qc
    try {
      string rgroup;
      if (!a.GetTag("RG",rgroup))
	cerr << "Failed to read rgroup" << endl;

      int this_isize = a.InsertSize;
      this_isize = (a.MateRefID != a.RefID || this_isize > 2000) ? 2000 : this_isize;

      assert(a.MapQuality <= 60);
      //assert(clipnum <= 101);
      //assert(as <= 101);
      //assert(xp <= 101);
      assert(a.Length <= 101 && a.Length  >= 0);
      //assert(phred <= 60 && phred  >= 0);
      //assert(nm <= 101);

      //qc.map[rgroup].nm[nm]++;
      qc.map[rgroup].mapq[a.MapQuality]++;
      //if (a.InsertSize > 0 && a.IsPaired() && (FR_f || FR_r) ) // only count "proper" reads
	//qc.map[rgroup].isize[this_isize]++;
      //qc.map[rgroup].xp[xp]++;
      //qc.map[rgroup].len[a.Length]++;
      //qc.map[rgroup].as[as]++;
      //qc.map[rgroup].clip[clipnum]++;
      //qc.map[rgroup].phred[phred]++;
      qc.map[rgroup].num_reads++;
      if (!a.IsMapped())
	qc.map[rgroup].unmap++;
      if (a.IsFailedQC()) 
	qc.map[rgroup].qcfail++;
      if (a.IsDuplicate())
	qc.map[rgroup].duplicate++;
      if (!a.IsPrimaryAlignment())
	qc.map[rgroup].supp++;
    } catch (...) {
      cerr << "Failed at adding to QC" << endl;
      //cerr << "Readgroup " << "NM " << nm << " mapq " << a.MapQuality << " xp " << xp << " len " << a.Length <<
      //	" as " << as << " phred " << phred << endl;
    }

    // counter
    //mapq0_counter += (a.MapQuality == 0 ) ? 1 : 0; 
    //n_counter += npass ? 0 : 1;
    //clipped_counter += clip ? 1 : 0;
    //discordant_counter += discr ? 1 : 0;

    //bool qual_read = (npass || !exclude_n) && nmpass && lenpass && !hardclip;
    //bool save_read = (discr && !hardclip) || (unmap && qual_read) || (clipp && qual_read);

    // check if read or pair is in full region
    /*
    bool in_full_region = false;
    if (a.RefID == m_region.chr && a.MateRefID == m_region.chr && havetree) {
      GenomicIntervalVector grv;
      tree.findOverlapping(a.Position, a.Position + a.Length, grv);
      tree.findOverlapping(a.MatePosition, a.MatePosition + a.Length, grv);
      
      if (grv.size() > 0) {
	in_full_region = true;
	a.AddTag("WR", "i", 0);
      }
      
    }
    */

    //bool save_read = true;
    if ( rule_pass != "" && !qc_only ) {

      mapq0_keep_counter += (a.MapQuality == 0 ) ? 1 : 0; 
      //n_keep_counter += npass ? 0 : 1;
      //clipped_keep_counter += clipp ? 1 : 0; // count AFTER filtering
      //discordant_keep_counter += discr ? 1 : 0;

      // keep track of pile
      if (a.MapQuality == 0) 
	pileup++;

      // add a tag to say which rule it pass
      a.AddTag("RL","i",rule_pass);

      //if (!a.AddTag("SS", "i", startpoint)) {
      //	cerr << "Can't add SS Tag" << endl;
      //	exit(EXIT_FAILURE);
      //}
      //a.AddTag("SE", "i", endpoint);
      //a.RemoveTag("MD");
      //a.RemoveTag("RD");
      //a.RemoveTag("SA");
      //a.RemoveTag("AS");
      //a.RemoveTag("XS");
      //a.RemoveTag("RG");
      //a.Qualities = "";
      //writer.SaveAlignment(a);

      // add a weird tag
      //if (save_read)
      //	a.AddTag("WR", "i", 0);
      
      bam_buffer.push_back(a);
      keep_counter++;
      keep_counter_MAIN++;
      
      int buffer_lim = 100;
      // deal with bam buff
      if (bam_buffer.size() >= buffer_lim/* && !in_full_region*/) {
	// check if bad region
	int32_t wr;
	if (pileup >= buffer_lim * 0.8 && (bam_buffer.back().Position - bam_buffer[0].Position <= 40)) {
	  for (auto it = bam_buffer.begin(); it != bam_buffer.end(); it++) 
	    if (it->MapQuality > 0)
	      writer.SaveAlignment(*it);
	  if (m_verbose > 2)
	    cout << "Detected mapq 0 pileup of " << pileup << " at " << a.RefID+1 << ":" << bam_buffer[0].Position << "-" << bam_buffer.back().Position << endl;
	} 
	// it's OK or its in full region
	else if (bam_buffer.size() >= buffer_lim) {
	  for (auto it = bam_buffer.begin(); it != bam_buffer.end(); it++) 
	    writer.SaveAlignment(*it);
	}

	// further filter on clusters
	//GenomicRegionVector clusters;
	//typedef unordered_map<string, string> RMap;
	//RMap clustered_reads;
	//ClusterReads::clusterReads(bam_buffer, clusters, clustered_reads,
	//			   "BB", m_isize, 1000);
	
	//debug
	//for (auto cc : clusters)
	//  cout << "cluster "<< cc.cluster << " count: " << cc.ncount << endl;
	
	bam_buffer.clear();
	pileup = 0;

      } // end buffer check
      
    } // end save read checking

  } // end read while loop

  // write the final buffer
  for (auto it = bam_buffer.begin(); it != bam_buffer.end(); it++)
    writer.SaveAlignment(*it);
  
  // print the final message
  if (m_verbose > 0) {
    char buffer[100];
    int perc  = SnowUtils::percentCalc<int>(keep_counter_MAIN, total_MAIN); 
    string posstring = SnowUtils::AddCommas<int>(a.Position);
    //sprintf (buffer, "Finished region at %20s. Kept %11s (%2d%%) [running count across whole BAM]",  m_region.toStringOffset().c_str(), SnowUtils::AddCommas<int>(keep_counter_MAIN).c_str(), perc);
    printf ("%s\n",buffer);
    char buffer2[100];
    sprintf(buffer2, "   Filter (%% of kept)  -- Reads with N (%2d%%), Mapq0 (%2d%%), Discordant (%2d%%), Clipped (%2d%%)", 
	    SnowUtils::percentCalc<int>(n_keep_counter, keep_counter), 
	    SnowUtils::percentCalc<int>(mapq0_keep_counter, keep_counter), 
	    SnowUtils::percentCalc<int>(discordant_keep_counter, keep_counter), 
	    SnowUtils::percentCalc<int>(clipped_keep_counter, keep_counter));
    char buffer3[100];
    sprintf(buffer3, "   Filter (%% of total) -- Reads with N (%2d%%), Mapq0 (%2d%%), Discordant (%2d%%), Clipped (%2d%%)", 
	    SnowUtils::percentCalc<int>(n_counter, total), 
	    SnowUtils::percentCalc<int>(mapq0_counter, total), 
	    SnowUtils::percentCalc<int>(discordant_counter, total), 
	    SnowUtils::percentCalc<int>(clipped_counter, total));
    if (m_verbose > 1) {
      printf("%s\n", buffer2);
      printf("%s\n", buffer3);
    }
  }
  
  
  return true;

}
