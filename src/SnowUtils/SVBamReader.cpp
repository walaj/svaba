#include "SVBamReader.h"
#include "unistd.h"
#include "api/algorithms/Sort.h"
#include <time.h> // for now
#include "SnowUtils.h"

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
    bool mapqr = (a.MapQuality >= m_mapq || unmap) && !a.IsFailedQC() && !a.IsDuplicate();
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
      softClip(m_qualthresh, trimmed_bases, trimmed_quals, clipnum, tooshort, not_real_clip, startpoint, endpoint);
      
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

/*
bool SVBamReader::bamToBAVecCluster(BamAlignmentVector &bav) { 

  // setup a map of reads to make sure we aren't duplicating
  //std::unordered_map<std::string, int> name_map;
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

  if (!m_reader.SetRegion(m_region[0], m_region[1], m_region[2], m_region[3])) {
    std::cerr << "FAILED TO SET REGION\n"; 
    return false;
  }

  unsigned counter = 0;

  BamTools::BamAlignment a;

  // loop through the reads
  while (m_reader.GetNextAlignmentCore(a)) {
 
    total_reads++;

    bool discr = (m_isize > 0) && (std::abs(a.InsertSize) > m_isize || a.RefID != a.MateRefID); 
    bool mapqr = (a.MapQuality >= m_mapq) && !a.IsFailedQC() && !a.IsDuplicate();
    bool supp  = a.IsPrimaryAlignment() || !m_skip_supp; // if m_skip_supp is true, only keep primary alignments
    mapqr = mapqr && supp;

    if ( discr && mapqr )   {

      a.BuildCharData(); //populate the rest of the string data

      // get the number of matches
      int nm;
      if (a.GetTag("NM", nm)) {} else { nm = 0; }

      // add a tag for position of pairmate
      uint32_t val = GenomicRegion::convertPos(a.MateRefID, a.MatePosition);
      a.AddTag("RP", "I", val);
      
      // want to index by position of anchor read
      std::string hp = string_numf(a.RefID, 2) + "_" + string_numf(a.Position, 9);
      a.AddTag("HP", "Z", hp);
      
      //set the prefix for this read name 
      std::ostringstream s_first;

      s_first << m_prefix << a.AlignmentFlag << "_" << a.Name;
      a.AddTag("JW", "Z", s_first.str());
	
      used_num++;
      a.RemoveTag("Q2"); 
      a.RemoveTag("R2"); 
      a.RemoveTag("XM"); 
      a.RemoveTag("SA"); 
      a.RemoveTag("XS"); 
      a.RemoveTag("RG"); 
      a.RemoveTag("AS"); 
      a.RemoveTag("XS");
      //a.Qualities=""; // extreme
      //a.QueryBases = ""; // extreme	
      bav.push_back(a);
      counter++;

      if (counter > limit) { // hit the read limit
	bav.clear();
	cout << "Limit of " << limit << " hit on region: " << 
	  m_region[0]+1 << ":" << m_region[1] << "-" << m_region[3] << 
	  " in " << m_bam << endl;
	return false;
      }
      
    }
  } // end while
	
  // ensure that it is sorted by HP
  std::sort( bav.begin(), bav.end(), BamTools::Algorithms::Sort::ByTag<std::string>("HP", BamTools::Algorithms::Sort::AscendingOrder) );
  
  return true;
    
}
*/

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
// pay attention to IR flag (for R2 reads) so that you preferentially add
// real reads before adding R2 created reads
void SVBamReader::deduplicateReads(const BamAlignmentVector &inbav, BamAlignmentVector &outbav) {
  
  unordered_map<string, bool> name_map;
  BamAlignmentVector r2vec;

  for (BamAlignmentVector::const_iterator it = inbav.begin(); it != inbav.end(); it++) {
    
    string uname = it->Name + "_" + to_string(it->IsFirstMate());
    bool not_found = name_map.find(uname) == name_map.end();
    bool isR2 = it->HasTag("IR");

    if (not_found && !isR2) {  // its not already added and is not R2. Definitely add
      name_map.insert(pair<string, int>(uname, true));
      outbav.push_back(*it);
    } else if (not_found && isR2) { // is new read, but is R2. Wait until done to see if should add
      r2vec.push_back(*it);
    } // else it is a total dupe, don't addd
  }
  
  // loop through the R2s and see what is not added
  for (BamAlignmentVector::const_iterator it = r2vec.begin(); it != r2vec.end(); it++) {
    string uname = it->Name + "_" + to_string(it->IsFirstMate());
    bool not_found = (name_map.find(uname) == name_map.end());
    if (not_found)
      outbav.push_back(*it);
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
    
    string uname = to_string(it->RefID) + "_" + to_string(it->Position) + "_" + to_string(it->MateRefID) + "_" + to_string(it->MatePosition);
    
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
void SVBamReader::getSamHeader(string bamfile, SamHeader &sam) {

  BamReader read;
  if (!read.Open(bamfile)) {
    cerr << "Failed to open contig BAM to get header on bam " << bamfile << endl;
    cerr << "   Setting default of 'none'" << endl;
    sam = SamHeader("none");
    read.Close();
  }

  sam = read.GetHeader();
  read.Close();
  return;

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

bool SVBamReader::preprocessBam(BamTools::BamWriter &writer) {

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
 
  int counter = 0;
  int total = 0;

  BamTools::BamAlignment a;

  // loop through the reads
  while (m_reader.GetNextAlignmentCore(a)) {
    total++;

    // get the number of soft-clipped bases BEFORE TRIMMING
    unsigned clipnum = getClipCount(a);
    
    if (m_verbose > 0 && total % 1000000 == 0) {
      char buffer[100];
      int perc  = static_cast<int>(floor((float)counter / (float)total * 100.0));
      string posstring = SnowUtils::AddCommas<int>(a.Position);
      sprintf (buffer, "Reading read %11s at position %2s:%-11s. Kept %11s (%2d%%)",  SnowUtils::AddCommas<int>(total).c_str(), GenomicRegion::chrToString(a.RefID).c_str(), posstring.c_str(),  SnowUtils::AddCommas<int>(counter).c_str(), perc);
      printf ("%s\n",buffer);

      // kill if seen 50m reads, and it's looking bad
      if (perc >= perclimit && total > 25000000) { 
	cerr << "This is a a really bad BAM after checking out 50m+ reads. Killing job. Percent weird reads: " << perc << " is above limit of " << perclimit << endl;
	cerr << "Reading in region" << m_region << endl;
	exit(EXIT_FAILURE);
      }
    }
      
    // get some other info about the alignment before filling string fields
    //bool single_align = a.MateRefID == -1;
    bool unmap = !a.IsMapped() || !a.IsMateMapped();
    bool discr = (m_isize > 0) && (std::abs(a.InsertSize) > m_isize || a.RefID != a.MateRefID); 
    bool mapqr = (a.MapQuality >= m_mapq || unmap) && !a.IsFailedQC() && !a.IsDuplicate();
    bool supp  = a.IsPrimaryAlignment() || !m_skip_supp; // if m_skip_supp is true, only keep primary alignments
    bool clip  = clipnum >= m_minclip || m_minclip == 0;
    mapqr = mapqr && supp;

    if ( (discr || unmap || clip) && (mapqr || discr)  && a.Length >= 60) {

      a.BuildCharData(); //populate the rest of the string data
      // clip bases with low optical quality
      std::string trimmed_bases = a.QueryBases;
      std::string trimmed_quals = a.Qualities;
      bool tooshort = false;
      bool not_real_clip = false;
      int startpoint;
      int endpoint;

      softClip(m_qualthresh, trimmed_bases, trimmed_quals, clipnum, tooshort, not_real_clip, startpoint, endpoint);
      
      // check that there are no N bases
      bool npass = passNtest(a.QueryBases);

      // get the number of matches
      uint32_t nm;
      if (a.GetTag("NM", nm)) {} else { nm = 0; }
      bool nmpass = nm <= m_nmlim;
      
      // make sure raw clip is sufficient, plus processed
      bool clipp = (!not_real_clip && !tooshort && clip) || m_minclip == 0;

      //      if ( (discr || unmap || clipp) && mapqr && !tooshort && npass && nmpass ) {
      if ( discr || unmap || (clipp && npass && nmpass) ) {
	//size_t readlen = a.Length;
	//if (trimmed_bases.length() != readlen) {
	if (!a.AddTag("SS", "i", startpoint)) {
	  cerr << "Can't add SS Tag" << endl;
	  exit(EXIT_FAILURE);
	}
	a.AddTag("SE", "i", endpoint);
	  //}
	a.RemoveTag("MD");
	a.RemoveTag("RD");
	a.RemoveTag("SA");
	a.RemoveTag("AS");
	a.RemoveTag("XS");
	//a.RemoveTag("RG");
        a.Qualities = "";
	writer.SaveAlignment(a);
	counter++;
      }
    }

  }

  return true;

}
