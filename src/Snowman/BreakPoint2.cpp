#include "BreakPoint2.h"

#include <getopt.h>
#include <iomanip>
#include <cassert>

#include "SnowTools/gzstream.h"

#include "SnowmanUtils.h"

//#define T_SPLIT_BUFF 15
//#define N_SPLIT_BUFF 8
#define T_SPLIT_BUFF 8
#define N_SPLIT_BUFF 8
// if the insertion is this big or larger, don't require splits to span both sides
#define INSERT_SIZE_TOO_BIG_SPAN_READS 16

// if homology is greater than homology / HOMOLOGY_FACTOR, then reject for assembly-only
#define HOMOLOGY_FACTOR 4
#define MIN_SOMATIC_RATIO 15

// define repeats
static std::vector<std::string> repr = {"AAAAAAAA", "TTTTTTTT", "CCCCCCCC", "GGGGGGGG",
				 "TATATATATATATA", "ATATATATATATAT", 
				 "GCGCGCGCGCGCGC", "CGCGCGCGCGCGCG", 
				 "TGTGTGTGTGTGTG", "GTGTGTGTGTGTGT", 
				 "TCTCTCTCTCTCTC", "CTCTCTCTCTCTCT", 
				 "CACACACACACACA", "ACACACACACACAC", 
				 "GAGAGAGAGAGAGA", "AGAGAGAGAGAGAG"};

// define large repeats
static std::vector<std::string> hirepr = {"AAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTT", "CCCCCCCCCCCCCCCC", "GGGGGGGGGGGGGGG"};

// check a sequence for a large homo/dinuc repeat
bool __check_homopolymer(const std::string& s) {

  for (auto& i : repr)
    if (s.find(i) != std::string::npos)
      return true;

  return false;

}

double scale_factor = 5.0;
static std::unordered_map<int, double> ERROR_RATES = {{0,scale_factor * 1e-4}, {1, scale_factor * 1e-4}, {2, scale_factor * 1e-4}, {3, scale_factor * 1e-4}, {4, scale_factor * 1e-4}, {5, scale_factor * 2e-4}, {6, scale_factor * 5e-4}, {7, scale_factor * 1e-3},
						      {scale_factor * 8,2e-3}, {9,scale_factor * 3e-3}, {10, scale_factor * 1e-2}, {11, scale_factor * 2e-2}, {12, scale_factor * 3e-5}};

namespace SnowTools {

  double __myround(double x) { return std:: floor(x * 10) / 10; }

  // make the file string
  std::string BreakPoint::toFileString(bool noreads) {
    
    // make sure we already ran scoring
    assert(evidence.length());
    assert(confidence.length());
    
    std::string sep = "\t";
    std::stringstream ss;
    
    // put the read names into a string
    if (!noreads)  {
      __format_readname_string();
    }
    double max_lod = 0;
    for (auto& s : allele) {
      max_lod = std::max(max_lod, s.second.LO);
    }

    ss << b1.chr_name << sep << b1.gr.pos1 << sep << b1.gr.strand << sep 
       << b2.chr_name << sep << b2.gr.pos1 << sep << b2.gr.strand << sep 
       << ref << sep << alt << sep 
       << getSpan() << sep
       << b1.mapq << sep << b2.mapq << sep 
       << b1.nm << sep << b2.nm << sep 
       << dc.mapq1 << sep << dc.mapq2 << sep
       //<< dc.ncount << sep << dc.tcount << sep
       << b1.sub_n << sep << b2.sub_n << sep      
       << (homology.length() ? homology : "x") << sep 
       << (insertion.length() ? insertion : "x") << sep 
       << cname << sep
       << num_align << sep 
       << confidence << sep << evidence << sep
       << quality << sep
       << secondary << sep << somatic_score << sep << somatic_lod << sep 
       << max_lod << sep 
       << pon << sep << (repeat_seq.length() ? repeat_seq : "x") << sep 
       << blacklist << sep << (rs.length() ? rs : "x") << sep 
       << (read_names.length() ? read_names : "x");

    for (auto& a : allele)
      ss << sep << a.second.toFileString();

    return ss.str();

  }
  
  BreakEnd::BreakEnd(const GenomicRegion& g, int mq, const std::string& chr_n) {
    gr = g;
    mapq = mq;
    chr_name = chr_n;
    cpos = -1; nm = -1; matchlen = -1; sub_n = 0; local = false;
  }
  
  ReducedBreakEnd::ReducedBreakEnd(const GenomicRegion& g, int mq, const std::string& chr_n) {
    
    gr = g;
    mapq = mq;
    chr_name = chr_n;
    nm = 0;
  }


  // make a breakpoint from a discordant cluster 
  BreakPoint::BreakPoint(DiscordantCluster& tdc, const BWAWrapper * bwa, SnowTools::DiscordantClusterMap& dmap) {
    
    num_align = 0;
    dc = tdc;
    
    std::string chr_name1, chr_name2;

    try {
       // if this throw error, it means that there are more chr in 
       // reads than in reference
       chr_name1 = bwa->ChrIDToName(dc.m_reg1.chr); //bwa->ChrIDToName(tdc.reads.begin()->second.ChrID());
       chr_name2 = bwa->ChrIDToName(dc.m_reg2.chr); //bwa->ChrIDToName(tdc.reads.begin()->second.ChrID());
    } catch (...) {
      std::cerr << "Warning: Found mismatch between reference genome and BAM genome." << std::endl;
      chr_name1 = "Unknown";
      chr_name2 = "Unknown";
    }

    assert(chr_name1.length());
    assert(chr_name2.length());

    int pos1 = (dc.m_reg1.strand == '+') ? dc.m_reg1.pos2 : dc.m_reg1.pos1;
    int pos2 = (dc.m_reg2.strand == '+') ? dc.m_reg2.pos2 : dc.m_reg2.pos1;
    b1 = BreakEnd(GenomicRegion(dc.m_reg1.chr, pos1, pos1), dc.mapq1, chr_name1);
    b2 = BreakEnd(GenomicRegion(dc.m_reg2.chr, pos2, pos2), dc.mapq2, chr_name2);
    b1.gr.strand = dc.m_reg1.strand;
    b2.gr.strand = dc.m_reg2.strand;

    // set the alt counts, counting only unique qnames
    std::unordered_map<std::string, std::unordered_set<std::string>> alt_counts;

    // add the supporting read info to allels
    for (auto& rr : dc.reads) {
      std::string sr = rr.second.GetZTag("SR");
      allele[sr.substr(0,4)].supporting_reads.insert(sr);
      alt_counts[sr.substr(0,4)].insert(rr.second.Qname());
    }
    for (auto& rr : dc.mates) {
      std::string sr = rr.second.GetZTag("SR");
      allele[sr.substr(0,4)].supporting_reads.insert(rr.second.Qname());
      alt_counts[sr.substr(0,4)].insert(rr.second.Qname());
    }
    
    // add the alt counts
    for (auto& i : allele) {
      i.second.indel = false;
      i.second.disc = alt_counts[i.first].size(); //allele[i.first].supporting_reads.size();
      i.second.__adjust_alt_counts();
      //i.second.alt = i.second.disc;
    }
      
    // give a unique id
    cname = dc.toRegionString();

    // check if another cluster overlaps, but different strands
    if (getSpan() > 800 || getSpan() == -1) // only check for large events.
    for (auto& d : dmap) {

      // don't overlap if on different chr, or same event
      if (dc.m_reg1.chr != d.second.m_reg1.chr || dc.m_reg2.chr != d.second.m_reg2.chr || d.second.ID() == dc.ID())
	continue;

      // isolate and pad
      SnowTools::GenomicRegion gr1 = d.second.m_reg1;
      SnowTools::GenomicRegion gr2 = d.second.m_reg2;
      gr1.pad(100);
      gr2.pad(100);
      
      if (dc.m_reg1.getOverlap(gr1) && dc.m_reg2.getOverlap(gr2))
	if (dc.m_reg1.strand != d.second.m_reg1.strand || dc.m_reg2.strand != d.second.m_reg2.strand) {
	  dc.m_id_competing = d.second.ID();
	  tdc.m_id_competing = d.second.ID();
	  d.second.m_id_competing = dc.ID();
	}
    }

  }
  
  bool BreakPoint::hasDiscordant() const {
    return (dc.ncount || dc.tcount);
  }
  
  bool BreakPoint::hasMinimal() const {
    int count;
    count = dc.tcount + dc.ncount + t.split + n.split;
    if (count >= 2)
      return true;
    return false;
  }
  
  bool BreakPoint::operator==(const BreakPoint &bp) const {
    return (b1.gr == bp.b1.gr && b2.gr == bp.b2.gr && bp.insertion == insertion); 
  }
    
  void BreakPoint::repeatFilter() {

    if (seq.length() == 0)
      return;

    // look at all fwd / rev repeats up to 5-mers
    repeat_seq = std::string();
    for (int fwd = 0; fwd <= 1; ++fwd)
      for (int i = 1; i <= 5; ++i) {
	std::string rr;
	__rep(i, rr, fwd > 0);
	repeat_seq = rr.length() > repeat_seq.length() ? rr : repeat_seq;
    }

    // look for massive repeats in contig
    /*for (auto& r : repr)
      if (seq.find(r) != std::string::npos && r.length() > repeat_seq.length())
	repeat_seq = r;
    */
  }


  BreakPoint::BreakPoint(const std::string &line, bam_hdr_t* h) {

    if (!h) {
      std::cerr << "BreakPoint::BreakPoint - Must supply non-empty header" << std::endl;
      exit(EXIT_FAILURE);
    }

    std::istringstream iss(line);
    std::string val;
    size_t count = 0;

    std::string chr1, pos1, chr2, pos2, chr_name1, chr_name2, id; 
    id = "";//dummmy
    char strand1 = '*', strand2 = '*';
    while (std::getline(iss, val, '\t')) {
      SampleInfo aaa;
      try{
	switch(++count) {
	case 1: chr1 = val; chr_name1 = val; break;
	case 2: pos1 = val; break; 
	case 3: assert(val.length()); strand1 = val.at(0); break;
	case 4: chr2 = val; chr_name2 = val; break;
	case 5: pos2 = val; break; 
	case 6: assert(val.length()); strand2 = val.at(0); break;
	case 7: ref = val; break;
	case 8: alt = val; break;
	case 9: break; //span = stoi(val); break; // automatically calculated
	case 10: 
	  b1 = BreakEnd(GenomicRegion(chr1, pos1, pos1, h), std::stoi(val), chr_name1); b1.gr.strand = strand1; break;
	case 11:
	  b2 = BreakEnd(GenomicRegion(chr2, pos2, pos2, h), std::stoi(val), chr_name2); b2.gr.strand = strand2; break;
	case 12: b1.nm = std::stoi(val); break;
	case 13: b2.nm = std::stoi(val); break;
	  //case 14: dc.ncount = std::stoi(val); break;
	  //case 15: dc.tcount = std::stoi(val); break;
	case 14: dc.mapq1 = std::stoi(val); break;  
	case 15: dc.mapq2 = std::stoi(val); break;  
	case 16: b1.sub_n = std::stoi(val); break;
	case 17: b2.sub_n = std::stoi(val); break;
	case 18: homology = (val == "x" ? "" : val); break;
	case 19: insertion = (val == "x" ? "" : val); break;
	case 20: cname = val; break;
	case 21: num_align = std::stoi(val); break;
	case 22: confidence = val; break;
	case 23: evidence = val; break;
	case 24: quality = std::stoi(val); break;
	case 25: secondary = val == "1";
	case 26: somatic_score = std::stod(val); break;
	case 27: somatic_lod = std::stod(val); break;
	case 28: a.LO = std::stod(val); break;
	case 29: pon = std::stoi(val); break;
	case 30: repeat_seq = val; break;
	case 31: blacklist = (val=="1"); break;
	case 32: rs = val; break;
	case 33: read_names = val; break;
        default: 
	  aaa.indel = evidence == "INDEL";
	  aaa.fromString(val);
	  id += "A";
	  allele[id] = aaa; //id is dummy. just keep in order as came in;

	  // fill in the discordant info
	  if (id == "A") // tumor
	    dc.tcount = aaa.disc;
	  else
	    dc.ncount += aaa.disc;
	    

	  break;
	}
      } catch(...) {
	std::cerr << "caught stoi error on: " << val << std::endl;
	std::cerr << line << std::endl;
	exit(1);
      }
    }

  }
  
  void BreakPoint::splitCoverage(BamReadVector &bav) {
    
    // keep track of if first and second mate covers same split. 
    // if so, this is fishy and remove them both
    std::unordered_map<std::string, bool> qname_and_num;

    // keep track of which ones are valid w/respect to unique qname / num combo (num = first, second mate)
    std::set<std::string> valid_reads; // rname
    
    // keep track of reads to reject
    std::set<std::string> reject_qnames;

    // get the homology length. this is useful because if read alignment ends in homologous region, it is not split
    int homlen = b1.cpos - b2.cpos;
    if (homlen < 0)
      homlen = 0;
   
    // loop all of the aligned reads
    for (auto& j : bav) {

      std::vector<std::string> cnvec = j.GetSmartStringTag("CN");
      
      // for indels, reads with ins alignments to del contig dont count, vice versa
      bool read_should_be_skipped = false;
      if (num_align == 1) {
	std::vector<std::string> cigvec = j.GetSmartStringTag("SC"); // read against contig CIGAR
	//if (cigvec.size() != cnvec.size())//debug
	//  std::cerr << cigvec.size() << " _ " << cnvec.size() << std::endl;
	assert(cigvec.size() == cnvec.size());
	size_t kk = 0;
	
	for (; kk < cnvec.size(); kk++)  // find the right contig r2c cigar (since can have one read to many contig)
	  if (cnvec[kk] == cname)  { // found the right cigar string

	    Cigar tcig = cigarFromString(cigvec[kk]);
	    std::vector<int> del_breaks;
	    std::vector<int> ins_breaks;
	    int pos = 0;
	    
	    // if this is a nasty repeat, don't trust non-perfect alignmentx on r2c alignment
	    if ( (repeat_seq.length() > 6 || __check_homopolymer(j.Sequence())) && tcig.size() > 1) 
	      read_should_be_skipped = true;
	    
	    // loop through r2c cigar and see positions
	    for(auto& i : tcig) {

	      if (i.Type() == 'D') 
		del_breaks.push_back(pos);
	      else if (i.Type() == 'I')
		ins_breaks.push_back(pos);

	      // update position on contig
	      if (i.ConsumesReference())
		pos += i.Length(); // update position on contig

	    }
	    
	    if (insertion.length()) { // for insertions, skip reads that have del to insertion at same pos
	      for (auto& i : del_breaks)
		if (i > b1.cpos - 2 || i < b1.cpos + 2) // if start of insertion is at start of a del of r2c
		  read_should_be_skipped = true;
	    } else { // for del, skip reads that have ins r2c to deletion at sam pos
	      for (auto& i : ins_breaks)
		if (i > b1.cpos - 2 || i < b1.cpos + 2) // if start of insertion is at start of a del of r2c
		  read_should_be_skipped = true;
	    }
	       
	  }
      }
      
      if (read_should_be_skipped) 
	continue;
      
      // get read ID
      std::string sr = j.GetZTag("SR");
      assert(sr.at(0) == 't' || sr.at(0) == 'n');
      bool tumor_read = sr.at(0) == 't';
      std::string sample_id = sr.substr(0,4);

      int rightbreak1 = b1.cpos + (tumor_read ? T_SPLIT_BUFF : N_SPLIT_BUFF); // read must extend this far right of break1
      int leftbreak1  = b1.cpos - (tumor_read ? T_SPLIT_BUFF : N_SPLIT_BUFF); // read must extend this far left of break1
      int rightbreak2 = b2.cpos + (tumor_read ? T_SPLIT_BUFF : N_SPLIT_BUFF);
      int leftbreak2  = b2.cpos - (tumor_read ? T_SPLIT_BUFF : N_SPLIT_BUFF);

      std::string contig_qname; // for sanity checking
      // get the alignment position on contig
      int pos = 0, te = 0;
      try {
	std::vector<int> posvec = j.GetSmartIntTag("SL");
	std::vector<int> tevec  = j.GetSmartIntTag("SE");
	
	for (size_t y = 0; y < cnvec.size(); ++y) {
	  if (cnvec[y] == cname) {
	    pos = posvec[y];
	    te  = tevec[y];
	    contig_qname = cname; // for sanity checking
	  }
	}
      } catch (...) {
	std::cerr << "error grabbing SE tag for tag " << j.GetZTag("SE") << std::endl;
      }
      assert(contig_qname == cname); // for sanity checking

      int rightend = te; 
      int leftend  = pos;
      bool issplit1 = (leftend <= leftbreak1) && (rightend >= rightbreak1);
      bool issplit2 = (leftend <= leftbreak2) && (rightend >= rightbreak2);
      //bool issplit1I = (leftend <= leftbreak1I) && (rightend >= rightbreak1I);
      //bool issplit2I = (leftend <= leftbreak2I) && (rightend >= rightbreak2I);

      bool both_split = issplit1 && issplit2;
      bool one_split = issplit1 || issplit2;
      //bool both_splitI = issplit1I && issplit2I;
      //bool one_splitI = issplit1I || issplit2I;
      // be more permissive for NORMAL, so keep out FPs
      bool valid  = (both_split && (tumor_read || homlen > 0)) || (one_split && !tumor_read && homlen == 0) || (one_split && insertion.length() >= INSERT_SIZE_TOO_BIG_SPAN_READS);
      //bool validI = (both_splitI && homlen > 0) || (one_splitI && homlen == 0) || (one_splitI && insertion.length() >= INSERT_SIZE_TOO_BIG_SPAN_READS);
      // requiring both break ends to be split for homlen > 0 is for situation beow
      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>A..........................
      // ............................B>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      // if we set break2 at pos B, then reads that end at A will count as splits on B,
      // when they should count as non-splits on A. For any read, if it is non-split on one 
      // segment, it should be non-split on all, for overlapping alignments like above. Not true of 
      // insertions at junctions, where one can split at one and not the other because of the intervening sequence buffer


      //debug
      /*if (sr == "t000_163_H01PEALXX140819:2:2202:14804:18907")
	std::cerr << " te " << te << " pos " << pos << " CIG " << j.GetZTag("SC") << " SL " << j.GetZTag("SL") << " SE " << j.GetZTag("SE") << 
	  " leftbreak " << leftbreak1 << " rightbreak " << rightbreak1 << " leftbreak2 " << leftbreak2 << " rightbreak2 " << rightbreak2 << 
	  " issplit1 " << issplit1 << " issplit2 " << issplit2 << " valid " << valid << std::endl;
      */

      // add the split reads for each end of the break
      // a read is split if it is spans both break ends for tumor, one break end for normal (to
      // be more sensitive to germline) and if it spans both ends for deletion (should be next to 
      // each other), or one end for insertions larger than 10, or this is a complex breakpoint

      if (valid) { 
	
	std::string qn = j.Qname();

	// if read seen and other read was other mate, then reject
	if (num_align > 1 && qname_and_num.count(qn) && qname_and_num[qn] != j.FirstFlag()) {
	  // need to reject all reads of this qname
	  // because we saw both first and second in pair hit same split
	  reject_qnames.insert(qn);
	  
	} else {

	  // if haven't seen this read, add here. If have, then dont because want to avoid re-setting first/second convention
	  // e.g. if we see a split read with first designation, we want to reject all reads with same qname if at any time
	  // we see one with second mate designation. If we don't have this conditional, we can get fooled if the order is 
	  // 1, 2, 1
	  if (!qname_and_num.count(qn)) 
	    qname_and_num[qn] = j.FirstFlag();
	  
	  // this is a valid read
	  if (valid)
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

      //if (issplit1I && validI)
      //	++b1.splitI[sample_id];
      //if (issplit2I && validI)
      //	++b2.splitI[sample_id];	

      
    } // end read loop

    // process valid reads
    for (auto& i : bav) {

      std::string sr = i.GetZTag("SR");
      if (valid_reads.count(sr)) {

	std::string qn = i.Qname();
	// i don't know why I wanted to keep this to large events...
	// changing back to only one split from a read pair
	if (qnames.count(qn)/* && span > 1000*/)
	  continue; // don't count support if already added and not a short event
	// check that it's not a bad 1, 2 split
	if (reject_qnames.count(qn)) {
	  continue; 
	}

	reads.push_back(i);
	split_reads.insert(sr);

	// get read ID
	//assert(sr.at(0) == 't' || sr.at(0) == 'n');
	//bool tumor_read = sr.at(0) == 't';
	std::string sample_id = sr.substr(0,4);

	// keep track of qnames of split reads
	qnames.insert(qn);
	allele[sample_id].supporting_reads.insert(sr);
	
      	++allele[sample_id].split;
      }
    }

      // adjust the alt count
    for (auto& i : allele) {
      i.second.indel = num_align == 1;
      i.second.__adjust_alt_counts();
      //i.second.alt = std::max((int)i.second.supporting_reads.size(), i.second.cigar);
    }

  }
  
  std::string BreakPoint::getHashString() const {
    
    bool isdel = insertion.length() == 0;
    //if (isdel) // del breaks are stored as last non-deleted base. CigarMap stores as THE deleted base
    //  pos1++;
    std::string st = std::to_string(b1.gr.chr) + "_" + std::to_string(b1.gr.pos1) + "_" + std::to_string(this->getSpan()) + (isdel ? "D" : "I");
    //std::string st = std::to_string(b1.gr.chr) + "_" + std::to_string(b1.gr.pos1); // + "_" + std::to_string(this->getSpan()) + (isdel ? "D" : "I");
    return st;
  }
  
  int BreakPoint::checkPon(const PONFilter * p) {
    
    // only built for indels for now
    if (num_align != 1)
      return 0;
    
    std::string chr = std::to_string(b1.gr.chr);
    
    std::string key1, key2, key3, key4, key5;
    std::string type = (insertion != "") ? "I" : "D";

    key1 = chr + "_" + std::to_string(b1.gr.pos1-1);// + "_" + std::to_string(span) + type;
    key2 = chr + "_" + std::to_string(b1.gr.pos1  );// + "_" + std::to_string(span) + type;
    key3 = chr + "_" + std::to_string(b1.gr.pos1+1);// + "_" + std::to_string(span) + type;
    key4 = chr + "_" + std::to_string(b1.gr.pos1-2);// + "_" + std::to_string(span) + type;
    key5 = chr + "_" + std::to_string(b1.gr.pos1+2);//+ "_" + std::to_string(span) + type;
    
    pon = std::max(p->NSamps(key1), pon);
    pon = std::max(p->NSamps(key2), pon);
    pon = std::max(p->NSamps(key3), pon);
    pon = std::max(p->NSamps(key4), pon);
    pon = std::max(p->NSamps(key5), pon);
    
    return pon;
    
  }
  
    std::ostream& operator<<(std::ostream& out, const BreakPoint& b) {
      
      if (b.isindel) {
	out << ">" << (b.insertion.size() ? "INS: " : "DEL: ") << b.getSpan() << " " << b.b1.gr;
	  //<< " T/N split: " << b.t.split << "/" << b.n.split << " T/N cigar: " 
          //  << b.t.cigar << "/" << b.n.cigar << " T/N Cov " << b.t.cov << "/" << b.n.cov << " DBSNP: " << rs_t;
	for (auto& i : b.allele)
	  out << " " << i.first << ":" << i.second.split;  
      } else {
	out << ": " << b.b1.gr.pointString() << " to " << b.b2.gr.pointString() << " SPAN " << b.getSpan();
	  //<< " T/N split: " << b.t.split << "/" << b.n.split << " T/N disc: " 
	  //  << b.dc.tcount << "/" << b.dc.ncount << " " << b.evidence;
	for (auto& i : b.allele)
	  out << " " << i.first << ":" << i.second.split;  
      }
      
      return out;
      
    }
  
  void BreakPoint::checkBlacklist(GRC &grv) {
    
    // only check for indels
    //if (num_align != 1)
    //  return;
    
    if (grv.findOverlapping(b1.gr) || grv.findOverlapping(b2.gr)) {
      blacklist = true;
    }
  }
  
  void BreakPoint::__set_homologies_insertions() {
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
  
  BreakEnd::BreakEnd(const BamRead& b) {
    gr.chr = b.ChrID(); 
    gr.pos1 = -1;
    gr.pos2 = -1;
    cpos = -1;
    mapq = b.MapQuality();
    chr_name = b.GetZTag("MC"); 
    assert(chr_name.length());
    assert(chr_name != "23");
    nm = std::max(b.GetIntTag("NM") - (int)b.MaxInsertionBases() - (int)b.MaxDeletionBases(), 0);
    as_frac = (double)b.GetIntTag("AS") / (double) b.NumMatchBases();
  }

  void BreakPoint::__combine_with_discordant_cluster(DiscordantClusterMap& dmap)
  {
    const int PAD = 50;
    GenomicRegion bp1 = b1.gr;
    GenomicRegion bp2 = b2.gr;
    bp1.pad(PAD);
    bp2.pad(PAD);

    for (auto& d : dmap)
      {
	if (!d.second.valid())
	  continue;

	bool bp1reg1 = bp1.getOverlap(d.second.m_reg1) > 0;
	bool bp2reg2 = bp2.getOverlap(d.second.m_reg2) > 0;

	bool s1 = bp1.strand == d.second.m_reg1.strand;
	bool s2 = bp2.strand == d.second.m_reg2.strand;
	
	int pos1 = d.second.m_reg1.strand == '+' ? d.second.m_reg1.pos2 : d.second.m_reg1.pos1; // get the edge of the cluster
	int pos2 = d.second.m_reg2.strand == '+' ? d.second.m_reg2.pos2 : d.second.m_reg2.pos1;

	bool pass = bp1reg1 && bp2reg2 && s1 && s2;

	// check that the ends are not way off
	if (std::abs(pos1 - b1.gr.pos1) > PAD || std::abs(pos2 - b2.gr.pos1) > PAD)
	  pass = false;
	
	//std::cerr << " cname " << cname << " pad " << PAD << " pass " << pass << " DC pos1 " << pos1 << " DC pos2 " << pos2 << 
	//  " s1 " << s1 << " s2 " << s2 << " bp1reg1 " << bp1reg1 << " bp2reg2 " << bp2reg2 << " diff1 " << std::abs(pos1 - b1.gr.pos1) <<
	//  " diff2 " << std::abs(pos2 - b2.gr.pos1) << " bp pos1 " << b1.gr.pos1 << " bp pos2 " << b2.gr.pos1 << std::endl;

	if (pass)
	  // check that we haven't already added a cluster to this breakpoint
	  // if so, chose the one with more tumor support
	  if (dc.isEmpty() || dc.tcount < d.second.tcount) {
	      dc = d.second;
	      d.second.m_contig = cname;

	      // add the counts
	      for (auto& c : d.second.counts) {
		allele[c.first].disc = c.second;
		allele[c.first].indel = num_align == 1;
	      }

	      // add the discordant reads names to supporting reads for each sampleinfo
	      for (auto& rr : d.second.reads) {
		std::string sr = rr.second.GetZTag("SR");
		allele[sr.substr(0,4)].supporting_reads.insert(sr);
	      }
	      for (auto& rr : d.second.mates) {
		std::string sr = rr.second.GetZTag("SR");
		allele[sr.substr(0,4)].supporting_reads.insert(sr);
	      }

	      // adjust the alt counts
	      for (auto& aa : allele)
		aa.second.__adjust_alt_counts();
	      //aa.second.alt = aa.second.supporting_reads.size();

	  } 
	
      }
    
  }

  void BreakPoint::__set_evidence() {

    bool isdisc = (dc.tcount + dc.ncount) != 0;
    //bool issplit = (t.split + n.split) != 0;

    if (num_align == 1)
      evidence = "INDEL";
    else if ( isdisc && !complex && num_align != 0) //num_align == 2 )
      evidence = "ASDIS";
    else if ( isdisc && num_align < 3)
      evidence = "DSCRD";
    else if (!complex) //num_align == 2)
      evidence = "ASSMB";
    else if (complex) //num_align > 2)
      evidence = "COMPL";
    else {
      //std::cerr << "cname " << cname << " num_align" << num_align << " isdisc " << " issplit "  << std::endl;
      //assert(false);
      
    }

    assert(evidence.length());

  }

  void BreakPoint::__score_assembly_only() {

    // get read length
    int readlen = 0;
    if (reads.size())
      for (auto& i : reads)
	readlen = std::max(i.Length(), readlen);

    int span = getSpan();
    int num_split = t.split + n.split;
    int cov_span = split_cov_bounds.second - split_cov_bounds.first ;

    // check for high repeats
    bool hi_rep = false;
    for (auto& rr : hirepr)
      if (seq.find(rr) != std::string::npos)
	hi_rep = true;

    if (!b1.local && !b2.local) // added this back in snowman71. 
      // issue is that if a read is secondary aligned, it could be 
      // aligned to way off region. Saw cases where this happend in tumor
      // and not normal, so false-called germline event as somatic.
      confidence = "NOLOCAL";
    if (has_local_alignment)
      confidence = "LOCALMATCH";
    else if ( num_split > 1 && ( (cov_span <= (readlen + 5 ) && cov_span > 0) || cov_span < 0) )
      confidence = "DUPREADS"; // the same sequences keep covering the split
    else if (homology.length() >= 20 && (span > 1500 || span == -1) && std::max(b1.mapq, b2.mapq) < 60)
      confidence = "NODISC";
    else if ((int)seq.length() < readlen + 30)
      confidence = "TOOSHORT";
    else if (blacklist)
      confidence = "BLACKLIST";
    else if (a.split < 7 && (span > 1500 || span == -1))  // large and inter chrom need 7+
      confidence = "NODISC";
    else if (std::max(b1.mapq, b2.mapq) <= 40 || std::min(b1.mapq, b2.mapq) <= 10) 
      confidence = "LOWMAPQ";
    else if ( std::min(b1.mapq, b2.mapq) <= 30 && a.split <= 8 ) 
      confidence = "LOWMAPQ";
    else if (std::max(b1.nm, b2.nm) >= 10 || std::min(b1.as_frac, b2.as_frac) < 0.8) 
      confidence = "LOWAS";
    else if ( (std::max(b1.nm, b2.nm) >= 3 || std::min(b1.as_frac, b2.as_frac) < 0.85) && getSpan() < 0 )
      confidence = "LOWAS";      
    else if (seq.length() - aligned_covered > 10)
      confidence = "LOWAS";        
    else if ( (b1.matchlen < 50 && b1.mapq < 60) || (b2.matchlen < 50 && b2.mapq < 60) )
      confidence = "LOWMAPQ";
    else if ( std::min(b1.nm, b2.nm) >= 10)
      confidence = "LOWMAPQ";
    else if (a.split <= 3 && span <= 1500 && span != -1) // small with little split
      confidence = "LOWSPLITSMALL";
    else if (num_align == 2 && b1.gr.chr != b2.gr.chr && std::min(b1.matchlen, b2.matchlen) < 60) // inter-chr, but no disc reads, weird alignment
      confidence = "LOWICSUPPORT";
    else if (num_align == 2 && b1.gr.chr != b2.gr.chr && std::max(b1.nm, b2.nm) >= 3 && std::min(b1.matchlen, b2.matchlen) < 150) // inter-chr, but no disc reads, and too many nm
      confidence = "LOWICSUPPORT";
    else if (std::min(b1.matchlen, b2.matchlen) < 0.6 * readlen)
      confidence = "LOWICSUPPORT";      
    else if (num_align == 2 && std::min(b1.mapq, b2.mapq) < 50 && b1.gr.chr != b2.gr.chr) // interchr need good mapq for assembly only
      confidence = "LOWMAPQ";
    else if (std::min(b1.matchlen, b2.matchlen) < 40 || (complex_local && std::min(b1.matchlen, b2.matchlen) < 100)) // not enough evidence
      confidence = "LOWMATCHLEN";    
    else if (std::min(b1.matchlen - homology.length(), b2.matchlen - homology.length()) < 40)
      confidence = "LOWMATCHLEN";          
    else if ((b1.sub_n && b1.mapq < 50) || (b2.sub_n && b2.mapq < 50)) // || std::max(b1.sub_n,b2.sub_n) >= 2)
      confidence = "MULTIMATCH";
    else if (secondary && std::min(b1.mapq, b2.mapq) < 30)
      confidence = "SECONDARY";
    else if ((repeat_seq.length() >= 10 && std::max(t.split, n.split) < 7) || hi_rep)
      confidence = "WEAKSUPPORTHIREP";
    else if (num_split < 6 && getSpan() < 300 && b1.gr.strand==b2.gr.strand) 
      confidence = "LOWQINVERSION";
    else if ( (b1.matchlen - b1.simple < 15 || b2.matchlen - b2.simple < 15) )
      confidence = "SIMPLESEQUENCE";
    else if ((int)homology.length() * HOMOLOGY_FACTOR > readlen) // if homology is too high, tough to tell from mis-assemly
      confidence = "HIGHHOMOLOGY";
	      
    else
      confidence = "PASS";

    assert(confidence.length());

  }
  
  void BreakPoint::__score_somatic(double NODBCUTOFF, double DBCUTOFF) {
    
    double ratio = n.alt > 0 ? (double)t.alt / (double)n.alt : 100;
    //double taf = t.cov > 0 ? (double)t.alt / (double)t.cov : 0;
    bool immediate_reject = (ratio <= 12 && n.cov > 10) || n.alt >= 2; // can't call somatic with 3+ normal reads or <5x more tum than norm ALT

    double somatic_ratio = 0;
    size_t ncount = 0;
  
  if (evidence == "INDEL") {
    
    somatic_lod = n.LO_n;

    if (rs.empty())
      somatic_score = somatic_lod > NODBCUTOFF && !immediate_reject;
    else
      somatic_score = somatic_lod > DBCUTOFF && !immediate_reject;	

    // reject if reall low AF and at DBSNP
    //if (taf < 0.2 && !rs.empty())
    //  somatic_score = 0;

  } else {
    
    somatic_lod = n.LO_n;

    // old model
    ncount = std::max(n.split, dc.ncount);
    somatic_ratio = ncount > 0 ? std::max(t.split,dc.tcount)/ncount : 100;
    
    somatic_score = somatic_ratio >= MIN_SOMATIC_RATIO && n.split < 2;
  }
  
  // kill all if too many normal support
  if (n.alt > 2)
    somatic_score = 0;
  
  // kill if single normal read in discordant clsuter
  if (evidence == "DSCRD" && n.alt > 0)
    somatic_score = 0;

  // kill if bad ratio
  double r = n.alt > 0 ? (double)t.alt / (double)n.alt : 100;
  if (r < 10)
    somatic_score = 0;

  //std::cout << " n.alt " << n.alt << " n.disc " << n.disc << " n.split " << n.split << " somatic_score " << somatic_score << " cname " << cname << " somatic ratio " << somatic_ratio << " ncount " << ncount << " r artio " r << std::endl;

}

  /*  void BreakPoint::__set_total_reads() {

    // total unique read support

    t_reads = 0;
    n_reads = 0;

    for (auto& a : allele) {
      if (a.first.at(0) == 't')
	t_reads += a.alt;
      else
	n_reads += a.alt;
    }

    //if (evidence == "ASDIS") {
      std::unordered_set<std::string> this_reads_t;
      std::unordered_set<std::string> this_reads_n;
      for (auto& i : split_reads) 
	if (!this_reads_t.count(i) && i.at(0) == 't')
	  this_reads_t.insert(i);
      for (auto& i : dc.reads)
	if (!this_reads_t.count(i.first) && i.first.at(0) == 't')
	  this_reads_t.insert(i.first);
      for (auto& i : split_reads)
	if (!this_reads_n.count(i) && i.at(0) == 'n')
	  this_reads_n.insert(i);
      for (auto& i : dc.reads)
	if (!this_reads_n.count(i.first) && i.first.at(0) == 'n')
	  this_reads_n.insert(i.first);
      
      t_reads = this_reads_t.size();
      n_reads = this_reads_n.size();

  }
  */

  void BreakPoint::__score_assembly_dscrd() {

    int this_mapq1 = b1.mapq;
    int this_mapq2 = b2.mapq;
    int span = getSpan();
    bool germ = dc.ncount > 0 || n.split > 0;

    int max_a_mapq = std::max(this_mapq1, dc.mapq1);
    int max_b_mapq = std::max(this_mapq2, dc.mapq2);

    //int min_assm_mapq = std::min(this_mapq1, this_mapq2);
    //int max_assm_mapq = std::max(this_mapq1, this_mapq2);
    
    // get read length
    int readlen = 0;
    if (reads.size())
      for (auto& i : reads)
	readlen = std::max(i.Length(), readlen);

    // how much of contig is covered by split reads
    int cov_span = split_cov_bounds.second - split_cov_bounds.first ;

    // set the total number of supporting reads for tumor / normal
    // these alt counts should already be one qname per alt (ie no dupes)
    int t_reads = 0;
    int n_reads = 0;
    for (auto& a : allele) {
      if (a.first.at(0) == 't')
	t_reads += a.second.alt;
      else
	n_reads += a.second.alt;
    }

    int total_count = t_reads + n_reads; //n.split + t.split + dc.ncount + dc.tcount;
    int disc_count = dc.tcount + dc.ncount;
    int hq = dc.tcount_hq + dc.ncount_hq;

    if ( (max_a_mapq < 30 && !b1.local && hq < 3) || (max_b_mapq < 30 && !b2.local && hq < 3) || (b1.sub_n > 7 && b1.mapq < 10 && !b1.local && hq < 3) || (b2.sub_n > 7 && b2.mapq < 10 && !b2.local && hq < 3) )
      confidence = "LOWMAPQ";
    else if ( std::min(b1.nm, b2.nm) >= 10)
      confidence = "LOWMAPQ";
    else if ( std::min(b1.mapq, b2.mapq) < 10/* && std::min(dc.mapq1, dc.mapq2) < 10 */&& hq < 2)
      confidence = "LOWMAPQ";      
    else if ( total_count < 4 || (std::max(t.split, n.split) <= 5 && cov_span < (readlen + 5) && disc_count < 7) )
      confidence = "LOWSUPPORT";
    else if ( total_count < 15 && germ && span == -1) // be super strict about germline interchrom
      confidence = "LOWSUPPORT";
    else if ( std::min(b1.matchlen, b2.matchlen) < 50 && b1.gr.chr != b2.gr.chr ) 
      confidence = "LOWICSUPPORT";
    else if (secondary && getSpan() < 1000) // local alignments are more likely to be false for alignemnts with secondary mappings
      confidence = "SECONDARY";	
    else if (dc.tcount_hq + dc.ncount_hq < 3) { // multimathces are bad if we don't have good disc support too
      if ( ((b1.sub_n && dc.mapq1 < 1) || (b2.sub_n && dc.mapq2 < 1))  )
	confidence = "MULTIMATCH";
      else if ( ( (secondary || b1.sub_n > 1) && !b1.local) && ( std::min(max_a_mapq, max_b_mapq) < 30 || std::max(dc.tcount, dc.ncount) < 10)) 
	confidence = "SECONDARY";
      else 
	confidence = "PASS";
    }
    else
      confidence = "PASS";
  }

  void BreakPoint::__score_dscrd(int min_dscrd_size) {

    t.alt = dc.tcount;
    n.alt = dc.ncount;

    int disc_count = dc.ncount + dc.tcount;
    int hq_disc_count = dc.ncount_hq+ dc.tcount_hq;
    int disc_cutoff = 6;
    int hq_disc_cutoff = disc_count >= 10 ? 3 : 5; // reads with both pair-mates have high MAPQ

    if (getSpan() >0 &&  (getSpan() < min_dscrd_size && b1.gr.strand == '+' && b2.gr.strand == '-')) // restrict span for del (FR) type 
      confidence = "LOWSPANDSCRD";
    else if (hq_disc_count < hq_disc_cutoff) // && getSpan() > 1e5 && (std::min(dc.mapq1, dc.mapq2) < 30 || std::max(dc.mapq1, dc.mapq2) <= 30)) // mapq here is READ mapq (37 std::max)
      confidence = "LOWMAPQDISC";
    else if (getSpan() <= 1e4 && (std::min(dc.mapq1, dc.mapq2) < 25))
      confidence = "LOWMAPQDISC";
    else if (!dc.m_id_competing.empty())
      confidence = "COMPETEDISC";
    else if (disc_count >= disc_cutoff || hq_disc_count >= hq_disc_cutoff) // && n.cov >= 10 && n.cov <= 200 && t.cov >= 10 && t.cov <= 1000)
      confidence = "PASS";
    else if ( disc_count < disc_cutoff) // || (dc.ncount > 0 && disc_count < 15) )  // be stricter about germline disc only
      confidence = "WEAKDISC";
    //else if ( getSpan() < min_span && disc_count < 12)
    //  confidence = "LOWSPAN";
    //else 
    //  confidence = "PASS";
    
    assert(confidence.length());
  }

  void BreakPoint::__score_indel(double LOD_CUTOFF) {

    assert(b1.mapq == b2.mapq);
    
    //double ratio = n.alt > 0 ? t.alt / n.alt : 100;

    double max_lod = 0;
    for (auto& s : allele) {
      max_lod = std::max(max_lod, s.second.LO);
    }

    double af_t = t.cov > 0 ? (double)t.alt / (double)t.cov : 0;
    double af_n = n.cov > 0 ? (double)n.alt / (double)n.cov : 0;
    double af = std::max(af_t, af_n);

    if (b1.mapq < 10) // || std::max(b1.nm, b2.nm) > 6)
      confidence="LOWMAPQ";
    else if (max_lod < LOD_CUTOFF) 
      confidence = "LOWLOD";
    else if (af < 0.05) // if really low AF, get rid of 
      confidence = "VLOWAF";
    else if (std::min(left_match, right_match) < 20) 
      confidence = "SHORTALIGNMENT"; // no conf in indel if match on either side is too small
    else
      confidence="PASS";
    
  }

  void BreakPoint::scoreBreakpoint(double LOD_CUTOFF, double DBCUTOFF, double NODBCUTOFF, double LRCUTOFF, int min_dscrd_size) {
    
    // set the evidence (INDEL, DSCRD, etc)
    __set_evidence();
    
    // ensure that breakpoint is oriented
    assert(valid()); 

    // 
    double er = repeat_seq.length() > 10 ? 0.04 : ERROR_RATES[repeat_seq.length()];
    if (evidence == "INDEL" || true)
      for (auto& i : allele) {
	i.second.modelSelection(er);
      }
    __combine_alleles();

    // kludge. make sure we have included the DC counts (should have done this arleady...)
    if (evidence == "DSCRD" || evidence == "ASDIS") {
      t.disc = dc.tcount;
      n.disc = dc.ncount;
    }

    // scale the LOD for MAPQ    
    int mapqr1 = b1.local ? std::max(30, b1.mapq) : b1.mapq; // if local, don't drop below 30
    int mapqr2 = b2.local ? std::max(30, b2.mapq) : b2.mapq; // if local (aligns to within window), don't drop below 30
    double scale = (double)( std::min(mapqr1, mapqr2) - 2 * b1.nm) / (double)60;
    for (auto& i : allele) 
      i.second.SLO = i.second.LO * scale;
    t.LO = t.SLO * scale;
    n.LO = n.SLO * scale;
    a.LO = a.SLO * scale;
   
    // sanity check
    int split =0;
    for (auto& i : allele) {
      split += i.second.split;
    }
    assert( (split == 0 && t.split == 0 && n.split==0) || (split > 0 && (t.split + n.split > 0)));


    //__set_total_reads();

    // do the scoring
    if (evidence == "ASSMB" || (evidence == "COMPL" && (dc.ncount + dc.tcount)==0))
      __score_assembly_only();
    if (evidence == "ASDIS" || (evidence == "COMPL" && (dc.ncount + dc.tcount))) 
      __score_assembly_dscrd();
    if (evidence == "DSCRD")
      __score_dscrd(min_dscrd_size);
    if (evidence == "ASDIS" && confidence != "PASS") {
      evidence = "DSCRD";
      __score_dscrd(min_dscrd_size);
    }

    if (evidence == "INDEL") 
      __score_indel(LOD_CUTOFF);
    // else {
    //   std::cerr << "evidence " << evidence << std::endl;
    //   std::cerr << "BreakPoint not classified. Exiting" << std::endl;
    //   exit(EXIT_FAILURE);
    // }

    __score_somatic(NODBCUTOFF, DBCUTOFF);

    if (confidence == "PASS")
      quality = 99;
    else
      quality = 0;
    
    double LR = -1000000;
    for (auto& i : allele) {
      LR = std::max(-i.second.LO_n, LR); 
      // neg because want to evaluate likelihood that is ALT in normal 
      // the original use was to get LL that is REF in normal (in order to call somatic)
      // but for next filter, we want to see if we are confident in the germline call
    }

    // for the germline hits, check now that they have a high LR score (that they are AF of 0.5+, as expected for germline)
    if (!somatic_score && num_align == 1) {
      if (LR < LRCUTOFF)
	confidence = "GERMLOWAF";
    }

    assert(getSpan() > -2);

  }

  void BreakPoint::order() {

    if (b1.gr < b2.gr)
      return;
    
    std::swap(b1, b2);
    
  }

  std::string BreakPoint::__format_readname_string() {
    
    std::string supporting_reads = "";
    std::unordered_map<std::string, bool> supp_reads;
    
    //add the discordant reads
    for (auto& r : dc.reads) {
      std::string tmp = r.second.GetZTag("SR");
      supp_reads[tmp] = true;
    }
    //for (auto& r : dc.mates) {
    //  std::string tmp = r.second.GetZTag("SR");
    //  supp_reads[tmp] = true;
    //}
    
    //add the reads from the breakpoint
    for (auto& r : reads) {
      std::string tmp = r.GetZTag("SR");
      supp_reads[tmp] = true;
    }
    
    // print reads to a string, delimit with a ,
    size_t lim = 0;
    for (auto& i : supp_reads) {
      if (++lim > 50)
	break;

      //boost::regex regr(".*?_[0-9]+_(.*)"); 
      //boost::cmatch pmatch;
      if (i.first.at(0) == 't') {
	size_t posr = i.first.find("_", 5) + 1;
	if (posr != std::string::npos)
	  supporting_reads = supporting_reads + "," + i.first.substr(posr, i.first.length() - posr);
	//if (boost::regex_match(i.first.c_str(), pmatch, regr))
	//  supporting_reads = supporting_reads + "," + pmatch[1].str();
      }
    }
    if (supporting_reads.size() > 0)
      supporting_reads = supporting_reads.substr(1, supporting_reads.size() - 1); // remove first _
    
    if (read_names.length() == 0)
      read_names = supporting_reads;
    
    return read_names;
  }

  bool BreakPoint::valid() const {

    // debug
    return true;
    
    if (!(b1.gr.strand == '+' || b1.gr.strand == '-') || !(b2.gr.strand == '+' || b2.gr.strand == '-')) {
      std::cerr << "b1.strand " << b1.gr.strand << " b2.strand " << b2.gr.strand << std::endl;
      return false;
    }

    // b1 is less than b2, or the same but signifies inverted connection
    if ((b1.gr < b2.gr) || (b1.gr.chr == b2.gr.chr && b1.gr.pos1 == b2.gr.pos1)) 
      return true;
    
    std::cerr << b1.gr << " " << b2.gr << std::endl;
    return false;
  }

  void BreakPoint::setRefAlt(faidx_t * main_findex, faidx_t * viral_findex) {
    
    int len;

    assert(main_findex);

    if (evidence != "INDEL") {
      
      // get the reference for BP1
      char * ref1 = faidx_fetch_seq(main_findex, const_cast<char*>(b1.chr_name.c_str()), b1.gr.pos1-1, b1.gr.pos1-1, &len);
      if (!ref1) {
	if (viral_findex)
	  ref1 = faidx_fetch_seq(viral_findex, const_cast<char*>(b1.chr_name.c_str()), b1.gr.pos1-1, b1.gr.pos1-1, &len);
      }
      if (!ref1) {
	std::cerr << "couldn't find reference on BP1 for ref " << b1.chr_name << " in either viral or human" << std::endl;
      }
      
      char * ref2 = faidx_fetch_seq(main_findex, const_cast<char*>(b2.chr_name.c_str()), b2.gr.pos1-1, b2.gr.pos1-1, &len);
      if (!ref2) {
	if (viral_findex)
	  ref2 = faidx_fetch_seq(viral_findex, const_cast<char*>(b2.chr_name.c_str()), b2.gr.pos1-1, b2.gr.pos1-1, &len);
      } 
      if (!ref2) {
	std::cerr << "couldn't find reference on BP2 for ref " << b2.chr_name << " in either viral or human" << std::endl;
      }

      // couldn't find something, so bail
      if (!ref1 || !ref2)
	return;

      // by convention, set ref to 1 and alt to 2. Gets sorted in VCF creation
      ref = std::string(ref1);
      alt = std::string(ref2);

      if (!ref.length())
	ref = "N";
      if (!alt.length())
	alt = "N";

      if (ref1)
	free(ref1);
      if (ref2)
	free(ref2);
      
    } else {

      if (insertion.length() && !insertion.empty()) {

	// reference
	char * refi = faidx_fetch_seq(main_findex, const_cast<char*>(b1.chr_name.c_str()), b1.gr.pos1-1, b1.gr.pos1-1, &len);
	if (!refi) {
	  std::cerr << "couldn't find reference sequence for ref " << b1.chr_name << " for indel, on human reference " << std::endl;
	  return;
	}
	ref = std::string(refi);
	if (!ref.length()) {
	  ref = "N";
	}

	// alt 
	alt = ref + insertion;
      
	if (refi)
	  free(refi);

      // deletion
      } else {	

	// reference
	assert(b2.gr.pos1 - b1.gr.pos1 - 1 >= 0);
	char * refi = faidx_fetch_seq(main_findex, const_cast<char*>(b1.chr_name.c_str()), b1.gr.pos1-1, b2.gr.pos1-2, &len);
	if (!refi) {
	  std::cerr << "couldn't find reference sequence for ref " << b1.chr_name << " for indel, on human reference " << std::endl;
	  return;
	}
	ref = std::string(refi);
	if (!ref.length())
	  ref = "N";
	alt = ref.substr(0,1);
	if (refi)
	  free(refi);
      }
    }
  }

  int BreakPoint::getSpan() const { 
    if (num_align == 1 && insertion.empty()) {// deletion
      return (abs((int)b1.gr.pos1 - (int)b2.gr.pos1) - 1);
    }
    if (num_align == 1) 
      return (insertion.length()); // insertion
    if (b1.gr.chr == b2.gr.chr)
      return abs((int)b1.gr.pos1-(int)b2.gr.pos1);
    else
      return -1;
  }

  ReducedBreakPoint::ReducedBreakPoint(const std::string &line, bam_hdr_t* h) {
    
    if (!h) {
      std::cerr << "ReducedBreakPoint::ReducedBreakPoint - Must supply non-empty header" << std::endl;
      exit(EXIT_FAILURE);
    }

    std::istringstream iss(line);
    std::string val;
    size_t count = 0;

    ref = nullptr;
    alt = nullptr;
    cname = nullptr;
    evidence = nullptr;
    confidence = nullptr;
    insertion = nullptr;
    homology = nullptr;

    //float afn, aft;
    std::string ref_s, alt_s, cname_s, insertion_s, homology_s, evidence_s, confidence_s, read_names_s;
    
    std::string chr1, pos1, chr2, pos2, chr_name1, chr_name2, repeat_s; 
    char strand1 = '*', strand2 = '*';
    while (std::getline(iss, val, '\t')) {
      try{
	switch(++count) {
	case 1: chr1 = val; chr_name1 = val; break;
	case 2: pos1 = val; break; 
	case 3: assert(val.length()); strand1 = val.at(0); break;
	case 4: chr2 = val; chr_name2 = val; break;
	case 5: pos2 = val; break; 
	case 6: assert(val.length()); strand2 = val.at(0); break;
	case 7: 
	  ref_s = val;
	  break; 
	case 8: 
	  alt_s = val;
	  break;
	case 9: break; //span = stoi(val); break; // automatically calculated
	case 10: 
	  b1 = ReducedBreakEnd(GenomicRegion(chr1, pos1, pos1, h), std::stoi(val), chr_name1); b1.gr.strand = strand1; break;
	case 11:
	  b2 = ReducedBreakEnd(GenomicRegion(chr2, pos2, pos2, h), std::stoi(val), chr_name2); b2.gr.strand = strand2; break;
	case 12: b1.nm = std::stoi(val); break;
	case 13: b2.nm = std::stoi(val); break;
	case 16: b1.sub_n = std::min((int)255, std::stoi(val)); break;
	case 17: b2.sub_n = std::min((int)255, std::stoi(val)); break;
	  //case 14: dc.ncount = std::min((int)255, std::stoi(val)); break;
	  //case 15: dc.tcount = std::min((int)255,std::stoi(val)); break;
	case 14: dc.mapq1 = std::stoi(val); break;  
	case 15: dc.mapq2 = std::stoi(val); break;  
	case 18: 
	  homology_s = val;
	  break; 
	case 19: 
	  //insertion_size = (val == "x") ? 0 : val.length();
	  insertion_s = val;
	  break; 
	case 20: cname_s = val; break;
	case 21: num_align = std::min((int)31, std::stoi(val)); break;
	case 22: 
	  pass = val == "PASS";
	  confidence_s = val;
	  break;
	case 23: 
	  evidence_s = val;
	  indel = val == "INDEL"; 
	  imprecise = val == "DSCRD"; 
	  break; 
	case 24: quality = std::min((int)255,std::stoi(val)); break;
	case 25: secondary = val == "1" ? 1 : 0;
	case 26: somatic_score = std::stod(val); break;
	case 27: somatic_lod = std::stod(val); break;
	case 28: true_lod = std::stod(val); break;
	case 29: pon = std::min(255,std::stoi(val)); break;
	case 30: repeat_s = val; break; // repeat_seq
	case 31: blacklist = (val=="1" ? 1 : 0); break;
	case 32: dbsnp = val != "x"; break;
	case 33: read_names_s = val; break; //reads
	default:
	  format_s.push_back(val);
	}

      } catch(...) {
	std::cerr << "caught stoi/stod/stof error on: " << val << " for count " << count << std::endl;
	std::cerr << line << std::endl;
	exit(1);
      }
    }
    
    confidence = __string_alloc2char(confidence_s, confidence);
    evidence   = __string_alloc2char(evidence_s, evidence);
    insertion  = __string_alloc2char(insertion_s, insertion);
    homology   = __string_alloc2char(homology_s, homology);
    cname      = __string_alloc2char(cname_s, cname);
    ref        = __string_alloc2char(ref_s, ref);
    alt        = __string_alloc2char(alt_s, alt);
    repeat     = repeat_s.empty() ? nullptr : __string_alloc2char(repeat_s, repeat);
    if (somatic_score && confidence_s == "PASS")
      read_names     = read_names_s;// == "x" ? nullptr : __string_alloc2char(read_names_s, read_names);

  }

  double SampleInfo::__log_likelihood(int ref, int alt, double f, double e) {
    
    // mutect log liklihood against error
    double ll = 0;
    const double back_mutate_chance = 0; // this should be zero? assume indel never accidentally back mutates
    ref = ref <= 0 ? 0 : ref;

    double arg1 = f * e * back_mutate_chance  /* p(alt mut to ref) */ + (1-f)  * (1-e); /* p(ref not mutated) */ 
    if (arg1 > 0)
      ll += ref * log10(arg1); // ref
    //else //if (ref == 0)// get rid of NaN issue
    //  ll = 0;

    arg1 = f * (1 - e) /* p(alt not mut) */ + (1-f) * e; /* p (ref mut to alt) */
    if (arg1 > 0)
      ll += alt * log10(arg1); // alt
    //else
    // ll = 0

    return ll;

  }

  void SampleInfo::modelSelection(double er) {

    alt = std::max(alt, cigar);

    if (alt >= cov) {
      cov = alt;
    }

    // make sure we get correct alt
    //__adjust_alt_counts();
    //alt = std::max(alt, std::max((int)supporting_reads.size(), cigar));

    af = cov > 0 ? (double)alt / (double)cov : 1;
    af = af > 1 ? 1 : af;

    // mutect log liklihood against error
    double ll_alt = __log_likelihood(cov - alt, alt, af    , er);
    double ll_err = __log_likelihood(cov - alt, alt, 0, er);
    LO = ll_alt - ll_err;

    //mutetct log likelihood normal
    //er = 0.0005; // make this low, so that ALT in REF is rare and NORM in TUM gives low somatic prob
    // actually, dont' worry about it too much. 3+ alt in ref excludes somatic anyways.
    double ll_alt_norm = __log_likelihood(cov - alt, alt, std::max(0.5, af), er); // likelihood that variant is 0.5 or above
    double ll_ref_norm = __log_likelihood(cov - alt, alt, 0  , er);
    LO_n = ll_ref_norm - ll_alt_norm; // higher number means more likely to be REF than ALT

    //std::cerr << LO_n << " ll_ref_norm " << ll_ref_norm << " ll_alt_norm " << ll_alt_norm << 
    //  " cov-al " << (cov-alt) << " alt " << alt << " er " << er << " test 100 " << 
    //  __log_likelihood(100, 100, 0, 0) << " " << __log_likelihood(100, 100, 0.5, 0) << std::endl;

    // genotype calculation as provided in 
    // http://bioinformatics.oxfordjournals.org/content/early/2011/09/08/bioinformatics.btr509.full.pdf+html
    int scaled_cov = std::floor((double)cov * 0.90);
    int this_alt = std::min(alt, scaled_cov);

    genotype_likelihoods[0] = __genotype_likelihoods(2, er, this_alt, scaled_cov); // ref/ref
    genotype_likelihoods[1] = __genotype_likelihoods(1, er, this_alt, scaled_cov);
    genotype_likelihoods[2] = __genotype_likelihoods(0, er, this_alt, scaled_cov);

    double max_likelihood = *std::max_element(genotype_likelihoods.begin(), genotype_likelihoods.end());
    if (max_likelihood == genotype_likelihoods[0])
      genotype = "0/0";
    else if (max_likelihood == genotype_likelihoods[1])
      genotype = "0/1";
    else 
      genotype = "1/1";

    // scale to max
    genotype_likelihoods[0] = max_likelihood - genotype_likelihoods[0];
    genotype_likelihoods[1] = max_likelihood - genotype_likelihoods[1];
    genotype_likelihoods[2] = max_likelihood - genotype_likelihoods[2];

    // get the genotype quality
    GQ = 99;
    for (auto& g : genotype_likelihoods)
      if (g != 0)
	GQ = std::min(GQ, __myround(g));

    // set the PL string
    std::stringstream sss;
    sss << std::setprecision(4) << __myround(genotype_likelihoods[0]) << "," << __myround(genotype_likelihoods[1]) << "," << __myround(genotype_likelihoods[2]);
    PL = sss.str();

  }

  void BreakPoint::addCovs(const std::unordered_map<std::string, STCoverage*>& covs, const std::unordered_map<std::string, STCoverage*>& clip_covs) {

    // setting to min coverage because we want to accurately genotype indels
    for (auto& i : covs) 
      allele[i.first].cov = std::min(i.second->getCoverageAtPosition(b1.gr.chr, b1.gr.pos1), i.second->getCoverageAtPosition(b2.gr.chr, b2.gr.pos1));
    for (auto& i : clip_covs)
      allele[i.first].clip_cov = std::min(i.second->getCoverageAtPosition(b1.gr.chr, b1.gr.pos1), i.second->getCoverageAtPosition(b2.gr.chr, b2.gr.pos1));
  }

  std::ostream& operator<<(std::ostream& out, const BreakEnd& b) {
    out << b.gr << " - " << b.id << " mapq " << b.mapq << " subn " << b.sub_n << std::endl;
    return out;
  }
  
  std::ostream& operator<<(std::ostream& out, const SampleInfo& a) {
    out << " split: " << a.split << " cigar " << a.cigar << " alt " << a.alt << " clip_cov " << a.clip_cov << " cov " << a.cov << " disc " << a.disc;
    return out;
  }

  SampleInfo operator+(const SampleInfo& a1, const SampleInfo& a2) {

    SampleInfo a;

    a.disc = a1.disc + a2.disc;
    a.split = a1.split + a2.split;
    a.cigar = a1.cigar + a2.cigar;
    a.clip_cov = a1.clip_cov + a2.clip_cov;
    a.cov = a1.cov + a2.cov;

    // add the reads
    for (auto& i : a1.supporting_reads)
      a.supporting_reads.insert(i);
    for (auto& i : a2.supporting_reads)
      a.supporting_reads.insert(i);
    
    //a.alt = std::max((int)a.supporting_reads.size(), a.cigar);
    a.__adjust_alt_counts();

    return a;
  }

  void BreakPoint::__combine_alleles() {

    for (auto& s : allele) {
      if (s.first.at(0) == 't') {
	t = t + s.second;
      } else {
	n = n + s.second;
      }
    }

    a = t + n;
    
    double er = repeat_seq.length() > 10 ? 0.04 : ERROR_RATES[repeat_seq.length()];
    //a.modelSelection(er);
    t.modelSelection(er);
    n.modelSelection(er);
    
  }


  std::string SampleInfo::toFileString() const {

    std::stringstream ss;

    if (indel)
      ss << std::setprecision(4) << genotype << ":" << 
	std::max(alt, cigar) << ":" << cov << ":" << GQ << ":" << PL << ":" << split << ":" << cigar 
	 << ":" << LO_n << ":" << LO << ":" << SLO;
    else
      ss << std::setprecision(4) << genotype << ":" << 
	alt << ":" << cov << ":" << GQ << ":" << PL << ":" << split
	 << ":" << disc << ":" << LO_n << ":" << LO << ":" << SLO;
          
    return ss.str();

  }

  void SampleInfo::fromString(const std::string& s) { 

    std::string val;
    int count = 0;
    std::istringstream input(s);
    
    while (std::getline(input, val, ':')) {
      
      // tmp fix
      if (val.find("nan") != std::string::npos)
	val = "0";

      //      std::cerr << val << std::endl;
      switch(++count) {
      case 6: split = std::stoi(val); break;
      case 8: LO_n =   std::stod(val); break;
      case 9: LO = std::stod(val); break;
      case 10: SLO = std::stod(val); break;
      case 7: 
	if (indel)
	  cigar = std::stoi(val);
	else
	  disc = std::stoi(val);
	break;
      case 2: alt = std::stoi(val); break;
      case 3: cov = std::stoi(val); break;      
      case 4: GQ = std::stod(val);  break;
      case 5: PL = val;  break;
      case 1: genotype = val; break;
      }
    }
    //assert(s == toFileString()); // tmp turn off
  }

  void BreakPoint::__rep(int rep_num, std::string& rseq, bool fwd) {
    
    // move left and right from breakend 1
    //int replen = 0;
    //int curr_replen = 1;
    rseq = "";

    // return if we are too close to the edge for some reason
    if (b1.cpos > (int)seq.length() - 5)
      return;

    const int REP_BUFF = 50;

    std::string sss = seq;
    if (!fwd)
      std::reverse(sss.begin(), sss.end());
    int cpos = b1.cpos + 1;
    if (!fwd)
      cpos = seq.length() - b1.cpos;

    int i = cpos; //b1.cpos + 1;
    int stop = std::min((int)seq.length() - 1, cpos + REP_BUFF); //fwd ? std::min((int)seq.length() - 1, b1.cpos + REP_BUFF) : std::max(0, b1.cpos - REP_BUFF);
    int end = -1;

    bool no_rep = true;

    assert(stop - i < 80);

    std::string c = sss.substr(cpos, std::min(rep_num, (int)seq.length() -  cpos));

    i += rep_num; //*fr;
    
    //for (; (i < stop && fwd)/* || (i > stop && !fwd)*/; i += rep_num*fr) {
    for (; i < stop; i+= rep_num) {

      if (i >= (int)seq.length() - 1 || i + rep_num > (int)seq.length()) // shouldn't happen, but ensures no substr errors
	break;

      // terminating
      if (c != sss.substr(i,rep_num)) { // || i >= (stop - rep_num)) {

	end = i;
	break;

      } else {

	no_rep = false;
	//curr_replen += rep_num;
      }

    }

    if (end == -1)
      end = stop;
	
    if (!no_rep)
      rseq = sss.substr(cpos, std::min(end - cpos, (int)seq.length() - cpos));
  }

  std::string BreakEnd::hash(int offset) const {
    
    return (std::to_string(gr.chr) + ":" + std::to_string(gr.pos1 + offset));
    
  }

  void BreakPoint::checkLocal(const GenomicRegion& window) {

    b1.checkLocal(window);
    b2.checkLocal(window);

  }

  void BreakEnd::checkLocal(const GenomicRegion& window) {
    
    if (gr.getOverlap(window))
      local = true;
  }

  void SampleInfo::__adjust_alt_counts() {
    
    // get unique qnames
    std::unordered_set<std::string> qn;
    for (auto& r : supporting_reads) {
      size_t posr = r.find("_", 5) + 1; // extract the qname from tXXX_XXX_QNAME
      qn.insert(r.substr(posr, r.length() - posr));
    }
    
    // alt count is max of cigar or unique qnames (includes split and discordant)
    alt = std::max((int)qn.size(), cigar);

  }

  // g is the number of reference alleles (e.g. g = 2 is homozygous reference)
  // assumes biallelic model
  double SampleInfo::__genotype_likelihoods(int g, double er, int alt, int cov) {
    double val =  - cov * log10(2) + (cov - alt) * log10( (2 - g) * er + g  * (1 - er) ) + alt * log10( (2 - g) * (1 - er) + g * er);
    return val;
  }

}
