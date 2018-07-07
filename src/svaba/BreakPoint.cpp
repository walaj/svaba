#include "BreakPoint.h"

#include <getopt.h>
#include <iomanip>
#include <cassert>

#include "gzstream.h"
#include "svabaUtils.h"

#include "svaba_params.h"

// define repeats
static std::vector<std::string> repr = {"AAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTT", 
					"CCCCCCCCCCCCCCCC", "GGGGGGGGGGGGGGGG",
				        "TATATATATATATATA", "ATATATATATATATAT", 
				        "GCGCGCGCGCGCGCGC", "CGCGCGCGCGCGCGCG", 
				        "TGTGTGTGTGTGTGTG", "GTGTGTGTGTGTGTGT", 
				        "TCTCTCTCTCTCTCTC", "CTCTCTCTCTCTCTCT", 
				        "CACACACACACACACA", "ACACACACACACACAC", 
				        "GAGAGAGAGAGAGAGA", "AGAGAGAGAGAGAGAG"};

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
static std::unordered_map<int, double> ERROR_RATES = {{0, scale_factor * 1e-4}, {1, scale_factor * 1e-4}, {2,  scale_factor * 1e-4}, {3,  scale_factor * 1e-4}, {4,  scale_factor * 1e-4}, {5, scale_factor * 2e-4}, {6, scale_factor * 5e-4}, {7, scale_factor * 1e-3},
						      {8, scale_factor * 2e-3}, {9, scale_factor * 3e-3}, {10, scale_factor * 1e-2}, {11, scale_factor * 2e-2}, {12, scale_factor * 3e-5}};


using namespace SeqLib;

  double __myround(double x) { return std:: floor(x * 10) / 10; }

  // make the file string
  std::string BreakPoint::toFileString(bool noreads) {
    
    // make sure we already ran scoring
    assert(evidence.length());
    assert(confidence.length());
    
    std::string sep = "\t";
    std::stringstream ss;
    
    // put the read names into a string
    if (!noreads)  
      format_readname_string();
    
    // make the BX table
    format_bx_string();

    double max_lod = 0;
    for (auto& s : allele) 
      max_lod = std::max(max_lod, s.second.LO);

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
       << (read_names.length() ? read_names : "x") << sep
       << (!bxtable.empty() ? bxtable : "x");

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
BreakPoint::BreakPoint(DiscordantCluster& tdc, const BWAWrapper * bwa, DiscordantClusterMap& dmap, 
		       const GenomicRegion& region) {
    
    num_align = 0;
    dc = tdc;
    
    std::string chr_name1, chr_name2;

    try {
       // if this throw error, it means that there are more chr in 
       // reads than in reference
       chr_name1 = bwa->ChrIDToName(dc.m_reg1.chr); //bwa->ChrIDToName(tdc.reads.begin()->second.ChrID());
       chr_name2 = bwa->ChrIDToName(dc.m_reg2.chr); //bwa->ChrIDToName(tdc.reads.begin()->second.ChrID());
    } catch (...) {
      std::cerr << "Warning: Found mismatch between reference genome and BAM genome for discordant cluster " << dc << std::endl;
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
      //std::string sr = SRTAG(rr.second);
      std::string sr = rr.second.SR();
      allele[rr.second.Prefix()].supporting_reads.insert(sr);
      alt_counts[rr.second.Prefix()].insert(rr.second.Qname());
    }
    for (auto& rr : dc.mates) {
      std::string sr = rr.second.SR();
      allele[rr.second.Prefix()].supporting_reads.insert(rr.second.Qname());
      alt_counts[rr.second.Prefix()].insert(rr.second.Qname());
    }
    
    // add the alt counts
    for (auto& i : allele) {
      i.second.indel = false;
      i.second.disc = alt_counts[i.first].size(); //allele[i.first].supporting_reads.size();
      i.second.adjust_alt_counts();
      //i.second.alt = i.second.disc;
    }
      
    // give a unique id
    cname = dc.toRegionString() + "__" + std::to_string(region.chr+1) + "_" + std::to_string(region.pos1) + 
      "_" + std::to_string(region.pos2) + "D";

    // check if another cluster overlaps, but different strands
    if (getSpan() > 800 || getSpan() == -1) // only check for large events.
    for (auto& d : dmap) {

      // don't overlap if on different chr, or same event
      if (dc.m_reg1.chr != d.second.m_reg1.chr || dc.m_reg2.chr != d.second.m_reg2.chr || d.second.ID() == dc.ID())
	continue;

      // isolate and pad
      GenomicRegion gr1 = d.second.m_reg1;
      GenomicRegion gr2 = d.second.m_reg2;
      gr1.Pad(100);
      gr2.Pad(100);
      
      if (dc.m_reg1.GetOverlap(gr1) && dc.m_reg2.GetOverlap(gr2))
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
  }

// used in the refilter module only
BreakPoint::BreakPoint(const std::string &line, const SeqLib::BamHeader& h) {

  if (h.isEmpty()) {
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
	case 34: bxtable = val; break;
        default: 
	  aaa.indel = evidence == "INDEL";
	  aaa.fromString(val);
	  id += "A";
	  allele[id] = aaa; //id is dummy. just keep in order as came in;
	  // fill in the discordant info
	  //if (id == "A") // tumor
	  //  dc.tcount = aaa.disc;
	  //else
	  //  dc.ncount += aaa.disc;
	  break;
	}
      } catch(...) {
	std::cerr << "caught stoi error on: " << val << std::endl;
	std::cerr << line << std::endl;
	exit(1);
      }
    }

  }
  
//void BreakPoint::splitCoverage(SeqLib::BamRecordVector &bav) {
  void BreakPoint::splitCoverage(svabaReadVector &bav) {
    
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
    for (auto& j : bav) {

      r2c& this_r2c = j.GetR2C(cname);

      bool read_should_be_skipped = false;
      if (num_align == 1) {

	std::vector<int> del_breaks;
	std::vector<int> ins_breaks;
	int pos = 0;
	
	// if this is a nasty repeat, don't trust non-perfect alignmentx on r2c alignment
	if (__check_homopolymer(j.Sequence())) 
	  read_should_be_skipped = true;
	
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
      std::string sample_id = j.Prefix(); //substr(0,4); // maybe just make this prefix
      std::string sr = j.SR(); 

      // need read to cover past variant by some buffer. If there is a repeat,
      // then this needs to be even longer to avoid ambiguity
      int this_tbuff = T_SPLIT_BUFF + repeat_seq.length();
      int this_nbuff = N_SPLIT_BUFF + repeat_seq.length();      

      int rightbreak1 = b1.cpos + (j.Tumor() ? this_tbuff : this_nbuff); // read must extend this far right of break1
      int leftbreak1  = b1.cpos - (j.Tumor() ? this_tbuff : this_nbuff); // read must extend this far left of break1
      int rightbreak2 = b2.cpos + (j.Tumor() ? this_tbuff : this_nbuff);
      int leftbreak2  = b2.cpos - (j.Tumor() ? this_tbuff : this_nbuff);

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
      bool valid  = (both_split && (j.Tumor() || homlen > 0)) || (one_split && !j.Tumor() && homlen == 0) || (one_split && insertion.length() >= INSERT_SIZE_TOO_BIG_SPAN_READS);
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
      
      r2c& this_r2c = i.GetR2C(cname);
      if (valid_reads.count(i.SR())) {

	std::string qn = i.Qname();
	if (qnames.count(qn))
	  continue; // don't count support if already added and not a short event
	// check that it's not a bad 1, 2 split
	if (reject_qnames.count(qn)) {
	  this_r2c.supports_var = false;
	  i.AddR2C(cname, this_r2c); // update that this actually does not support
	  continue; 
	}

	reads.push_back(i);

	// keep track of qnames of split reads
	qnames.insert(qn);
	allele[i.Prefix()].supporting_reads.insert(i.SR());
	
      	++allele[i.Prefix()].split;
      }
    }

      // adjust the alt count
    for (auto& i : allele) {
      i.second.indel = num_align == 1;
      i.second.adjust_alt_counts();
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
  
/*  int BreakPoint::checkPon(const PONFilter * p) {
    
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
*/

std::ostream& operator<<(std::ostream& out, const BreakPoint& b) {
  
      if (b.isindel) {
	out << ">" << (b.insertion.size() ? "INS: " : "DEL: ") << b.getSpan() << " " << b.b1.gr << " " << b.cname << " " << b.evidence;
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
      
    }
  
  void BreakPoint::checkBlacklist(GRC &grv) {
    if (grv.CountOverlaps(b1.gr) || grv.CountOverlaps(b2.gr)) 
      blacklist = true;
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
  
BreakEnd::BreakEnd(const SeqLib::BamRecord& b) {
  b.GetIntTag("SQ", sub_n);
  gr.chr = b.ChrID(); 
  gr.pos1 = -1;
  gr.pos2 = -1;
  cpos = -1;
  mapq = b.MapQuality();
  b.GetZTag("MC", chr_name); 
  assert(chr_name.length());

  int tnm=0;
  b.GetIntTag("NM",tnm);
  nm = std::max(tnm - (int)b.MaxInsertionBases() - (int)b.MaxDeletionBases(), 0);
  int thisas=0;
  b.GetIntTag("AS",thisas);
  as_frac = (double)thisas / (double) b.NumMatchBases();
}

  void BreakPoint::__combine_with_discordant_cluster(DiscordantClusterMap& dmap)
  {
    const int PAD = 50;
    GenomicRegion bp1 = b1.gr;
    GenomicRegion bp2 = b2.gr;
    bp1.Pad(PAD);
    bp2.Pad(PAD);

    for (auto& d : dmap) {

	if (!d.second.valid())
	  continue;

	bool bp1reg1 = bp1.GetOverlap(d.second.m_reg1) > 0;
	bool bp2reg2 = bp2.GetOverlap(d.second.m_reg2) > 0;

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

	      // add the read counts
	      for (auto& c : d.second.counts) {
		allele[c.first].disc = c.second;
		allele[c.first].indel = num_align == 1;
	      }

	      // add the discordant reads names to supporting reads for each sampleinfo
	      for (auto& rr : d.second.reads) {
		allele[rr.second.Prefix()].supporting_reads.insert(rr.second.SR());

	      }
	      for (auto& rr : d.second.mates) {
		allele[rr.second.Prefix()].supporting_reads.insert(rr.second.SR());
	      }

	      // adjust the alt counts
	      for (auto& aa : allele)
		aa.second.adjust_alt_counts();
	      //aa.second.alt = aa.second.supporting_reads.size();
	  } 
	
      }
    
  }

  void BreakPoint::set_evidence() {

    // if we are in refilter, then this is already set
    if (!evidence.empty())
      return;

    bool isdisc = (dc.tcount + dc.ncount) != 0;

    if (num_align == 1)
      evidence = "INDEL";
    else if ( isdisc && !complex && num_align > 0)
      evidence = "ASDIS";
    else if ( isdisc && num_align < 3)
      evidence = "DSCRD";
    else if (!complex) 
      evidence = "ASSMB";
    else if (complex && !complex_local) // is A-C of an ABC
      evidence = "TSI_G";
    else if (complex && complex_local) // is AB or BC of an ABC 
      evidence = "TSI_L";

    assert(evidence.length());

  }

  void BreakPoint::score_assembly_only() {

    int span = getSpan();
    int num_split = t.split + n.split;
    int cov_span = split_cov_bounds.second - split_cov_bounds.first ;

    // check for high repeats
    bool hi_rep = false;
    for (auto& rr : hirepr)
      if (seq.find(rr) != std::string::npos)
	hi_rep = true;

    if (!b1.local && !b2.local && !complex_local) // added this back in v71
      // issue is that if a read is secondary aligned, it could be 
      // aligned to way off region. Saw cases where this happend in tumor
      // and not normal, so false-called germline event as somatic.
      confidence = "NOLOCAL";
    else if (has_local_alignment)
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
    else if ((double)aligned_covered / (double)seq.length() < 0.80) // less than 80% of read is covered by some alignment
      confidence = "LOWAS";        
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
    else if (std::min(b1.matchlen, b2.matchlen) < 0.6 * readlen)
      confidence = "LOWICSUPPORT";      
    else if (std::min(b1.mapq, b2.mapq) < 50 && b1.gr.chr != b2.gr.chr) // interchr need good mapq for assembly only
      confidence = "LOWMAPQ";
    else if (std::min(b1.matchlen, b2.matchlen) < 40 || (complex_local && std::min(b1.matchlen, b2.matchlen) < 100)) // not enough evidence
      confidence = "LOWMATCHLEN";    
    else if (std::min(b1.matchlen - homology.length(), b2.matchlen - homology.length()) < 40)
      confidence = "LOWMATCHLEN";          
    else if ((b1.sub_n && b1.mapq < 50) || (b2.sub_n && b2.mapq < 50)) 
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
  
  void BreakPoint::score_somatic(double NODBCUTOFF, double DBCUTOFF) {

    // this is LOD of normal being REF vs AF = 0.5+
    // We want this to be high for a somatic call
    somatic_lod = n.LO_n;

    // find the somatic to normal ratio
    double ratio = n.alt > 0 ? (double)t.alt / (double)n.alt : 100;    
    
    if (evidence == "INDEL") {
      
      // somatic score is just true or false for now
      // use the specified cutoff for indels, taking into account whether at dbsnp site
      somatic_score = somatic_lod > ( (rs.empty() || rs=="x") ? NODBCUTOFF : DBCUTOFF);

    // can't call somatic with 5+ normal reads or <5x more tum than norm ALT
    //if ((ratio <= 12 && n.cov > 10) || n.alt > 5)
    //if (n.alt > 5)
    //  somatic_score = 0;

  // for SVs, just use a hard cutoff for gauging somatic vs germline
  } else {

    // require no reads in normal or 1 read and tons of tumor reads

    somatic_score = ratio >= MIN_SOMATIC_RATIO && n.split < 2 && dc.ncount < 2;
  }
    
  // set germline if single normal read in discordant clsuter
  if (evidence == "DSCRD" && n.alt > 0)
    somatic_score = 0;
  
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

  void BreakPoint::score_dscrd(int min_dscrd_size) {

    t.alt = dc.tcount;
    n.alt = dc.ncount;

    int disc_count = dc.ncount + dc.tcount;
    int hq_disc_count = dc.ncount_hq + dc.tcount_hq;
    int disc_cutoff = 8;
    int hq_disc_cutoff = disc_count >= 10 ? 3 : 5; // reads with both pair-mates have high MAPQ

    if (getSpan() > 0 && (getSpan() < min_dscrd_size && b1.gr.strand == '+' && b2.gr.strand == '-')) // restrict span for del (FR) type 
      confidence = "LOWSPANDSCRD";
    else if (hq_disc_count < hq_disc_cutoff && (disc_count < disc_cutoff || std::min(dc.mapq1, dc.mapq2) < 15))
      confidence = "LOWMAPQDISC";
    else if (!dc.m_id_competing.empty())
      confidence = "COMPETEDISC";
    else if ( disc_count < disc_cutoff)
      confidence = "WEAKDISC";
    else 
      confidence = "PASS";
    
    assert(confidence.length());
  }

void BreakPoint::score_indel(double LOD_CUTOFF, double LOD_CUTOFF_DBSNP) {

    assert(b1.mapq == b2.mapq);

    bool is_refilter = !confidence.empty(); // act differently if this is refilter run
    
    // for refilter, only consider ones that were low lod or PASS
    // ie ones that with a different lod threshold may be changed
    // if confidence is empty, this is original run so keep going
    if (confidence != "LOWLOD" && confidence != "PASS" && is_refilter)
      return;
    
    double max_lod = 0;
    for (auto& s : allele) 
      max_lod = std::max(max_lod, s.second.LO);

    // check if homozygous reference is most likely GT
    bool homozygous_ref = true;
    for (auto& s : allele) {
      if (s.second.genotype_likelihoods[0] > s.second.genotype_likelihoods[1] ||
	  s.second.genotype_likelihoods[0] > s.second.genotype_likelihoods[2])
	homozygous_ref = false;
    }
    
    // get the allelic ractions, just for VLOWAF filter
    double af_t = t.cov > 0 ? (double)t.alt / (double)t.cov : 0;
    double af_n = n.cov > 0 ? (double)n.alt / (double)n.cov : 0;
    double af = std::max(af_t, af_n);

    if (b1.mapq < 10) 
      confidence="LOWMAPQ";
    else if (!is_refilter && (double)aligned_covered / (double)seq.length() < 0.80) // less than 80% of read is covered by some alignment
      confidence = "LOWAS";  
    else if ((b1.sub_n && b1.mapq < 50) || (b2.sub_n && b2.mapq < 50)) 
      confidence = "MULTIMATCH";      
    else if (max_lod < LOD_CUTOFF && rs.empty())        // non db snp site
      confidence = "LOWLOD";
    else if (max_lod < LOD_CUTOFF_DBSNP && !rs.empty()) // be more permissive for dbsnp site
      confidence = "LOWLOD";
    else if (af < 0.05) // if really low AF, get rid of 
      confidence = "VLOWAF";
    else if (!is_refilter && std::min(left_match, right_match) < 20) 
      confidence = "SHORTALIGNMENT"; // no conf in indel if match on either side is too small
    else if (homozygous_ref)
      confidence = "NONVAR";
    else
      confidence="PASS";

  }

  // LOD cutoff is just whether this is ref / alt
  // LOD cutoff dbsnp is same, but at DBSNP site (should be lower threshold since we have prior)
  // LOD SOM cutoff is cutoff for whether NORMAL BAM(s) are AF = 0
  // LOD SOM DBSNP cutoff same as above, but at DBSNP site (should be HIGHER threshold,
  //   since we have prior that it's NOT somatic
void BreakPoint::scoreBreakpoint(double LOD_CUTOFF, double LOD_CUTOFF_DBSNP, double LOD_CUTOFF_SOMATIC, double LOD_CUTOFF_SOMATIC_DBSNP, double scale_errors, int min_dscrd_size) {
    
    // set the evidence (INDEL, DSCRD, etc)
    set_evidence();

    __combine_alleles();

    // 
    error_rate = repeat_seq.length() > 10 ? MAX_ERROR : ERROR_RATES[repeat_seq.length()];
    for (auto& i : allele) {
      i.second.readlen = readlen;
      i.second.modelSelection(error_rate);
    }

    // kludge. make sure we have included the DC counts (should have done this arleady...)
    if (evidence == "DSCRD" || evidence == "ASDIS") {
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
    
    // sanity check
    int split =0;
    for (auto& i : allele) 
      split += i.second.split;
    assert( (split == 0 && t.split == 0 && n.split==0) || (split > 0 && (t.split + n.split > 0)));

    // do the scoring
    bool iscomplex = evidence.find("TSI") != std::string::npos;
    if (confidence.empty() && evidence != "INDEL") {
      if (evidence == "ASSMB" || (iscomplex  && (dc.ncount + dc.tcount)==0))
	score_assembly_only();
      if (evidence == "ASDIS" || (iscomplex && (dc.ncount + dc.tcount))) 
	score_assembly_dscrd();
      if (evidence == "DSCRD")
	score_dscrd(min_dscrd_size);
      // it failed assembly filters, but might pass discordant filters
      if (evidence == "ASDIS" && confidence != "PASS") { 
	evidence = "DSCRD";
	score_dscrd(min_dscrd_size);
      }
    }
    else if (evidence == "INDEL") {
      score_indel(LOD_CUTOFF, LOD_CUTOFF_DBSNP);
    }

    // filter out SVs with only 1 BX tag supporting
    // if we have bx tags
    format_bx_string();
    if (bx_count == 1 && confidence == "PASS")
      confidence = "SINGLEBX";

    // set the somatic_score field to true or false
    score_somatic(LOD_CUTOFF_SOMATIC, LOD_CUTOFF_SOMATIC_DBSNP);

    // quality score is odds that read is non-homozygous reference (max 99)
    quality = 0;
    for (auto& a : allele)
      quality = std::max(a.second.NH_GQ, (double)quality); 
    
  }

  void BreakPoint::order() {

    if (b1.gr < b2.gr)
      return;
    
    std::swap(b1, b2);
    
  }

// format the text BX tag table (for 10X reads)x
void BreakPoint::format_bx_string() {

  // only make if alraedy empty
  if (!bxtable.empty())
    return;

  std::unordered_set<std::string> qn; // only add tag once per qname
    std::unordered_map<std::string, size_t> supp_tags;

    //add the discordant reads
    for (const auto& r : dc.reads) {
      std::string qname = r.second.Qname();
      if (qn.count(qname))
	continue;
      std::string tmp;
      r.second.GetZTag("BX", tmp);
      if (!tmp.empty()) 
	++supp_tags[tmp];
      qn.insert(qname);
    }
    
    //add the reads from the breakpoint
    for (const auto& r : reads) {
      std::string qname = r.Qname();
      if (qn.count(qname))
	continue;
      std::string tmp;
      r.GetZTag("BX", tmp);
      if (!tmp.empty())
	++supp_tags[tmp];
      qn.insert(qname);
    }

    bx_count = supp_tags.size();
    
    for (const auto& s : supp_tags)
      bxtable = bxtable += s.first + "_" + std::to_string(s.second) + ",";
    if (bxtable.length())
      bxtable.pop_back(); // remove last comma
}

  void BreakPoint::format_readname_string() {
    
    // only operate it not already formated
    if (!read_names.empty())
      return;
    
    std::unordered_set<std::string> supp_reads;
    
    //add the discordant reads
    for (auto& r : dc.reads) 
      supp_reads.insert(r.second.SR());
    //supp_reads.insert(SRTAG(r.second));

    //add the reads from the breakpoint
    for (auto& r : reads) 
      //supp_reads.insert(SRTAG(r)); // BamRecords (not svabaBamRead)
    supp_reads.insert(r.SR());
    
    // print reads to a string, delimit with a ,
    size_t lim = 0;
    for (auto& i : supp_reads) {
      if (++lim > 50)
	break;

      if (i.at(0) == 't') {
	size_t posr = i.find("_", 5) + 1;
	if (posr != std::string::npos)
	  read_names = read_names + "," + i.substr(posr, i.length() - posr);
      }
    }
    if (!read_names.empty())
      read_names = read_names.substr(1, read_names.size() - 1); // remove first _

  }

  void BreakPoint::setRefAlt(const RefGenome * main_rg, const RefGenome * viral) {

    assert(!main_rg->IsEmpty());
    assert(ref.empty());
    assert(alt.empty());
    
    if (evidence != "INDEL") {

      try {
	// get the reference for BP1
	if (b1.chr_name.find("gi|") == std::string::npos) {
	  ref = main_rg->QueryRegion(b1.chr_name, b1.gr.pos1-1, b1.gr.pos1-1);
	} else {
	  throw std::invalid_argument("dummy to get to viral");
	}
      } catch (const std::invalid_argument& ia) {}
      
      // try viral approach
      if (viral && ref.empty() && !viral->IsEmpty())
	try {
	  ref = viral->QueryRegion(b1.chr_name, b1.gr.pos1-1, b1.gr.pos1-1); 
	} catch (const std::invalid_argument& ia) {
	  ref = "N";
	  std::cerr << "Caught exception in BreakPoint:setRefAlt for SV Ref: " << ia.what() << std::endl;
	}
      
      try {
	if (b2.chr_name.find("gi|") == std::string::npos)
	  alt = main_rg->QueryRegion(b2.chr_name, b2.gr.pos1-1, b2.gr.pos1-1);
	else
	  throw std::invalid_argument("dummy to get to viral");
      } catch (const std::invalid_argument& ia) {}

      // try viral alt
      if (viral && alt.empty() && !viral->IsEmpty())
	try {
	  alt = viral->QueryRegion(b2.chr_name, b2.gr.pos1-1, b2.gr.pos1-1);
	} catch (const std::invalid_argument& ia) {
	  alt = "N";
	  std::cerr << "Caught exception in BreakPoint:setRefAlt for SV Alt: " << ia.what() << std::endl;
	}
      
    } else {

      if (insertion.length() && !insertion.empty()) {
	
	try {
	  ref = main_rg->QueryRegion(b1.chr_name, b1.gr.pos1-1, b1.gr.pos1-1);
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
	  ref = main_rg->QueryRegion(b1.chr_name, b1.gr.pos1-1, b2.gr.pos1-2);
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

ReducedBreakPoint::ReducedBreakPoint(const std::string &line, const SeqLib::BamHeader& h) {
    
  if (h.isEmpty()) {
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
    std::string ref_s, alt_s, cname_s, insertion_s, homology_s, evidence_s, confidence_s, read_names_s, bxtable_s;
    
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
	case 24: quality = std::stod(val); break; //std::min((int)255,std::stoi(val)); break;
	case 25: secondary = val == "1" ? 1 : 0;
	case 26: somatic_score = std::stod(val); break;
	case 27: somatic_lod = std::stod(val); break;
	case 28: true_lod = std::stod(val); break;
	case 29: pon = std::min(255,std::stoi(val)); break;
	case 30: repeat_s = val; break; // repeat_seq
	case 31: blacklist = (val=="1" ? 1 : 0); break;
	case 32: dbsnp = val != "x"; break;
	case 33: read_names_s = val; break; //reads
	case 34: bxtable_s = val; break; //bx tags
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
      read_names     = read_names_s; // == "x" ? nullptr : __string_alloc2char(read_names_s, read_names);
    bxtable = bxtable_s;

  }

  double SampleInfo::__log_likelihood(double ref, double alt, double f, double e) {
    
    // less negative log-likelihoods means more likely
    // eg for low error rate, odds that you see 5 ALT and 5 REF
    // if you are testing for AF = 0 is going to be very low (eg -40)
    // To test if something is true, we want to test the log-likelihood that
    // it's AF is != 0, so we test LL(ref, alt, AF=0, er). If this is 
    // a large negative number, it means that AF = 0 is very unlikely.
    // If we use a larger error rate, then it is more likely that we will
    // see ALT reads even if true AF = 0, so as ER goes up, then 
    // LL(ref, alt, AF=0, er) becomes less negative.

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
    else if (alt > 0)
      ll = -10000000 ; // error is zero, but have artifact so impossible

    return ll;

  }

  void SampleInfo::modelSelection(double er) {

    // can't have more alt reads than total reads
    // well you can for SVs...
    int thiscov = cov;
    if (alt >= cov)  
      thiscov = alt;

    // adjust the alt count 
    if (alt < cigar)
      alt = cigar;
    if (alt < split)
      alt = split;

    // adjust the coverage to be more in line with restrictions on ALT.
    // namely that ALT reads must overlap the variant site by more than T_SPLIT_BUFF
    // bases, but the raw cov calc does not take this into account. Therefore, adjust here
    double a_cov;
    if (readlen) {
      a_cov = (double)thiscov * (double)(readlen - 2 * T_SPLIT_BUFF)/readlen;
      a_cov = a_cov < 0 ? 0 : a_cov;
    } else {
      a_cov = thiscov;
    }
    af = a_cov > 0 ? (double)alt / (double)a_cov : 1;
    af = af > 1 ? 1 : af;

    int scaled_alt = std::min(alt, (int)a_cov);

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
    double ll_alt = __log_likelihood(a_cov - scaled_alt, scaled_alt, af, er);
    double ll_err = __log_likelihood(a_cov - scaled_alt, scaled_alt, 0 , er);
    LO = ll_alt - ll_err; 

    //mutetct log likelihood normal
    // er = 0.0005; // make this low, so that ALT in REF is rare and NORM in TUM gives low somatic prob
    // actually, dont' worry about it too much. 3+ alt in ref excludes somatic anyways.
    // so a high LO_n for the normal means that we are very confident that site is REF only in 
    // normal sample. This is why LO_n from the normal can be used as a discriminant for 
    // germline vs somatic. eg if somatic_lod = normal.LO_n is above a threshold X
    // or above a larger threshold XX if at DBSNP site, then we accept as somatic
    // LO_n should not be used for setting the confidence that something is real, just the 
    // confidence that it is somatic
    double ll_ref_norm = __log_likelihood(a_cov - scaled_alt, scaled_alt, 0 , er); // likelihood that varaint is actually true REF
    double ll_alt_norm = __log_likelihood(a_cov - scaled_alt, scaled_alt, std::max(af, 0.5), er); // likelihood that variant is 0.5
    //std::cerr << " COV " << a_cov << " ALT " << scaled_alt << " LL ALT " << ll_alt_norm << " LL REF " << ll_ref_norm << " ER " << er << " LL TOTAL " << (ll_ref_norm - ll_alt_norm) << std::endl;
    LO_n = ll_ref_norm - ll_alt_norm; // higher number means more likely to be AF = 0 (ref) than AF = 0.5 (alt). 

    // genotype calculation as provided in 
    // http://bioinformatics.oxfordjournals.org/content/early/2011/09/08/bioinformatics.btr509.full.pdf+html
    //int scaled_cov = std::floor((double)cov * 0.90);
    //int this_alt = std::min(alt, scaled_cov);
    genotype_likelihoods[0] = __genotype_likelihoods(2, er, scaled_alt, a_cov); // 0/0
    genotype_likelihoods[1] = __genotype_likelihoods(1, er, scaled_alt, a_cov); // 0/1
    genotype_likelihoods[2] = __genotype_likelihoods(0, er, scaled_alt, a_cov); // 1/1

    //debug
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

    // scale GT likelihoods to max
    genotype_likelihoods[0] = max_likelihood - genotype_likelihoods[0];
    genotype_likelihoods[1] = max_likelihood - genotype_likelihoods[1];
    genotype_likelihoods[2] = max_likelihood - genotype_likelihoods[2];

    // get the genotype quality
    GQ = 99;
    for (auto& g : genotype_likelihoods)
      if (g != 0)
	GQ = std::min(GQ, __myround(g));

    // get the genotype quality that it is not hom ref
    NH_GQ = std::min((double)99, __myround(genotype_likelihoods[0]));

    // set the PL string
    std::stringstream sss;
    sss << std::setprecision(4) << __myround(genotype_likelihoods[0]) << "," << __myround(genotype_likelihoods[1]) << "," << __myround(genotype_likelihoods[2]);
    PL = sss.str();

  }

  void BreakPoint::addCovs(const std::unordered_map<std::string, STCoverage*>& covs) {

    for (auto& i : covs)  {
      int c = 0;
      for (int j = -COVERAGE_AVG_BUFF; j <= COVERAGE_AVG_BUFF; ++j) {
	c +=  i.second->getCoverageAtPosition(b1.gr.chr, b1.gr.pos1 + j);
	c +=  i.second->getCoverageAtPosition(b2.gr.chr, b2.gr.pos1 + j);
      }
      allele[i.first].cov = c / 2 / (COVERAGE_AVG_BUFF*2 + 1); // std::max(i.second->getCoverageAtPosition(b1.gr.chr, b1.gr.pos1), i.second->getCoverageAtPosition(b2.gr.chr, b2.gr.pos1));

      if (cname=="c_22_16905001_16930001_191C") //debug
	std::cerr << i.first << " COV " << allele[i.first].cov << std::endl;
    }
    
  }

  std::ostream& operator<<(std::ostream& out, const BreakEnd& b) {
    out << b.gr << " - " << b.id << " mapq " << b.mapq << " subn " << b.sub_n << std::endl;
    return out;
  }
  
  std::ostream& operator<<(std::ostream& out, const SampleInfo& a) {
    out << " split: " << a.split << " cigar " << a.cigar << " alt " << a.alt << " cov " << a.cov << " disc " << a.disc;
    return out;
  }

  SampleInfo operator+(const SampleInfo& a1, const SampleInfo& a2) {

    SampleInfo a;

    a.disc = a1.disc + a2.disc;
    a.split = a1.split + a2.split;
    a.cigar = a1.cigar + a2.cigar;
    a.cov = a1.cov + a2.cov;

    // add the reads
    for (auto& i : a1.supporting_reads)
      a.supporting_reads.insert(i);
    for (auto& i : a2.supporting_reads)
      a.supporting_reads.insert(i);
    
    //a.alt = std::max((int)a.supporting_reads.size(), a.cigar);
    
    if (a.supporting_reads.size()) // we have the read names, so do that (non-refilter run)
      a.adjust_alt_counts();
    else // no read names stored, so just get directly
      a.alt = a1.alt + a2.alt;

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

    a.readlen = readlen;
    t.readlen = readlen;
    n.readlen = readlen;
    error_rate = (repeat_seq.length() > 10) ? MAX_ERROR : ERROR_RATES[repeat_seq.length()];
    t.modelSelection(error_rate);
    n.modelSelection(error_rate);
  }


  std::string SampleInfo::toFileString() const {

    std::stringstream ss;

    if (indel)
      ss << std::setprecision(4) << genotype << ":" << 
	std::max(alt, cigar) << ":" << cov << ":" << GQ << ":" << PL << ":" << split << ":" << cigar 
	 << ":" << LO_n << ":" << LO;// << ":" << SLO;
    else
      ss << std::setprecision(4) << genotype << ":" << 
	alt << ":" << cov << ":" << GQ << ":" << PL << ":" << split
	 << ":" << disc << ":" << LO_n << ":" << LO;// << ":" << SLO;
          
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
    
    if (gr.GetOverlap(window))
      local = true;
  }

  void SampleInfo::adjust_alt_counts() {
    
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
  // http://bioinformatics.oxfordjournals.org/content/early/2011/09/08/bioinformatics.btr509.full.pdf+html
  double SampleInfo::__genotype_likelihoods(int g, double er, int alt, int cov) {
    double val =  - cov * log10(2) + (cov - alt) * log10( (2 - g) * er + g  * (1 - er) ) + alt * log10( (2 - g) * (1 - er) + g * er);
    return val;
  }

bool ReducedBreakPoint::operator<(const ReducedBreakPoint& bp) const { 

  //ASDIS > ASSMB > COMP > DSCRD
  if (std::strcmp(evidence,bp.evidence) < 0) // <
    return true;
  else if (std::strcmp(evidence, bp.evidence) > 0) // >
    return false;
  
  if (nsplit > bp.nsplit) 
    return true;
  else if (nsplit < bp.nsplit)
    return false;
  
  if (tsplit > bp.tsplit)
    return true;
  else if (tsplit < bp.tsplit)
    return false;
  
  if (dc.ncount > bp.dc.ncount)
    return true;
  else if (dc.ncount < bp.dc.ncount)
    return false;
  
  if (dc.tcount > bp.dc.tcount)
    return true;
  else if (dc.tcount < bp.dc.tcount)
    return false;

  // break the tie somehow
  if (cname > bp.cname)
    return true;
  else if (cname < bp.cname)
    return false;
  
  return false;
}
