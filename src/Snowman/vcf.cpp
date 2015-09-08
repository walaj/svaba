#include "vcf.h"

#include <string> 
#include <iostream>
#include <getopt.h>
#include <sstream>
#include <regex>

#include "SnowTools/SnowUtils.h"
#include "SnowTools/SnowToolsCommon.h"
#include "SnowTools/gzstream.h"
#include "SnowTools/GenomicRegionCollection.h"

#include "htslib/tbx.h"
#include "htslib/bgzf.h"

using namespace std;

static faidx_t * findex;
static InfoMap flag_map;
static InfoMap union_info_fields;

static std::unordered_map<size_t, std::string> chr_number_to_string;

// forward declare
void tabixVcf(const string &fn);

// comparator for info fields
// lhs < rhs
// want READ_ID to be > than everything
bool compareInfoFields(const pair<string,string> &lhs, const pair<string,string> &rhs) {
  return ( (rhs.first == "READ_ID" && lhs.first != "READ_ID") || ( (rhs.first != "READ_ID" && lhs.first != "READ_ID") && lhs.first < rhs.first));
}

// print out the VCF header
std::ostream& operator<<(std::ostream& out, const VCFHeader& v) {

  out << "##fileformat=" << v.fileformat << endl;
  out << "##fileDate="   << v.filedate << endl;
  out << "##source="     << v.source << endl;
  out << "##reference="  << v.reference << endl;

  //output the contig information
  for (ContigFieldMap::const_iterator it = v.contigfieldmap.begin(); it != v.contigfieldmap.end(); it++)
    out << "##contig=<ID=" << it->first << "," << it->second << ">" << endl;
  
  //output the info information
  for (InfoMap::const_iterator it = v.infomap.begin(); it != v.infomap.end(); it++)
    out << "##INFO=<ID=" << it->first << "," << it->second << ">" << endl;
										
  //output the filter information
  for (FilterMap::const_iterator it = v.filtermap.begin(); it != v.filtermap.end(); it++) 
    out << "##FILTER=<ID=" << it->first << "," << it->second << ">" << endl;  

  //output the format information
  for (FormatMap::const_iterator it = v.formatmap.begin(); it != v.formatmap.end(); it++) 
    out << "##FORMAT=<ID=" << it->first << "," << it->second << ">" << endl;  

  //output the sample information
  for (SampleMap::const_iterator it = v.samplemap.begin(); it != v.samplemap.end(); it++) 
    out << "##SAMPLE=<ID=" << it->first << ">" << endl;  


  // output the colnames
  out << v.colnames;

  return out;
}

//add an info field
void VCFHeader::addInfoField(string field, string number, string type, string description) {

  if (infomap.find(field) != infomap.end()) {
    cerr << "Warning: Info field already exists: " << field << endl;
    return;
  }
    
  if (type == "Flag")
    flag_map.insert(pair<string,string>(field, type));

  string net = "Number=" + number + ",Type=" + type + ",Description=\"" + description + "\"";
  infomap[field] = net;
  return;

}

//add a filter field
void VCFHeader::addFilterField(string field, string description) {

  if (filtermap.find(field) != filtermap.end()) {
    cerr << "Warning: Filter field already exists" << endl;
    return;
  }
    
  string net = "Description=\"" + description + "\"";
  filtermap[field] = net;
  return;

}

//add a format field
void VCFHeader::addFormatField(string field, string number, string type, string description) {

  if (formatmap.find(field) != formatmap.end()) {
    cerr << "Warning: Format field already exists" << endl;
    return;
  }

  string net = "Number=" + number + ",Type=" + type + ",Description=\"" + description + "\"";    
  formatmap[field] = net;
  return;

}

//add a sample field
void VCFHeader::addSampleField(string field) {

  if (samplemap.find(field) != samplemap.end()) {
    cerr << "Warning: Sample field already exists" << endl;
    return;
  }
    
  samplemap[field] = field;
  return;

}

// print out the VCF Entry
std::ostream& operator<<(std::ostream& out, const VCFEntry& v) {

  // move to a vector to be sorted
  bool imprecise = false;
  vector<pair<string, string> > tmpvec; // id, evertythign else
  for (InfoMap::const_iterator it = v.info_fields.begin(); it != v.info_fields.end(); it++) {
    //if (it->first == "IMPRECISE")
    //  imprecise = true;
    //if (it->second != "")
    tmpvec.push_back(pair<string,string>(it->first, it->second)); // id, d
  }
  sort(tmpvec.begin(), tmpvec.end(), compareInfoFields); // sort it

  string info;
  string equals = "=";
  for (vector<pair<string, string> >::const_iterator it = tmpvec.begin(); it != tmpvec.end(); it++)
    if (!(it->first == "HOMSEQ" && imprecise) && !(it->first=="HOMLEN" && imprecise) && !(it->first=="INSERTION" && imprecise))// dont print some fields if imprecise
      info = info + it->first + ( (flag_map.count(it->first) == 0) ? "=" : "") + it->second + ";"; // dont print = for flags

  // trim the last semicolon from info
  if (info.length() > 0)
    info = info.substr(0, info.length() - 1);

  string sep = "\t";
  string qual_string = (v.qual >= 0) ? std::to_string(v.qual) : ".";
  //out << SnowTools::GenomicRegion::chrToString(v.chr) << sep 
  out << v.be.chr_name << sep //chr_number_to_string[v.chr] << sep 
      << v.be.gr.pos1 << sep << v.id << sep << v.ref << sep << v.alt << sep << v.qual << sep
      << v.filter << sep << info << sep << v.format << sep << v.samp1 << sep << v.samp2;

  return out;
}

// sort the VCFEntry by genomic position
bool VCFEntry::operator<(const VCFEntry &v) const {
  return be.gr < v.be.gr;
}

// create a VCFFile from a snowman breakpoints file
VCFFile::VCFFile(string file, const char* index, string id, bam_hdr_t * h) {

  analysis_id = id;

  //open the file
  igzstream infile(file.c_str(), ios::in);
  
  // confirm that it is open
  if (!infile) {
    cerr << "Can't read file " << file << " for parsing VCF" << endl;
    exit(EXIT_FAILURE);
  }

  // read the reference if not open
  std::cerr << "...vcf - reading in the reference" << endl;
  if (!findex)
    findex = fai_load(index);  // load the reference

  // read in the header of the csv
  string line;

  string sample_id_tum = analysis_id + "T";
  string sample_id_norm= analysis_id + "N";

  // make the VCFHeader
  sv_header.filedate = "";
  sv_header.source = "snowmanSV";
  indel_header.filedate = "";
  indel_header.source = "snowmanSV";

  // add the filters that apply to SVs
  sv_header.addFilterField("NODISC","Rearrangement was not detected independently by assembly");
  sv_header.addFilterField("LOWMAPQ","Assembly contig has non 60/60 mapq and no discordant support");
  sv_header.addFilterField("WEAKASSEMBLY","4 or fewer split reads and no discordant support and span > 1500bp");
  sv_header.addFilterField("WEAKDISC","Fewer than 7 supporting discordant reads and no assembly support");
  sv_header.addFilterField("TOOSHORT","Contig alignment for part of this rearrangement has <= 25bp match to reference");
  sv_header.addFilterField("PASS", "Strong assembly support, strong discordant support, or combined support. Strong MAPQ"); //3+ split reads, 0 normal split reads, 60/60 contig MAPQ OR 3+ discordant reads or 60/60 MAPQ with 4+ split reads");
  sv_header.addSampleField(sample_id_norm);
  sv_header.addSampleField(sample_id_tum);
  sv_header.colnames = sv_header.colnames + "\t" + sample_id_norm + "\t" + sample_id_tum;

  // add the filters that apply to indels
  indel_header.addFilterField("LOWMAPQ","Assembly contig has less than MAPQ 60");
  indel_header.addFilterField("WEAKASSEMBLY","4 or fewer split reads");
  indel_header.addFilterField("WEAKCIGARMATCH","For indels <= 5 bp, require 8+ split reads or 4+ and 3+ cigar matches");
  indel_header.addFilterField("REPEAT", "3+ split reads, 60 contig MAPQ");
  indel_header.addFilterField("PASS", "3+ split reads, 60 contig MAPQ");
  indel_header.addSampleField(sample_id_norm);
  indel_header.addSampleField(sample_id_tum);
  indel_header.colnames = indel_header.colnames + "\t" + sample_id_norm + "\t" + sample_id_tum;

  // set the time string
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  stringstream month;
  stringstream mdate;
  if ( (now->tm_mon+1) < 10)
    month << "0" << now->tm_mon+1;
  else 
    month << now->tm_mon+1;
  mdate << (now->tm_year + 1900) << month.str() <<  now->tm_mday;

  indel_header.filedate = mdate.str();
  sv_header.filedate = mdate.str();

  //add the SV info fields
  sv_header.addInfoField("REPSEQ","1","String","Repeat sequence near the event");
  sv_header.addInfoField("SVTYPE","1","String","Type of structural variant");
  sv_header.addInfoField("HOMSEQ","1","String","Sequence of base pair identical micro-homology at event breakpoints. Plus strand sequence displayed.");
  sv_header.addInfoField("IMPRECISE","0","Flag", "Imprecise structural variation");
  sv_header.addInfoField("SECONDARY","0","Flag", "SV calls comes from a secondary alignment");
  sv_header.addInfoField("HOMLEN","1","Integer","Length of base pair identical micro-homology at event breakpoints");
  sv_header.addInfoField("BKDIST","1","Integer","Distance between breakpoints (-1 if difference chromosomes");
  sv_header.addInfoField("MAPQ","1","Integer","Mapping quality (BWA-MEM) of this fragement of the contig (-1 if discordant only)");
  sv_header.addInfoField("MATEMAPQ","1","Integer","Mapping quality of the partner fragment of the contig");
  sv_header.addInfoField("NSPLIT","1","Integer","Number of split reads from the normal BAM");
  sv_header.addInfoField("TSPLIT","1","Integer","Number of split reads from the tumor BAM");
  sv_header.addInfoField("TDISC","1","Integer","Number of discordant read pairs from the tumor BAM");
  sv_header.addInfoField("NDISC","1","Integer","Number of discordant read pairs from the normal BAM");
  sv_header.addInfoField("MATEID","1","String","ID of mate breakends");

  sv_header.addFormatField("READ_ID",".","String","ALT supporting Read IDs");
  sv_header.addFormatField("NALT_SR","1","Integer","Number of ALT support Split Reads");           
  sv_header.addFormatField("NALT_RP","1","Integer","Number of ALT support aberrant Read Pairs");
  sv_header.addFormatField("NREF","1","Integer","Number of REF support Reads");
  sv_header.addFormatField("NALT","1","Integer","Number of ALT support reads or pairs");

  sv_header.addInfoField("NUMPARTS","1","Integer","If detected with assembly, number of parts the contig maps to. Otherwise 0");
  sv_header.addInfoField("EVDNC","1","String","Provides type of evidence for read. ASSMB is assembly only, ASDIS is assembly+discordant. DSCRD is discordant only.");
  sv_header.addInfoField("SCTG","1","String","Identifier for the contig assembled by SnowmanSV to make the SV call");
  sv_header.addInfoField("INSERTION","1","String","Sequence insertion at the breakpoint.");
  sv_header.addInfoField("SPAN","1","Integer","Distance between the breakpoints. -1 for interchromosomal");

  // add the indel header fields
  indel_header.addInfoField("MAPQ","1","Integer","Mapping quality (BWA-MEM) of the assembled contig");
  indel_header.addInfoField("SPAN","1","Integer","Size of the indel");
  indel_header.addInfoField("NCIGAR","1","Integer","Number of normal reads with cigar strings supporting this indel");
  indel_header.addInfoField("TCIGAR","1","Integer","Number of tumor reads with cigar strings supporting this indel");
  indel_header.addInfoField("NSPLIT","1","Integer","Number of normal reads with read-to-contig alignments supporting this indel");
  indel_header.addInfoField("TSPLIT","1","Integer","Number of tumor reads with read-to-contig alignments supporting this indel");
  indel_header.addFormatField("READ_ID",".","String","ALT supporting Read IDs");
  indel_header.addFormatField("NALT_SR","1","Integer","Number of ALT support Split Reads");           
  indel_header.addFormatField("NALT_RP","1","Integer","Number of ALT support aberrant Read Pairs");
  indel_header.addFormatField("NREF","1","Integer","Number of REF support Reads");
  indel_header.addFormatField("NALT","1","Integer","Number of ALT support reads or pairs");
  indel_header.addInfoField("REPSEQ","1","String","Repeat sequence near the event");
  indel_header.addInfoField("TCOV","1","Integer","Tumor coverage at break");
  indel_header.addInfoField("NCOV","1","Integer","Normal coverage at break");
  indel_header.addInfoField("TFRAC","1","String","Tumor allelic fraction at break. -1 for undefined");
  indel_header.addInfoField("NFRAC","1","String","Normal allelic fraction at break. -1 for undefined");
  indel_header.addInfoField("GRAYLIST","1","String","Normal allelic fraction at break. -1 for undefined");

  // keep track of exact positions to keep from duplicating
  unordered_map<string,bool> double_check;

  // read the reference if not open
  cerr << "...vcf - reading in the breakpoints file" << endl;

  // read it in line by line
  getline(infile, line, '\n'); // skip first line
  size_t line_count = 0;
  while (getline(infile, line, '\n')) {

    // parse the breakpoint from the file
    SnowTools::BreakPoint bp(line, h);
    //assert(bp.valid());

    // add the VCFentry Pair
    std::string idcommon = std::to_string(++line_count) + ":" + analysis_id;
    VCFEntryPair vpair(bp, idcommon);
    
    if (vpair.e1.info_fields["EVDNC"] == "INDEL")
      indels.insert(pair<string, VCFEntry>(idcommon, vpair.e1));
    else
      entry_pairs.insert(pair<string, VCFEntryPair>(idcommon, vpair));
      
  }

  // dedupe
  std::cerr << "...vcf - deduplicating " << entry_pairs.size() << " events" << std::endl;
  deduplicate();
  std::cerr << "...vcf - deduplicated down to " << (entry_pairs.size() - dups.size()) << " break pairs" << std::endl;
  
}

// return a sequence from the reference
string getRefSequence(const std::string& chr_name, const SnowTools::GenomicRegion& gr, faidx_t * fi) {

  int len;
  char * seq = faidx_fetch_seq(fi, const_cast<char*>(chr_name.c_str()), gr.pos1-1, gr.pos2-1, &len);
  
  if (seq) {
    return string(seq);
  } else {
    return "N";
  }

}

// make a class to hold break end + id
class GenomicRegionWithID : public SnowTools::GenomicRegion 
{
  public: 
  GenomicRegionWithID(int32_t c, uint32_t p1, uint32_t p2, const std::string& i, const std::string& q) : SnowTools::GenomicRegion(c,p1,p2) { id = i; qual = q; }
  std::string id, qual;
};

// deduplicate
void VCFFile::deduplicate() {

  std::cerr << "...vcf - deduping events" << endl;

  // create the interval tree maps
  // grv1 are left entries, grv2 are right
  // keep it sorted so grv1 always has left most
  SnowTools::GenomicRegionCollection<GenomicRegionWithID> grv1;
  SnowTools::GenomicRegionCollection<GenomicRegionWithID> grv2;
  for (auto& i : entry_pairs) {
    grv1.add(GenomicRegionWithID(i.second.e1.be.gr.chr, i.second.e1.be.gr.pos1, i.second.e1.be.gr.pos2, i.first, i.second.e1.filter)); 
    grv2.add(GenomicRegionWithID(i.second.e2.be.gr.chr, i.second.e2.be.gr.pos1, i.second.e2.be.gr.pos2, i.first, i.second.e2.filter)); 
  }
  grv1.createTreeMap();
  grv2.createTreeMap();
  assert(grv1.size() == grv2.size());

  int pad = 1;

  size_t count = 0;
  
  for (auto& i : entry_pairs) {

    // if it's already de-duped, dont do it again
    if (dups.count(i.first))
      continue;

    ++count;
    SnowTools::GenomicIntervalVector giv1, giv2;
    grv1.m_tree[i.second.e1.be.gr.chr].findContained(i.second.e1.be.gr.pos1-pad, i.second.e1.be.gr.pos1+pad, giv1);
    grv2.m_tree[i.second.e2.be.gr.chr].findContained(i.second.e2.be.gr.pos1-pad, i.second.e2.be.gr.pos1+pad, giv2);
    
    // loop through hits and only add if not the current site.
    // If key_count is 2, then it hit on each side. This is a dup
    // If this is a dup, remove all the things it dups to
    // then when you come across one that was marked as dup, just skip
    // the intersection step

    bool is_pass = i.second.e1.filter == "PASS";

    std::unordered_map<std::string, size_t> key_count;
    // loop hits to left end
    for (auto& j : giv1)
      if (grv1.at(j.value).id != i.first && ( is_pass == (grv1.at(j.value).qual == "PASS") )) //j is hit. Make sure have same pass status
	++key_count[grv1.at(j.value).id];
    // loop hits to right end
    for (auto& j : giv2)
      if (grv2.at(j.value).id != i.first && ( is_pass == (grv2.at(j.value).qual == "PASS") )) //j is hit, grv2.at(j.value).id is key of hit
	++key_count[grv2.at(j.value).id];

    //loop through hit keys and if key is hit twice (1 left, 1 right), it is an overlap
    for (auto& j : key_count) {
      if (j.second == 2) { // left and right hit for this key. add 
	dups.insert(j.first); 
      }
    }
  }
  
}

// print a breakpoint pair
ostream& operator<<(ostream& out, const VCFEntryPair& v) {

  out << v.e1 << endl;
  out << v.e2 << endl;
  return (out);

}

bool VCFEntry::operator==(const VCFEntry &v) const {

  return (v.be.gr == be.gr) ; //chr == v.chr && pos == v.pos);

}

// write out somatic and germline INDEL vcfs
void VCFFile::writeIndels(string basename, bool zip) const {

  string gname = basename + "germline.indel.vcf.gz";
  string sname = basename + "somatic.indel.vcf.gz";
  string gname_nz = basename + "germline.indel.vcf";
  string sname_nz = basename + "somatic.indel.vcf";

  ofstream out_g, out_s;

  BGZF* g_bg = NULL;
  BGZF* s_bg = NULL;

  if (zip) {
    g_bg = bgzf_open(gname.c_str(), "w");
    s_bg = bgzf_open(sname.c_str(), "w");
    stringstream indel_h;
    indel_h << indel_header << endl;
    if (!bgzf_write(g_bg, indel_h.str().c_str(), indel_h.str().length())) {
      cerr << "Could not write bzipped vcf" << endl;
    }
    if (!bgzf_write(s_bg, indel_h.str().c_str(), indel_h.str().length())) {
      cerr << "Could not write bzipped vcf" << endl;
    }
  } else {
    out_g.open(gname_nz.c_str());
    out_s.open(sname_nz.c_str());
    out_g << indel_header << endl;
    out_s << indel_header << endl;
  }

  VCFEntryVec tmpvec;

  // put the indels into a sorted vector
  for (auto& i : indels) {
    tmpvec.push_back(i.second);
  }

  // sort the temp entry vec
  sort(tmpvec.begin(), tmpvec.end());  

  // print out the entries
  for (auto& i : tmpvec) { 
    if (i.filter == "PASS" || include_nonpass) {

      if (i.info_fields.count("SOMATIC")) {
	if (zip) {
	  stringstream ss;
	  ss << i << endl;
	  if (!bgzf_write(s_bg, ss.str().c_str(), ss.str().length())) {
	    cerr << "Could not write bzipped vcf for line " << ss.str() << endl;
	  }
	} else {
	  out_s << i << endl;
	}
	
      } else {
	if (zip) {
	  stringstream ss;
	  ss << i << endl;
	  if (!bgzf_write(g_bg, ss.str().c_str(), ss.str().length())) {
	    cerr << "Could not write bzipped vcf for line " << ss.str() << endl;
	  }
	} else {
	  out_g << i << endl;
	}


      }
    }
  }
 
  if (zip) {
    bgzf_close(g_bg);
    bgzf_close(s_bg);
  } else {
    out_g.close();
    out_s.close();
  }

  if (zip) {
    // tabix it
    tabixVcf(gname);
    tabixVcf(sname);
  }

}

// write out somatic and germline SV vcfs
void VCFFile::writeSVs(string basename, bool zip) const {
  
  string gname, sname, gname_nz, sname_nz; 
  gname = basename + "germline.sv.vcf.gz";
  sname = basename + "somatic.sv.vcf.gz";
  gname_nz = basename + "germline.sv.vcf";
  sname_nz = basename + "somatic.sv.vcf";

  ofstream out_g, out_s;

  BGZF* s_bg = NULL;
  BGZF* g_bg = NULL;

  if (zip) {
    g_bg = bgzf_open(gname.c_str(), "w");
    s_bg = bgzf_open(sname.c_str(), "w");
    stringstream sv_h;
    sv_h << sv_header << endl;
    if (!bgzf_write(g_bg, sv_h.str().c_str(), sv_h.str().length())) {
      cerr << "Could not write bzipped vcf" << endl;
    }
    if (!bgzf_write(s_bg, sv_h.str().c_str(), sv_h.str().length())) {
      cerr << "Could not write bzipped vcf" << endl;
    }
  } else {
    out_g.open(gname_nz.c_str());
    out_s.open(sname_nz.c_str());

    out_g << sv_header << endl;
    out_s << sv_header << endl;
  }
    
  VCFEntryVec tmpvec;

  // put the pair maps into a vector
  for (VCFEntryPairMap::const_iterator it = entry_pairs.begin(); it != entry_pairs.end(); it++) {
    
    if (!dups.count(it->first)) { // dont include duplicate entries

      // renumber the ids
      //++id_counter;
      //VCFEntryPair tmppair = it->second;

      //tmppair.e1.idcommon = to_string(id_counter) + ":" + analysis_id;
      //tmppair.e2.idcommon = to_string(id_counter) + ":" + analysis_id;
      //tmppair.e1.id = tmppair.e1.idcommon + ":1";
      //tmppair.e2.id = tmppair.e2.idcommon + ":2";
      //tmppair.e1.info_fields["MATEID"] = tmppair.e2.id;
      //tmppair.e2.info_fields["MATEID"] = tmppair.e1.id;
    
      tmpvec.push_back(it->second.e1);
      tmpvec.push_back(it->second.e2);
      //tmpvec.push_back(tmppair.e1);
      //tmpvec.push_back(tmppair.e2);
    }

  }

  // sort the temp entry vec
  sort(tmpvec.begin(), tmpvec.end());  

  // print out the entries
  for (auto& i : tmpvec) { 
    
    if (i.filter == "PASS" || include_nonpass) {

      // somatic
      if ( i.info_fields.count("SOMATIC") ) { 

	if (zip) {
	  stringstream ss;
	  ss << i << endl;
	  if (!bgzf_write(s_bg, ss.str().c_str(), ss.str().length())) {
	    cerr << "Could not write bzipped vcf for line " << ss.str() << endl;
	  }
	} else {
	  out_s << i << endl;
	}
	
      // germline
      } else {
	if (zip) {
	  stringstream ss;
	  ss << i << endl;
	  if (!bgzf_write(g_bg, ss.str().c_str(), ss.str().length())) {
	    cerr << "Could not write bzipped vcf for line " << ss.str() << endl;
	  }
	} else {
	  out_g << i << endl;
	}
      }
    }
  }

  if (zip) {
    bgzf_close(g_bg);
    bgzf_close(s_bg);
  } else {
    out_s.close();
    out_g.close();
  }

  // tabix it
  if (zip) 
    tabixVcf(sname);

}


// tabix the vcf
void tabixVcf(const string &fn) {

  // tabix it
  tbx_conf_t conf = tbx_conf_gff;
  tbx_conf_t * conf_ptr = &tbx_conf_vcf;
  conf = *conf_ptr;
  if ( tbx_index_build(fn.c_str(), 0, &conf) ) 
    cerr << "tbx_index_build failed: " << fn << endl;

}

VCFEntryPair::VCFEntryPair(const SnowTools::BreakPoint& b, const std::string& id) {

  e1.be = b.b1;
  e2.be = b.b2;

  e1.qual = b.quality;
  e2.qual = b.quality;
  
  idcommon = id;
  e1.id = id;

  std::unordered_map<std::string, std::string> common_fields;

  // put all the common fields in 
  common_fields["SPAN"] = std::to_string(b.getSpan());
  common_fields["NSPLIT"] = std::to_string(b.nsplit);
  common_fields["TSPLIT"] = std::to_string(b.tsplit);
  common_fields["SCTG"] = b.cname;
  common_fields["EVDNC"] = b.evidence;
  if (b.somatic_score >= 4)
    common_fields["SOMATIC"] = "";
  //common_fields["PONCOUNT"]  = std::to_string(pon); 
  //if (b.repeat_seq)
  //  common_fields["REPSEQ"]    = b.repeat_seq; 

  assert( (b.num_align == 1) == (b.evidence == "INDEL") );

  // put all the common fields for SVs
  if (b.num_align != 1) {

    common_fields["NDISC"] = std::to_string(b.dc.ncount);
    common_fields["TDISC"] = std::to_string(b.dc.tcount);
    if (b.homology != "x") common_fields["HOMSEQ"] = b.homology;
    if (b.insertion != "x") common_fields["INSERTION"] = b.insertion;
    common_fields["NUMPARTS"] = std::to_string(b.num_align);
    if (b.evidence=="DSCRD") 
      common_fields["IMPRECISE"] = ""; 
    if (b.secondary) 
      common_fields["SECONDARY"] = "";

    //
    e1.info_fields["SUBN"] = b.b1.sub_n;
    e2.info_fields["SUBN"] = b.b2.sub_n;
  
  }
  // put all the common fields in for indels
  else {
    e1.info_fields["DISC_MAPQ"] = b.dc.mapq1;
    e2.info_fields["DISC_MAPQ"] = b.dc.mapq2;
  
    common_fields["NCIGAR"] = std::to_string(b.ncigar);
    common_fields["TCIGAR"] = std::to_string(b.tcigar);
    common_fields["NCOV"]      = std::to_string(b.ncov); 
    common_fields["TCOV"]      = std::to_string(b.tcov);
    common_fields["NFRAC"]     = std::to_string(b.af_n); 
    common_fields["TFRAC"]     = std::to_string(b.af_t); 
    if (b.blacklist)
      common_fields["GRAYLIST"]  = "";
    if (!b.rs.empty() && b.rs != "x")
      common_fields["DBSNP"] = b.rs; 
  }

  e1.filter = b.confidence;
  e2.filter = b.confidence;

  // add the common fields
  e1.info_fields = common_fields;
  e2.info_fields = common_fields;
  
  // unique fields
  e1.info_fields["MAPQ"] = std::to_string(b.b1.mapq);
  e2.info_fields["MAPQ"] = std::to_string(b.b2.mapq);

  // add mate info for SVs
  if (b.num_align != 1) {

    e1.id = idcommon + ":1";
    e2.id = idcommon + ":2";
    e1.info_fields["MATEID"] = e2.id;
    e2.info_fields["MATEID"] = e1.id;
    
    e1.info_fields["SVTYPE"] = "BND";
    e2.info_fields["SVTYPE"] = "BND";

    e1.format = "NALT:NALT_RP:NALT_SR:READ_ID";
    e2.format = "NALT:NALT_RP:NALT_SR:READ_ID";
    
    // put the reads into the format string
    std::string new_readid_t = formatReadString(b.read_names, 't');
    std::string new_readid_n = formatReadString(b.read_names, 'n');
   
    int numt = std::stoi(common_fields["TDISC"]) + std::stoi(common_fields["TSPLIT"]);
    int numn = std::stoi(common_fields["NDISC"]) + std::stoi(common_fields["NSPLIT"]);
    e1.samp2 = to_string(numt) + ":" + common_fields["TDISC"] + ":" + common_fields["TSPLIT"] + ":" + new_readid_t;
    e2.samp2 = e1.samp2;
    e1.samp1 = to_string(numn) + ":" + common_fields["NDISC"] + ":" + common_fields["NDISC"] + ":" + new_readid_n;
    e2.samp1 = e1.samp1;
    samples = {e1.samp1, e1.samp2};

    // grab the reference sequence
    e1.ref = b.ref;
    e2.ref = b.alt;
    //if (e1.be.gr.chr < 24 && e2.be.gr.chr < 24) { // TODO FIX 
    //  e1.ref = getRefSequence(e1.be.chr_name, e1.be.gr, findex);
    //  e2.ref = getRefSequence(e2.be.chr_name, e1.be.gr, findex);
    //}
      
    // set the reference position for making ALT tag
    stringstream ptag1, ptag2;
    ptag1 << e2.be.chr_name << ":" << e2.be.gr.pos1;
    ptag2 << e1.be.chr_name << ":" << e1.be.gr.pos1;
    
    string insertion = common_fields["INSERTION"];
    string ttag1, ttag2;
    // TODO insertion + vcf1.ref not right
    //ttag1 = (insertion != "" && e1.be.gr.strand == '+' && false) ? (insertion + e1.ref) : e1.ref;
    //ttag1 = (insertion != "" && e1.be.gr.strand == '-' && false) ? (e1.ref + insertion) : e1.ref;
    //ttag2 = (insertion != "" && e2.be.gr.strand == '+' && false) ? (insertion + e2.ref) : e2.ref;
    //ttag2 = (insertion != "" && e2.be.gr.strand == '-' && false) ? (e2.ref + insertion) : e2.ref;
    ttag1 = e1.ref;
    ttag2 = e2.ref;

    // set the alternate
    stringstream alt1, alt2;
    if (e1.be.gr.strand == '+' && e2.be.gr.strand == '+') {
      alt1 << ttag1 << "]" << ptag1.str() << "]";
      alt2 << ttag2 << "]" << ptag2.str() << "]";
    } else if (e1.be.gr.strand =='+' && e2.be.gr.strand == '-') {
      alt1 << ttag1 << "[" << ptag1.str() << "[";
      alt2 << "]" << ptag2.str() << "]" << ttag2;
    } else if (e1.be.gr.strand == '-' && e2.be.gr.strand == '+') {
      alt1 << "]" << ptag1.str() << "]" << ttag1;
      alt2 << ttag2 << "[" << ptag2.str() << "[";
    } else {
      alt1 << "[" << ptag1.str() << "[" << ttag1;      
      alt2 << "[" << ptag2.str() << "[" << ttag2;      
    }
    
    e1.alt = alt1.str();
    e2.alt = alt2.str();

  }
  
  // indel
  if (common_fields["EVDNC"] == "INDEL") {

    e1.idcommon = id;

    // get the alt sequence
    //SnowTools::GenomicRegion ref(e1.be.gr.chr, e1.be.gr.pos1, e2.be.gr.pos1-1);
    //if (e1.be.gr.chr < 24)
    //  e1.ref = getRefSequence(e1.be.chr_name, ref, findex);
    e1.ref = b.ref;
    e1.alt = b.alt;
    //e1.alt = e1.ref.substr(0,1);    

    // add the sequences if insertion
    //if (common_fields["INSERTION"].length()) 
    //  e1.alt += b.insertion;
  }

}

std::string formatReadString(const std::string& readid, char type) {

  if (readid == "x" || readid.empty())
    return std::string();
  
  // parse the readids
  SupportingReadsMap suppr;
  istringstream iss_r(readid);
  string thisread;
  string new_readid = "";
  set<string> dup;

  // regex to clean out the t, n identifer
  regex regc("[a-z][0-9]+_(.*?)$");
  smatch smatchr;

  while (getline(iss_r, thisread, ',')) {
    
    string thisread_clean; // get only the read name
    if (!regex_search(thisread, smatchr, regc))
      cerr << "FAILED TO MATCH ON "<< thisread << endl;
    else 
      thisread_clean = smatchr[1].str();

    suppr.insert(pair<string,bool>(thisread_clean, false));
    
    if (thisread.at(0) == type && !dup.count(thisread_clean)) {
      new_readid += thisread_clean + ",";
      dup.insert(thisread_clean);
    }
  }

  // remove the last comma
  if (new_readid.length() > 0)
    new_readid = new_readid.substr(0, new_readid.length() - 1);
  
  // ok, so we cant have any colons in the read name, so swtich : for -
  std::replace( new_readid.begin(), new_readid.end(), ':', '-');  

  return new_readid;

}

