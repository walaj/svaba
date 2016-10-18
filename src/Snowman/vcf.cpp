#include "vcf.h"

#include <getopt.h>
#include <regex>
#include <string> 
#include <iomanip>
#include <sstream>
#include <iostream>
#include <unordered_set> 

#include "htslib/tbx.h"
#include "htslib/bgzf.h"

#include "gzstream.h"
#include "SeqLib/GenomicRegionCollection.h"

#define VCF_SECONDARY_CAP 200
#define SOMATIC_LOD 1


using namespace std;

static std::string sv_format = "GT:AD:DP:GQ:PL:SR:DR:LR:LO:SL"; //"NALT:NALT_RP:NALT_SR:READ_ID";
static std::string indel_format = "GT:AD:DP:GQ:PL:SR:CR:LR:LO:SL";
static InfoMap flag_map;
static int global_id = 0;
static std::stringstream lod_ss;

static std::unordered_map<std::string, int> cname_count;

void __write_to_zip_vcf(const VCFEntry& v, BGZF * f) {
  std::stringstream ss;
  ss << v << endl;
  if (!bgzf_write(f, ss.str().c_str(), ss.str().length())) 
    cerr << "Could not write bzipped vcf for line " << ss.str() << endl;
}

// forward declare
void tabixVcf(const std::string &fn);

// comparator for info fields
// lhs < rhs
// want READ_ID to be > than everything
bool compareInfoFields(const pair<std::string,std::string> &lhs, const pair<std::string,std::string> &rhs) {
  return ( (rhs.first == "READ_ID" && lhs.first != "READ_ID") || ( (rhs.first != "READ_ID" && lhs.first != "READ_ID") && lhs.first < rhs.first));
}


bool pairCompare(const std::pair<int, std::pair<std::string, std::string>>& firstElem, const std::pair<int, std::pair<std::string, std::string>>& secondElem) {
  return firstElem.first > secondElem.first;
}

// print out the VCF header
std::ostream& operator<<(std::ostream& out, const VCFHeader& v) {

  out << "##fileformat=" << v.fileformat << endl;
  out << "##fileDate="   << v.filedate << endl;
  out << "##source="     << v.source << endl;
  out << "##reference="  << v.reference << endl;
									
  // order the contigs by length
  typedef std::pair<std::string, std::string> CC; // contig struct
  typedef std::pair<int, CC> CI; // contig with integer length
  std::vector<CI> contig_vec;
  for (auto& i : v.contigfieldmap)
    contig_vec.push_back(CI(std::stoi(i.second), CC(i.first, i.second)));
  std::sort(contig_vec.begin(), contig_vec.end(), pairCompare);
  
  //for (ContigFieldMap::const_iterator it = v.contigfieldmap.begin(); it != v.contigfieldmap.end(); it++)
  //  out << "##contig=<ID=" << it->first << ",length=" << it->second << ">" << endl;
  for (auto& i : contig_vec)
    out << "##contig=<ID=" << i.second.first << ",length=" << i.second.second << ">" << std::endl;
  for (InfoMap::const_iterator it = v.infomap.begin(); it != v.infomap.end(); it++)
    out << "##INFO=<ID=" << it->first << "," << it->second << ">" << endl;
  for (FilterMap::const_iterator it = v.filtermap.begin(); it != v.filtermap.end(); it++) 
    out << "##FILTER=<ID=" << it->first << "," << it->second << ">" << endl;  
  for (FormatMap::const_iterator it = v.formatmap.begin(); it != v.formatmap.end(); it++) 
    out << "##FORMAT=<ID=" << it->first << "," << it->second << ">" << endl;  
  for (SampleMap::const_iterator it = v.samplemap.begin(); it != v.samplemap.end(); it++) 
    out << "##SAMPLE=<ID=" << it->first << ">" << endl;  

  // output the colnames
  out << v.colnames;

  return out;
}

//add an info field
void VCFHeader::addInfoField(std::string field, std::string number, std::string type, std::string description) {

  if (infomap.find(field) != infomap.end()) {
    cerr << "Warning: Info field already exists: " << field << endl;
    return;
  }
    
  if (type == "Flag")
    flag_map.insert(pair<std::string,std::string>(field, type));

  std::string net = "Number=" + number + ",Type=" + type + ",Description=\"" + description + "\"";
  infomap[field] = net;
  return;

}

//add a filter field
void VCFHeader::addFilterField(std::string field, std::string description) {

  if (filtermap.find(field) != filtermap.end()) {
    cerr << "Warning: Filter field already exists" << endl;
    return;
  }
    
  std::string net = "Description=\"" + description + "\"";
  filtermap[field] = net;
  return;

}

//add a format field
void VCFHeader::addFormatField(std::string field, std::string number, std::string type, std::string description) {

  if (formatmap.find(field) != formatmap.end()) {
    cerr << "Warning: Format field already exists" << endl;
    return;
  }

  std::string net = "Number=" + number + ",Type=" + type + ",Description=\"" + description + "\"";    
  formatmap[field] = net;
  return;

}

//add a sample field
void VCFHeader::addSampleField(std::string field) {

  if (samplemap.find(field) != samplemap.end()) {
    cerr << "Warning: Sample field already exists" << endl;
    return;
  }
    
  samplemap[field] = field;
  return;

}

// print out the VCF Entry
std::ostream& operator<<(std::ostream& out, const VCFEntry& v) {

  std::unordered_map<std::string, std::string> info_fields = v.fillInfoFields();

  // move to a vector to be sorted
  vector<pair<string, std::string> > tmpvec; // id, evertythign else
  for (InfoMap::const_iterator it = info_fields.begin(); it != info_fields.end(); it++) 
    tmpvec.push_back(pair<std::string,std::string>(it->first, it->second)); 
  sort(tmpvec.begin(), tmpvec.end(), compareInfoFields); // sort it

  std::string info;
  std::string equals = "=";
  for (vector<pair<std::string, std::string> >::const_iterator it = tmpvec.begin(); it != tmpvec.end(); it++) {
    if (!(it->first == "HOMSEQ" && v.bp->imprecise) && !(it->first=="HOMLEN" && v.bp->imprecise) && !(it->first=="INSERTION" && v.bp->imprecise))// dont print some fields if imprecise
      info = info + it->first + ( (flag_map.count(it->first) == 0) ? "=" : "") + it->second + ";"; // dont print = for flags
  }

  // trim the last semicolon from info
  if (info.length() > 0)
    info = info.substr(0, info.length() - 1);

  std::string sep = "\t";
  ReducedBreakEnd * be = v.id_num == 1 ? &v.bp->b1 : &v.bp->b2;

  //std::pair<std::string, std::string> samps = v.getSampStrings();
  out << be->chr_name << sep  
      << be->gr.pos1 << sep << v.getIdString() << sep << v.getRefString() << sep << v.getAltString() << sep 
      << v.bp->quality << sep
      << v.bp->confidence << sep << info << sep 
      << (v.bp->indel ? indel_format : sv_format); // << sep << samps.first << sep << samps.second;
  for (auto& i : v.bp->format_s)
    out << sep << i;
  return out;
}

// sort the VCFEntry by genomic position
bool VCFEntry::operator<(const VCFEntry &v) const {
  ReducedBreakEnd * be = id_num == 1 ? &bp->b1 : &bp->b2;
  ReducedBreakEnd * vbe = v.id_num == 1 ? &v.bp->b1 : &v.bp->b2;
  return be->gr < vbe->gr;    
}

// create a VCFFile from a snowman breakpoints file
VCFFile::VCFFile(std::string file, std::string id, const SeqLib::BamHeader& h, const VCFHeader& vheader, bool nopass) {

  analysis_id = id;

  //open the file
  igzstream infile(file.c_str(), ios::in);
  
  // confirm that it is open
  if (!infile) {
    cerr << "Can't read file " << file << " for parsing VCF" << endl;
    exit(EXIT_FAILURE);
  }

  // read in the header of the csv
  std::string line;

  //string sample_id_tum = analysis_id + "T";
  //string sample_id_norm= analysis_id + "N";

  sv_header    = vheader;
  indel_header = vheader;

  // add the filters that apply to SVs
  sv_header.addFilterField("NOLOCAL","Contig realigned to region outside of local assembly region, and no disc support.");
  sv_header.addFilterField("LOCALMATCH","Contig realigned to assembly region without clipping");
  sv_header.addFilterField("HIGHHOMOLOGY","Contig realigns with > 25% of readlength of homology. High probaility of assembly/mapping artifact");
  sv_header.addFilterField("DUPREADS","Contig built from what appear to be duplicate reads (split reads all same contig cov))");
  sv_header.addFilterField("NODISC","Rearrangement was not detected independently by assembly");
  sv_header.addFilterField("LOWSUPPORT","Fewer than 2 split reads or < 4 total alt reads for ASDISC");
  sv_header.addFilterField("COMPETEDISC","Discordant cluster found with nearly same breakpoints, but different strands for DSCRD event");
  sv_header.addFilterField("LOWSPAN","Discordant read cluster (no split read support), and less than 10kb span and < 12 reads");
  //sv_header.addFilterField("FOLDBACK","Rearrangement is inversion type with span < 80. Very likely fold-back Illumina error");
  sv_header.addFilterField("LOWMAPQ","Assembly contig has non 60/60 mapq and no discordant support");
  sv_header.addFilterField("LOWQINVERSION","Assembly-only inversion of span < 300 and < 6 split reads. Common artifact in Illumina data");
  sv_header.addFilterField("LOWMAPQDISC","Both clusters of reads failed to achieve mean mapq of > 30 for DSCRD");
  sv_header.addFilterField("LOWSPLITSMALL","Fewer than 4 split reads for small events ( < 1500 bp)");
  sv_header.addFilterField("LOWICSUPPORT","Less than 60bp of contig match on one end of an inter-chromosomal break");
  sv_header.addFilterField("LOWAS","Alignment score of one end is less than 80% of contig length, or number of mismatch bases (NM) on one end is >= 10");
  sv_header.addFilterField("WEAKSUPPORTHIREP","Fewer then 7 split reads for variant with >= 10 bases of repeat sequence (need to be more strict)");
  sv_header.addFilterField("WEAKDISC","Fewer than 7 supporting discordant reads and no assembly support");
  sv_header.addFilterField("TOOSHORT","Contig alignment for part of this rearrangement has <= 25bp match to reference");
  sv_header.addFilterField("PASS", "Strong assembly support, strong discordant support, or combined support. Strong MAPQ"); 
  sv_header.addFilterField("MULTIMATCH", "Low MAPQ and this contig fragment maps well to multiple locations");
  sv_header.addFilterField("LOWSPANDSCRD", "Discordant-only cluster is too small given isize distribution to call confidently"); 
  sv_header.addFilterField("SIMPLESEQUENCE", "Major portion of one contig mapping falls in a simple sequence, as given by -R flag. Assembly-only filter"); 
  //sv_header.addSampleField(sample_id_norm);
  //sv_header.addSampleField(sample_id_tum);
  //sv_header.colnames = sv_header.colnames + "\t" + sample_id_norm + "\t" + sample_id_tum;

  // add the filters that apply to indels
  indel_header.addFilterField("LOWMAPQ","Assembly contig has less than MAPQ 10");
  indel_header.addFilterField("SHORTALIGNMENT","Matched (M) contig frag to left or right of indel < 20 bp");
  indel_header.addFilterField("LOWLOD","LOD score is less than the cutoff");
  //indel_header.addFilterField("WEAKASSEMBLY","4 or fewer split reads");
  //indel_header.addFilterField("WEAKCIGARMATCH","For indels <= 5 bp, require 8+ split reads or 4+ and 3+ cigar matches");
  indel_header.addFilterField("PASS", "LOD score pass");
  //indel_header.addFilterField("GRAYLISTANDPON", "Indel overlaps with panel of normals, and has overlap with tricky genomic region");
  indel_header.addFilterField("VLOWAF", "allelic fraction < 0.05");
  //indel_header.addFilterField("LOWNORMCOV", "Fewer than 5 normal reads at this site");
  //indel_header.addFilterField("GERMLOWAF", "Germline variant with support for being AF of 0.5+ (LR) less than cutoff");

  //indel_header.addSampleField(sample_id_norm);
  //indel_header.addSampleField(sample_id_tum);
  //indel_header.colnames = indel_header.colnames + "\t" + sample_id_norm + "\t" + sample_id_tum;

  //add the SV info fields
  sv_header.addInfoField("REPSEQ","1","String","Repeat sequence near the event");
  sv_header.addInfoField("READNAMES",".","String","IDs of ALT reads");
  sv_header.addInfoField("BX",".","String","Table of BX tag counts for supporting reads");
  sv_header.addInfoField("NM","1","Integer","Number of mismatches of this alignment fragment to reference");
  sv_header.addInfoField("MATENM","1","Integer","Number of mismatches of partner alignment fragment to reference");
  sv_header.addInfoField("SVTYPE","1","String","Type of structural variant");
  sv_header.addInfoField("HOMSEQ","1","String","Sequence of base pair identical micro-homology at event breakpoints. Plus strand sequence displayed.");
  sv_header.addInfoField("IMPRECISE","0","Flag", "Imprecise structural variation");
  sv_header.addInfoField("SECONDARY","0","Flag", "SV calls comes from a secondary alignment");
  sv_header.addInfoField("HOMLEN","1","Integer","Length of base pair identical micro-homology at event breakpoints");
  sv_header.addInfoField("DBSNP","1","String","TRUE if variant overlaps a dbSNP site");
  //sv_header.addInfoField("BKDIST","1","Integer","Distance between breakpoints (-1 if difference chromosomes");
  sv_header.addInfoField("MAPQ","1","Integer","Mapping quality (BWA-MEM) of this fragement of the contig (-1 if discordant only)");
  sv_header.addInfoField("MATEMAPQ","1","Integer","Mapping quality of the partner fragment of the contig");
  //sv_header.addInfoField("NSPLIT","1","Integer","Number of split reads from the normal BAM");
  //sv_header.addInfoField("TSPLIT","1","Integer","Number of split reads from the tumor BAM");
  //sv_header.addInfoField("TDISC","1","Integer","Number of discordant read pairs from the tumor BAM");
  //sv_header.addInfoField("NDISC","1","Integer","Number of discordant read pairs from the normal BAM");
  sv_header.addInfoField("MATEID","1","String","ID of mate breakends");
  //sv_header.addInfoField("SOMATIC","0","Flag","Variant is somatic");
  sv_header.addInfoField("SUBN","1","Integer","Number of secondary alignments associated with this contig fragment");
  //sv_header.addInfoField("TCOV","1","Integer","Max tumor coverage at break");
  //sv_header.addInfoField("NCOV","1","Integer","Max normal coverage at break");
  //sv_header.addInfoField("TFRAC","1","String","Tumor allelic fraction at break. -1 for undefined");
  //sv_header.addInfoField("NFRAC","1","String","Normal allelic fraction at break. -1 for undefined");

  sv_header.addInfoField("NUMPARTS","1","Integer","If detected with assembly, number of parts the contig maps to. Otherwise 0");
  sv_header.addInfoField("EVDNC","1","String","Provides type of evidence for read. ASSMB is assembly only, ASDIS is assembly+discordant. DSCRD is discordant only.");
  sv_header.addInfoField("SCTG","1","String","Identifier for the contig assembled by SnowmanSV to make the SV call");
  sv_header.addInfoField("INSERTION","1","String","Sequence insertion at the breakpoint.");
  sv_header.addInfoField("SPAN","1","Integer","Distance between the breakpoints. -1 for interchromosomal");
  sv_header.addInfoField("DISC_MAPQ","1","Integer","Mean mapping quality of discordant reads mapped here");

  // add the indel header fields
  indel_header.addInfoField("SCTG","1","String","Identifier for the contig assembled by SnowmanSV to make the indel call");
  indel_header.addInfoField("MAPQ","1","Integer","Mapping quality (BWA-MEM) of the assembled contig");
  indel_header.addInfoField("SPAN","1","Integer","Size of the indel");

  indel_header.addFormatField("GT", "1","String", "Genotype (currently not supported. Always 0/1)");
  indel_header.addFormatField("AD", "1","Integer", "Allele depth: Number of reads supporting the variant");
  indel_header.addFormatField("DP","1","Integer","Depth of coverage: Number of reads covering site.");
  indel_header.addFormatField("GQ", "1","String", "Genotype quality (currently not supported. Always 0)");
  indel_header.addFormatField("PL","1","Integer","Normalized likelihood of the current genotype (currently not supported, always 0)");
  indel_header.addFormatField("SR","1","Integer","Number of spanning reads for this variants");
  indel_header.addFormatField("CR","1","Integer","Number of cigar-supported reads for this variant");
  indel_header.addFormatField("LR","1","Float","Log-odds that this variant is AF=0 vs AF>=0.5");
  indel_header.addFormatField("LO","1","Float","Log-odds that this variant is real vs artifact");
  indel_header.addFormatField("SL","1","Float","Alignment-quality Scaled log-odds, where LO is LO * (MAPQ - 2*NM)/60");

  sv_header.addFormatField("GT", "1","String", "Genotype (currently not supported. Always 0/1)");
  sv_header.addFormatField("AD", "1","Integer", "Allele depth: Number of reads supporting the variant");
  sv_header.addFormatField("DP","1","Integer","Depth of coverage: Number of reads covering site.");
  sv_header.addFormatField("GQ", "1","String", "Genotype quality (currently not supported. Always 0)");
  sv_header.addFormatField("PL","1","Integer","Normalized likelihood of the current genotype (currently not supported, always 0)");
  sv_header.addFormatField("SR","1","Integer","Number of spanning reads for this variants");
  sv_header.addFormatField("DR","1","Integer","Number of discordant-supported reads for this variant");
  sv_header.addFormatField("LR","1","Float","Log-odds that this variant is REF vs AF=0.5");
  sv_header.addFormatField("LO","1","Float","Log-odds that this variant is real vs artifact");
  sv_header.addFormatField("SL","1","Float","Alignment-quality Scaled log-odds, where LO is LO * (MAPQ - 2*NM)/60");

  indel_header.addInfoField("REPSEQ","1","String","Repeat sequence near the event");
  indel_header.addInfoField("GRAYLIST","0","Flag","Indel is low quality and cross a difficult region of genome");
  //indel_header.addInfoField("SOMATIC","0","Flag","Variant is somatic");
  indel_header.addInfoField("PON","1","Integer","Number of normal samples that have this indel present");
  indel_header.addInfoField("NM","1","Integer","Number of mismatches of this alignment fragment to reference");
  indel_header.addInfoField("READNAMES",".","String","IDs of ALT reads");
  indel_header.addInfoField("BX",".","String","Table of BX tag counts for supporting reads");
  indel_header.addInfoField("LOD","1","Float","Log of the odds that variant is real vs artifact");

  // keep track of exact positions to keep from duplicating
  // read the reference if not open
  cerr << "...vcf - reading in the breakpoints file" << endl;
  
  include_nonpass = nopass;

  // read it in line by line
  getline(infile, line, '\n'); // skip first line
  size_t line_count = 0;
  while (getline(infile, line, '\n')) {

    if (line.find("mapq") != std::string::npos)
      continue;

    if (line.find("Unknown") != std::string::npos)
      continue;

    // parse the breakpoint from the file
    std::shared_ptr<ReducedBreakPoint> bp(new ReducedBreakPoint(line, h));

    // add the VCFentry Pair
    ++line_count;
    std::shared_ptr<VCFEntryPair> vpair(new VCFEntryPair(bp));

    // skip non pass if not emitting unfiltered
    if (!include_nonpass && !bp->pass)
      continue;

    ++cname_count[std::string(bp->cname)];
    if (cname_count[std::string(bp->cname)] >= VCF_SECONDARY_CAP)
      {
	//delete bp;
	//delete vpair;
	continue;
      }

    if (bp->indel) {
      indels.insert(pair<int, std::shared_ptr<VCFEntryPair>>(line_count, vpair));
    }
    else  {
      entry_pairs.insert(pair<int, std::shared_ptr<VCFEntryPair>>(line_count, vpair));
    }
   
  }
  
  cname_count.clear();
  std::cerr << "...vcf sizeof empty VCFEntryPair " << sizeof(VCFEntryPair) << " bytes " << std::endl;
  std::cerr << "...read in " << SeqLib::AddCommas(indels.size()) << " indels and " << SeqLib::AddCommas(entry_pairs.size()) << " SVs " << std::endl;
  
  std::cerr << "...vcf - deduplicating " << SeqLib::AddCommas(entry_pairs.size()) << " events" << std::endl;
  deduplicate();
  std::cerr << "...vcf - deduplicated down to " << SeqLib::AddCommas((entry_pairs.size() - dups.size())) << " break pairs" << std::endl;
  
}

// make a class to hold break end + id
class GenomicRegionWithID : public SeqLib::GenomicRegion 
{
  public: 
  GenomicRegionWithID(int32_t c, uint32_t p1, uint32_t p2, int i, int p) : SeqLib::GenomicRegion(c,p1,p2), id(i), pass(p) {}
  uint32_t id:30, pass:2;
};

// deduplicate
void VCFFile::deduplicate() {

  std::cerr << "...vcf - deduping events" << endl;

  // create the interval tree maps
  // grv1 are left entries, grv2 are right
  // keep it sorted so grv1 always has left most
  SeqLib::GenomicRegionCollection<GenomicRegionWithID> grv1;
  SeqLib::GenomicRegionCollection<GenomicRegionWithID> grv2;
  for (auto& i : entry_pairs) {
    grv1.add(GenomicRegionWithID(i.second->bp->b1.gr.chr, i.second->bp->b1.gr.pos1, i.second->bp->b1.gr.pos2, i.first, i.second->e1.bp->pass)); 
    grv2.add(GenomicRegionWithID(i.second->bp->b2.gr.chr, i.second->bp->b2.gr.pos1, i.second->bp->b2.gr.pos2, i.first, i.second->e2.bp->pass)); 
  }
  grv1.CreateTreeMap();
  grv2.CreateTreeMap();
  assert(grv1.size() == grv2.size());

  int pad = 1;
  size_t count = 0;
  
  for (auto& i : entry_pairs) {

    // if it's already de-duped, dont do it again
    if (dups.count(i.first))
      continue;
    
    // if both ends are close (within 10) then they match
    pad = (i.second->bp->b1.gr.chr != i.second->bp->b2.gr.chr) || std::abs(i.second->bp->b1.gr.pos1 - i.second->bp->b2.gr.pos1) > 100 ? 10 : 1;

    ++count;
    SeqLib::GenomicIntervalVector giv1, giv2;
    SeqLib::GenomicIntervalTreeMap::const_iterator ff1 = grv1.GetTree()->find(i.second->bp->b1.gr.chr);
    SeqLib::GenomicIntervalTreeMap::const_iterator ff2 = grv2.GetTree()->find(i.second->bp->b2.gr.chr);
    assert(ff1 != grv1.GetTree()->end());
    assert(ff2 != grv2.GetTree()->end());
    ff1->second.findContained(i.second->bp->b1.gr.pos1-pad, i.second->bp->b1.gr.pos1+pad, giv1);
    ff2->second.findContained(i.second->bp->b2.gr.pos1-pad, i.second->bp->b2.gr.pos1+pad, giv2);
      //grv1.m_tree[i.second->bp->b1.gr.chr].findContained(i.second->bp->b1.gr.pos1-pad, i.second->bp->b1.gr.pos1+pad, giv1);
      //grv2.m_tree[i.second->bp->b2.gr.chr].findContained(i.second->bp->b2.gr.pos1-pad, i.second->bp->b2.gr.pos1+pad, giv2);
    
    // loop through hits and only add if not the current site.
    // If key_count is 2, then it hit on each side. This is a dup
    // If this is a dup, remove all the things it dups to
    // then when you come across one that was marked as dup, just skip
    // the intersection step

    bool is_pass = i.second->e1.bp->pass; 

    std::unordered_map<int, size_t> key_count;
    // loop hits to left end
    for (auto& j : giv1)
      if (grv1.at(j.value).id != i.first && ( is_pass == (grv1.at(j.value).pass) )) //j is hit. Make sure have same pass status
	++key_count[grv1.at(j.value).id];
    // loop hits to right end
    for (auto& j : giv2)
      if (grv2.at(j.value).id != i.first && ( is_pass == (grv2.at(j.value).pass) )) //j is hit, grv2.at(j.value).id is key of hit
	++key_count[grv2.at(j.value).id];

    //loop through hit keys and if key is hit twice (1 left, 1 right), it is an overlap
    for (auto& j : key_count) {
      if (j.second == 2) { // left and right hit for this key. add 
	dups.insert(j.first); 
      }
    }
  }

  // dedupe the indels
  std::cerr << "...hashing " << SeqLib::AddCommas(indels.size()) << " indels for dedupe" << std::endl;
  std::unordered_set<std::string> hashr;
  VCFEntryPairMap tmp_indels;

  for (auto& i : indels) {
    
    std::string hh;
    try {
      hh = std::to_string(i.second->e1.bp->b1.gr.chr) + ":" + std::to_string(i.second->e1.bp->b1.gr.pos1) + 
	     "_" + i.second->e1.getRefString() + "_" + i.second->e1.getAltString();
      } catch (...) {
      	std::cerr << " error " << std::endl;
   }
    if (!hashr.count(hh)) {
      hashr.insert(hh);
      tmp_indels.insert(pair<int, std::shared_ptr<VCFEntryPair>>(i.first, i.second));
    }
  }
  std::cerr << "...done deduping indels" << std::endl;
  indels = tmp_indels;
}

// print a breakpoint pair
ostream& operator<<(ostream& out, const VCFEntryPair& v) {

  out << v.e1 << endl;
  out << v.e2 << endl;
  return (out);

}

bool VCFEntry::operator==(const VCFEntry &v) const {
  ReducedBreakEnd * be = id_num == 1 ? &bp->b1 : &bp->b2;
  ReducedBreakEnd * vbe = v.id_num == 1 ? &v.bp->b1 : &v.bp->b2;
  return (vbe->gr == be->gr) ; //chr == v.chr && pos == v.pos);
}

// write out somatic and germline INDEL vcfs
void VCFFile::writeIndels(string basename, bool zip, bool onefile) const {

  std::string gname = basename + "germline.indel.vcf.gz";
  std::string sname = basename + "somatic.indel.vcf.gz";
  std::string gname_nz = basename + "germline.indel.vcf";
  std::string sname_nz = basename + "somatic.indel.vcf";

  if (onefile) {
    gname_nz = basename + "indel.vcf";
    gname    = basename + "indel.vcf.gz";
  }
    
  ofstream out_g, out_s;

  BGZF* g_bg = NULL;
  BGZF* s_bg = NULL;

  if (zip) {
    g_bg = bgzf_open(gname.c_str(), "w");
    if (!onefile) s_bg = bgzf_open(sname.c_str(), "w");
    std::stringstream indel_h;
    indel_h << indel_header << endl;
    if (!bgzf_write(g_bg, indel_h.str().c_str(), indel_h.str().length())) {
      cerr << "Could not write bzipped vcf" << endl;
    }
    if (!onefile)
      if (!bgzf_write(s_bg, indel_h.str().c_str(), indel_h.str().length())) {
	cerr << "Could not write bzipped vcf" << endl;
    }
  } else {
    out_g.open(gname_nz.c_str());
    if (!onefile)
      out_s.open(sname_nz.c_str());
    out_g << indel_header << endl;
    if (!onefile)
      out_s << indel_header << endl;
  }

  VCFEntryVec tmpvec;

  // put the indels into a sorted vector
  for (auto& i : indels) {
    tmpvec.push_back(i.second->e1);
  }

  // sort the temp entry vec
  sort(tmpvec.begin(), tmpvec.end());  

  // print out the entries
  for (auto& i : tmpvec) { 

    if (!i.bp->pass && !include_nonpass)
      continue;

    std::stringstream ss;
    if (!onefile && i.bp->somatic_score >= SOMATIC_LOD) {
      if (zip) 
	__write_to_zip_vcf(i, s_bg);
      else 
	out_s << i << endl;
      
    } else {
      if (zip) 
	__write_to_zip_vcf(i, g_bg);
      else
	out_g << i << endl;
    }
    
  }

  if (zip) {
    bgzf_close(g_bg);
    if (!onefile)
      bgzf_close(s_bg);
  } else {
    out_g.close();
    if (!onefile)
      out_s.close();
  }
  
  if (zip) {
    // tabix it
    tabixVcf(gname);
    if (!onefile)
      tabixVcf(sname);
  }
  
}

// write out somatic and germline SV vcfs
void VCFFile::writeSVs(std::string basename, bool zip, bool onefile) const {

  std::string gname, sname, gname_nz, sname_nz; 
  gname = basename + "germline.sv.vcf.gz";
  sname = basename + "somatic.sv.vcf.gz";
  gname_nz = basename + "germline.sv.vcf";
  sname_nz = basename + "somatic.sv.vcf";

  if (onefile) {
    gname    = basename + "sv.vcf.gz";
    gname_nz = basename + "sv.vcf";

  }

  ofstream out_g, out_s;

  BGZF* s_bg = NULL;
  BGZF* g_bg = NULL;

  if (zip) {
    g_bg = bgzf_open(gname.c_str(), "w");
    if (!onefile)
      s_bg = bgzf_open(sname.c_str(), "w");
    std::stringstream sv_h;
    sv_h << sv_header << endl;
    if (!bgzf_write(g_bg, sv_h.str().c_str(), sv_h.str().length())) {
      cerr << "Could not write bzipped vcf" << endl;
    }
    if (!onefile)
      if (!bgzf_write(s_bg, sv_h.str().c_str(), sv_h.str().length())) {
	cerr << "Could not write bzipped vcf" << endl;
      }
  } else {
    out_g.open(gname_nz.c_str());
    if (!onefile)
      out_s.open(sname_nz.c_str());

    out_g << sv_header << endl;
    if (!onefile)
      out_s << sv_header << endl;
  }
    
  VCFEntryVec tmpvec;

  // put the pair maps into a vector
  for (VCFEntryPairMap::const_iterator it = entry_pairs.begin(); it != entry_pairs.end(); it++) {
    
    if (!dups.count( it->first)) { // dont include duplicate entries

      // renumber the ids
      //++id_counter;
      //VCFEntryPair tmppair = it->second;

      //tmppair.e1.idcommon = std::to_string(id_counter) + ":" + analysis_id;
      //tmppair.e2.idcommon = std::to_string(id_counter) + ":" + analysis_id;
      //tmppair.e1.id = tmppair.e1.idcommon + ":1";
      //tmppair.e2.id = tmppair.e2.idcommon + ":2";
      //tmppair.e1.info_fields["MATEID"] = tmppair.e2.id;
      //tmppair.e2.info_fields["MATEID"] = tmppair.e1.id;
    
      tmpvec.push_back(it->second->e1);
      tmpvec.push_back(it->second->e2);
      //tmpvec.push_back(tmppair.e1);
      //tmpvec.push_back(tmppair.e2);
    }

  }

  // sort the temp entry vec
  sort(tmpvec.begin(), tmpvec.end());  

  // print out the entries
  for (auto& i : tmpvec) { 
    
    if (!i.bp->pass && !include_nonpass)
      continue;
    
    // somatic
    if (!onefile &&  i.bp->somatic_score >= SOMATIC_LOD) { 
      if (zip) 
	__write_to_zip_vcf(i, s_bg);
      else
	out_s << i << endl;
      // germline
    } else {
      if (zip)
	__write_to_zip_vcf(i, g_bg);
      else 
	out_g << i << endl;
    }

  }
  
  if (zip) {
    bgzf_close(g_bg);
    if (!onefile)
      bgzf_close(s_bg);
  } else {
    out_s.close();
    if (!onefile)
      out_g.close();
  }

  // tabix it
  if (zip) {
    if (!onefile)
      tabixVcf(sname); 
    tabixVcf(gname);
  }

}


// tabix the vcf
void tabixVcf(const std::string &fn) {

  // tabix it
  tbx_conf_t conf = tbx_conf_gff;
  tbx_conf_t * conf_ptr = &tbx_conf_vcf;
  conf = *conf_ptr;
  if ( tbx_index_build(fn.c_str(), 0, &conf) ) 
    cerr << "tbx_index_build failed: " << fn << endl;

}

VCFEntryPair::VCFEntryPair(std::shared_ptr<ReducedBreakPoint>& b) {

  bp = b;
  e1.bp = bp;
  e2.bp = bp;
  ++global_id;
  e1.id = global_id;
  e2.id = global_id;
  e1.id_num = 1;
  e2.id_num = 2;

}

std::string formatReadString(const std::string& readid, char type) {

  if (readid == "x" || readid.empty())
    return std::string();
  
  // parse the readids
  SupportingReadsMap suppr;
  std::istringstream iss_r(readid);
  std::string thisread;
  std::string new_readid = "";
  set<std::string> dup;

  // regex to clean out the t, n identifer
  regex regc("[a-z][0-9]+_(.*?)$");
  smatch smatchr;

  while (getline(iss_r, thisread, ',')) {
    
    std::string thisread_clean; // get only the read name
    if (!regex_search(thisread, smatchr, regc))
      cerr << "FAILED TO MATCH ON "<< thisread << endl;
    else 
      thisread_clean = smatchr[1].str();

    suppr.insert(pair<std::string,bool>(thisread_clean, false));
    
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

std::unordered_map<std::string, std::string> VCFEntry::fillInfoFields() const {
  
  std::unordered_map<std::string, std::string> info_fields;

  // put all the common fields in 
  info_fields["SPAN"] = std::to_string(bp->getSpan());
  info_fields["SCTG"] = bp->cname;
  if (!bp->indel) {
    info_fields["EVDNC"] = bp->evidence;
    info_fields["SVTYPE"] = "BND";
  }

  if (!bp->read_names.empty() && bp->read_names != "x")
    info_fields["READNAMES"] = bp->read_names;

  if (!bp->bxtable.empty() && bp->bxtable != "x")
    info_fields["BX"] = bp->bxtable;

  if (bp->repeat)
    info_fields["REPSEQ"] = std::string(bp->repeat);

  if (bp->pon)
    info_fields["PON"] = std::to_string(bp->pon);

  if (bp->num_align != 1) {
    info_fields["MATEID"] = std::to_string(id) + ":" + std::to_string(id_num == 1 ? 2 : 1);
    if (id_num == 1) {
      info_fields["NM"] = std::to_string(bp->b1.nm);
      info_fields["MATENM"] = std::to_string(bp->b2.nm);
    } else {
      info_fields["NM"] = std::to_string(bp->b2.nm);
      info_fields["MATENM"] = std::to_string(bp->b1.nm);
    }
  }

  if (bp->num_align == 1)
    info_fields["NM"] = std::to_string(bp->b1.nm);

  if (id_num == 1)
    info_fields["MAPQ"] = std::to_string(bp->b1.mapq);
  else
    info_fields["MAPQ"] = std::to_string(bp->b2.mapq);

  //if (bp->somatic_score >= SOMATIC_LOD)
  //  info_fields["SOMATIC"] = "";
  
  // put all the info fields for SVs
  if (bp->num_align != 1) {

    if (id_num == 1) {
      if (bp->b1.sub_n)
	info_fields["SUBN"] = std::to_string(bp->b1.sub_n);
      else if (bp->b2.sub_n)
	info_fields["SUBN"] = std::to_string(bp->b2.sub_n);
    }

    if (bp->homology) info_fields["HOMSEQ"] = std::string(bp->homology);
    if (bp->insertion) info_fields["INSERTION"] = std::string(bp->insertion);
    info_fields["NUMPARTS"] = std::to_string(bp->num_align);

    if (bp->imprecise)
      info_fields["IMPRECISE"] = ""; 
    if (bp->secondary) 
      info_fields["SECONDARY"] = "";

    if (info_fields["EVDNC"] != "ASSMB") {
      if (id_num == 1)
	info_fields["DISC_MAPQ"] = std::to_string(bp->dc.mapq1);
      else
	info_fields["DISC_MAPQ"] = std::to_string(bp->dc.mapq2);
    }

  }

  else {

    lod_ss << std::setprecision(4) << bp->true_lod;
    info_fields["LOD"] = lod_ss.str();
    lod_ss.str(std::string());
    if (bp->blacklist)
      info_fields["GRAYLIST"]  = "";
    if (bp->dbsnp)
      info_fields["DBSNP"] = "TRUE"; //bp->rs; 
  }

  return info_fields;

}

std::string VCFEntry::getRefString() const {

  char* p;
  if (bp->indel || id_num == 1) 
  	p = bp->ref;
  else
    p = bp->alt;
    
   if (!p) {
   	  std::cerr << "WARNING: Empty ref/alt field for bp " << std::endl;
   	  return (std::string("N"));
    }
   
   return (std::string(p));
}

std::string VCFEntry::getAltString() const {


  if (bp->indel) {
  
    char* p;
    p = bp->alt;

 	if (!p) {
  	 	 std::cerr << "WARNING: Empty ref/alt field for bp " << std::endl;
   	 	 return (std::string("N"));
 	 }
 	 
 	 return (std::string(p));
  
  }
  
  std::string ref = getRefString();

  std::stringstream ptag;
  if (id_num == 1) {
    ptag << bp->b2.chr_name << ":" << bp->b2.gr.pos1;
  } else {
    ptag << bp->b1.chr_name << ":" << bp->b1.gr.pos1;
  }
  
  std::stringstream alt;
  if (bp->b1.gr.strand == '+' && bp->b2.gr.strand == '+') {
    alt << ref << "]" << ptag.str() << "]";
  } else if (bp->b1.gr.strand =='+' && bp->b2.gr.strand == '-') {
    if (id_num == 1)
      alt << ref << "[" << ptag.str() << "[";
    else
      alt << "]" << ptag.str() << "]" << ref;
  } else if (bp->b1.gr.strand == '-' && bp->b2.gr.strand == '+') {
    if (id_num == 1)
      alt << "]" << ptag.str() << "]" << ref;
    else
      alt << ref << "[" << ptag.str() << "[";
  } else {
    alt << "[" << ptag.str() << "[" << ref;      
  }
  
  return alt.str();
}

std::string VCFEntry::getIdString() const {

  if (!bp->indel)
    return(std::to_string(id) + ":" + std::to_string(id_num));

  return(std::to_string(id));

}

std::pair<std::string, std::string> VCFEntry::getSampStrings() const {

  // put the reads into the format string
  //std::string new_readid_t = formatReadString(b.read_names, 't');
  //std::string new_readid_n = formatReadString(b.read_names, 'n');
  std::string new_readid_t = "";
  std::string new_readid_n = "";
  
  int numt, numn;
  numt = bp->dc.tcount + bp->tsplit; 
  numn = bp->dc.ncount + bp->nsplit; 
  
  std::string samp1, samp2;
  if (!bp->indel) {
    samp2 = std::to_string(numt) + ":" + std::to_string(bp->dc.tcount) + ":" + std::to_string(bp->tsplit) + ":" + new_readid_t;
    samp1 = std::to_string(numn) + ":" + std::to_string(bp->dc.ncount) + ":" + std::to_string(bp->nsplit) + ":" + new_readid_n;
  } else {
    samp2 = std::to_string(bp->tsplit) + ":" + new_readid_t;
    samp1 = std::to_string(bp->nsplit) + ":" + new_readid_n;
  }

  return std::pair<std::string, std::string>(samp1, samp2);


}

void VCFHeader::addContigField(std::string id, int len) {
  contigfieldmap[id] = std::to_string(len);
}
