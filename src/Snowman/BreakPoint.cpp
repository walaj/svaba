#include "BreakPoint.h"
#include <getopt.h>
#include "gzstream.h"
#include "vcf.h"

#define SPLIT_BUFF 12

namespace bopt {

  static string input_file = "";
  static string output_file = "";
  static string pon = "";
  static string analysis_id = "noid";
  static bool noreads = false;

  static string ref_index = REFHG19;

  static int verbose = 1;
}

static const char* shortopts = "hxi:a:v:q:g:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "input-bps",               required_argument, NULL, 'i'},
  { "panel-of-normals",        required_argument, NULL, 'q'},
  { "reference-genome",        required_argument, NULL, 'g'},
  { "analysis-id",             required_argument, NULL, 'a'},
  { "no-reads",                no_argument, NULL, 'x'},
  { "verbose",                 required_argument, NULL, 'v' },
  { NULL, 0, NULL, 0 }
};

static const char *BP_USAGE_MESSAGE =
"Usage: snowman refilter [OPTION] -i bps.txt.gz -o bps.new.txt.gz\n\n"
"  Description: \n"
"\n"
"  General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h, --help                           Display this help and exit\n"
"  -g, --reference-genome               Path to indexed reference genome to be used by BWA-MEM. Default is Broad hg19 (/seq/reference/...)\n"
"  Required input\n"
"  -i, --input-bps                      Original bps.txt.gz file\n"
"  Optional input\n"                       
"  -q, --panel-of-normals               Panel of normals files generated from snowman pon\n"                       
"  -a, --id-string                      String specifying the analysis ID to be used as part of ID common.\n"
"  -x, --no-reads                       Flag to turn off recording of supporting reads. Setting this flag greatly reduces file size.\n"

"\n";

// send breakpoint to a string
string BreakPoint::toString() const {
    stringstream out;
    out << cname << " " << gr1.chr+1 << ":" << gr1.pos1 << " to " << gr2.chr+1 << ":" << gr2.pos1 <<
      " Span: " << getSpan() << " MAPQ: " << 
       gr1.mapq << "/" << gr2.mapq << " homology: " << 
       homology << " insertion: " << insertion << " NumDups: " << num_dups << " Nsplit: " << 
       nsplit << " Tsplit: " << tsplit << " Tdisc: " << dc.tcount << " Ndisc: " << dc.ncount 
	<< " Ncigar " << ncigar << " Tcigar " << tcigar<< " cname " << cname; // << " isBest: " << isBest;
    return out.str();
}

// test whether are same break 
/*bool BreakPoint::sameBreak(BreakPoint &bp) const {
    return bp.refID1 == refID1 && bp.pos1 == pos1 && bp.refID2 == refID2 && bp.pos2 == pos2;
    }*/

// order them
void BreakPoint::order() {

    if (gr1 < gr2)
      return;

    GenomicRegion tmp = gr1;
    gr1 = gr2;
    gr2 = tmp;
    unsigned tmptsplit1 = tsplit1;  tsplit1 = tsplit2;  tsplit2 = tmptsplit1;
    unsigned tmpnsplit1 = nsplit1;  nsplit1 = nsplit2;  nsplit2 = tmpnsplit1;

}

// return whether a bp is good to move on
bool BreakPoint::isGoodSomatic(int mapq, size_t tsplit_cutoff, size_t nsplit_cutoff) const {

  mapq = 0;

  // set whether there is enoguh somatic read support
  bool issom = min(tsplit1, tsplit2) >= tsplit_cutoff; // split
  issom = issom || (dc.tcount >= tsplit_cutoff && dc.getMeanMapq() >= 0); //discordant

  // make sure that there are no germline reads
  bool nogerm = min(nsplit1, nsplit2) <= nsplit_cutoff; //split
  bool good_enough = static_cast<double>(min(tsplit1, tsplit2)) / static_cast<double>(max(nsplit1, nsplit2)) >= 5 && max(nsplit1, nsplit2) < 2;
  nogerm = nogerm && (dc.ncount == 0 || (dc.ncount == 1 && dc.tcount >= 10) );
  nogerm = nogerm || good_enough;

  // if it was aligned to the right spot, forget mapq, it's OK
  //int tmp_mapq1 = local1 ? 100 : mapq1; 
  //int tmp_mapq2 = local2 ? 100 : mapq2;
  
  int tmp_mapq1 = gr1.mapq;
  int tmp_mapq2 = gr2.mapq;

  // if the matching length is good enough, keep it
  tmp_mapq1 = matchlen1 > 250 ? 100 : tmp_mapq1;
  tmp_mapq2 = matchlen2 > 250 ? 100 : tmp_mapq2;

  // make sure the mapq quality is high enough
  bool qual = (min(tmp_mapq1, tmp_mapq2) >= mapq) || dc.tcount >= 1;

  // if it has split and discordant support, forget about mapq if it's high enoguh
  qual = qual || (min(tsplit1, tsplit2) >= 1 && dc.tcount >= 3 && max(gr1.mapq, gr2.mapq) >= 30);

  //bool enough_tum_reads = static_cast<double>(tall)/static_cast<double>(nall) >= 2 && max(tsplit1, tsplit2) > 0 && max(nsplit1, nsplit2) == 0 && tall >= 5;
  
  bool passall = issom && nogerm && qual; 

  return passall;
}

// return whether a bp is good to move on
bool BreakPoint::isGoodGermline(int mapq, size_t allsplit) const {

  mapq = 0;

  // set whether there is enoguh somatic read support
  bool isger = min(nsplit1, nsplit2) >= allsplit; // split
  int total_disc_count = dc.ncount + dc.tcount;
  isger = isger || (dc.ncount >= 1 && total_disc_count >= 12 && dc.getMeanMapq() > 20); //discordant

  // is this close enoguht to be somatic?
  bool good_enough = static_cast<double>(min(tsplit1, tsplit2)) / static_cast<double>(max(nsplit1, nsplit2)) >= 5 && max(nsplit1, nsplit2) < 2;

  // be extra extra careful for germline inter-chrom
  bool inter_ok = (gr1.chr == gr2.chr) || 
    ( (gr1.mapq == 60 && gr2.mapq == 60) && (gr1.nm < 5 && gr2.nm < 5)  && (matchlen1 > 70 && matchlen2 > 70)  && total_disc_count >= 6 && dc.getMeanMapq() > 20);

  int tmp_mapq1 = gr1.mapq;
  int tmp_mapq2 = gr2.mapq;

  // if the matching length is good enough, keep it
  tmp_mapq1 = matchlen1 > 250 ? 100 : tmp_mapq1;
  tmp_mapq2 = matchlen2 > 250 ? 100 : tmp_mapq2;

  // make sure the mapq quality is high enough
  bool qual = min(tmp_mapq1, tmp_mapq2) >= mapq || dc.ncount >= 1;

  // if it has split and discordant support, forget about mapq, it's OK
  qual = qual || (min(nsplit1+tsplit, nsplit2+tsplit) >= 2 && total_disc_count >= 3);

  //bool enough_tum_reads = static_cast<double>(tall)/static_cast<double>(nall) >= 2 && max(tsplit1, tsplit2) > 0 && max(nsplit1, nsplit2) == 0 && tall >= 5;
  
 bool passall = isger && qual && inter_ok && !good_enough; 

  return passall;

}

// make the file string
string BreakPoint::toFileString(bool noreads) {

  string sep = "\t";
  stringstream ss;

  int discordant_tum = dc.tcount;
  int discordant_norm = dc.ncount;
  string contig_name = cname;
  if (discovar) {
    discordant_tum = disco_tum;
    discordant_norm = disco_norm;
    contig_name = "discovar_" + cname;
  } 
  
  // set the evidence string
  bool isdisc = (dc.tcount + dc.ncount) != 0;
  //bool issplit = (tsplit + nsplit) != 0;
  bool issplit = dc.contig != ""; 
  if (discovar)
    evidence = "DSCVR";
  else if (num_align == 1)
    evidence = "INDEL";
  else if ( isdisc && issplit )
    evidence = "ASDIS";
  else if ( isdisc )
    evidence = "DSCRD";
  else if (num_align == 2)
    evidence = "ASSMB";
  else if (num_align > 2)
    evidence = "COMPL";
  assert(evidence.length());
  
  confidence = "";
  // low split, high span, no discordant. Suspicous

  int disc_count = dc.tcount + dc.ncount;
  int split_count = tsplit + nsplit;
  bool germ = dc.ncount > 0 || nsplit > 0;
  int total_count = disc_count + split_count;

  int span = getSpan();

  // check assembly -only ones
  if (num_align > 1 && !hasDiscordant()) {

    if (split_count < 6 && span > 1500) 
      confidence = "NODISC";
    else if (max(gr1.mapq, gr2.mapq) != 60 || min(gr1.mapq, gr2.mapq) <= 50) 
      confidence = "LOWMAPQ";
    else if ( split_count <= 3 && span <= 1500)
      confidence = "WEAKASSEMBLY";
    else if (min_end_align_length <= 40 || (germ && span == -1) || (germ && span > 1000000)) // super short alignemtns are not to be trusted. Also big germline events
      confidence = "WEAKASSEMBLY";
    else
      confidence = "PASS";

    // score ones with both assembly and discordant
  } else if (num_align > 1 && hasDiscordant()) {

    double min_disc_mapq = min(dc.getMeanMapq(true), dc.getMeanMapq(false));
    int min_assm_mapq = min(gr1.mapq, gr2.mapq);
    //double max_disc_mapq = max(dc.getMeanMapq(true), dc.getMeanMapq(false));
    int max_assm_mapq = max(gr1.mapq, gr2.mapq);
    
    if ( (min_disc_mapq < 10 && min_assm_mapq < 30) || (max_assm_mapq < 40))
      confidence = "LOWMAPQ";
    else if ( total_count < 4 || (germ && (total_count <= 6) )) // stricter about germline
      confidence = "WEAKASSEMBLY";
    else if ( total_count < 15 && germ && span == -1) // be super strict about germline interchrom
      confidence = "WEAKASSEMBLY";
    else
      confidence = "PASS";

    // score ones with discordant only
  } else if (num_align == 0) {

    if (min(gr1.mapq, gr2.mapq) < 30 || max(gr1.mapq, gr2.mapq) <= 36) // mapq here is READ mapq (37 max)
      confidence = "LOWMAPQ";
    else if ( disc_count < 8 || (dc.ncount > 0 && disc_count < 12) )  // be stricter about germline disc only
      confidence = "WEAKDISC";
    else 
      confidence = "PASS";
    
    // indels
  } else if (num_align == 1) {

    if ( split_count < 4)
      confidence="WEAKASSEMBLY";
    else if ( split_count < 8 && (tcigar+ncigar) < 3 && getSpan() < 6)
      confidence="WEAKCIGARMATCH";
    else if (gr1.mapq != 60)
      confidence="LOWMAPQ";
    else if (seq.find("AAAAAAAAAAA") != string::npos || seq.find("TTTTTTTTTTT") != string::npos || seq.find("TGTGTGTGTGTGTGTGTGTGTGTGTG") != string::npos)
      confidence="REPEAT";
    else
      confidence="PASS";
  } 

  assert(confidence.length());

  assert(getSpan() > -2);

  if (!noreads) {
    string supporting_reads = "";
    unordered_map<string, bool> supp_reads;
    
    //add the discordant reads
    for (auto& r : dc.reads) {
      string tmp;
      r_get_SR(r.second,tmp);
      supp_reads[tmp] = true;
    }
    for (auto& r : dc.mates) {
      string tmp;
      r_get_SR(r.second, tmp);
      supp_reads[tmp] = true;
    }
    
    //add the reads from the breakpoint
    for (auto& r : reads) {
      string tmp;
      r_get_SR(r,tmp);
      supp_reads[tmp];
    }
    
    // print reads to a string
    // delimit with a ,
    size_t lim = 0;
    for (auto& i : supp_reads) {
      if (++lim > 50)
	break;
      supporting_reads = supporting_reads + "," + i.first;
    }
    if (supporting_reads.size() > 0)
      supporting_reads = supporting_reads.substr(1, supporting_reads.size() - 1); // remove first _
    
    if (read_names.length() == 0)
      read_names = supporting_reads;
  } else {
    read_names = "";
  }

  // TODO convert chr to string with treader
  ss << gr1.chr+1 << sep << gr1.pos1 << sep << gr1.strand << sep 
     << gr2.chr+1 << sep << gr2.pos1 << sep << gr2.strand << sep 
     << getSpan() << sep
     << gr1.mapq <<  sep << gr2.mapq << sep 
     << nsplit << sep << tsplit << sep
     << discordant_norm << sep << discordant_tum << sep
     << ncigar << sep << tcigar << sep
    //<< nall << sep << tall << sep 
     << (homology.length() ? homology : "x") << sep 
     << (insertion.length() ? insertion : "x") << sep 
     << contig_name << sep
     << num_align << sep 
     << confidence << sep << evidence << sep << (read_names.length() ? read_names : "x") << sep 
     << pon;

  return ss.str();
  
}

// print to file string
/*void BreakPoint::printToFile(ofstream &of, const BamAlignmentVector &bamreads) {
  string sep = ",";
  
  int discordant_tum = dc.tcount;
  int discordant_norm = dc.ncount;
  string contig_name = cname;
  if (discovar) {
    discordant_tum = disco_tum;
    discordant_norm = disco_norm;
    contig_name = "discovar_" + cname;
  } 
  
  // set the evidence string
  bool isdisc = (dc.tcount + dc.ncount) != 0;
  //bool issplit = (tsplit + nsplit) != 0;
  bool issplit = dc.contig != ""; 
  if (discovar)
    evidence = "DSCVR";
  else if ( isdisc && issplit )
    evidence = "ASDIS";
  else if ( isdisc )
    evidence = "DSCRD";
  else if (num_align == 2)
    evidence = "ASSMB";
  else if (num_align > 2)
    evidence = "COMPL";
  
  confidence = "";
  // low split, high span, no discordant. Suspicous
  if ((tsplit + nsplit) < 6 && span > 1500 && !hasDiscordant()) 
    confidence = "NODISC";
  // assembly break, low mapq, no discordant
  else if ((tsplit + nsplit) != 0 && gr1.mapq != 60 && gr2.mapq != 60 && !hasDiscordant())
    confidence = "LOWMAPQ";
  // discordant only
  else if ((tsplit + nsplit) == 0 && hasDiscordant() && (dc.tcount + dc.ncount) < 6)
    confidence = "WEAKDISC";
  // assembly only, low evidence
  else if ( (tsplit + nsplit) <= 3 && !hasDiscordant() && span > 1500)
    confidence = "WEAKASSEMBLY";
  // disc is 6+, has assembly + disc, assembly is 3+ with 60/60 mapq
  else 
    confidence = "PASS";

  assert(span > -2);

  string supporting_reads = "";
  unordered_map<string, bool> reads;
  //add the discordant reads
  for (unordered_map<string, bool>::const_iterator it = dc.qnames.begin(); it != dc.qnames.end(); it++) {
    if (reads.find(it->first) == reads.end())
      if (it->first.size() > 0)
	reads.insert(pair<string, bool>(it->first.substr(1, it->first.size()-1), true));
  }
  //add the contig reads
  //BamAlignmentVector bamreads = (*contigs)[cname].m_bamreads;
  for (BamAlignmentVector::const_iterator it = bamreads.begin(); it != bamreads.end(); it++) {
    if (reads.find(it->Name) == reads.end())
      reads.insert(pair<string, bool>(it->Name, true));
  }
   
  // print reads to a string
  for (unordered_map<string, bool>::const_iterator it = reads.begin(); it != reads.end(); it++) 
    supporting_reads = supporting_reads + "_" + it->first;
  if (supporting_reads.size() > 0)
    supporting_reads = supporting_reads.substr(1, supporting_reads.size() - 1); // remove first _

  of << evidence << sep << refID1+1 << sep << pos1 << sep << strand1 << sep 
     << refID2+1 << sep << pos2 << sep << strand2 << sep
     << mapq1 <<  sep << mapq2 << sep 
     << nsplit << sep << tsplit << sep
     << discordant_norm << sep << discordant_tum << sep
    //<< nall << sep << tall << sep 
     << homology << sep << insertion << sep << contig_name << sep
     << span << sep << num_dups << sep << num_align << sep 
     << confidence << sep << supporting_reads << endl;
}
*/

// make a breakpoint from a discordant cluster 
BreakPoint::BreakPoint(DiscordantCluster tdc) {

  dc = tdc;
  gr1.pos1 = (tdc.reg1.strand == '+') ? tdc.reg1.pos2 : tdc.reg1.pos1;
  gr1.pos2 = gr1.pos1;
  gr2.pos1 = (tdc.reg2.strand == '+') ? tdc.reg2.pos2 : tdc.reg2.pos1;
  gr2.pos2 = gr2.pos1;
  gr1.chr = tdc.reg1.chr;
  gr2.chr = tdc.reg2.chr;
  cname = tdc.cluster;
  gr1.strand = tdc.reg1.strand;
  gr2.strand = tdc.reg2.strand;
  
  gr1.mapq = tdc.getMeanMapq(false); 
  gr2.mapq = tdc.getMeanMapq(true); // mate

}

bool BreakPoint::hasDiscordant() const {
  return !dc.reg1.isEmpty();
}


/** 
 * Has at least two supporting reads
 */
bool BreakPoint::hasMinimal() const {

  int total = tsplit + nsplit + dc.tcount + dc.ncount;

  if (total >= 2)
    return true;
  else
    return false;

}

bool BreakPoint::operator==(const BreakPoint &bp) const {

  return (gr1.chr == bp.gr1.chr && gr1.pos1 == bp.gr1.pos1 && gr2.pos1 == bp.gr2.pos1);

}


void runRefilterBreakpoints(int argc, char** argv) {
  
  parseBreakOptions(argc, argv);

  bopt::output_file = bopt::analysis_id + ".filtered.bps.txt.gz";
  if (bopt::verbose > 0) {
    cout << "Input bps file:  " << bopt::input_file << endl;
    cout << "Output bps file: " << bopt::output_file << endl;
    cout << "Panel of normals file: " << bopt::pon << endl;
    cout << "Analysis id: " << bopt::analysis_id << endl;
  }

  igzstream iz(bopt::input_file.c_str());
  if (!iz) {
    cerr << "Can't read file " << bopt::input_file << endl;
    exit(EXIT_FAILURE);
  }

  unique_ptr<PON> pmap;
  BreakPoint::readPON(bopt::pon, pmap);
  
  ogzstream oz(bopt::output_file.c_str(), ios::out);
  if (!oz) {
    cerr << "Can't write to output file " << bopt::output_file << endl;
    exit(EXIT_FAILURE);
  }
  
  if (bopt::verbose)
    cout << "...refiltering variants" << endl;
  
  // set the header
  oz << BreakPoint::header() << endl;

  string line;
  //skip the header
  getline(iz, line, '\n');

  while (getline(iz, line, '\n')) {
    BreakPoint bp(line);

    // check if in panel of normals
    bp.checkPon(pmap);

    if (bp.hasMinimal() && bp.gr1.chr < 24 && bp.gr2.chr < 24)
      oz << bp.toFileString(bopt::noreads) << endl;
  }
  
  oz.close();

  // make the VCF files
  if (bopt::verbose)
    cout << "...converting filtered.bps.txt.gz to vcf files" << endl;
  VCFFile filtered_vcf(bopt::output_file, bopt::ref_index.c_str(), '\t', bopt::analysis_id);

  // output the vcfs
  string basename = bopt::analysis_id + ".broad-snowman.DATECODE.filtered.";
  filtered_vcf.writeIndels(basename, true); // zip them -> true
  filtered_vcf.writeSVs(basename, true);
  
}

BreakPoint::BreakPoint(string &line) {

  istringstream iss(line);
  string val;
  size_t count = 0;
  while (getline(iss, val, '\t')) {
    switch(++count) {
    case 1: gr1.chr = stoi(val) - 1; break;
    case 2: gr1.pos1 = stoi(val); gr1.pos2 = gr1.pos1; break;
    case 3: gr1.strand = val.at(0); break;
    case 4: gr2.chr = stoi(val) - 1; break;
    case 5: gr2.pos1 = stoi(val); gr2.pos2 = gr2.pos1; break;
    case 6: gr2.strand = val.at(0); break;
      //case 7: span = stoi(val); break;
    case 8: gr1.mapq = stoi(val); break;
    case 9: gr2.mapq = stoi(val); break;
    case 10: nsplit = stoi(val); break;
    case 11: tsplit = stoi(val); break;
    case 12: dc.ncount = stoi(val); break;
    case 13: dc.tcount = stoi(val); break;
    case 14: ncigar = stoi(val); break;
    case 15: tcigar = stoi(val); break;
    case 16: homology = val; break;
    case 17: insertion = val; break;
    case 18: cname = val; break;
    case 19: num_align = stoi(val); break;
    case 20: confidence = val; break;
    case 21: evidence = val; break;
    case 22: read_names = val; break;
    }
  }

  if (evidence == "INDEL")
    isindel = true;

  //debug
  if (num_align == 1) {
    dc.tcount = 0;
    dc.ncount = 0;
  }


}


// parse the command line options
void parseBreakOptions(int argc, char** argv) {
  bool die = false;

  if (argc <= 2) 
    die = true;

  string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
      case 'h': die = true; break;
      case 'g': arg >> bopt::ref_index; break;
      case 'i': arg >> bopt::input_file; break;
      case 'v': arg >> bopt::verbose; break;
      case 'q': arg >> bopt::pon; break;
      case 'a': arg >> bopt::analysis_id; break;
      case 'x': bopt::noreads = true; break;
    }
  }

  if (bopt::input_file.length() == 0)
    die = true;
  if (bopt::output_file.length() == 0)
    die = true;

  if (die) {
      cout << "\n" << BP_USAGE_MESSAGE;
      exit(1);
    }
}

void BreakPoint::splitCoverage(ReadVec &bav) {

  assert(cname.length());

  // make sure we're not double counting
  assert(tsplit == 0 && nsplit == 0 && nsplit1 == 0 && nsplit2 == 0 && tsplit1 == 0 && tsplit2 == 0);

  // total number of valid split reads
  tsplit = 0;
  nsplit = 0;
  
  int rightbreak1= cpos1 + SPLIT_BUFF; // read must extend this far right of break1
  int leftbreak1 = cpos1 - SPLIT_BUFF; // read must extend this far left break2
  int rightbreak2= cpos2 + SPLIT_BUFF;
  int leftbreak2 = cpos2 - SPLIT_BUFF;
  
  for (auto& j : bav) {
    
    string sr;
    r_get_Z_tag(j, "SR", sr);
    assert(sr.at(0) == 't' || sr.at(0) == 'n');

    bool tumor_read = (sr.at(0) == 't');
    string seq;
    r_get_trimmed_seq(j, seq);
    assert(seq.length() > 0);
    
    string qname = SnowUtils::GetStringTag(j, "CN").back();
    int pos = SnowUtils::GetIntTag(j, "AL").back();
    
    if (qname != cname)
      cerr << "qname " << qname << "cname " << cname << endl;

    assert(qname == cname);

    if (qname == cname) { // make sure we're comparing the right alignment
      int rightend = pos + seq.length();
      int leftend  = pos;
      bool issplit1 = (leftend <= leftbreak1) && (rightend >= rightbreak1);
      bool issplit2 = (leftend <= leftbreak2) && (rightend >= rightbreak2);
      
      // add the split reads for each end of the break
      if (issplit1 || issplit2)
	reads.push_back(j);
      
      // update the counters for each end
      if (issplit1 && tumor_read)
	tsplit1++;
      if (issplit1 && !tumor_read)
	nsplit1++;
      if (issplit2 && tumor_read)
	tsplit2++;
      if (issplit2 && !tumor_read)
	nsplit2++;
      
      // read spans both ends
      if ((issplit1 || issplit2) && tumor_read)
	tsplit++;
      if ((issplit1 || issplit2) && !tumor_read)
	nsplit++;

      /*if (issplit1 || issplit2) {
	if (tumor_read)
	  tunique++;
	else
	  nunique++;
	  }*/
      
    }
  }
  //cout << "tsplit1 " << tsplit1 << " nsplit1 " << nsplit1 << " tsplit2 " << tsplit2 << " nsplit2 " << nsplit2 << " mapq1 " << gr1.mapq << " gr2.mapq " << gr2.mapq << " cname " << cname << " bp " << *this << endl;
}

string BreakPoint::getHashString() const {
  
  string st = to_string(gr1.chr) + "_" + to_string(gr1.pos1) + "_" + to_string(this->getSpan()) + (insertion.length() ? "I" : "D");
  return st;
}

int BreakPoint::checkPon(unique_ptr<PON> &p) {

  // only built for indels for now
  if (!isindel || !p)
    return 0;

  /*string chr = to_string(gr1.chr+1);
  if (chr == "23")
    chr = "X";
  if (chr == "24")
    chr = "Y";
  if (chr == "25")
    chr = "M";
  */
  string chr = to_string(gr1.chr);

  string key1, key2, key3, key4, key5;
  int span = getSpan();
  string type = (insertion != "") ? "I" : "D";
  key1 = chr + "_" + to_string(gr1.pos1-1) + "_" + to_string(span) + type;
  key2 = chr + "_" + to_string(gr1.pos1  ) + "_" + to_string(span) + type;
  key3 = chr + "_" + to_string(gr1.pos1+1) + "_" + to_string(span) + type;
  key4 = chr + "_" + to_string(gr1.pos1-2) + "_" + to_string(span) + type;
  key5 = chr + "_" + to_string(gr1.pos1+2) + "_" + to_string(span) + type;
  
  if (p->count(key1))
    pon = max((*p)[key1], pon);
  if (p->count(key2))
    pon = max((*p)[key2], pon);
  if (p->count(key3))
    pon = max((*p)[key3], pon);
  if (p->count(key4))
    pon = max((*p)[key4], pon);
  if (p->count(key5))
    pon = max((*p)[key5], pon);

  if (pon > 0)
    cout << "pon key " << key1 << " val " << pon << endl;
  
  return pon;
  
}


void BreakPoint::readPON(string &file, unique_ptr<PON> &pmap) {

  // import the pon
  igzstream izp(file.c_str());
  if (!izp) {
    cerr << "Can't read file " << file << endl;
    exit(EXIT_FAILURE);
  }

  if (bopt::verbose)
    cout << "...importing PON data" << endl;
  pmap = unique_ptr<PON>(new PON());
  string pval;
  while (getline(izp, pval, '\n')) {
    istringstream gg(pval);
    string tval;
    string key;
    //vector<int> scount;
    int scount_total = 0;
    size_t c = 0;
    while (getline(gg, tval, '\t')) {
      if (++c == 1 && tval.length()) {
	key = tval.substr(1,tval.length() - 1);
	if (key.at(0) == 'T') // only accumulate normal
	  break;
      }
      else if (tval.length())
	try { scount_total += stoi(tval); } catch(...) { cerr << "stoi error in PON read with val " << tval << " on line " << pval << endl; }
    }

    (*pmap)[key] = scount_total;
  }

  return;

}
