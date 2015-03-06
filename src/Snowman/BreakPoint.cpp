#include "BreakPoint.h"
#include <getopt.h>
#include "gzstream.h"

namespace bopt {

  static string input_file = "";
  static string output_file = "";

  static int verbose = 1;
}

static const char* shortopts = "hi:v:o:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "input-bps",               required_argument, NULL, 'i'},
  { "output-bps",              required_argument, NULL, 'o'},
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
"  Required input\n"
"  -i, --tumor-bam                      Tumor BAM file\n"
"  -o, --tumor-bam                      Tumor BAM file\n"
"  Optional input\n"                       
"\n";

// send breakpoint to a string
string BreakPoint::toString() const {
    stringstream out;
    out << cname << " " << gr1.chr+1 << ":" << gr1.pos1 << " to " << gr2.chr+1 << ":" << gr2.pos1 <<
       " Span: " << span << " MAPQ: " << 
       gr1.mapq << "/" << gr2.mapq << " homology: " << 
       homology << " insertion: " << insertion << " NumDups: " << num_dups << " Nsplit: " << 
       nsplit << " Tsplit: " << tsplit << " Tdisc: " << dc.tcount << " Ndisc: " << dc.ncount; // << " isBest: " << isBest;
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
string BreakPoint::toFileString() {

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

  // check assembly -only ones
  if (num_align > 1 && !hasDiscordant()) {

    if ((tsplit + nsplit) < 6 && span > 1500) 
      confidence = "NODISC";
    else if (max(gr1.mapq, gr2.mapq) != 60 || min(gr1.mapq, gr2.mapq) <= 50) 
      confidence = "LOWMAPQ";
    else if ( (tsplit + nsplit) <= 3 && span <= 1500)
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
    else if ( (dc.tcount + dc.ncount + nsplit + tsplit) < 4)
      confidence = "WEAKASSEMBLY";
    else
      confidence = "PASS";

    // score ones with discordant only
  } else if (num_align == 0) {

    if (min(gr1.mapq, gr2.mapq) < 30 || max(gr1.mapq, gr2.mapq) <= 36) // mapq here is READ mapq (37 max)
      confidence = "LOWMAPQ";
    else if (dc.tcount + dc.ncount < 8)
      confidence = "WEAKDISC";
    else
      confidence = "PASS";
    
    // indels
  } else if (num_align == 1) {

    if ( (tsplit + nsplit) < 4)
      confidence="WEAKASSEMBLY";
    else if (gr1.mapq != 60)
      confidence="LOWMAPQ";
    else if (seq.find("AAAAAAA") != string::npos || seq.find("TTTTTTT") != string::npos || seq.find("TGTGTGTG") != string::npos)
      confidence="REPEAT";
    else
      confidence="PASS";
  } 

  assert(confidence.length());

  assert(span > -2);

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
  for (auto& i : supp_reads) 
    supporting_reads = supporting_reads + "," + i.first;
  if (supporting_reads.size() > 0)
    supporting_reads = supporting_reads.substr(1, supporting_reads.size() - 1); // remove first _

  if (read_names.length() == 0)
    read_names = supporting_reads;

  // TODO convert chr to string with treader
  ss << gr1.chr+1 << sep << gr1.pos1 << sep << gr1.strand << sep 
     << gr2.chr+1 << sep << gr2.pos1 << sep << gr2.strand << sep 
     << span << sep
     << gr1.mapq <<  sep << gr2.mapq << sep 
     << nsplit << sep << tsplit << sep
     << discordant_norm << sep << discordant_tum << sep
     << ncigar << sep << tcigar << sep
    //<< nall << sep << tall << sep 
     << (homology.length() ? homology : "x") << sep 
     << (insertion.length() ? insertion : "x") << sep 
     << contig_name << sep
     << num_align << sep 
     << confidence << sep << evidence << sep << (read_names.length() ? read_names : "x");

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

  
  if (gr1.chr != gr2.chr)
    span = -1;
  else
    span = abs(gr1.pos1 - gr2.pos1);

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

  if (bopt::verbose > 0) {
    cout << "Input bps file:  " << bopt::input_file << endl;
    cout << "Output bps file: " << bopt::output_file << endl;
  }

  igzstream iz(bopt::input_file.c_str());
  if (!iz) {
    cerr << "Can't read file " << bopt::input_file << endl;
    exit(EXIT_FAILURE);
  }

  ogzstream oz(bopt::output_file.c_str(), ios::out);
  if (!oz) {
    cerr << "Can't write to output file " << bopt::output_file << endl;
    exit(EXIT_FAILURE);
  }
  // set the header
  oz << BreakPoint::header() << endl;


  string line;
  //skip the header
  getline(iz, line, '\n');

  while (getline(iz, line, '\n')) {
    BreakPoint bp(line);
    if (bp.hasMinimal())
      oz << bp.toFileString() << endl;
  }
  
  oz.close();
}

/**
 *
 */
BreakPoint::BreakPoint(string &line) {

  istringstream iss(line);
  string val;
  size_t count = 0;
  while (getline(iss, val, '\t')) {
    switch(++count) {
    case 1: gr1.chr = stoi(val); break;
    case 2: gr1.pos1 = stoi(val); gr1.pos2 = gr1.pos1; break;
    case 3: gr1.strand = val.at(0); break;
    case 4: gr2.chr = stoi(val); break;
    case 5: gr2.pos1 = stoi(val); gr2.pos2 = gr2.pos1; break;
    case 6: gr2.strand = val.at(0); break;
    case 7: span = stoi(val); break;
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
      case 'i': arg >> bopt::input_file; break;
      case 'o': arg >> bopt::output_file; break;
      case 'v': arg >> bopt::verbose; break;
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
