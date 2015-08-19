#include "vcf.h"

#include <string> 
#include <iostream>
#include <getopt.h>
#include <time.h>
#include <sstream>
#include <regex>

#include "SnowTools/SnowUtils.h"
#include "SnowTools/SnowToolsCommon.h"
#include "SnowTools/gzstream.h"
#include "SnowTools/GenomicRegionCollection.h"

//#include "SnowTools/reads.h"
#include "htslib/tbx.h"
#include "htslib/bgzf.h"

using namespace std;

static const char *VCF_USAGE_MESSAGE =
"Usage: snowman vcf [OPTION] \n\n"
"  Description: Utility for handling VCF files from Snowman (Broad), dRanger (Broad), DELLY (EMBL) and Brass (WTSI)\n"
"\n"
" General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h, --help                           Display this help and exit\n"
"  -o, --output-vcf                     File to out the VCF to.\n"
"  -c, --write-csv                      File to output the VCF as a comma separated file. Default: not output\n"
"  -g, --reference-genome               Path to indexed reference genome to be used by BWA-MEM. Default is Broad hg19 (/seq/reference/...)\n"
"      --no-zip                         Dont bgzip or tabix the outputs.\n"
  //"      --stats                          Write a simple stats file as specified.\n"
" Optional Input\n"
"  -i, --snowman-csv                    A csv breakpoint file (usually breakpoints.somatic.txt) that needs to be converted to VCF.\n"
"  -a, --id-string                      String specifying the analysis ID to be used as part of ID common.\n"
"  -t, --tumor-bam-string               Tumor BAM\n"
"  -n, --normal-bam-string              Normal BAM\n"
"  -s, --snowman-vcf                    A VCF produced from SnowmanSV.\n"
"  -r, --dranger-vcf                    A VCF produced from dRanger.\n"
"  -b, --brass-vcf                      A VCF produced from Brass.\n"
"  -d, --delly-vcf                      A VCF produced from DELLY.\n"
" Optional\n"
"  -y, --padding                        Amount of padding to place around breakpoint to consider overlap (BPs off by less than this are overlaps). Default: 100\n"
"      --old-snowman                    \n"
"      --include-nonpass                Include rearrangements that don't PASS in the final output.\n"
"  -c, --write-csv                      Provide filename to Output the VCF as a comma separated file. Default: not output\n"
"\n";

namespace vopt {
  static size_t verbose = 1;
  
  static bool old_snowman = false;

  static string snowman = "";
  static string dranger = "";
  static string brass = "";
  static string delly = "";

  static bool filtered = false;

  static bool zip = true;

  static string csv = "";
  static string normal = "normal";
  static string tumor = "tumor";
  static string analysis_id = "no_id";

  static string ref_index = SnowTools::REFHG19;

  static int pad = 100;
  //static string outvcf = "vars.vcf";

  static bool include_nonpass = false;
  static string as_csv = "";;

  static string stats = "stats.txt";
}

static const char* shortopts = "hv:a:t:n:r:b:s:d:i:o:y:m:g:xzjc:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "snowman-csv",             required_argument, NULL, 'i' },
  { "verbose",                 required_argument, NULL, 'v' },
  { "tumor-bam-string",        required_argument, NULL, 't' },
  { "normal-bam-string",       required_argument, NULL, 'n' },
  { "dranger-vcf",             required_argument, NULL, 'r' },
  { "reference-genome",        required_argument, NULL, 'g' },
  { "snowman-vcf",             required_argument, NULL, 's' },
  { "delly-vcf",               required_argument, NULL, 'd' },
  { "brass-vcf",               required_argument, NULL, 'b' },
  { "id-string",               required_argument, NULL, 'a' },
  { "padding",                 required_argument, NULL, 'y' },
  { "output-vcf",              required_argument, NULL, 'o' },
  { "old-snowman",             no_argument, NULL, 'x'},
  { "include-nonpass",         no_argument, NULL, 'j'},
  { "no-zip",                  no_argument, NULL, 'z'},
  { "write-csv",               required_argument, NULL, 'c'},
  { NULL, 0, NULL, 0 }
};

static faidx_t * findex;
static InfoMap flag_map;
static InfoMap union_info_fields;

// comparator for info fields
// lhs < rhs
// want READ_ID to be > than everything
bool compareInfoFields(const pair<string,string> &lhs, const pair<string,string> &rhs) {
  return ( (rhs.first == "READ_ID" && lhs.first != "READ_ID") || ( (rhs.first != "READ_ID" && lhs.first != "READ_ID") && lhs.first < rhs.first));
}

void tabixVcf(const string &fn);

void runVCF(int argc, char** argv) {

  /*
  //debug
  char* fnameb = "tmp.bcf";
  write_bcf(fnameb);
  char* fnamev = "tmp.vcf";
  bcf_to_vcf(fnameb);
  exit(1);
  ///////////
  */

  cout << " running snowman vcf" << endl;
  parseVCFOptions(argc, argv);

  if (vopt::verbose > 0) {
    cout << "-- Reference genome " << vopt::ref_index << endl;
    cout << "-- Input VCFs:  " << endl;
    cout << "   Snowman:           " << vopt::snowman << endl;
    cout << "   dRanger:           " << vopt::dranger << endl;
    cout << "   DELLY:             " << vopt::delly << endl;
    cout << "   Brass:             " << vopt::brass << endl;
    cout << "-- Input CSV:         " << endl;
    cout << "   Snowman CSV:       " << vopt::csv << endl;
    //    cout << "-- Output VCF:        " << vopt::outvcf << endl;
    cout << "-- Annotation strings:" << endl;
    cout << "   Analysis ID        " << vopt::analysis_id << endl;
    cout << "   Tumor:             " << vopt::tumor << endl;
    cout << "   Normal:            " << vopt::normal << endl;
    cout << "-- Other options:     " << endl;
    cout << "   Padding:           " << vopt::pad << endl;
    cout << "   Old Snowman:       " << vopt::old_snowman << endl;
    cout << "   Include non-PASS:  " << (vopt::include_nonpass ? "ON" : "OFF") << endl;
    cout << "   Write as CSV:      " << vopt::as_csv << endl;
    cout << "   Is from a filtered CSV?:      " << (vopt::filtered ? "TRUE" : "FALSE") << endl;
    cout << endl;
  }

  // if dranger only exists, move it over and stop
  /*if (vopt::snowman == "" && vopt::dranger != "") {
    cerr << "Warning: SnowmanSV VCF does not exist. Copying input dRanger VCF to output" << endl;
    string runthis = "cp " + vopt::dranger;
    system(runthis.c_str());
    return;
    } */
  
  // load the reference genome
  if (vopt::csv != "") {
    if (vopt::verbose > 0)
      cout << "attempting to load: " << vopt::ref_index << endl;
    findex = fai_load(vopt::ref_index.c_str());  // load the reference
  } else {
    findex = NULL;
  }

  VCFFile snowvcf, dranvcf;
  if (vopt::snowman != "")
    snowvcf = VCFFile(vopt::snowman, "snowman");
  if (vopt::dranger != "")
    dranvcf = VCFFile(vopt::dranger, "dranger");

  // if only one file specified, outpt
  if (vopt::dranger == "" && vopt::snowman != "") {
    //if (vopt::outvcf != "")
    //  snowvcf.write(vopt::analysis_id + ".broad-snowman.DATECODE.");
    if (vopt::as_csv != "")
      snowvcf.writeCSV();
  } else if (vopt::snowman == "" && vopt::dranger != "") {
    //if (vopt::outvcf != "")
    //  dranvcf.write(vopt::analysis_id + ".broad-dRanger.DATECODE.");
    if (vopt::as_csv != "")
      dranvcf.writeCSV();
  } else if (vopt::dranger != "" && vopt::snowman != "") {
    VCFFile merged_vcf = mergeVCFFiles(dranvcf, snowvcf);
    merged_vcf.write(vopt::analysis_id + ".broad-merged.DATECODE.");
    //if (vopt::outvcf != "") {
    //  merged_vcf.write(vopt::analysis_id + ".broad-merged.DATECODE.");
    // }
    if (vopt::as_csv != "")
      merged_vcf.writeCSV();
  }
  
  if (vopt::csv != "") {
    cout << "...converting Snowman csv to vcf" << endl;
    VCFFile snowcsv(vopt::csv, vopt::ref_index.c_str(), '\t', vopt::analysis_id);
    //snowcsv.filename = vopt::outvcf;
    snowcsv.write(vopt::analysis_id + ".broad-snowman.DATECODE.");
   }

  return;

}

//parse the options supplied to snowman vcf
void parseVCFOptions(int argc, char** argv) {

  bool die = false;

  if (argc < 2) 
    die = true;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'a': arg >> vopt::analysis_id; break;
    case 'g': arg >> vopt::ref_index; break;
    case 'h': die = true; break;
    case 'z': vopt::zip = false; break;
    case 't': arg >> vopt::tumor; break;
    case 'n': arg >> vopt::normal; break;
    case 's': arg >> vopt::snowman; break;
    case 'd': arg >> vopt::delly; break;
    case 'v': arg >> vopt::verbose; break;
    case 'r': arg >> vopt::dranger; break;
    case 'b': arg >> vopt::brass; break;
    case 'y': arg >> vopt::pad; break;
      //case 'o': arg >> vopt::outvcf; break;
    case 'i': arg >> vopt::csv; break;
    case 'x': vopt::old_snowman = true; break;
    case 'j': vopt::include_nonpass = true; break;
    case 'c': arg >> vopt::as_csv; break;
    }
  }

  /*if (vopt::as_csv == "" && vopt::outvcf == "") {
    die = true;
    cerr << "Must output either csv or vcf" << endl;
    }*/
  
  // make sure that the csv exists
  if (vopt::csv != "" && !SnowTools::read_access_test(vopt::csv)) {
    die = true;
    cerr << "CSV file " << vopt::csv << " does not exist " << endl;
  }

  // make sure that the inputs exist
  if (!SnowTools::read_access_test(vopt::dranger) && vopt::dranger != "") {
    die = true;
    cerr << "############Supplied dRanger VCF does not exist: " << vopt::dranger << endl;
  }
  if (!SnowTools::read_access_test(vopt::brass) && vopt::brass != "") {
    die = true;
    cerr << "############Supplied Brass VCF does not exist: " << vopt::brass << endl;
  }
  if (!SnowTools::read_access_test(vopt::snowman) && vopt::snowman != "") {
    die = true;
    cerr << "############Supplied snowman VCF does not exist: " << vopt::snowman << endl;
  }
  if (!SnowTools::read_access_test(vopt::delly) && vopt::delly != "") {
    die = true;
    cerr << "############Supplied DELLY VCF does not exist: " << vopt::delly << endl;
  }

  // check if this is a filtered bps
  if (vopt::csv.find("filter") != string::npos)
    vopt::filtered = true;

  // something went wrong, kill
  if (die) {
    cout << "\n" << VCF_USAGE_MESSAGE;
    exit(1);
  }

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

//set the filedate string to the current date
void VCFHeader::setCurrentDate() {

  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  stringstream month;
  stringstream mdate;
  if ( (now->tm_mon+1) < 10)
    month << "0" << now->tm_mon+1;
  else 
    month << now->tm_mon+1;
  mdate << (now->tm_year + 1900) << month.str() <<  now->tm_mday;

  filedate = mdate.str();
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

//create a VCF header record from a file
VCFHeader::VCFHeader(string file) {
  
  //open the file
  igzstream inFile(file.c_str());
  
  if (!inFile) {
    cerr << "Can't read file " << file << " for parsing VCF" << endl;
    return;
  }

  // file format
  regex reg_ff("^##fileformat=(.*+)");
  // file date
  regex reg_date("^##fileDate=(.*)");
  // source
  regex reg_source("^##source=(.*)");
  // reference
  regex reg_ref("^##reference=(.*)");
  // info
  regex reg_info("^##INFO=<ID=(.*?),Number=([0-9.]),Type=(.*?),Description=\"(.*?)\">");
  // filter
  regex reg_filter("^##FILTER=<ID=(.*?),Description=\"(.*?)\">");
  // format
  //regex reg_format("^##FORMAT=<ID=(.*?),Description=\"(.*?)\">");
  regex reg_format("^##FORMAT=<ID=(.*?),Number=([0-9.]),Type=(.*?),Description=\"(.*?)\">");
  // sample
  regex reg_sample("^##SAMPLE=<ID=(.*?)>");
  // colnames
  regex reg_colnam("^#CHR(.*)");
  // contig
  regex reg_contig("contig=<ID=(.*?),assembly=(.*?),length=([0-9]+),species=([A-Z]+)>");

  // read it in
  string line;
  while (getline(inFile, line, '\n')) {
    
    if (line.size() > 0)
      if (line.at(0) != '#')
	break;
    smatch match;
    // perform the regex matching
    if (regex_search(line, match, reg_ff))
      fileformat = match[1].str();
    else if (regex_search(line, match, reg_date))
      filedate = match[1].str();
    else if (regex_search(line, match, reg_source))
      source = match[1].str();
    else if (regex_search(line, match, reg_ref))
      reference = match[1].str();
    else if (regex_search(line, match, reg_info))
      addInfoField(match[1].str(), match[2].str(), match[3].str(), match[4].str());
    else if (regex_search(line, match, reg_filter))
      addFilterField(match[1].str(), match[2].str());
    else if (regex_search(line, match, reg_format))
      addFormatField(match[1].str(), match[2].str(), match[3].str(), match[4].str());
    else if (regex_search(line, match, reg_sample))
      addSampleField(match[1].str());
    else if (regex_search(line, match, reg_colnam))
      colnames = "#CHR" + match[1].str();
    else if (regex_search(line, match, reg_contig))
      addContigField(match[1].str(), match[2].str(), match[3].str(), match[4].str());
  }

  inFile.close();
}

// add a contig field
void VCFHeader::addContigField(string id, string assembly, string length, string species) {
  
  if (contigfieldmap.find(id) != contigfieldmap.end()) {
    cerr << "Warning: Contig with this ID already exists" << endl;
    return;
  }
  
  string net = "assembly=" + assembly + ",length=" + length + ",species=" + species; 
  contigfieldmap[id] = net;
  return;
}

// create a VCFEntry out of the string
VCFEntry::VCFEntry(string line, string method) {

  string info;
  regex reg;
  if (method == "dranger" || method == "snowman") 
    reg = regex("(.*?)\t([0-9]+)\t(.*?)\t(.*?)\t(.*?)\t(.*)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)(\t|$).*");    
  else 
    cerr << "Method not recognized: " << method << endl;

  smatch match;
  if (regex_search(line, match, reg)) {
    chr = SnowTools::GenomicRegion::chrToNumber(match[1].str());
    pos = stoi(match[2].str());
    id  = match[3].str();
    ref = match[4].str();
    alt = match[5].str();
    qual= match[6].str();
    filter = match[7].str();
    info = match[8].str();
    format = match[9].str();
    samp1 = match[10].str();
    samp2 = match[11].str();

    // scrub the chr from alt
    alt = SnowTools::scrubString(alt, "chr");

  } else {
    cerr << "VCF line " << line << " not matched" << endl;
  }

  // parse the INFO field
  string val;
  regex inforeg("(.*?)=(.*)");  
  istringstream iss(info);
  while (getline(iss, val, ';')) {
    smatch infomatch;
    if (regex_match(val, infomatch, inforeg)) {

      /* tmp hack to remove discordant "contigs" */
      if (infomatch[1].str() != "SCTG" || (infomatch[1].str() == "SCTG" && infomatch[2].str().at(0) == 'c'))
	info_fields[infomatch[1].str()] = infomatch[2].str();
    }
  }

  // parse the FORMAT fields
  format_fields = FormatStringToFormatRecordMap(format, samp1, samp2);
  
  // parse the ID field
  //regex reg_sanger("^([0-9]+)_[0-9]+");
  regex reg_id("(.*?):[1-2]$");
  smatch match_id;
  if (regex_match(id, match_id, reg_id))
    idcommon = match_id[1];
  else
    cerr << "Warning: commond pair ID not properly parsed from ID string " << id << endl;

}

// print out the VCF Entry
std::ostream& operator<<(std::ostream& out, const VCFEntry& v) {

  /* print the INFO field */
  // move to a vector to be sorted
  bool imprecise = false;
  vector<pair<string, string> > tmpvec; // id, evertythign else
  for (InfoMap::const_iterator it = v.info_fields.begin(); it != v.info_fields.end(); it++) {
    if (it->first == "IMPRECISE")
      imprecise = true;
    if (it->second != "")
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
  out << SnowTools::GenomicRegion::chrToString(v.chr) << sep 
      << v.pos << sep << v.id << sep << v.ref << sep << v.alt << sep << v.qual << sep
      << v.filter << sep << info << sep << v.format << sep << v.samp1 << sep << v.samp2;

  return out;
}

// create a VCFFile from a VCF file
VCFFile::VCFFile(string file, string tmethod) {

  if (vopt::verbose > 0) 
    cout << "Reading/parsing VCF file from method " << tmethod << endl;

  filename = file;
  method = tmethod;

  //open the file
  igzstream inFile(file.c_str());
  if (!inFile) {
    cerr << "Can't read file " << file << " for parsing VCF" << endl;
    return;
  }

  //grab the header
  sv_header = VCFHeader(file);
  if (sv_header.source == "")
    sv_header.source = tmethod;

  VCFEntryVec entries;

  // parse and create the VCFEntry lines
  string line;
  while (getline(inFile, line, '\n')) {
    if (line.size() > 0)
      if (line.at(0) != '#')
	entries.push_back(VCFEntry(line, tmethod));
  }

  sort(entries.begin(), entries.end());

  VCFEntryPairMap tmp_map;
  // populate the pairs
  for (VCFEntryVec::iterator it = entries.begin(); it != entries.end(); it++) {
    string entry_pair_key = it->idcommon + tmethod;
    VCFEntryPairMap::iterator ff = tmp_map.find(entry_pair_key);
    if (ff == tmp_map.end()) {

      VCFEntryPair thispair(*it);
      thispair.method = tmethod;
      thispair.idcommon = it->idcommon;

      // add the split/discordant reads
      try {
	thispair.tdisc  = stoi(it->format_fields["NALT_RP"].second);
	thispair.tsplit = stoi(it->format_fields["NALT_SR"].second);
      } catch (...) {
	cerr << "Caught error with stoi in VCFFile(filename, method) for method " << tmethod << endl;
	cerr << "Tumor disc "  << it->format_fields["NALT_RP"].second << endl;
	cerr << "Tumor split " << it->format_fields["NALT_SR"].second << endl;
	exit(EXIT_FAILURE);
      }
      tmp_map.insert(pair<string, VCFEntryPair>(entry_pair_key, thispair));
    } else {
      ff->second.e2 = *it; // add the second pair
    }
  }

  // remove bad entries
  for (VCFEntryPairMap::const_iterator it = tmp_map.begin(); it != tmp_map.end(); it++) {  
    if (it->second.e1 < it->second.e2)
      entry_pairs.insert(pair<string, VCFEntryPair>(it->first, it->second));
  }

}

// print out the VCFFile
std::ostream& operator<<(std::ostream& out, const VCFFile& v) {

  out << v.sv_header << endl;

  VCFEntryVec tmpvec;
  size_t id_counter = 0;

  // put the pair maps into a vector
  for (VCFEntryPairMap::const_iterator it = v.entry_pairs.begin(); it != v.entry_pairs.end(); it++) {

    if (v.dups.count(it->first) == 0) { // dont include duplicate entries
      // renumber the ids
      id_counter++;
      VCFEntryPair tmppair = it->second;
      if (vopt::dranger != "" && vopt::snowman != "") {
	
	tmppair.e1.idcommon = to_string(id_counter) + ":" + vopt::analysis_id;
	tmppair.e2.idcommon = to_string(id_counter) + ":" + vopt::analysis_id;
	tmppair.e1.id = tmppair.e1.idcommon + ":1";
	tmppair.e2.id = tmppair.e2.idcommon + ":2";
	tmppair.e1.info_fields["MATEID"] = tmppair.e2.id;
	tmppair.e2.info_fields["MATEID"] = tmppair.e1.id;
      }
      
      tmpvec.push_back(tmppair.e1);
      if (!tmppair.indel)
	tmpvec.push_back(tmppair.e2);
    }
  }

  // put the indels into a sorted vector
  for (auto& i : v.indels) {
    tmpvec.push_back(i.second);
  }

  // sort the temp entry vec
  sort(tmpvec.begin(), tmpvec.end());  
  
  // print out the entries
  for (VCFEntryVec::const_iterator it = tmpvec.begin(); it != tmpvec.end(); it++) {
    if (it->filter == "PASS" || vopt::include_nonpass)
      out << *it << endl;
  }

  return out;
}


// sort the VCFEntry by genomic position
bool VCFEntry::operator<(const VCFEntry &v) const {
  return chr < v.chr || (chr == v.chr && pos < v.pos);
}

// merge two VCF headers. Right now, h1 takes some priority (in description field)
VCFHeader mergeVCFHeaders(VCFHeader const &h1, VCFHeader const &h2) {

  if (vopt::verbose > 0)
    cout << "...merging VCF headers" << endl;

  VCFHeader output;

  assert(h1.fileformat == h2.fileformat);
  assert(h1.reference == h2.reference);
  //assert(h1.colnames == h2.colnames);

  output.fileformat = h1.fileformat;
  output.reference = h1.reference;
  output.colnames = (h1.source != "snowmanSV") ? h1.colnames : h2.colnames; // take dRanger col
  
  output.filedate = (h1.filedate > h2.filedate) ? h1.filedate : h2.filedate;
  output.source = h1.source + "_" + h2.source + "__merged";
  
  output.infomap =   mergeHeaderMaps<InfoMap>(h1.infomap, h2.infomap);
  output.filtermap = mergeHeaderMaps<FilterMap>(h1.filtermap, h2.filtermap);
  output.formatmap = mergeHeaderMaps<FormatMap>(h1.formatmap, h2.formatmap);
  output.addInfoField("CALLER","1","String","SnowmanSV/dRanger agree and SnowmanSV shown (SD) or dRanger shown (DS). Private SnowmanSV (S) or dRanger (D)");

  // keep only dRanger sample
  output.samplemap = (h1.source != "snowmanSV") ? h1.samplemap : h2.samplemap;

  return output;

}

// merge two VCF files. 
VCFFile mergeVCFFiles(VCFFile const &v1, VCFFile const &v2) {

  if (vopt::verbose > 0)
    cout << "...merging VCF files" << endl;

  // merge the VCFHeaders
  VCFHeader header = mergeVCFHeaders(v1.sv_header, v2.sv_header);

  if (vopt::verbose > 0)
    cout << "...merging VCF entries" << endl;

  // set the final map
  VCFEntryPairMap final_map;

  // make the combined VCFEntryPair
  VCFEntryPairMap combined_map = v1.entry_pairs;
  combined_map.insert(v2.entry_pairs.begin(), v2.entry_pairs.end());

  string ospan;
  string dspan;
  string sspan;
  int ocount = 0, total = 0, sprivate = 0, dprivate = 0, dpass = 0, spass = 0, dpass_ovl = 0, spass_ovl = 0;
  for (VCFEntryPairMap::iterator it = combined_map.begin(); it != combined_map.end(); it++) {
    bool dranger_pass = it->second.method == "dranger" && it->second.e1.filter == "PASS";
    bool snowman_pass = it->second.method == "snowman" && it->second.e1.filter == "PASS";
    dpass += dranger_pass ? 1 : 0;
    spass += snowman_pass ? 1 : 0;
    for (VCFEntryPairMap::iterator jt = combined_map.begin(); jt != combined_map.end(); jt++) {
      if (it->second.hasOverlap()) // it overlap already found
	break;
      if (it->second.method != jt->second.method && !jt->second.hasOverlap()) { // not part of same VCF
	if (it->second.getOverlaps(vopt::pad, jt->second)) {
	  VCFEntryPair merged(it->second, jt->second);
	  final_map.insert(pair<string,VCFEntryPair>(it->first + "_" + jt->first, merged));
	  jt->second.overlap_partner = it->second.idcommon;
	  it->second.overlap_partner = jt->second.idcommon;
	  assert(jt->second.hasOverlap());
	  assert(it->second.hasOverlap());
	  dpass_ovl += ( (it->second.e1.filter == "PASS" && it->second.method == "dranger") || (jt->second.e1.filter == "PASS" && jt->second.method == "dranger") ) ? 1 : 0;
	  spass_ovl += ( (it->second.e1.filter == "PASS" && it->second.method == "snowman") || (jt->second.e1.filter == "PASS" && jt->second.method == "snowman") ) ? 1 : 0;
	  ospan = ospan + it->second.e1.info_fields["SPAN"] + ";";
	  break;
	}
      }
    }
  }

  ocount = final_map.size();

  // store the private ones
  for (VCFEntryPairMap::iterator it = combined_map.begin(); it != combined_map.end(); it++) {
    if (!it->second.hasOverlap()) {
      bool scaller = it->second.method == "snowman";
      sprivate += scaller ? 1 : 0;
      dprivate += scaller ? 0 : 1;
      dspan = dspan + (!scaller ? (it->second.e1.info_fields["SPAN"] + ";") : "");
      sspan = sspan + (scaller  ? (it->second.e1.info_fields["SPAN"] + ";") : "");
      it->second.addCommonInfoTag("CALLER", scaller ? "S" : "D");
      final_map.insert(pair<string, VCFEntryPair>(it->first, it->second));
    }
  }

  VCFFile out = v1;
  out.sv_header = header;
  out.entry_pairs = final_map;
  //out.filename = vopt::outvcf;

  total = final_map.size();
  int dperc = SnowTools::percentCalc<int>(dprivate, total);
  int sperc = SnowTools::percentCalc<int>(sprivate, total);
  int operc = SnowTools::percentCalc<int>(ocount, total);

  // print some statistics
  if (vopt::verbose > 0) {

    char buffer[100];
  
    sprintf(buffer, " [%5d] - [%4d](%2d%%) : [%4d](%2d%%) - [%4d](%d%%)", 
	    total, ocount, operc, sprivate, sperc, dprivate, dperc);
    cout << "     dRanger total PASS: " << dpass << endl;
    cout << "     Snowman total PASS: " << spass << endl;
    cout << " [TOTAL] - [OVLP](xx%) : [SPRI](xx%) - [DPRI](xx%)" << endl;
    cout << buffer << endl;
  }

  ospan = SnowTools::cutLastChar(ospan);
  sspan = SnowTools::cutLastChar(sspan);
  dspan = SnowTools::cutLastChar(dspan);

  // write the stats file
  ofstream out_stat;
  out_stat.open(vopt::stats);
  out_stat << "Total,Overlaps,SnowPrivate,DranPrivate,OverlapsSpan,SnowPrivateSpans,DranPrivateSpans" << endl;
  out_stat << total << "," << ocount << "," << sprivate << "," << dprivate << "," << ospan << "," << sspan << "," << dspan << endl;
  
  return out;
}
		
// write the VCF to the output file								  
bool VCFFile::write(string basename) const {

  if (vopt::verbose > 0)
    cout << "...writing final VCFs: " << endl; // << vopt::outvcf << endl;

  //ogzstream out(vopt::outvcf.c_str());
  /*ofstream out;
  out.open(vopt::outvcf);
  out << *this;
  */
  // write the indel one
  if (vopt::dranger == "") //basename.find("merged") == string::npos)
    writeIndels(basename, vopt::zip);
  writeSVs(basename, vopt::zip);
  
  return true;

}

// create a VCFFile from a SNOWMANCSV
VCFFile::VCFFile(string file, const char* index, char sep, string analysis_id) {

  //open the file
  igzstream infile(file.c_str(), ios::in);
  
  // confirm that it is open
  if (!infile) {
    cerr << "Can't read file " << file << " for parsing VCF" << endl;
    exit(EXIT_FAILURE);
  }

  // read the reference if not open
  if (vopt::verbose)
    cout << "...vcf - reading in the reference" << endl;
  if (!findex)
    findex = fai_load(index);  // load the reference

  // read in the header of the csv
  string line;
  /*if (!getline(infile, line, sep)) {
    cerr << "CSV file " << file << " is empty" << endl;
    exit(EXIT_FAILURE);
    }*/

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
  indel_header.addInfoField("BLACKLIST","1","String","Normal allelic fraction at break. -1 for undefined");

  // keep track of exact positions to keep from duplicating
  unordered_map<string,bool> double_check;

  // read the reference if not open
  if (vopt::verbose)
    cout << "...vcf - reading in the breakpoints file" << endl;

  // read it in line by line
  getline(infile, line, '\n'); // skip first line
  size_t line_counter = 0;
  while (getline(infile, line, '\n')) {
    line_counter++;
    istringstream iss(line);

    VCFEntry vcf1;
    VCFEntry vcf2;
    vcf1.idcommon = to_string(line_counter) + ":" + vopt::analysis_id;
    vcf2.idcommon = to_string(line_counter) + ":" + vopt::analysis_id;

    string strand1;
    string strand2;
    
    unordered_map<string, string> info_fields;

    string readid;

    size_t val_counter = 0;
    string val;
    while (getline(iss, val, sep)) {
      val_counter++;
      //cout << val_counter << " " << val << endl;
      switch (val_counter) {
      case 1 : vcf1.chr = SnowTools::GenomicRegion::chrToNumber(val); break;
      case 2 : vcf1.pos = stoi(val); break;
      case 3 : strand1 = val; break;
      case 4 : vcf2.chr = SnowTools::GenomicRegion::chrToNumber(val); break;
      case 5 : vcf2.pos = stoi(val); break;
      case 6 : strand2 = val; break;
      case 7: info_fields["SPAN"] = val; break;
      case 8 : vcf1.info_fields["MAPQ"] = val; vcf2.info_fields["MATEMAPQ"] = val; break;
      case 9 : vcf2.info_fields["MAPQ"] = val; vcf1.info_fields["MATEMAPQ"] = val; break;
      case 10: info_fields["NSPLIT"] = val; break;
      case 11: info_fields["TSPLIT"] = val; break;
      case 12: info_fields["NDISC"]  = val; break;
      case 13: info_fields["TDISC"]  = val; break;
      case 14: info_fields["NCIGAR"]  = val; break;
      case 15: info_fields["TCIGAR"]  = val; break;
      case 16: info_fields["HOMSEQ"] = val; break;
      case 17: info_fields["INSERTION"] = val; break;
      case 18: info_fields["SCTG"] = val; break;
      case 19: info_fields["NUMPARTS"] = val; break;
      case 20: vcf1.filter = val; vcf2.filter = val; break; // CONFIDENCE
      case 21: info_fields["EVDNC"] = val; if (val=="DSCRD") info_fields["IMPRECISE"] = ""; break;
      case 22: info_fields["PONCOUNT"] = val; break;
      case 23: info_fields["REPSEQ"] = val; break;
      case 24: info_fields["NCOV"] = val; break;
      case 25: info_fields["TCOV"] = val; break;
      case 26: info_fields["NFRAC"] = val; break;
      case 27: info_fields["TFRAC"] = val; break;
      case 28: info_fields["BLACKLIST"] = val; break;
      case 29: info_fields["DBSNP"] = val; break;
      case 30: readid = val; break;
      }
    }

    // fix the x
    if (info_fields["HOMSEQ"] == "x")
      info_fields["HOMSEQ"] = "";
    if (info_fields["INSERTION"] == "x")
      info_fields["INSERTION"] = "";
    if (info_fields["REPSEQ"] == "x")
      info_fields["REPSEQ"] = "";

    //string nalt = to_string(stoi(info_fields["TSPLIT"]) + stoi(info_fields["TDISC"]));
    string nalt_rp = info_fields["TDISC"];
    string nalt_sp = info_fields["TSPLIT"];

    // remove empty fields
    if (info_fields["DBSNP"] == "x")
      info_fields.erase("DBSNP");

    // treak indels separatley
    if (info_fields["EVDNC"] != "INDEL" && (vcf1.filter == "PASS" || vopt::include_nonpass)) {

      // remove some fields not relevant to SVs
      info_fields.erase("TCOV");
      info_fields.erase("NCOV");
      info_fields.erase("TFRAC");
      info_fields.erase("NFRAC");
      info_fields.erase("TCIGAR");
      info_fields.erase("NCIGAR");
      info_fields.erase("PONCOUNT");

      vcf1.info_fields["MATEID"] = vcf2.id;
      vcf2.info_fields["MATEID"] = vcf1.id;
      
      vcf1.info_fields["SVTYPE"] = "BND";
      vcf2.info_fields["SVTYPE"] = "BND";
      
      vcf1.id = vcf1.idcommon + ":1";
      vcf2.id = vcf2.idcommon + ":2";

      vcf1.info_fields["MATEID"] = vcf2.id;
      vcf2.info_fields["MATEID"] = vcf1.id;
      
      // add the info_fields
      vcf1.info_fields.insert(info_fields.begin(), info_fields.end());
      vcf2.info_fields.insert(info_fields.begin(), info_fields.end());
      
      vcf1.format = "NALT:NALT_RP:NALT_SR:READ_ID";
      vcf2.format = "NALT:NALT_RP:NALT_SR:READ_ID";
      
      // parse the readids
      SupportingReadsMap suppr;
      istringstream iss_r(readid);
      string thisread;
      string new_readid_t = "";
      string new_readid_n = "";
      
      // regex to clean out the t, n identifer
      regex regc("[a-z][0-9]+_(.*?)$");
      smatch smatchr;

      size_t numt = 0;
      size_t numn = 0;
      unordered_map<string, bool> tdup, ndup;
      
      if (readid != "x")
      while (getline(iss_r, thisread, ',')) {

	string thisread_clean; // get only the read name
	if (!regex_search(thisread, smatchr, regc))
	  cerr << "FAILED TO MATCH ON "<< thisread << endl;
	else 
	  thisread_clean = smatchr[1].str();
	//	cout << "readid : " << thisread_clean << endl;
	suppr.insert(pair<string,bool>(thisread_clean, false));

	if (thisread.at(0) == 't' && tdup.count(thisread_clean) == 0) {
	  numt++;
	  new_readid_t += thisread_clean + ",";
	  tdup[thisread_clean]= true;
	} else if (thisread.at(0) == 'n' && ndup.count(thisread_clean) == 0) {
	  numn++;
	  new_readid_n += thisread_clean + ",";
	  ndup[thisread_clean]= true;
	}
      }
      
      // remove the last comma
      if (new_readid_t.length() > 0)
	new_readid_t = new_readid_t.substr(0, new_readid_t.length() - 1);
      if (new_readid_n.length() > 0)
	new_readid_n = new_readid_n.substr(0, new_readid_n.length() - 1);	
      
      assert(readid == "x" || (readid != "x" && suppr.size() > 0) );
      //assert(new_readid_t.length() == readid.length());
      
      
      // make reads a . if empty
      new_readid_t = (new_readid_t == "") ? "." : new_readid_t;
      new_readid_n = (new_readid_n == "") ? "." : new_readid_n;
      
      // ok, so we cant have any colons in the read name, so swtich : for -
      std::replace( new_readid_n.begin(), new_readid_n.end(), ':', '-'); // replace all 'x' to 'y'
      std::replace( new_readid_t.begin(), new_readid_t.end(), ':', '-'); // replace all 'x' to 'y'
     
      string t_sample_string = to_string(numt) + ":" + nalt_rp + ":" + nalt_sp + ":" + new_readid_t; 
      string n_sample_string = to_string(numn) + ":" + info_fields["NDISC"] + ":" + info_fields["NSPLIT"] + ":" + new_readid_n; 
  
      //debug
      vcf1.samp2 = t_sample_string;
      vcf2.samp2 = t_sample_string;
      vcf1.samp1 = n_sample_string;
      vcf2.samp1 = n_sample_string;
      // set the tumor NALT
      //vcf1.samp2 = to_string(numt) + ":" + nalt_rp + ":" + nalt_sp + ":" + new_readid_t; 
      //vcf2.samp2 = to_string(numt) + ":" + nalt_rp + ":" + nalt_sp + ":" + new_readid_t; 

      // set the normal NALT
      //vcf1.samp1 = to_string(numn) + ":" + info_fields["NDISC"] + ":" + info_fields["NSPLIT"] + ":" + new_readid_n; 
      //vcf2.samp1 = to_string(numn) + ":" + info_fields["NDISC"] + ":" + info_fields["NSPLIT"] + ":" + new_readid_n; 
      
      // grab the reference sequence
      SnowTools::GenomicRegion gr1(vcf1.chr, vcf1.pos, vcf1.pos);
      SnowTools::GenomicRegion gr2(vcf2.chr, vcf2.pos, vcf2.pos);
      if (vcf1.chr < 24 && vcf2.chr < 24) { // TODO FIX 
	vcf1.ref = getRefSequence(gr1, findex);
	vcf2.ref = getRefSequence(gr2, findex);
      }
      
      // set the reference position for making ALT tag
      stringstream ptag1, ptag2;
      ptag1 << SnowTools::GenomicRegion::chrToString(vcf2.chr) << ":" << vcf2.pos;
      ptag2 << SnowTools::GenomicRegion::chrToString(vcf1.chr) << ":" << vcf1.pos;
      
      //
      string insertion = vcf1.info_fields["INSERTION"];
      string ttag1, ttag2;
      // TODO insertion + vcf1.ref not right
      ttag1 = (insertion != "" && strand1 == "+" && false) ? (insertion + vcf1.ref) : vcf1.ref;
      ttag1 = (insertion != "" && strand1 == "-" && false) ? (vcf1.ref + insertion) : vcf1.ref;
      ttag2 = (insertion != "" && strand2 == "+" && false) ? (insertion + vcf2.ref) : vcf2.ref;
      ttag2 = (insertion != "" && strand2 == "-" && false) ? (vcf2.ref + insertion) : vcf2.ref;
      
      // set the alternate
      stringstream alt1, alt2;
      if (strand1 == "+" && strand2 == "+") {
	alt1 << ttag1 << "]" << ptag1.str() << "]";
	alt2 << ttag2 << "]" << ptag2.str() << "]";
      } else if (strand1 =="+" && strand2 == "-") {
	alt1 << ttag1 << "[" << ptag1.str() << "[";
	alt2 << "]" << ptag2.str() << "]" << ttag2;
      } else if (strand1 == "-" && strand2 == "+") {
	alt1 << "]" << ptag1.str() << "]" << ttag1;
	alt2 << ttag2 << "[" << ptag2.str() << "[";
      } else {
	alt1 << "[" << ptag1.str() << "[" << ttag1;      
	alt2 << "[" << ptag2.str() << "[" << ttag2;      
      }
      
      vcf1.alt = alt1.str();
      vcf2.alt = alt2.str();
      
      /*
      // do additional filtering
      int tsplit = stoi(vcf1.info_fields["TSPLIT"]);
      int tdisc = stoi(vcf1.info_fields["TDISC"]);
      int mapq1 = stoi(vcf1.info_fields["MAPQ"]);
      int mapq2 = stoi(vcf1.info_fields["MATEMAPQ"]);
      //bool asdisc = vcf1.info_fields["EVDNC"] == "ASDIS";
      int span = stoi(vcf1.info_fields["SPAN"]);
      
      if (span < 500 && span > 0) { // weird bug fix
      vcf1.info_fields["TDISC"] = "0";
      vcf2.info_fields["TDISC"] = "0";
      tdisc = 0;
      vcf1.info_fields["EVDNC"] == "ASSMB";
      vcf2.info_fields["EVDNC"] == "ASSMB";
      }
      
      if (tsplit <= 4 && vcf1.info_fields["EVDNC"] == "ASSMB") {
      vcf1.filter = "WEAKASSEMBLY";
      vcf2.filter = "WEAKASSEMBLY";
      } else if ((vcf1.info_fields["EVDNC"] == "ASSMB") && (mapq1 != 60 || mapq2 != 60) || 
      (span < 1000 && vcf1.info_fields["NUMPARTS"] == "2" && (vcf1.info_fields["HOMSEQ"].length() > 8 || vcf1.info_fields["INSERTION"].length() > 8) ) ) {
      vcf1.filter = "LOWMAPQ";
      vcf2.filter = "LOWMAPQ";
      } else if ((vcf1.info_fields["EVDNC"] == "DSCRD") && ((mapq1 <= 40 && tdisc < 8) || (mapq1 <= 30 && tdisc < 12)) ) {
      vcf1.filter = "WEAKDISC";
      vcf2.filter = "WEAKDISC";
      } else if ( vcf1.info_fields["EVDNC"] == "ASDIS" && tdisc < 3 && tsplit < 3 ) {
      vcf1.filter = "WEAKASSEMBLY";
      vcf2.filter = "WEAKASSEMBLY";
      }
      
      
      // reject all "DSCRD"
      if (vcf1.info_fields["EVDNC"] == "DSCRD") {
      vcf1.filter = "NOASSMB";
      vcf2.filter = "NOASSMB";
      }
      
      // reject anything with too big of homoseq
      if (vcf1.info_fields["HOMSEQ"].length() > 9) {
      vcf1.filter = "LOWMAPQ";
      vcf2.filter = "LOWMAPQ";
      }
      */

      VCFEntryPair vpair;
      if (vcf1 < vcf2)
	vpair = VCFEntryPair(vcf1, vcf2);
      else
	vpair = VCFEntryPair(vcf2, vcf1);
      
      //vpair.supp_reads = ReadIDToReads(readid);
      vpair.samples = {n_sample_string, t_sample_string};



      // add the supporting reads
      
      // double check that it's not doubled. Bug with doubled events
      /*string dc = to_string(vcf1.chr)+"_"+to_string(vcf1.pos)+"_"+strand1+
	to_string(vcf2.chr)+"_"+to_string(vcf2.pos)+"_"+strand2;
      if (double_check.count(dc) > 0) { // already found, don't continue with this entry
	//cerr << "found double with " << dc << endl;
      } 
      else { */
      if (vcf1.chr < 24 && vcf2.chr < 24)
	entry_pairs.insert(pair<string, VCFEntryPair>(vcf1.idcommon, vpair));
	//	double_check[dc] = true;
	//}
      
    }  else if (info_fields["EVDNC"] == "INDEL"  && (vcf1.filter == "PASS" || vopt::include_nonpass)) { // END THE INDEL CONDITIONAL

      // set the ID
      vcf1.id = vcf1.idcommon;
      
      // fix the info field
      vcf1.info_fields.insert(info_fields.begin(), info_fields.end());
      for(auto it = begin(vcf1.info_fields); it != end(vcf1.info_fields);) {
	if (it->first == "MATEMAPQ" || it->first == "NUMPARTS")
	  it = vcf1.info_fields.erase(it); // previously this was something like info_fields.erase(it++);
	else
	  ++it;
      }

      // get the reference for insertions
      string ins = info_fields["INSERTION"];
      //cout << ins << endl;
      if (ins != "") {

	// get the reference
	//vcf1.pos--; 
	SnowTools::GenomicRegion ref(vcf1.chr, vcf1.pos, vcf2.pos-1);
	if (vcf1.chr < 24) {
	  vcf1.ref = getRefSequence(ref, findex);
	}
	
	vcf1.alt = vcf1.ref.substr(0,1) + ins; //.substr(1,ins.length());
	//vcf1.alt = vcf1.ref.substr(0,2) + ins; //.substr(1,ins.length());

      } else { // get the reference for deletions

	// fix the span
	//vcf1.info_fields["SPAN"] = to_string(abs(vcf1.pos - vcf2.pos) - 1);

	// get the reference
	SnowTools::GenomicRegion ref(vcf1.chr, vcf1.pos, vcf2.pos-1); //-1 so we don't include the final unaltered base
	if (vcf1.chr < 24) 
	  vcf1.ref = getRefSequence(ref, findex);

	vcf1.alt = vcf1.ref.substr(0,1);
	
      }
      
      stringstream id_string;
      id_string << vcf1.chr << "_" << vcf1.pos;

      // this indel is already present
      if (indels.count(id_string.str())) {

	VCFEntry orig = indels[id_string.str()];

	// new one passes, old one doesnt, so replace
	if (orig.filter != "PASS" && vcf1.filter == "PASS") {
	  indels[id_string.str()] = vcf1;
	// else if both good, take one with better normal, then better tumor
	} else if (orig.filter == "PASS" && vcf1.filter == "PASS") {
	  int v_n = stoi(vcf1.info_fields["NCIGAR"]) + stoi(vcf1.info_fields["NSPLIT"]);
	  int o_n = stoi(orig.info_fields["NCIGAR"]) + stoi(orig.info_fields["NSPLIT"]);
	  int v_t = stoi(vcf1.info_fields["TCIGAR"]) + stoi(vcf1.info_fields["TSPLIT"]);
	  int o_t = stoi(orig.info_fields["TCIGAR"]) + stoi(orig.info_fields["TSPLIT"]);
	  // new one has more normal, or same number but more tumor
	  if ( (v_n > o_n) || (v_n == o_n && v_t > o_t))
	    indels[id_string.str()] = vcf1;
	}
      } else {
	indels[id_string.str()] = vcf1;
      }
      
    }
  } // END THE WHILE LOOP

  // dedupe
  deduplicate();

  
}

// return a sequence from the reference
string getRefSequence(const SnowTools::GenomicRegion &gr, faidx_t * fi) {

  int len;
  string chrstring = SnowTools::GenomicRegion::chrToString(gr.chr);
  //char * seq = faidx_fetch_seq(fi, const_cast<char*>(chrstring.c_str()), gr.pos1-1, gr.pos2-1, &len);
  char * seq = faidx_fetch_seq(fi, const_cast<char*>(chrstring.c_str()), gr.pos1, gr.pos2, &len);
  
  if (seq) {
    return string(seq);
  } else {
    //cout << "Failed to get reference sequence at " << gr << endl;
    return "LOAD_FAIL";
  }

}

// parse the READ_ID tag and store as separate reads
SupportingReadsMap ReadIDToReads(string readid) {
  
  SupportingReadsMap supp_reads;

  if (readid != "") {
    // parse the supporting reads
    istringstream issr(readid);
    string val;
    while (getline(issr, val, ';')) 
      supp_reads.insert(pair<string, bool>(val,false));
  } 
  return supp_reads;

}

// merge string maps from the header (e.g. INFO, FORMAT, SAMPLE)
template<typename T> T mergeHeaderMaps(T const &m1, T const &m2) {

  T shared_fields;

  // merge the info maps
  for (typename T::const_iterator it = m1.begin(); it != m1.end(); it++) {
    if (m2.count(it->first) > 0)
      shared_fields.insert(pair<string,string>(it->first, it->second));
  }
  for (typename T::const_iterator it = m2.begin(); it != m2.end(); it++) {
    if (m1.count(it->first) > 0 && shared_fields.count(it->first) == 0)
      shared_fields.insert(pair<string,string>(it->first, it->second));
  }  
  
  T final_map;
  // add the private infos
  for (typename T::const_iterator it = m1.begin(); it != m1.end(); it++) {
    if (shared_fields.count(it->first) == 0)
      final_map.insert(pair<string,string>(it->first, it->second));
  }

  // add the private infos 1
  for (typename T::const_iterator it = m1.begin(); it != m1.end(); it++) {
    if (shared_fields.count(it->first) == 0)
      final_map.insert(pair<string,string>(it->first, it->second));
  }
  // add the private infos 2
  for (typename T::const_iterator it = m2.begin(); it != m2.end(); it++) {
    if (shared_fields.count(it->first) == 0)
      final_map.insert(pair<string,string>(it->first, it->second));
  }
  // add the shared
  for (typename T::const_iterator it = shared_fields.begin(); it != shared_fields.end(); it++) {
    final_map.insert(pair<string,string>(it->first, it->second));
  }

  return final_map;

}

// find if one VCFEntryPair (breakpoint pair) overlaps with another
bool VCFEntryPair::getOverlaps(int pad, VCFEntryPair &v) {

  SnowTools::GenomicRegion gr1(e1.chr, e1.pos, e1.pos);
  gr1.pad(pad);
  SnowTools::GenomicRegion gr2(e2.chr, e2.pos, e2.pos);
  gr2.pad(pad);

  SnowTools::GenomicRegion t1(v.e1.chr, v.e1.pos, v.e1.pos);
  t1.pad(pad);
  SnowTools::GenomicRegion t2(v.e2.chr, v.e2.pos, v.e2.pos);
  t2.pad(pad);
  
  if (!(t1 <= t2 && gr1 <= gr2)) {
    cerr << "T1 " << t1 << " T2 " << t2 << " GR1 " << gr1 << " GR2 " << gr2 << endl;
    assert(t1 <= t2 && gr1 <= gr2);
  }

  if (gr1.getOverlap(t1) > 0 && gr2.getOverlap(t2) > 0) {
    return true;
  }

  return false;
}

// to each of the two VCFEntry in a pair, add the info tag
void VCFEntryPair::addCommonInfoTag(string tag, string value) {
  e1.info_fields[tag] = value;
  e2.info_fields[tag] = value;
  return;
}

// make a new VCF entry pair
VCFEntryPair::VCFEntryPair(VCFEntryPair &v1, VCFEntryPair &v2) {

  VCFEntryPair dran = (v1.method == "dranger") ? v1 : v2;
  VCFEntryPair snow = (v1.method == "dranger") ? v2 : v1;

  // merge the info fields
  //InfoMap merged_fields = mergeInfoFields(dran.e1.info_fields, snow.e1.info_fields);

  // dranger wins if snow+dran && BPRESLULT
  if (dran.tsplit > 0) {

    method = "dranger";
    idcommon = dran.idcommon;
    string overlap_partner = snow.idcommon;
    dran.addCommonInfoTag("CALLER", "DS");
    e1 = dran.e1;
    e2 = dran.e2;

  } else if (dran.tsplit == 0 && snow.tsplit > 0) {

    method = "snowman";
    idcommon = snow.idcommon;
    string overlap_partner = dran.idcommon;
    snow.addCommonInfoTag("CALLER", "SD");
    e1 = snow.e1;
    e2 = snow.e2;
  }

  // merge the supporting reads
  supp_reads = dran.supp_reads;
  for (SupportingReadsMap::const_iterator it = snow.supp_reads.begin(); it != snow.supp_reads.end(); it++) {
    if (supp_reads.count(it->first) == 0) // read is new
      supp_reads.insert(pair<string, bool>(it->first, true));
  }
}


// merge the info fields
/*InfoMap mergeInfoFields(InfoMap const &m1, InfoMap const &m2) {

  InfoMap shared_fields;

  // merge the info maps
  for (InfoMap::const_iterator it = m1.begin(); it != m1.end(); it++) {
    if (m2.count(it->first) > 0) {
      shared_fields.insert(pair<string,string>(it->first, it->second));
    }
  }
  for (InfoMap::const_iterator it = m2.begin(); it != m2.end(); it++) {
    if (m1.count(it->first) > 0 && shared_fields.count(it->first) == 0)
      shared_fields.insert(pair<string,string>(it->first, it->second));
  }  
    
  return shared_fields;
  }*/

// write the VCF as a csv file
bool VCFFile::writeCSV() const {

  if (vopt::verbose > 0)
    cout << "...writing final VCF as CSV: " << vopt::as_csv << endl;

  ofstream out; 
  out.open(vopt::as_csv);

  // find out what are all the info fields to print
  for (VCFEntryPairMap::const_iterator it = entry_pairs.begin(); it != entry_pairs.end(); it++) 
    for (InfoMap::const_iterator jt = it->second.e1.info_fields.begin(); jt != it->second.e1.info_fields.end(); jt++)
      union_info_fields[jt->first] = jt->second;
  
  // print out the header
  out << "chr1,pos1,strand1,chr2,pos2,strand2";
  for (InfoMap::const_iterator it = union_info_fields.begin(); it != union_info_fields.end(); it++) 
    out << "," << it->first;
  out << endl;
  
  // print out the entries
  for (VCFEntryPairMap::const_iterator it = entry_pairs.begin(); it != entry_pairs.end(); it++) 
    out << it->second.toCSVString() << endl;
  
  return true;

}

// 
string VCFEntryPair::toCSVString() const {
  
  string sep = ",";

  stringstream ss;

  regex regPP(".*?\\].*?\\]$");
  regex regPM(".*?\\[.*?\\[$");
  regex regMP("^\\].*?\\].*");
  regex regMM("^\\[.*?\\[.*");

  bool pp = regex_match(e1.alt, regPP);
  bool pm = regex_match(e1.alt, regPM);
  bool mp = regex_match(e1.alt, regMP);
  bool mm = regex_match(e1.alt, regMM);

  // for sanity
  bool pp2 = regex_match(e2.alt, regPP);
  bool pm2 = regex_match(e2.alt, regPM);
  bool mp2 = regex_match(e2.alt, regMP);
  bool mm2 = regex_match(e2.alt, regMM);
  assert(pp == pp2);
  assert(mm == mm2);
  assert(mp == pm2);
  assert(pm == mp2);
  assert(pm || pp || mp || mm);

  string strand1 = (pp || pm) ? "+" : "-";
  string strand2 = (pp || mp) ? "+" : "-";

  string readid;
  for (auto it = supp_reads.begin(); it != supp_reads.end(); it++) {
    readid = it->first + ";";
  }
  if (readid.length() > 0)
    readid = readid.substr(0, readid.length() - 1);
  
  ss << SnowTools::GenomicRegion::chrToString(e1.chr) << sep << e1.pos << sep << strand1 << sep 
     << SnowTools::GenomicRegion::chrToString(e2.chr) << sep << e2.pos << sep << strand2;
  for (InfoMap::const_iterator it = union_info_fields.begin(); it != union_info_fields.end(); it++) {
    InfoMap::const_iterator ff = e1.info_fields.find(it->first);
    ss << sep;
    if (ff != e1.info_fields.end()) {
      if (ff->first != "READ_ID")
	ss << ff->second;
      else 
	ss << readid;
    }
    
  }

  return ss.str();
}


// turn a series of FORMAT string (identifier, samp1, samp2) into a set of string records
FormatRecordMap FormatStringToFormatRecordMap(string format, string samp1, string samp2) {

  string val;

  FormatRecordMap format_fields;

  // parse the FORMAT field
  istringstream issf(format);
  vector<string> format_names;
  vector<FormatPair> format_vals;
  while (getline(issf, val, ':')) {
    format_names.push_back(val);
  }
  reverse(format_names.begin(), format_names.end());
  assert(format_names.size() == 4);
  // add the FORMAT data for samp1
  istringstream iss1(samp1);
  while (getline(iss1, val, ':')) {
    FormatPair form;
    form.first = val;
    format_vals.push_back(form);
  }
  assert(format_vals.size() == 4);
  reverse(format_vals.begin(), format_vals.end());
  // add the FORMAT data for samp2
  istringstream iss2(samp2);
  while (getline(iss2, val, ':')) {
    FormatPair form = format_vals.back();
    format_vals.pop_back();
    form.second = val;
    string field_name = format_names.back();
    format_names.pop_back();
    format_fields.insert(pair<string,FormatPair>(field_name, form));
  }

  return format_fields;

}

// deduplicate
void VCFFile::deduplicate() {

  if (vopt::verbose)
    cout << "deduping events" << endl;

  SnowTools::GenomicRegionCollection<SnowTools::GenomicRegion> grv;

  // create the interval tree
  for (auto i : entry_pairs) { // todo add 
    grv.add(SnowTools::GenomicRegion(i.second.e1.chr, i.second.e1.pos, i.second.e1.pos));
    grv.add(SnowTools::GenomicRegion(i.second.e2.chr, i.second.e2.pos, i.second.e2.pos));
  }
  grv.createTreeMap();

  // interval tree
  SnowTools::GenomicIntervalTreeMap tree = grv.m_tree;

  size_t count = 0;

  for (auto i : entry_pairs) {
    count++;
    SnowTools::GenomicIntervalVector giv;
    tree[i.second.e1.chr].findContained(i.second.e1.pos-vopt::pad, i.second.e1.pos+vopt::pad, giv);
    tree[i.second.e2.chr].findContained(i.second.e2.pos-vopt::pad, i.second.e2.pos+vopt::pad, giv);

    if (giv.size() >= 3)
      for (auto j : entry_pairs) {
	if (i.first != j.first && dups.count(j.first) == 0) {
	  if (i.second.getOverlaps(vopt::pad, j.second)) 
	    dups[i.first] = true;
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

  return (chr == v.chr && pos == v.pos);

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
  // try BGZIP
  //gname += ".bz";
  //sname += ".bz";
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
    
  //out_g << indel_header << endl;
  //out_s << indel_header << endl;

  VCFEntryVec tmpvec;

  // put the indels into a sorted vector
  for (auto& i : indels) {
    tmpvec.push_back(i.second);
  }

  // sort the temp entry vec
  sort(tmpvec.begin(), tmpvec.end());  

  // print out the entries
  for (auto& i : tmpvec) { 
    if (i.filter == "PASS" || vopt::include_nonpass) {

      // get the counts
      size_t pon = i.info_fields.count("PONCOUNT") ? stoi(i.info_fields["PONCOUNT"]) : 0;
      size_t nsplit = i.info_fields.count("NSPLIT") ? stoi(i.info_fields["NSPLIT"]) : 0;
      size_t ncigar = i.info_fields.count("NCIGAR") ? stoi(i.info_fields["NCIGAR"]) : 0;
      size_t tcigar = i.info_fields.count("TCIGAR") ? stoi(i.info_fields["TCIGAR"]) : 0;
      size_t tsplit = i.info_fields.count("TSPLIT") ? stoi(i.info_fields["TSPLIT"]) : 0;
      //size_t repseq = i.info_fields.count("REPSEQ") ? i.info_fields["REPSEQ"].length() : 0;
      double naf = i.info_fields.count("NFRAC") ? stod(i.info_fields["NFRAC"]) : 0;
      //double taf = i.info_fields.count("TFRAC") ? stod(i.info_fields["TFRAC"]) : 0;

      double somatic_ratio = 100;
      size_t ncount = max(nsplit,ncigar);
      if (ncount > 0)
	somatic_ratio = (max(tsplit,tcigar)) / ncount;

      // if it's DB SNp, it's germline
      if (i.info_fields.count("DBSNP"))
	somatic_ratio = 0;

      if (somatic_ratio >= 20 && ncount < 2 && pon <= 1 && naf < 0.05) { // ok if its just one...
	//out_s << i << endl;
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
	//out_g << i << endl;

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
  string filt = vopt::filtered ? ".filtered" : "";
  gname = basename + "germline.sv" + filt + ".vcf.gz";
  sname = basename + "somatic.sv" + filt + ".vcf.gz";
  gname_nz = basename + "germline.sv" + filt + ".vcf";
  sname_nz = basename + "somatic.sv" + filt + ".vcf";

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
  size_t id_counter = 0;

  // put the pair maps into a vector
  for (VCFEntryPairMap::const_iterator it = entry_pairs.begin(); it != entry_pairs.end(); it++) {

    if (dups.count(it->first) == 0) { // dont include duplicate entries
      // renumber the ids
      id_counter++;
      VCFEntryPair tmppair = it->second;
      if (vopt::dranger != "" && vopt::snowman != "") {
	
	tmppair.e1.idcommon = to_string(id_counter) + ":" + vopt::analysis_id;
	tmppair.e2.idcommon = to_string(id_counter) + ":" + vopt::analysis_id;
	tmppair.e1.id = tmppair.e1.idcommon + ":1";
	tmppair.e2.id = tmppair.e2.idcommon + ":2";
	tmppair.e1.info_fields["MATEID"] = tmppair.e2.id;
	tmppair.e2.info_fields["MATEID"] = tmppair.e1.id;
      } 
      
      tmpvec.push_back(tmppair.e1);
      tmpvec.push_back(tmppair.e2);
    }

  }

  // sort the temp entry vec
  sort(tmpvec.begin(), tmpvec.end());  
  
  // for now all merged calls are somatic
  //bool ismerged = basename.find("merged") != string::npos;
  bool ismerged = vopt::dranger != "";

  // print out the entries
  for (auto& i : tmpvec) { 
    if (i.filter == "PASS" || vopt::include_nonpass) {

      size_t tsplit = i.info_fields.count("TSPLIT") ? stoi(i.info_fields["TSPLIT"]) : 0;
      size_t nsplit = i.info_fields.count("NSPLIT") ? stoi(i.info_fields["NSPLIT"]) : 0;
      size_t tdisc = i.info_fields.count("TDISC") ? stoi(i.info_fields["TDISC"]) : 0;
      size_t ndisc = i.info_fields.count("NDISC") ? stoi(i.info_fields["NDISC"]) : 0;

      double somatic_ratio = 100;
      size_t ncount = max(nsplit, ndisc);
      if (ncount > 0)
	somatic_ratio = (max(tsplit,tdisc)) / ncount;

      if ( (somatic_ratio >= 12 && ncount < 2) || (ismerged)) {
	if (zip) {
	  stringstream ss;
	  ss << i << endl;
	  if (!bgzf_write(s_bg, ss.str().c_str(), ss.str().length())) {
	    cerr << "Could not write bzipped vcf for line " << ss.str() << endl;
	  }
	} else {
	  out_s << i << endl;
	}

	//out_s << i << endl;
      } else if (!ismerged) {
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

  // dont need germline if merged
  if (ismerged) {
    string cmd;
    if (zip) {
      cmd = "rm " + gname;
    } else {
      cmd = "rm " + gname_nz;
    }
    system(cmd.c_str());
  }
  
  // tabix it
  if (zip) {
    if (!ismerged)
      tabixVcf(gname);
    tabixVcf(sname);
  }


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
