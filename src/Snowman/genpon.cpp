#include "genpon.h"
#include <getopt.h>
#include <string>
#include <unordered_map>
#include <memory>
#include "gzstream.h"
#include <sstream>
//#include "reads.h"
#include "SnowUtils.h"

using namespace std;

namespace popt {

  static string input_file = "";
  static string output_file = "";
  static string file_list = "";

  static int verbose = 1;
}


void runGeneratePONfromVCF(int argc, char** argv) {

  parsePONOptions(argc, argv);

  typedef unordered_map<string, size_t> pmap;
  pmap * pon = new pmap();

  // read in the file
  igzstream iz(popt::input_file.c_str());
  if (!iz) {
    cerr << "Can't read file " << popt::input_file << endl;
    exit(EXIT_FAILURE);
  }
  
  // count the files
  string file;
  size_t file_count = 0; 
  size_t fc = 0;
  while (getline(iz, file))
    file_count++;
  // close and reopen
  iz.close();

  igzstream izl(popt::input_file.c_str());

  while (getline(izl, file, '\n')) {
    fc++;
    // open this VCF file
    igzstream zz(file.c_str());
    if (!zz) {
      cerr << "Can't read file " << file << endl;
      exit(EXIT_FAILURE);
    }

    if (popt::verbose) 
      cout << "working on vcf "  << file << " (" << fc << " of " << file_count << ")" << endl;

    // loop through and add 
    string line;
    while (getline(zz, line, '\n')) {
      if (!line.length())
	break;
      size_t c = 0;
      if (line.at(0) != '#') {
	istringstream is(line);
	string val;
	string chr;
	string pos;
	
	while(getline(is, val, '\t')) {
	  switch (++c) {
	  case 1: chr = val; break;
	  case 2: pos = val;
	  }
	  if (c > 1)
	    break;
	}
	
	string key = chr + "_" + pos;
	if (!pon->count(key))
	  (*pon)[key] = 1;
	else
	  (*pon)[key]++;
      }
    }
  }

  // write the pon
  ogzstream oz(popt::output_file.c_str(), ios::out);
  if (!oz) {
    cerr << "Can't write to output file " << popt::output_file << endl;
    exit(EXIT_FAILURE);
  }
  for (auto& i : (*pon)) {
    oz << i.first << "\t" << i.second << endl;
  }
  oz.close();
  
  delete pon;
}

static const char* shortopts = "hi:v:o:l:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "input-vcfs",              required_argument, NULL, 'i'},
  { "output-pon",              required_argument, NULL, 'o'},
  { "file-list",               required_argument, NULL, 'l'},
  { "verbose",                 required_argument, NULL, 'v' },
  { NULL, 0, NULL, 0 }
};

static const char *PON_USAGE_MESSAGE =
"Usage: snowman pon [OPTION] -i vcf_files.txt -o pon.txt.gz\n\n"
"  Description: \n"
"\n"
"  General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h, --help                           Display this help and exit\n"
"  -l, --file-list                      File containing list of cigar.gz files for each sample to be panel\n"
"  -o, --output-pon                     Output panel of normals file for use with snowman run -q or snowman refilter -q \n"
"  Required input\n"
"  Optional input\n"                       
"\n";

void parsePONOptions(int argc, char** argv) {
  bool die = false;

  if (argc <= 2) 
    die = true;

  string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
      case 'h': die = true; break;
      case 'i': arg >> popt::input_file; break;
      case 'o': arg >> popt::output_file; break;
      case 'l': arg >> popt::file_list; break;
      case 'v': arg >> popt::verbose; break;
    }
  }

  //if (popt::input_file.length() == 0)
  //  die = true;
  if (popt::output_file.length() == 0)
    die = true;

  if (die) {
      cout << "\n" << PON_USAGE_MESSAGE;
      exit(1);
    }
}

void runGeneratePON(int argc, char** argv) {

  parsePONOptions(argc, argv);

  bool list_exist = SnowUtils::read_access_test(popt::file_list);
  if (list_exist) {
    if (popt::verbose)
      cout << "...combining PON" << endl;
    combinePON();
  } else if (!list_exist && popt::file_list.length()) {
    cerr << "!!!!!!!!!!" << endl << "!!!!!!!!!!" << endl << "!!!!!!!!!!" << endl;
    cerr << "List of cigar file does not exist or is not readable -- " << popt::file_list << endl;
    cerr << "!!!!!!!!!!" << endl << "!!!!!!!!!!" << endl << "!!!!!!!!!!" << endl;
  }

  return;

  /*
  ogzstream ogz;
  ogz.open(popt::output_file.c_str(), ios::out);

  typedef unordered_map<string, size_t> pmap;
  pmap cigmap;
  
#ifdef HAVE_HTSLIB
  // hts
  BGZF * fp = 0;
  //hts_idx_t * idx = 0; // hts_idx_load(bamfile.c_str(), HTS_FMT_BAI);
  hts_itr_t * hts_itr = 0; // sam_itr_queryi(idx, 3, 60000, 80000);
  bam_hdr_t * br = 0; // header
#endif  


  // open the HTS reader
#ifdef HAVE_HTSLIB
  fp = bgzf_open(popt::input_file.c_str(), "r"); 
  
  if (!fp) {
    cerr << "Error using HTS reader on opening " << popt::input_file << endl;
    exit(EXIT_FAILURE);
  }
  br = bam_hdr_read(fp);
  // open the header with HTS
  if (!br) {
    cerr << "Error using HTS reader on opening " << popt::input_file << endl;
    exit(EXIT_FAILURE);
  }

  //HTS set region
  //if (!idx)
  //  idx = hts_idx_load(bam.c_str(), HTS_FMT_BAI);
  //hts_itr = sam_itr_queryi(idx, interval.chr, interval.pos1, interval.pos2);
  //if (!hts_itr) {
  //  std::cerr << "Error: Failed to set region: " << interval << endl; 
  //  exit(EXIT_FAILURE);
  //}
#endif

#ifdef HAVE_HTSLIB
  void *dum; // needed by hts_itr_next, for some reason...
#endif 

  /// loop through all the reads
  size_t lastpos = -1;
  for (;;) {

    Read r;
#ifdef HAVE_HTSLIB
    bam1_t* b = bam_init1(); 
    if (hts_itr_next(fp, hts_itr, b, dum) <= 0) {
      bam_destroy1(b);
      break; 
    }
    r = std::shared_ptr<bam1_t> (b, free_delete());
#endif

    // immediately add cigar
    if (r_cig_size(r) > 1) {

      // dump the map
      if (r_pos(r) - lastpos > 500) {
	for (auto& i : cigmap) {
	  ogz << i.first << "\t" << i.second << endl;
	}
	cigmap.clear();
      }
      lastpos = r_pos(r);    
      
      stringstream ss; 
      int pos = r_pos(r); 
      
      for (int i = 0; i < r_cig_size(r); i++) {
	if (r_cig_type(r,i) == 'D' || r_cig_type(r,i) == 'I') {	
	  //ss << r_id(r) << "_" << pos << "_" <<  r_cig_len(r,i) <<  r_cig_type(r, i);
	  ss << r_id(r) << "_" << pos << "_" <<  r_cig_type(r, i);
	  cigmap[ss.str()]++;
	  ss.str("");
	  
	}
	if (!(r_cig_type(r, i) == 'I') && !(r_cig_type(r, i) == 'S') && !(r_cig_type(r,i) == 'H'))
	  pos += r_cig_len(r, i);
      }
      
    }

  }
  

  */
}


void combinePON() {

  typedef unordered_map<string, vector<size_t>> pmap;
  pmap * pon = new pmap();

  // read in the file
  igzstream iz(popt::file_list.c_str());
  if (!iz) {
    cerr << "Can't read file " << popt::file_list << endl;
    exit(EXIT_FAILURE);
  }
  
  // count the files
  string file;
  size_t file_count = 0; 
  size_t fc = 0;
  while (getline(iz, file))
    file_count++;
  // close and reopen
  iz.close();

  igzstream izl(popt::file_list.c_str());

  while (getline(izl, file, '\n')) {
    fc++;

    // open this PON file
    igzstream zz(file.c_str());
    if (!zz) {
      cerr << "Can't read file " << file << endl;
      exit(EXIT_FAILURE);
    }

    if (popt::verbose) 
      cout << "working on file "  << file << " (" << fc << " of " << file_count << ")" << endl;

    // loop through and add 
    string line;
    while (getline(zz, line, '\n')) {
      if (!line.length())
	break;
      size_t c = 0;
      if (line.at(0) != '#') {
	istringstream is(line);
	string val;
	string key;
	string num;
	string tum;
	
	while(getline(is, val, '\t')) {
	  switch (++c) {
	  case 1: key = val; break;
	  case 2: num = val;
	  case 3: tum = val; 
	  }
	  if (c > 2)
	    break;
	}
	
	key = tum + key;
	//string key = chr + "_" + pos;
	//if (!pon->count(key))
	  (*pon)[key].push_back(stoi(num));
	  //else
	  //(*pon)[key]++;
      }
    }
  }

  // write the pon
  ogzstream oz(popt::output_file.c_str(), ios::out);
  if (!oz) {
    cerr << "Can't write to output file " << popt::output_file << endl;
    exit(EXIT_FAILURE);
  }
  for (auto& i : (*pon)) {
    oz << i.first;
    for (auto& j : i.second) 
      oz << "\t" << j;
    oz << endl;
  }
  oz.close();
  
  delete pon;
}
