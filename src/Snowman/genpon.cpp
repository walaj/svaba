#include "genpon.h"

#include <getopt.h>
#include <string>
#include <sstream>
#include <unordered_map>

#include "SnowTools/gzstream.h"

#include "PonWalker.h"

static pthread_mutex_t gp_lock;

using namespace std;

namespace popt {

  static string input_file = "";
  static string output_file = "";
  static string file_list = "";
  static bool collapse = false;

  static int verbose = 1;
  static int numThreads = 1;
  static std::string regionFile = "";
  static std::string verify = "";
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

static const char* shortopts = "hci:v:o:l:p:k:V:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "input-vcfs",              required_argument, NULL, 'i'},
  { "output-pon",              required_argument, NULL, 'o'},
  { "file-list",               required_argument, NULL, 'l'},
  { "region-file",             required_argument, NULL, 'k' },
  { "num-threads",             required_argument, NULL, 'p'},
  { "verify",                  required_argument, NULL, 'V'},
  { "collapse",                no_argument, NULL, 'c'},
  { "verbose",                 required_argument, NULL, 'v' },
  { NULL, 0, NULL, 0 }
};

static const char *PON_USAGE_MESSAGE =
"Usage: snowman pon [OPTION] -i vcf_files.txt -o pon.txt.gz\n\n"
"  Description: \n"
"\n"
"  Required input\n"
"  -l, --file-list                      File containing list of cigar.gz files for each sample to be panel\n"
"  General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h, --help                           Display this help and exit\n"
"  -o, --output-pon                     Output panel of normals file for use with snowman run -q or snowman refilter -q \n"
"  -c, --collapse                       Collapse all of the counts from individual files to one point. Makes PON file much smaller\n"
"  -p, --num-threads                    Number of threads to use. Runs one BAM per thread\n"
"  -k, --region-file                    Set a region txt file. Format: one region per line, Ex: 1,10000000,11000000\n"
"  -V, --verify                         Read in a binary PON file generated from snowman pon and profile/verify it\n"
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
      case 'c': popt::collapse = true; break;
      case 'i': arg >> popt::input_file; break;
      case 'p': arg >> popt::numThreads; break;
      case 'o': arg >> popt::output_file; break;
      case 'l': arg >> popt::file_list; break;
      case 'k': arg >> popt::regionFile; break;
      case 'v': arg >> popt::verbose; break;
      case 'V': arg >> popt::verify; break;
    }
  }

  if (popt::output_file.length() == 0 && popt::verify.length() == 0)
    die = true;

  if (die) {
      cout << "\n" << PON_USAGE_MESSAGE;
      exit(1);
    }
}

void verifyPON() {
  
  if (!SnowTools::read_access_test(popt::verify)) {
    std::cerr << "!!!!!!!!!!!!!!!!" << std::endl;
    std::cerr << "Can't access PON for verification: " << popt::verify << std::endl;
    std::cerr << "!!!!!!!!!!!!!!!!" << std::endl;
  }

  std::vector<DiscRead> drv;

  std::cerr << "...reading in PON file: " << popt::verify << std::endl;
  std::ifstream ifs;
  ifs.open(popt::verify.c_str(), std::ios::in);

  std::unordered_map<uint16_t, size_t> id_count;

  while(!ifs.eof()) {
    DiscRead dr;
    dr.load(ifs);
    if (dr.pos1 > 0) 
      drv.push_back(dr);
    ++id_count[dr.id];
  }

  for (auto& i : id_count)
    std::cerr << "ID: " << i.first << " count " << i.second << std::endl;

  //debug
  if (popt::output_file.length()) {
    ogzstream outr(popt::output_file.c_str(), std::ios::out);
    if (!outr) {
      std::cerr << "ERROR: Could not open " << popt::output_file << " for writing" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    for (auto& i : drv)
      outr << i.toTabString() << std::endl;
    
    outr.close();
  }


}

void runGeneratePON(int argc, char** argv) {

  parsePONOptions(argc, argv);

  if (popt::verify.length()) {
    verifyPON();
    return;
  }

  cout << "PON file list:   " << popt::file_list << endl;
  cout << "PON output file: " << popt::output_file << endl;
  cout << "Collapse PON:    " << (popt::collapse ? "ON" : "OFF") << endl;

  std::vector<std::string> bams;
  if (SnowTools::read_access_test(popt::file_list)) {
    std::ifstream iss_file(popt::file_list.c_str());
    std::string line;
    while(std::getline(iss_file, line)) {
      if (line.length() && line.find("#") == std::string::npos)
	bams.push_back(line);
    }
  }
  else {
    cerr << "!!!!!!!!!!" << endl << "!!!!!!!!!!" << endl << "!!!!!!!!!!" << endl;
    cerr << "List of BAM files does not exist or is not readable -- " << popt::file_list << endl;
    cerr << "!!!!!!!!!!" << endl << "!!!!!!!!!!" << endl << "!!!!!!!!!!" << endl;
    exit(EXIT_FAILURE);
  }

  // parse the regions
  SnowTools::GRC file_regions;
  if (SnowTools::read_access_test(popt::regionFile))
    file_regions.regionFileToGRV(popt::regionFile, 0);
  else if (popt::regionFile.find(":") != std::string::npos && popt::regionFile.find("-") != std::string::npos
	   && bams.size() && SnowTools::read_access_test(*bams.begin())) {
    SnowTools::BamWalker bwalker(*bams.begin());
    file_regions.add(SnowTools::GenomicRegion(popt::regionFile, bwalker.header()));    
  } 

  popt::numThreads = std::min((int)bams.size(), popt::numThreads);
  
  std::ofstream outr;
  outr.open("dum", std::ios::binary | std::ios::out);

  std::string rl = "global@!hardclip;!supplementary;phred[4,100];%region@WG%discordant[0,800];mapq[1,1000]";
  SnowTools::MiniRulesCollection * mr = new SnowTools::MiniRulesCollection(rl);
  
  // open the mutex
  if (pthread_mutex_init(&gp_lock, NULL) != 0) {
      printf("\n mutex init failed\n");
      return;
  }

  // Create the queue and consumer (worker) threads
  wqueue<PONWorkItem*>  queue;
  std::vector<ConsumerThread<PONWorkItem>*> threadqueue;
  for (int i = 0; i < popt::numThreads; i++) {
    ConsumerThread<PONWorkItem>* threadr = new ConsumerThread<PONWorkItem>(queue, true);
    threadr->start();
    threadqueue.push_back(threadr);
  }
  
  SnowTools::BamWalker out_bw;
  SnowTools::BamWalker bwalker(*bams.begin());
  bam_hdr_t * r2c_hdr = bam_hdr_dup(bwalker.header());
  out_bw.SetWriteHeader(r2c_hdr);
  out_bw.OpenWriteBam(popt::output_file);
  
  // send the jobs
  size_t id = 0;
  for (auto& b : bams) {
    if (!SnowTools::read_access_test(b)) {
      std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      std::cerr << "Cannot read bam " << b << std::endl;\
      std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;      
    }
    PONWorkItem * item = new PONWorkItem(b, &outr, &gp_lock, id++, mr, &file_regions, &out_bw) ;
    queue.add(item);
  }

  // wait for the threads to finish
  for (int i = 0; i < popt::numThreads; i++) 
    threadqueue[i]->join();

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
	try {
	  (*pon)[key].push_back(stoi(num));
	} catch (...) {
	  cerr << "stoi exception on num " << num << " in line " << line << endl;
	  //exit(EXIT_FAILURE);
	}
      }
    }
  }

  // write the pon
  ogzstream oz(popt::output_file.c_str(), ios::out);
  if (!oz) {
    cerr << "Can't write to output file " << popt::output_file << endl;
    exit(EXIT_FAILURE);
  }

  size_t total_reads = 0;
  size_t total_samples = 0;
  for (auto& i : (*pon)) {
    if (i.second.size() > 1) {
      oz << i.first;
      for (auto& j : i.second) {
	if (!popt::collapse) {
	  oz << "\t" << j;
	} else {
	  total_reads += j;
	  total_samples += (j > 0 ? 1 : 0);
	}
      }
      if (popt::collapse)
	oz << "\t" << total_reads << "\t" << total_samples << endl;
      oz << endl;
    }
  }
  oz.close();
  
  delete pon;
}
