#include "extract.h"

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include "api/BamWriter.h"
#include "api/BamReader.h"
#include <getopt.h>
#include <unordered_map>
#include <unistd.h>
#include <time.h>
#include "seqan_tools.h"
#include "SVBamReader.h"

const int MAX_CHARS_PER_LINE = 512;
const int MAX_TOKENS_PER_LINE = 20;
const char* const DELIMITER = ",";

using namespace std;
using namespace BamTools;

static struct timespec start;

static const char *EXTRACT_USAGE_MESSAGE = 
"Usage: snowtools extract [OPTION] <bam1> <bam2> etc\n\n"
"  Description: Grab reads from each bam within a certain region or with certain properties, write new BAM.\n"
"\n"
"  General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h, --help                           Display this help and exit\n"
  //"  -p, --threads                        Use NUM threads to run snowtools. Default: 1\n"
"  Required input\n"
"  -i, --input-bam                      Input BAM file\n"
"  -o, --output-bam                     Output BAM file\n"
"  Optional input\n"                       
"  -r, --region-file                    Set a region txt file. Format: one region per line, Ex: 1,10000000,11000000\n"
"  -n, --nm-limit                       Only extract reads with n mismatches or greater (uses NM tag). Default: 1\n"
"      --ignore-skip-centromere         Set this flag to assemble in centromeric/telomeric regions. Default is to skip.\n"
"\n";

namespace opt {
  static unsigned verbose = 1;
  static string regionFile = "";
  static std::string outdir = "./";
  static bool ignore_skip_cent = false;
  static size_t numThreads = 1;
  static string inbam = "";
  static string outbam = "";
  static size_t nmlim = 1;
}

enum {
  OPT_IGNORE_SKIP_CENT = 1
};

static const char* shortopts = "i:p:r:ho:r:v:n:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "input-bam",               required_argument, NULL, 'i' },
  { "output-bam",        required_argument, NULL, 'o' },
  { "threads (\"processes\")", required_argument, NULL, 'p' },
  { "region-file",             required_argument, NULL, 'r' },
  { "nm-limit",             required_argument, NULL, 'n' },
  { "ignore-skip-centromere" , no_argument, NULL, OPT_IGNORE_SKIP_CENT   },
  { "verbose",                 required_argument, NULL, 'v' },
  { NULL, 0, NULL, 0 }
};

bool runExtractor(int argc, char** argv) {
  
  parseExtractOptions(argc, argv);

  // print out information at runtime
  if (opt::verbose > 0) {
    cout << "Files:" << endl;
    cout << "  Input file: " << opt::inbam       << endl;
    cout << "  Output file: " << opt::outbam     << endl;
    cout << "  Region file: " << opt::regionFile << endl;
    cout << "Read filters:\n"; 
    cout << "   Min NM:                    " << opt::nmlim  << endl;
    cout << "   Ignore skip centromere:    " << (opt::ignore_skip_cent ? "ON" : "OFF") << endl; 
  }

  GenomicRegionVector regions;
  parseRegionFile(regions);

  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);

  // print the jobs
  if (opt::verbose > 0)
    for (GenomicRegionVector::const_iterator it = regions.begin(); it != regions.end(); it++)
      cout << "Input Regions: " << *it << endl;

  BamAlignment a;
  
  // setup the reader
  BamReader reader;

  // find the bam index
  string bai = opt::inbam;
  bai = bai.substr(0, opt::inbam.size()-3).append("bai");
  string bam_bai = opt::inbam;
  bam_bai = bam_bai.append(".bai");

  // open the bam file
  if (!reader.Open(opt::inbam)) {
    std::cerr << "FAILED TO OPEN BAM: " << opt::inbam << std::endl;
    exit(EXIT_FAILURE);
  }
  
  // open the BAM index
  if (!reader.OpenIndex(bai)) {
    if (!reader.OpenIndex(bam_bai) && opt::regionFile != "") {
      cerr << "FAILED TO OPEN INDEX: " << bam_bai << endl;
      exit(EXIT_FAILURE);
    }
  }

  // open the bam for writing
  BamWriter writer;

  //get the header
  SamHeader sam;
  string header = SVBamReader::getSamHeader(opt::inbam, sam);

  // set the reference for the BAM
  RefVector ref;  
  SVBamReader::getRefVector(opt::inbam, ref);

  if (!writer.Open(opt::outbam, header, ref)) {
    cerr << "Error initializing the BAM for: " << opt::outbam << endl;
    exit(EXIT_FAILURE);
  }

  // loop through the regions, extract the reads, and write out
  int count = 0;
  int keep_count = 0;
  int count_print = (opt::verbose > 2) ? 100 : 50000; // set how often to print update

  for (GenomicRegionVector::const_iterator it = regions.begin(); it != regions.end(); it++) {

    // set the BAM region
    if (!reader.SetRegion(it->chr, it->pos1, it->chr, it->pos2)) {
      cerr << "Failed to set BAM region " << *it << endl;
    } else {
     
      // loop through the reads in this region
      while (reader.GetNextAlignment(a)) {
	count++;
	uint32_t nm;

	if (!a.GetTag("NM", nm))
	  nm = 0;

	// print the progress
	if ( (count % count_print == 0) && opt::verbose > 1) {
	  char buffer[100];
	  sprintf(buffer, " Checking read %27s at position %2d:%-9ds. Kept %-7d of %-7d. NM tag: %2d", a.Name.c_str(), a.RefID+1, a.Position,  keep_count, count, nm);
	  cout << buffer << endl;
	  cout << buffer << SnowUtils::displayRuntime(start) << endl;
	}
	
	if (nm >= opt::nmlim && a.RefID < 24) {
	  keep_count++;
	  writer.SaveAlignment(a);
	}
      }

    } // end bam-region OK
  } // end region loop
  
  writer.Close();
  
  sleep(1);
  // try to open the BAM just made to create index
  BamReader final_reader;
  if (!final_reader.Open(opt::outbam)) {
    cerr << "Failed to open the BAM that was just created: " << opt::outbam;
    exit(EXIT_FAILURE);
  }
  

  // create the index file. The input bam should be sorted.
  if (!final_reader.CreateIndex()) {
    cerr << "Failed to created index for "  << opt::outbam << endl;
    cerr << " This can happen if the input bam is not sorted. Use samtools sort inbam.bam inbam.sort; samtools index inbam.sort.bam" << endl;
  }

  if (opt::verbose > 0) {
    char buffer[150];
    sprintf(buffer, " Finished. Kept %7d of %-7d", keep_count, count);
    cout << buffer << endl;
  }

  return true;
}

// Parse the region file specified by the -r flag. Must be comma delimited
bool parseRegionFile(GenomicRegionVector &gr) {

  ifstream fin;
  fin.open(opt::regionFile); // open a file
  if (!fin.good()) 
    return false; // exit if file not found

  // establish valid chr to run on
  unordered_map<string, int> VALID_CHR;
  for (unsigned i = 0; i < CHR_NAME_NUM.size(); i++)
    VALID_CHR.insert(pair<string, int>(CHR_NAME_NUM[i], 0));

  // get the centromeric positions
  //GenomicRegion gen;
  //GenomicRegionVector cent = gen.getCentromeres();

  // read each line of the file
  while (!fin.eof()) {

    // read an entire line into memory
    char buf[MAX_CHARS_PER_LINE];
    fin.getline(buf, MAX_CHARS_PER_LINE);
  
    // parse the line into blank-delimited tokens
    int n = 0; // a for-loop index
  
    // array to store memory addresses of the tokens in buf
    const char* token[MAX_TOKENS_PER_LINE] = {}; // initialize to 0
  
    // parse the line
    token[0] = strtok(buf, DELIMITER); // first token
    if (token[0]) // zero if line is blank
      for (n = 1; n < MAX_TOKENS_PER_LINE; n++) {
        token[n] = strtok(0, DELIMITER); // subsequent tokens
  	if (!token[n]) 
	  break; // no more tokens
       }
  	
    // loop through the lines
    for (int i = 0; i < n; i += 3) {
	string mchr  = token[i];
  	string mpos1 = token[i+1];
  	string mpos2 = token[i+2];
  
  	// deal with x and y
  	bool isx = mchr.compare("X")==0;
  	bool isy = mchr.compare("Y")==0;
	mchr = isx ? "23" : mchr;
	mchr = isy ? "24" : mchr;

  	string error = "";
  	bool valid_chr = VALID_CHR.find(mchr) != VALID_CHR.end();
  	if (!valid_chr)
  	  error = "Invalid chromosome";

  	// create the genomic region
  	GenomicRegion this_gr;
  	try {
  	  this_gr.chr = atoi(mchr.c_str())-1;
  	  this_gr.pos1 = min(atoi(mpos1.c_str()), CHR_LEN[this_gr.chr]);
  	  this_gr.pos2 = min(atoi(mpos2.c_str()), CHR_LEN[this_gr.chr]);
  	} catch (...) {
  	  cerr << "Caught error with parsing region file for region " << mchr << ":" << mpos1 << "-" << mpos2 << endl;
  	}

  	// deteremine if this region overlaps with a centromere
  	int overlap_result = 0;
	if (!opt::ignore_skip_cent) {
	  overlap_result = this_gr.centromereOverlap();
	  if (overlap_result != 0) 
	    error = (overlap_result == 1) ? "Partial overlap with centromere" : "Full overlap with centromere";
	}
  	    
  	// keep only those regions that are valid
  	if (valid_chr && overlap_result != 2) 
  	  gr.push_back(this_gr);

  	// print out if verbose
	if (opt::verbose > 0 && error.length() > 0) 
  	  cout << "   " << this_gr.toStringOffset() << " " << error << endl;

    } // end i+=3 for
  } // end while

  return true;
  
}

void parseExtractOptions(int argc, char** argv) {

  bool die = false;

  if (argc < 2) 
    die = true;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
      case 'p': arg >> opt::numThreads; break;
      case 'h': die = true; break;
      case 'i': arg >> opt::inbam; break;
      case 'v': arg >> opt::verbose; break;
      case 'r': arg >> opt::regionFile; break;
      case 'o': arg >> opt::outbam; break;
      case 'n': arg >> opt::nmlim; break;
      case OPT_IGNORE_SKIP_CENT: opt::ignore_skip_cent = true; break;
    }
  }

  if (opt::inbam == "" || opt::outbam == "")
    die = true;
  
  // clean the outdir
  //opt::outdir = getDirPath(opt::outdir);

  if(opt::numThreads <= 0) {
      cout << "invalid number of threads: " << opt::numThreads << "\n";
      die = true;
    }

  if (die) {
      cout << "\n" << EXTRACT_USAGE_MESSAGE;
      exit(1);
    }
}
