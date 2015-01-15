#include "bamprocess.h"
#include <string>
#include <vector>
#include <getopt.h>
#include <iostream>
#include "GenomicRegion.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "SVBamReader.h"

using namespace std;
using namespace BamTools;

static const char *PREP_USAGE_MESSAGE =
"Usage: snowman preprocess [OPTION] \n\n"
"  Description: Process a BAM file for use with SnowmanSV by removing proper pairs and bad regions\n"
"\n"
" General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h, --help                           Display this help and exit\n"
"  -i, --input-bam                      BAM file to preprocess\n"
"  -o, --output-bam                     File to write final BAM to\n"
" Read Filters\n"
"  -w, --min-mapq                       Minimum mapping quality a read must have to be included. Default 0\n"
"  -n  --nm-limit                       Skip reads that have more than nm mismatches (NM tag). Default: 5\n"
"  -p  --perc-limit                     Exit with failure if more than this percent of reads are weird (Must see over 50m). Default: 30\n"
"  -s, --isize                          Insert size threshold to consider discordant. Default: 800\n"
" Optional Input\n"
"\n";

namespace opt {
  static string bam;
  static string out;
  static int isize = 800;
  static size_t verbose = 1;

  static int mapq = 0;
  static int qualthresh = 4;
  static int minOverlap = 30;
  static int min_clip = 5;
  static bool skip_supp = true;
  static bool skip_cent = true;

  static int nmlim = 5;
  static int perc_limit = 30;
  
}

static const char* shortopts = "hv:i:o:w:n:p:s:";
static const struct option longopts[] = {
  { "help",                       no_argument, NULL, 'h' },
  { "verbose",                    required_argument, NULL, 'v' },
  { "input-bam",                  required_argument, NULL, 'i' },
  { "output-bam",                 required_argument, NULL, 'o' },
  { "min-mapq",                 required_argument, NULL, 'w' },
  { "nm-limit",                 required_argument, NULL, 'n' },
  { "perc-limit",                 required_argument, NULL, 'p' },
  { "isize",                 required_argument, NULL, 's' },
  { NULL, 0, NULL, 0 }
};

void runPrep(int argc, char** argv) {

  parsePrepOptions(argc, argv);

  if (opt::verbose > 0) {
    cout << "Input BAM:  " << opt::bam << endl;
    cout << "Output BAM: " << opt::out << endl;
    cout << "   Read Filters:   " << endl;
    cout << "   Min insert size:           " << opt::isize << endl;
    cout << "   Min clip length:           " << opt::min_clip << endl;
    cout << "   Percent weird limit:       " << opt::perc_limit << endl;
    cout << "   Min MAPQ:                  " << opt::mapq  << endl;
    cout << "   Max NM:                    " << opt::nmlim  << endl;
    cout << "   Insert size:               " << opt::isize << endl;
    cout << "   Min Overlap:               " << opt::minOverlap << endl;
    cout << "   Min quality score:         " << opt::qualthresh << endl;
    cout << "   Skip supplementary reads:  " << (opt::skip_supp ? "ON" : "OFF") << endl;
    cout << "   Skip centromere:           " << (opt::skip_cent ? "ON" : "OFF") << endl; 
  }

  // open the BAM Reader
  BamReader reader;
  if (!reader.Open(opt::bam)) {
    cerr << "Could not open BAM: " << opt::bam << endl;
    exit(EXIT_FAILURE);
  }

  // get the BAM header
  RefVector ref;
  SamHeader sam;
  SVBamReader::getRefVector(opt::bam, ref);
  SVBamReader::getSamHeader(opt::bam, sam);

  // open the BAM Writer
  BamWriter writer;
  if (!writer.Open(opt::out, sam, ref)) {
    cerr << "Could not open output BAM: " << opt::out << endl;
    exit(EXIT_FAILURE);
  }

  BamAlignment a;

  GenomicRegion gr(0, 10000, 100000);

  SVBamReader sv_reader(opt::bam, "", opt::isize, opt::mapq, opt::qualthresh, opt::minOverlap, opt::skip_supp, opt::min_clip, opt::verbose);
  sv_reader.perclimit = opt::perc_limit;
  if (!sv_reader.findBamIndex()) {
    cerr  << "Failed to open BAM index" << endl;
    exit(EXIT_FAILURE);
  }

  GenomicRegionVector grv = GenomicRegion::non_centromeres;

  for (auto it = grv.begin(); it != grv.end(); it++) {
    if (!sv_reader.setBamRegion(*it)) {
      cerr << "Failed to set BAM region " << (*it) << endl;
      exit(EXIT_FAILURE);
    }
    if (opt::verbose > 0)
      cout << "Running region: "  << (*it) << endl;
    sv_reader.preprocessBam(writer);
  }

  writer.Close();

  // index it
  if (!reader.Open(opt::out)) {
    cerr << "Cant read output BAM" << endl;
    exit(EXIT_FAILURE);
  }

  // create the index file. The input bam should be sorted.                                                                                                                                                                                   
  if (!reader.CreateIndex()) {                                                                                                                                                                                                          
    cerr << "Failed to created index for "  << opt::out << endl;
    cerr << " This can happen if the input bam is not sorted. Use samtools sort inbam.bam inbam.sort; samtools index inbam.sort.bam" << endl;
    exit(EXIT_FAILURE);
  }
}

void parsePrepOptions(int argc, char** argv) {

  bool die = false;

  if (argc < 2) 
    die = true;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'h': die = true; break;
    case 'v': arg >> opt::verbose; break;
    case 'i': arg >> opt::bam; break;
    case 'o': arg >> opt::out; break;
    case 'w': arg >> opt::mapq; break;
    case 'n': arg >> opt::nmlim; break;
    case 'p': arg >> opt::perc_limit; break;
    case 's': arg >> opt::isize; break;
    }
  }

  // something went wrong, kill
  if (die) {
    cout << "\n" << PREP_USAGE_MESSAGE;
    exit(1);
  }

}
