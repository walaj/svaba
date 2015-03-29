#include "nozzle.h"
#include <getopt.h>
#include <cstdlib>
#include <iostream>
#include <sstream>

using namespace std;

static const char *NOZZLE_USAGE_MESSAGE =
"Usage: snowman nozzle [OPTION] \n\n"
"  Description: Wrapper to an Rscript, which generates the nozzle report\n"
"\n"
" General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h, --help                           Display this help and exit\n"
"  -p, --threads                        Use NUM threads to run Snowman. Default: 1\n"
" Required Input\n"
"  -i, --input-dir                      A reads2contig BAM produced from snowman run (or -> snowman clean).\n"
"  -o, --output-dir                     A bam file containing aligned contigs. Produced at end of snowman run\n"
"\n";

namespace opt {
  static string rpath = "/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_3.1.0-bioconductor-2.14/bin/Rscript";
  static string npath = "/home/unix/jwala/GIT/isva/Snowman/SnowNozzleSimple.R";
  static string indir = "./";
  static string outdir = "./";
  static int verbose = 1;
  static int numThreads = 1;
}

static const char* shortopts = "hv:p:i:o:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "input-dir",               required_argument, NULL, 'i' },
  { "output-dir",              required_argument, NULL, 'o' },
  { "threads (\"processes\")", required_argument, NULL, 'p' },
  { "verbose",                 required_argument, NULL, 'v' },
  { NULL, 0, NULL, 0 }
};

void runNozzle(int argc, char** argv) {

  parseNozzleOptions(argc, argv);

  string cmd = opt::rpath + " " + opt::npath + " --indir=" + opt::indir + " --outdir=" + opt::outdir + " --cores=" + to_string(opt::numThreads);
  if (opt::verbose > 0)
    cout << cmd << endl;
  system(cmd.c_str());

}

void parseNozzleOptions(int argc, char** argv) {
  bool die = false;

  if (argc < 2) 
    die = true;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
      case 'p': arg >> opt::numThreads; break;
      case 'h': die = true; break;
      case 'i': arg >> opt::indir; break;
      case 'o': arg >> opt::outdir; break;
      case 'v': arg >> opt::verbose; break;
    }
  }

  // if filename, convert to directory
  //opt::indir     = getDirPath(opt::indir);
  //opt::outdir    = getDirPath(opt::outdir);

  if(opt::numThreads <= 0) {
    cout << "Invalid number of threads: " << opt::numThreads << "\n";
    die = true;
  }

  if (die) {
    cout << "\n" << NOZZLE_USAGE_MESSAGE;
    exit(1);
  }

}



