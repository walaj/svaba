#include "benchmark.h"

#include <getopt.h>
#include <iostream>
#include <string>
#include <sstream>

#include "SnowTools/SnowToolsCommon.h"

namespace opt {

  static std::string refgenome = SnowTools::REFHG19;  
  static int mode = 0;
}

enum { 
  OPT_ASSEMBLY
};


static const char *BENCHMARK_USAGE_MESSAGE =
"Usage: snowman benchmark\n\n"
"  Description: Various benchmarking tests for Snowman\n"
"\n"
"  General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -G, --reference-genome               Path to indexed reference genome to be used by BWA-MEM. Default is Broad hg19 (/seq/reference/...)\n"
"      --test-assembly                  dd\n"
"\n";


static const char* shortopts = "haG:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "reference-genome",        required_argument, NULL, 'G' },
  { "assembly",                no_argument, NULL, OPT_ASSEMBLY},
  { NULL, 0, NULL, 0 }
};

void runBenchmark(int argc, char** argv) {

  parseBenchmarkOptions(argc, argv);

  std::cerr << 
    "-----------------------------------------" << std::endl << 
    "--- Running Snowman Benchmarking Test ---" << std::endl <<
    "-----------------------------------------" << std::endl;
  if (opt::mode == OPT_ASSEMBLY)
    std::cerr << "********* RUNNING ASSEMBLY TEST ***********" << std::endl;
}


void parseBenchmarkOptions(int argc, char** argv) {

  bool die = false;

  if (argc <= 2) 
    die = true;

  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'h': die = true; break;
    case OPT_ASSEMBLY: opt::mode = OPT_ASSEMBLY; break;
    case 'G': arg >> opt::refgenome; break;
    default: die= true; 
    }
  }

  if (die) {
    std::cerr << "\n" << BENCHMARK_USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }


}
