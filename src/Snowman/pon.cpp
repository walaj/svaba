#include "pon.h"

using namespace std;

static const char *PON_USAGE_MESSAGE =
"Usage: snowman pon [OPTION] \n\n"
"  Description: Generate a panel of normals file\n"
"\n"
" General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h, --help                           Display this help and exit\n"
" Optional Input\n"
"\n";

namespace popt {

  static size_t verbose = 1;
  static string input = "":
  static string output = "";
  static string input_bams = "";

}

static const char* shortopts = "hv:i:o:b:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "verbose",                 required_argument, NULL, 'v' },
  { "input-bam",              required_argument, NULL, 'i'},
  { "output-bam",              required_argument, NULL, 'o'},
  { "bam-list",              required_argument, NULL, 'b'},
  { NULL, 0, NULL, 0 }
};

static unordered_map<string, size_t> pon;

void runPON(int argc, char** argv) {

  parsePONOptions(argc, argv);

  if (popt::verbose > 0) {
    cout << "-- Input Normal BAM " << popt::input << endl;
    cout << "-- Output PON file:  " << popt::output << endl;
  }


  
  
}

void parsePONOptions(int argc, char** argv) {

  bool die = false;

  if (argc < 2) 
    die = true;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'i': arg >> popt::input; break;
    case 'v': arg >> popt::verbose; break;
    case 'b': arg >> popt::input_bams; break;
    }
  }

  // something went wrong, kill
  if (die) {
    cout << "\n" << PON_USAGE_MESSAGE;
    exit(1);
  }

}

