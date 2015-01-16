//-----------------------------------------------
// Copyright 2014 Broad Institute 
// Modified by Jeremiah Wala (jwala@broadinstitute.org)
// Released under the GPL
//-----------------------------------------------

//g++ -std=c++11 -I $bamt/include -L $bamt/lib snowtools0.cpp GenomicRegion.h GenomicRegion.cpp -o snowtools -lbamtools

#include <string>
#include <iostream>
#include "extract.h"
#include "deduper.h"
#include <unistd.h>

#define PROGRAM_BIN "snowtools"
#define AUTHOR "Jeremiah Wala <jwala@broadinstitute.org>"

using namespace std;
using namespace BamTools;

static const char *SNOWTOOLS_USAGE_MESSAGE =
"Contact: " AUTHOR "\n"
"Usage: " PROGRAM_BIN " <program> [options]\n\n"
"  Description: Collection of C++ tools for maniuplating NGS data. Developed with SnowmanSV.\n"
"\n"
"  Programs:\n"
"    extract             Extracts a regions of interest from a BAM file and writes a shorted BAM file.\n"
"    dedupe              Removes duplicate reads from a FASTA file, taking into account reverse complementarity.\n"
"\n";

int main(int argc, char** argv) {

  if (argc <= 1) {
    cerr << SNOWTOOLS_USAGE_MESSAGE; 
    return 0;
  } else {
    string command(argv[1]);
    if (command == "help" || command == "--help") {
      cerr << SNOWTOOLS_USAGE_MESSAGE;
      return 0;
    } else if (command == "extract") {
      runExtractor(argc -1, argv + 1);
    } else if (command == "dedupe") {
      runDeduper(argc -1, argv + 1);
    } else {
      cerr << SNOWTOOLS_USAGE_MESSAGE;
      return 0;
    }
  }
  
  return 0;

}
