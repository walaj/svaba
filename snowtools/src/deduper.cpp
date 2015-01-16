#include "deduper.h"
#include <iostream>
#include <vector>
#include <fstream>
#include "GenomicRegion.h"
#include <getopt.h>
#include <unordered_map>
#include <unistd.h>

using namespace std;

namespace opt {
  static unsigned verbose = 1;
  static string infasta = "";
  static string outfasta = "";
}

static const char* shortopts = "i:ho:v:";
static const struct option longopts[] = {
  { "help",                no_argument, NULL, 'h' },
  { "input-fasta",         required_argument, NULL, 'i' },
  { "output-fasta",        required_argument, NULL, 'o' },
  { "verbose",             required_argument, NULL, 'v' },
  { NULL, 0, NULL, 0 }
};


static const char *DEDUPER_USAGE_MESSAGE = 
"Usage: snowtools deduper [OPTION] -i in.fasta -o out.fasta\n\n"
"  Description: Deduplicate a FASTA by sequences, taking into account reverse-complement duplicates..\n"
"\n"
"  General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h, --help                           Display this help and exit\n"
"  Required input\n"
"  -i, --input-fasta                    Input fasta file\n"
"  -o, --output-fasta                   Output fasta file\n"
"\n";

void runDeduper(int argc, char** argv) {

  parseDeduperOptions(argc, argv);

  if (opt::verbose > 0)
    cerr << "...running fasta de-duplicator. Outputing to " << opt::outfasta << endl;
  size_t count = 0;
  size_t kept = 0;
  
  unordered_map<string, string> * map = new unordered_map<string, string>();
  string line;
  string qname;
  while (getline(cin, line)) {
    if (line.at(0) == '>') {
      qname = line;
    } else {

      count++;
      if (count % 100000 == 0)
	cerr << "kept " << kept << " of " << count << endl;

      unordered_map<string, string>::const_iterator ff = map->find(line);
      if (ff == map->end()) {
	rcomplement(line);
	ff = map->find(line);
	if (ff == map->end()) {
	  map->insert(pair<string, string>(line, qname));
	  kept++;
	}
	
      }
    }
  }

  if (opt::verbose > 0)
    cout << "kept " << kept << " of " << count << endl;

  // writing the contig file
  ofstream ostream;
  string ofile;
  ostream.open(opt::outfasta);
  
  // write the contigs to the .tmp.fa file
  for(unordered_map<string, string>::const_iterator it = map->begin(); it != map->end(); ++it) {
    ostream << it->second << "\n";
    ostream << it->first << "\n";
  }
  ostream.close();

}

void rcomplement(std::string &a) {

  reverse(&a[0], &a[a.size()]);
  string::iterator it = a.begin();
  for (; it != a.end(); it++)
    if (*it == 'A')
      *it = 'T';
    else if (*it == 'T')
      *it = 'A';
    else if (*it == 'C')
      *it = 'G';
    else
      *it = 'C';
}

void parseDeduperOptions(int argc, char** argv) {

  bool die = false;

  if (argc < 2) 
    die = true;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
      case 'h': die = true; break;
      case 'i': arg >> opt::infasta; break;
      case 'v': arg >> opt::verbose; break;
      case 'o': arg >> opt::outfasta; break;
    }
  }

  if (die) {
      cout << "\n" << DEDUPER_USAGE_MESSAGE;
      exit(1);
    }
}
