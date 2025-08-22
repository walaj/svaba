#include <getopt.h>
#include <sstream>
#include <iostream>

#include "gzstream.h"
#include "SeqLib/BamReader.h"

#include "vcf.h"
#include "BreakPoint.h"
#include "svabaUtils.h"

void parseToVCFOptions(int argc, char** argv);

namespace opt {

  static std::string input_file;
  static std::string bam;
  static std::string analysis_id;
  static int verbose = 0;
}

static const char* shortopts = "hi:b:va:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "input-bps",               required_argument, NULL, 'i'},
  { "bam",                     required_argument, NULL, 'b'},
  { "verbose",                 no_argument, NULL, 'v' },
  { "id-string",               required_argument, NULL, 'a'},  
  { NULL, 0, NULL, 0 }
};

static const char *TOVCF_USAGE_MESSAGE =
"Usage: svaba tovcf [OPTION] -i bps.txt.gz -b <bam>\n\n"
"  Description: Convert a *bps.txt.gz file to a *vcf file\n"
"\n"
"  General options\n"
"  -v, --verbose                        Flag to make more verbose\n" 
"  -h, --help                           Display this help and exit\n"
"  Required input\n"
"  -a, --id-string                      String specifying the analysis ID to be used as part of ID common.\n"  
"  -i, --input-bps                      Original bps.txt.gz file\n"
"  -b, --bam                            BAM file used to grab header from\n"
"\n";

// parse the command line options
void parseToVCFOptions(int argc, char** argv) {
  bool die = false;
  
  if (argc <= 2) 
    die = true;
  
  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'h': die = true; break;
    case 'i': arg >> opt::input_file; break;
    case 'v': opt::verbose = 1; break;
    case 'a': arg >> opt::analysis_id; break;
    case 'b': arg >> opt::bam; break; 
    }
  }

  if (opt::analysis_id.empty()) {
    std::cerr << "Error: svaba tovcf -- -a analysis id required, this names the VCF file" << std::endl;
    die = true;
  }

  // just check if input here. Check later if readable
  if (opt::input_file.length() == 0) 
    die = true;
  
  if (opt::bam.length() == 0) {
    std::cerr << "BAM is required (for the header)" << std::endl;
    die = true;
  }

  if (die) {
    std::cerr << "\n" << TOVCF_USAGE_MESSAGE;
    exit(1);
  }
}

void runToVCF(int argc, char** argv) {
  
  parseToVCFOptions(argc, argv);
  
  if (opt::verbose > 0) {
    std::cerr << "Input bps file:  " << opt::input_file << std::endl;
    std::cerr << "Analysis id: " << opt::analysis_id << std::endl;
  }

  if (!SeqLib::read_access_test(opt::input_file)) {
    std::cerr << "ERROR: Cannot read file " << opt::input_file  << std::endl;
    exit(EXIT_FAILURE);
  }
  
  SeqLib::BamReader bwalker;
  assert(bwalker.Open(opt::bam));
  
  // start a new VCF file
  VCFHeader header;
  header.filedate = svabaUtils::fileDateString();
  header.source = opt::input_file;
  header.reference = "";//opt::refgenome;
  
  // open BAM
  SeqLib::BamHeader hdr = bwalker.Header();
  
  // read in the bps.txt.gz file
  std::vector<std::string> allele_names; // store with real name
  std::map<std::string, SampleInfo> tmp_alleles;
  igzstream infile(opt::input_file.c_str(), std::ios::in);
  size_t line_count = 0;
  
  // Read the header line first
  std::string headerLine;
  if (std::getline(infile, headerLine)) {
    std::vector<std::string> headerv = svabaUtils::tokenize_delimited(headerLine, '\t');

    // assume a certain format for bps.txt
    assert(headerv.size() >= 39);
    // everything at (0-based) 38 and above is a sample id
    for (size_t i = 38; i < headerv.size(); i++) {
      assert(headerv[i].at(0) == 't' || headerv[i].at(0) == 'n');
      allele_names.push_back(headerv[i]);
      header.colnames += "\t" + headerv[i].substr(5);
    }
  }

  
  // convert to VCF
  VCFFile snowvcf(opt::input_file, opt::analysis_id, bwalker.Header(), header, true,
		  opt::verbose > 0);
  std::string basename = opt::analysis_id + ".svaba.unfiltered.";
  snowvcf.include_nonpass = true;
  snowvcf.writeIndels(basename, false, allele_names.size() == 1);
  snowvcf.writeSVs(basename, false, allele_names.size() == 1);
  basename = opt::analysis_id + ".svaba.";
  snowvcf.include_nonpass = false;
  snowvcf.writeIndels(basename, false, allele_names.size() == 1);
  snowvcf.writeSVs(basename, false, allele_names.size() == 1);
  
  return;
}
