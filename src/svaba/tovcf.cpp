#include <getopt.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdio>

#include "gzstream.h"
#include "SeqLib/BamReader.h"

#include "vcf.h"
#include "BreakPoint.h"
#include "svabaUtils.h"
#include "SvabaSharedConfig.h"
#include "svabaOptions.h"
#include "svabaLogger.h"
#include "svabaOutputWriter.h"

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
  
  // Create minimal configuration objects for VCF conversion
  SvabaLogger logger;  // Use default constructor
  SvabaOptions opts;
  
  // Set default values for minimal functionality
  opts.lod = 8.0;
  opts.lodDb = 5.0;
  opts.lodSomatic = 2.5;
  opts.lodSomaticDb = 2.0;
  
  SvabaOutputWriter writer(logger, opts);  // Pass logger and opts as required
  
  SvabaSharedConfig config(logger, opts, writer);
  config.readlen = 150; // default read length
  config.header = bwalker.Header();
  
  // start a new VCF file
  VCFHeader header;
  header.filedate = svabaUtils::fileDateString();
  header.source = opt::input_file;
  header.reference = "";//opt::refgenome;
  
  // open BAM
  SeqLib::BamHeader hdr = bwalker.Header();
  
  // read in the bps.txt.gz file (which may or may not actually be gzipped)
  std::vector<std::string> allele_names; // store with real name
  
  // Try to determine if file is actually gzipped by checking magic bytes
  std::ifstream test_file(opt::input_file.c_str(), std::ios::binary);
  char magic[2] = {0};
  if (test_file.is_open()) {
    test_file.read(magic, 2);
    test_file.close();
  }
  
  bool is_gzipped = (magic[0] == '\x1f' && magic[1] == '\x8b');
  
  std::string headerLine;
  if (is_gzipped) {
    igzstream infile(opt::input_file.c_str(), std::ios::in);
    if (!infile) {
      std::cerr << "Can't read gzipped file " << opt::input_file << std::endl;
      exit(EXIT_FAILURE);
    }
    std::getline(infile, headerLine);
    infile.close();
  } else {
    std::ifstream infile(opt::input_file.c_str());
    if (!infile) {
      std::cerr << "Can't read file " << opt::input_file << std::endl;
      exit(EXIT_FAILURE);
    }
    std::getline(infile, headerLine);
    infile.close();
  }
  
  // Parse header to get sample names
  if (!headerLine.empty()) {
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

  
  // Handle the file format issue: VCFFile constructor expects .gz files to be gzipped
  // but svaba writes .bps.txt.gz files as plain text. Create a temporary file if needed.
  std::string vcf_input_file = opt::input_file;
  bool need_temp_file = false;
  
  if (!is_gzipped && opt::input_file.substr(opt::input_file.length() - 3) == ".gz") {
    // File has .gz extension but isn't gzipped - create temp file
    vcf_input_file = opt::analysis_id + ".temp.bps.txt";
    need_temp_file = true;
    
    std::ifstream src(opt::input_file.c_str());
    std::ofstream dst(vcf_input_file.c_str());
    dst << src.rdbuf();
    src.close();
    dst.close();
  }
  
  // convert to VCF
  VCFFile snowvcf(vcf_input_file, opt::analysis_id, bwalker.Header(), header, true,
		  opt::verbose > 0, &config);
  std::string basename = opt::analysis_id + ".svaba.unfiltered.";
  snowvcf.include_nonpass = true;
  snowvcf.writeIndels(basename, false, allele_names.size() == 1, bwalker.Header());
  snowvcf.writeSVs(basename, false, allele_names.size() == 1, bwalker.Header());
  basename = opt::analysis_id + ".svaba.";
  snowvcf.include_nonpass = false;
  snowvcf.writeIndels(basename, false, allele_names.size() == 1, bwalker.Header());
  snowvcf.writeSVs(basename, false, allele_names.size() == 1, bwalker.Header());
  
  // Clean up temp file if created
  if (need_temp_file) {
    std::remove(vcf_input_file.c_str());
  }
  
  return;
}
