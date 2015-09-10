#include "refilter.h"

#include "vcf.h"
#include "SnowTools/BreakPoint.h"
#include <iostream>
#include <getopt.h>

#include "SnowTools/SnowUtils.h"
#include "SnowTools/gzstream.h"

#include <sstream>

namespace opt {

  static std::string input_file = "";
  static std::string output_file = "";
  static std::string pon = "";
  static std::string analysis_id = "noid";
  static bool noreads = false;
  static std::string indel_mask = "";

  static std::string ref_index = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
  static std::string bam;

  static int verbose = 1;
}

static const char* shortopts = "hxi:a:v:q:g:m:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "input-bps",               required_argument, NULL, 'i'},
  { "panel-of-normals",        required_argument, NULL, 'q'},
  { "indel-mask",              required_argument, NULL, 'm'},
  { "reference-genome",        required_argument, NULL, 'g'},
  { "analysis-id",             required_argument, NULL, 'a'},
  { "no-reads",                no_argument, NULL, 'x'},
  { "verbose",                 required_argument, NULL, 'v' },
  { NULL, 0, NULL, 0 }
};


static const char *BP_USAGE_MESSAGE =
"Usage: snowman refilter [OPTION] -i bps.txt.gz -o bps.new.txt.gz\n\n"
"  Description: \n"
"\n"
"  General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h, --help                           Display this help and exit\n"
"  -g, --reference-genome               Path to indexed reference genome to be used by BWA-MEM. Default is Broad hg19 (/seq/reference/...)\n"
"  Required input\n"
"  -i, --input-bps                      Original bps.txt.gz file\n"
"  Optional input\n"                       
"  -q, --panel-of-normals               Panel of normals files generated from snowman pon\n"                       
"  -a, --id-string                      String specifying the analysis ID to be used as part of ID common.\n"
"  -x, --no-reads                       Flag to turn off recording of supporting reads. Setting this flag greatly reduces file size.\n"
"  -m, --indel-mask                     BED-file with blacklisted regions for indel calling. Default none\n"
"\n";

// parse the command line options
void parseBreakOptions(int argc, char** argv) {
  bool die = false;
  
  if (argc <= 2) 
    die = true;
  
  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'h': die = true; break;
    case 'b': arg >> opt::bam; break;
    case 'g': arg >> opt::ref_index; break;
    case 'm': arg >> opt::indel_mask; break;
    case 'i': arg >> opt::input_file; break;
    case 'v': arg >> opt::verbose; break;
    case 'q': arg >> opt::pon; break;
    case 'a': arg >> opt::analysis_id; break;
    case 'x': opt::noreads = true; break;
    }
  }
  
  if (opt::input_file.length() == 0)
    die = true;
  
  if (opt::bam.length() == 0) {
    std::cerr << "BAM is required (for the header)" << std::endl;
    die = true;
  }

  if (die) {
    std::cout << "\n" << BP_USAGE_MESSAGE;
    exit(1);
  }
}

void runRefilterBreakpoints(int argc, char** argv) {
  
  parseBreakOptions(argc, argv);
  
  opt::output_file = opt::analysis_id + ".filtered.bps.txt.gz";
  if (opt::verbose > 0) {
    std::cout << "Input bps file:  " << opt::input_file << std::endl;
    std::cout << "Output bps file: " << opt::output_file << std::endl;
    std::cout << "Panel of normals file: " << opt::pon << std::endl;
    std::cout << "Indel mask BED:      " << opt::indel_mask << std::endl;
    std::cout << "Analysis id: " << opt::analysis_id << std::endl;
  }
    
  if (!SnowTools::read_access_test(opt::input_file)) {
    std::cerr << "ERROR: Cannot read file " << opt::input_file  << std::endl;
    exit(EXIT_FAILURE);
  }
  
  SnowTools::BamWalker bwalker(opt::bam);

  // load the reference index
  //if (opt::verbose > 0)
  //  std::cout << "attempting to load: " << opt::ref_index << std::endl;
  //faidx_t * findex = fai_load(opt::ref_index.c_str());  // load the reference
  
  // open the output file
  igzstream iz(opt::input_file.c_str());
  if (!iz) {
    std::cerr << "Can't read file " << opt::input_file << std::endl;
    exit(EXIT_FAILURE);
  }
  
  // read in the PON
  //std::unique_ptr<SnowTools::PON> pmap;
  //SnowTools::BreakPoint::readPON(opt::pon, pmap);
  
  // read the indel mask
  //GenomicIntervalTreeMap grm_mask;
  if (!SnowTools::read_access_test(opt::indel_mask) && opt::indel_mask.length()) {
    std::cerr << "indel mask " << opt::indel_mask << " does not exist / is not readable. Skipping indel masking."  << std::endl;
    opt::indel_mask = "";
  }
  
  SnowTools::GRC grv_mask;
  if (opt::indel_mask.length()) 
    grv_mask.regionFileToGRV(opt::indel_mask, 0, NULL);
  
  ogzstream oz(opt::output_file.c_str(), std::ios::out);
  if (!oz) {
    std::cerr << "Can't write to output file " << opt::output_file << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (opt::verbose)
    std::cout << "...refiltering variants" << std::endl;
  
  // set the header
  oz << SnowTools::BreakPoint::header() << std::endl;
  
  std::string line;

  //skip the header
  std::getline(iz, line, '\n');
  
  size_t count = 0;
  while (std::getline(iz, line, '\n')) {
    
    if (++count % 10000 == 1 && opt::verbose > 0)
      std::cerr << "filtering breakpoint " << count << std::endl;
    
    SnowTools::BreakPoint bp(line, bwalker.header());
    
    // check if in panel of normals
    //bp.checkPon(pmap);
    
    // check for blacklist
    bp.checkBlacklist(grv_mask);
    
    // check for repeat sequence
    if (bp.b1.gr.chr < 24) {
      SnowTools::GenomicRegion gr = bp.b1.gr;
      gr.pad(20);
      
      // get the reference seqeuence for this piece
      //int len;
      //std::string chrstring = SnowTools::GenomicRegion::chrToString(gr.chr);
      //char * seq = faidx_fetch_seq(findex, const_cast<char*>(chrstring.c_str()), gr.pos1-1, gr.pos2-1, &len);
      //if (!seq) {
      //std::cerr <<  "faidx_fetch_seq fail" << std::endl;
      //}
      //std::string seqr = std::string(seq);

      //for (auto& i : repr)
      //if (seqr.find(i) != std::string::npos)
      //	  bp.repeat_seq = i;
      
    }
    
    //if (bp.hasMinimal() && bp.b1.gr.chr < 24 && bp.gr2.chr < 24)
    oz << bp.toFileString(opt::noreads) << std::endl;
  }
  
  oz.close();
  
  // make the VCF files
  /*
    if (opt::verbose)
    std::cout << "...converting filtered.bps.txt.gz to vcf files" << std::endl;
    VCFFile filtered_vcf(opt::output_file, opt::ref_index.c_str(), '\t', opt::analysis_id);
    
    // output the vcfs
      string basename = opt::analysis_id + ".broad-snowman.DATECODE.filtered.";
      filtered_vcf.writeIndels(basename, true); // zip them -> true
      filtered_vcf.writeSVs(basename, true);
  */
}
