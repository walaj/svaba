#include "assembly2vcf.h"

#include <getopt.h>
#include <string>
#include <sstream>
#include <iostream>

#include "vcf.h"
#include "AssemblyBamWalker.h"
#include "SnowmanBamWalker.h"
#include "SnowmanUtils.h"

#define MIN_ISIZE_FOR_DISC 700 

faidx_t * findex;

enum { 
  OPT_NO_READS
};

struct bidx_delete {
  void operator()(void* x) { hts_idx_destroy((hts_idx_t*)x); }
};

static const char* shortopts = "hi:m:q:p:v:k:z:x:a:r:t:n:G:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "assembly-bam",            required_argument, NULL, 'i' },
  { "tumor-reads-bam",         required_argument, NULL, 't' },
  { "normal-reads-bam",        required_argument, NULL, 'n' },
  { "indel-mask",              required_argument, NULL, 'm' },
  { "panel-of-normals",        required_argument, NULL, 'q' },
  { "id-string",               required_argument, NULL, 'a' },
  { "normal-bam",              required_argument, NULL, 'n' },
  { "threads",                 required_argument, NULL, 'p' },
  { "region-file",             required_argument, NULL, 'k' },
  { "rules",                   required_argument, NULL, 'r' },
  { "no-zip",                  no_argument, NULL, 'z' },
  { "no-read-tracking",        no_argument, NULL, OPT_NO_READS },
  { "no-r2c-bam",              no_argument, NULL, 'x' },
  { "threads",                 required_argument, NULL, 'p' },
  { "verbose",                 required_argument, NULL, 'v' },
  { "reference-genome",        required_argument, NULL, 'G' },
  { NULL, 0, NULL, 0 }
};


static const char *ASSEMBLY2VCF_USAGE_MESSAGE =
"Usage: snowman assembly2vcf -i <assembly_bam> -r <reads_bam>\n\n"
"  Description: Take a BAM (from de novo assembly) and aligned reads --> Call variants\n"
"\n"
"  General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h, --help                           Display this help and exit\n"
"  -p, --threads                        Use NUM threads to run snowman. Default: 1\n"
"  Required input\n"
"  -i, --assembly-bam                   BAM from aligned de-novo assembly\n"
"  -t, --tumor-reads-bam                BAM from tumor reads\n"
"  -n, --normal-reads-bam               BAM from normal reads\n"
"  -G, --reference-genome               Path to indexed reference genome to be used by BWA-MEM. Default is Broad hg19 (/seq/reference/...)\n"
"  Optional input\n"                       
"  -a, --id-string                      String specifying the analysis ID to be used as part of ID common.\n"
"  -r, --rules                          VariantBam style rules string to determine which reads to do assembly on. See documentation for default.\n"
"  -k, --region-file                    Set a region txt file. Format: one region per line, Ex: 1,10000000,11000000\n"
"  -q, --panel-of-normals               Panel of normals gzipped txt file generated from snowman pon\n"
"  -m, --indel-mask                     BED-file with blacklisted regions for indel calling. Default /xchip/gistic/Jeremiah/Projects/HengLiMask/um75-hs37d5.bed.gz\n"
"      --no-read-tracking               Don't track supporting reads. Reduces file size.\n"
"      --no-zip                         Don't tabix and gzip the output vcf files.\n"
"\n";

namespace opt {

  static std::string args = "snowman ";
  static std::string assembly_bam;
  static std::string tumor_reads_bam;
  static std::string normal_reads_bam;
  static std::string analysis_id = "assembly2vcf_noid";
  static size_t verbose = 1;
  static std::string rules = "{\"global\" : {\"duplicate\" : false, \"qcfail\" : false}, \"\" : { \"rules\" : [{\"isize\" : [MINS,0]},{\"rr\" : true},{\"ff\" : true}, {\"rf\" : true}, {\"ic\" : true}, {\"clip\" : 5, \"phred\" : 4}, {\"ins\" : true}, {\"del\" : true}, {\"mapped\": true , \"mate_mapped\" : false}, {\"mate_mapped\" : true, \"mapped\" : false}]}}"; 
  static std::string indel_mask = "";
  static std::string pon = "";
  static bool no_reads = false;
  static bool zip = true;
  static size_t numThreads = 1;
  static bool no_r2c = false;
  static std::string regionFile = "";
  static std::string refgenome = SnowTools::REFHG19;  
}



void runAssembly2VCF(int argc, char** argv)
{

  // parse the options
  parseAssembly2VCFOptions(argc, argv);

  // load the reference
  findex = fai_load(opt::refgenome.c_str());  // load the reference

  // set the min isize for FR discordants
  opt::rules.replace(opt::rules.find("MINS"), 4, std::to_string(MIN_ISIZE_FOR_DISC));

  SnowTools::MiniRulesCollection * mr = new SnowTools::MiniRulesCollection(opt::rules);

  SnowmanBamWalker twalk, nwalk;
  if (SnowTools::read_access_test(opt::tumor_reads_bam)) {
    twalk = SnowmanBamWalker(opt::tumor_reads_bam);
    twalk.prefix= "t000";  
    twalk.max_cov = 500;
    twalk.SetMiniRulesCollection(*mr);
  } 
  if (SnowTools::read_access_test(opt::normal_reads_bam)) {
    nwalk = SnowmanBamWalker(opt::normal_reads_bam);
    nwalk.prefix= "n000";
    nwalk.max_cov = 500;
    nwalk.SetMiniRulesCollection(*mr);    
  }
    
  // read in the assembly bam file
  AssemblyBamWalker awalk(opt::assembly_bam);
  awalk.twalk = twalk;
  awalk.nwalk = nwalk;
  awalk.findex = findex;
  awalk.id = opt::analysis_id;
  if (SnowTools::read_access_test(opt::tumor_reads_bam))
    awalk.tindex = std::shared_ptr<hts_idx_t>(hts_idx_load(opt::tumor_reads_bam.c_str(), HTS_FMT_BAI), bidx_delete());
  if (SnowTools::read_access_test(opt::normal_reads_bam))
    awalk.nindex = std::shared_ptr<hts_idx_t>(hts_idx_load(opt::normal_reads_bam.c_str(), HTS_FMT_BAI), bidx_delete());
  else
    opt::normal_reads_bam = std::string();

  std::cerr << awalk << std::endl;

  awalk.tbam = opt::tumor_reads_bam;
  awalk.nbam = opt::normal_reads_bam;
  awalk.numThreads = opt::numThreads;
  awalk.refGenome = opt::refgenome;
  awalk.walkDiscovar();

  // make the VCFs
  std::cerr << "...loading the bps files for conversion to VCF" << std::endl;
  std::string file = opt::analysis_id + ".bps.txt.gz";
  if (!SnowTools::read_access_test(file))
    file = opt::analysis_id + ".bps.txt";

  // put args into string for VCF later
  for (int i = 0; i < argc; ++i)
    opt::args += std::string(argv[i]) + " ";

  // make the header
  VCFHeader header;
  header.filedate = SnowmanUtils::fileDateString();
  header.source = opt::args;
  header.reference = opt::refgenome;
  
  //boost::filesystem::path tfp(opt::tumor_reads_bam);
  //boost::filesystem::path nfp(opt::normal_reads_bam);

  // TODO fix this
  header.addSampleField(opt::tumor_reads_bam);
  header.colnames += "\t" + opt::tumor_reads_bam;
  header.addSampleField(opt::normal_reads_bam);
  header.colnames += "\t" + opt::normal_reads_bam;
  //header.addSampleField(nfp.filename().string());
  //header.colnames += "\t" + nfp.filename().string(); 
  //header.addSampleField(tfp.filename().string());
  //header.colnames += "\t" + tfp.filename().string(); 
  
  bool zip = false;
  VCFFile snowvcf(file, opt::analysis_id, twalk.header(), header);
  std::string basename = opt::analysis_id + ".unfiltered.";
  snowvcf.include_nonpass = true;
  snowvcf.writeIndels(basename, zip, false);
  snowvcf.writeSVs(basename, zip, false);
  
  basename = opt::analysis_id + ".assembly."; //"assembly.";
  snowvcf.include_nonpass = false;
  snowvcf.writeIndels(basename, zip, false);
  snowvcf.writeSVs(basename, zip, false);

}

// parse the command line options
void parseAssembly2VCFOptions(int argc, char** argv) {

  bool die = false;

  if (argc <= 2) 
    die = true;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
      case 'p': arg >> opt::numThreads; break;
      case 'a': arg >> opt::analysis_id; break;
      case 'm': arg >> opt::indel_mask; break;
      case 'q': arg >> opt::pon; break;
      case 'z': opt::zip = false; break;
      case 'h': die = true; break;
      case 'x': opt::no_r2c = true; break;
      case OPT_NO_READS: opt::no_reads = true; break;
      case 'i': arg >> opt::assembly_bam; break;
      case 't': arg >> opt::tumor_reads_bam; break;
      case 'n': arg >> opt::normal_reads_bam; break;
      case 'k': arg >> opt::regionFile; break;
      case 'r': arg >> opt::rules; break;
      case 'v': arg >> opt::verbose; break;
      case 'G': arg >> opt::refgenome; break;
      default: die= true; 
    }
  }

  if (die) 
    {
      std::cout << "\n" << ASSEMBLY2VCF_USAGE_MESSAGE;
      exit(EXIT_FAILURE);
    }
}
