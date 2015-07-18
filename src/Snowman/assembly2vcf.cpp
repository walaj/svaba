#include "assembly2vcf.h"

#include <getopt.h>
#include <string>
#include <iostream>
#include <sstream>

#include "AssemblyBamWalker.h"

enum { 
  OPT_NO_READS
};

static const char* shortopts = "hi:m:q:p:v:k:z:x:a:r:t:n:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "assembly-bam",            required_argument, NULL, 'i' },
  { "tumor-reads-bam",            required_argument, NULL, 't' },
  { "normal-reads-bam",            required_argument, NULL, 'n' },
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
  { "verbose",                 required_argument, NULL, 'v' },
  { NULL, 0, NULL, 0 }
};


static const char *ASSEMBLY2VCF_USAGE_MESSAGE =
"Usage: snowman assembly2vcf -i <assembly_bam> -r <reads_bam>\n\n"
"  Description: Take a BAM (from de novo assembly) and aligned reads --> Call variants\n"
"\n"
"  General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h, --help                           Display this help and exit\n"
"  Required input\n"
"  -i, --assembly-bam                   BAM from aligned de-novo assembly\n"
"  -t, --tumor-reads-bam                BAM from tumor reads\n"
"  -n, --normal-reads-bam               BAM from normal reads\n"
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

  static std::string assembly_bam;
  static std::string tumor_reads_bam;
  static std::string normal_reads_bam;
  static std::string analysis_id = "assembly2vcf_noid";
  static size_t verbose = 1;
  static std::string rules = "global@nbases[0,0];!hardclip;!supplementary;!duplicate;!qcfail;phred[4,100];%region@WG%discordant[0,1000];mapq[1,1000]%mapq[1,1000];clip[5,1000]%ins[1,1000];mapq[1,100]%del[1,1000];mapq[1,1000]";
  static std::string indel_mask = "";
  static std::string pon = "";
  static bool no_reads = false;
  static bool zip = true;
  static size_t numThreads = 1;
  static bool no_r2c = false;
  static std::string regionFile = "";
}


void runAssembly2VCF(int argc, char** argv)
{

  // parse the options
  parseAssembly2VCFOptions(argc, argv);

  // read in the assembly bam file
  AssemblyBamWalker awalk(opt::assembly_bam);

  std::cout << awalk << std::endl;
  
  awalk.walkDiscovar();

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
      default: die= true; 
    }
  }

  if (die) 
    {
      std::cout << "\n" << ASSEMBLY2VCF_USAGE_MESSAGE;
      exit(EXIT_FAILURE);
    }
}
