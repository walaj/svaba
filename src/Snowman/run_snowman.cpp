#include "run_snowman.h"

#include <getopt.h>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>

#include "vcf.h"
#include "bwa/bwa.h"

#include "SnowTools/SnowUtils.h"
#include "SnowTools/MiniRules.h"
#include "SnowTools/SnowToolsCommon.h"
#include "SnowTools/DBSnpFilter.h"

//#define NO_LOAD_INDEX 1

#define LITTLECHUNK 20000 
#define WINDOW_PAD 500
#define MICROBE_MATCH_MIN 50
//#define MATE_READ_LIMIT 3000
#define GET_MATES 1
#define MICROBE 1

typedef std::unordered_map<std::string, std::string> BamMap;

static faidx_t * findex;
static faidx_t * findex_viral;

static std::string POLYA = "AAAAAAAAAAAAAAAAAAA";
static std::string POLYT = "TTTTTTTTTTTTTTTTTTT";
static std::string POLYC = "CCCCCCCCCCCCCCCCCCC";
static std::string POLYG = "GGGGGGGGGGGGGGGGGGG";
static std::string POLYAT = "ATATATATATATATATATATATAT";
static std::string POLYCG = "CGCGCGCGCGCGCGCGCGCGCGCG";
static std::string POLYTG = "TGTGTGTGTGTGTGTGTGTGTGTG";
static std::string POLYCA = "CACACACACACACACACACACACA";
static std::string POLYAG = "AGAGAGAGAGAGAGAGAGAGAGAG";
static std::string POLYTC = "TCTCTCTCTCTCTCTCTCTCTCTC";

static SnowTools::BamWalker bwalker, r2c_writer, er_writer, b_microbe_writer, b_allwriter;
static SnowTools::BWAWrapper * microbe_bwa = nullptr;
static SnowTools::BWAWrapper * main_bwa = nullptr;
static SnowTools::MiniRulesCollection * mr;
static SnowTools::GRC blacklist, indel_blacklist_mask;
static SnowTools::DBSnpFilter * dbsnp_filter;

// output files
static ogzstream all_align, os_allbps, all_disc_stream, os_cigmap, r2c_c;
static ofstream log_file;

static SnowTools::GRC file_regions, regions_torun;

static pthread_mutex_t snow_lock;
static struct timespec start;

namespace opt {

  namespace assemb {
    static int minOverlap = -1;
    static float error_rate = 0.00; 
    static bool writeASQG = false;
  }

  static bool write_extracted_reads = false;

  //static int isize = 1000;
  static bool no_assemble_normal = false;
  static std::string indel_mask = ""; //"/xchip/gistic/Jeremiah/Projects/HengLiMask/um75-hs37d5.bed.gz";
  //static bool output_cov = false;
  static bool no_reads = true;
  static int32_t readlen;
  static bool r2c = false;
  static bool zip = false;
  //  static std::string pon = "";

  // parameters for filtering reads
  static std::string rules = "global@!hardclip;!supplementary;!duplicate;!qcfail;phred[4,100];length[50,1000]%region@WG%discordant[0,1200];mapq[0,1000]%mapq[0,1000];clip[5,1000];!isize[0,200]%ins[1,1000];mapq[1,100]%del[1,1000];mapq[1,1000]%mapped;!mate_mapped;mapq[1,1000]%mate_mapped;!mapped";

  static int max_cov = 100;
  static int chunk = LITTLECHUNK; // 1000000;

  // runtime parameters
  static int verbose = 1;
  static int numThreads = 1;

  // data
  static BamMap bam;
  static std::string refgenome = SnowTools::REFHG19;  
  static std::string microbegenome = "/xchip/gistic/Jeremiah/Projects/SnowmanFilters/viral.1.1.genomic_ns.fna";  
  static std::string analysis_id = "no_id";

  //subsample
  float subsample = 1.0;

  static std::string regionFile = "";
  static std::string blacklist; // = "/xchip/gistic/Jeremiah/Projects/HengLiMask/um75-hs37d5.bed.gz";

  static std::string dbsnp = ""; // = "/xchip/gistic/Jeremiah/SnowmanFilters/dbsnp_138.b37_indel.vcf";

  // filters on when / how to assemble
  static bool disc_cluster_only = false;

  static bool refilter = false;
  static bool adapter_trim = true;

  static std::string gemcode;
}

enum { 
  OPT_ASQG,
  OPT_DISC_CLUSTER_ONLY,
  OPT_READ_TRACK,
  OPT_NO_ASSEMBLE_NORMAL,
  OPT_MAX_COV,
  OPT_DISCORDANT_ONLY,
  OPT_REFILTER,
  OPT_WRITE_EXTRACTED_READS,
  OPT_ADAPTER_TRIM,
  OPT_GEMCODE_AWARE
};

static const char* shortopts = "hzxt:n:p:v:r:G:r:e:g:k:c:a:m:B:M:D:Y:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "tumor-bam",               required_argument, NULL, 't' },
  { "indel-mask",              required_argument, NULL, 'M' },
  //{ "panel-of-normals",        required_argument, NULL, 'q' },
  { "id-string",               required_argument, NULL, 'a' },
  { "normal-bam",              required_argument, NULL, 'n' },
  { "threads",                 required_argument, NULL, 'p' },
  { "chunk-size",              required_argument, NULL, 'c' },
  { "region-file",             required_argument, NULL, 'k' },
  { "rules",                   required_argument, NULL, 'r' },
  { "reference-genome",        required_argument, NULL, 'G' },
  { "microbial-genome",        required_argument, NULL, 'Y' },
  { "min-overlap",             required_argument, NULL, 'm' },
  { "dbsnp-vcf",               required_argument, NULL, 'D' },
  { "g-zip",                 no_argument, NULL, 'z' },
  { "read-tracking",         no_argument, NULL, OPT_READ_TRACK },
  { "write-extracted-reads", no_argument, NULL, OPT_WRITE_EXTRACTED_READS },
  { "no-assemble-normal",    no_argument, NULL, OPT_NO_ASSEMBLE_NORMAL },
  { "discordant-only",       no_argument, NULL, OPT_DISCORDANT_ONLY },
  { "refilter",              no_argument, NULL, OPT_REFILTER },
  { "r2c-bam",               no_argument, NULL, 'x' },
  { "write-asqg",            no_argument, NULL, OPT_ASQG   },
  { "error-rate",            required_argument, NULL, 'e'},
  { "verbose",               required_argument, NULL, 'v' },
  { "blacklist",             required_argument, NULL, 'B' },
  { "max-coverage",          required_argument, NULL, OPT_MAX_COV },
  { "no-adapter-trim",       no_argument, NULL, OPT_ADAPTER_TRIM },
  { "assemble-by-tag",       required_argument, NULL, OPT_GEMCODE_AWARE },
  { NULL, 0, NULL, 0 }
};

static const char *RUN_USAGE_MESSAGE =
"Usage: snowman run [OPTION] -t Tumor BAM\n\n"
"  Description: Grab weird reads from the BAM and perform assembly with SGA\n"
"\n"
"  General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h, --help                           Display this help and exit\n"
"  -p, --threads                        Use NUM threads to run snowman. Default: 1\n"
"  -a, --id-string                      String specifying the analysis ID to be used as part of ID common.\n"
"      --refilter                      Re-create the VCF files from the bps.txt.gz file. Need to supply correct id-string.\n"
"  Required input\n"
"  -G, --reference-genome               Path to indexed reference genome to be used by BWA-MEM. Default is Broad hg19 (/seq/reference/...)\n"
"  -t, --tumor-bam                      Tumor BAM file\n"
"  Optional input\n"                       
"  -n, --normal-bam                     Normal BAM file\n"
"  -r, --rules                          VariantBam style rules string to determine which reads to do assembly on. See documentation for default.\n"
"  -m, --min-overlap                    Minimum read overlap, an SGA parameter. Default: 0.4* readlength\n"
"  -k, --region-file                    Set a region txt file. Format: one region per line, Ex: 1,10000000,11000000\n"
  //"  -q, --panel-of-normals               Panel of normals gzipped txt file generated from snowman pon\n"
"  -M, --indel-mask                     BED-file with graylisted regions for stricter indel calling.\n"
"  -D, --dbsnp-vcf                      DBsnp database (VCF) to compare indels against\n"
"  -B, --blacklist                      BED-file with blacklisted regions to not extract any reads from.\n"
"  -z, --g-zip                          Gzip and tabix the output VCF files. Default: off\n"
"      --r2c-bam                        Output a BAM of reads that aligned to a contig, and fasta of kmer corrected sequences\n"
"      --discordant-only                Only run the discordant read clustering module, skip assembly. Default: off\n"
"      --read-tracking                  Track supporting reads. Increases file sizes.\n"
"      --max-coverage                   Maximum weird read coverage to send to assembler (per BAM). Subsample reads to achieve max if overflow. Defualt 500\n"
"  -V, --microbial-genome               Path to indexed reference genome of microbial sequences to be used by BWA-MEM to filter reads.\n"
"      --write-extracted-reads          For the tumor BAM, write the extracted reads (those sent to assembly) to a BAM file. Good for debugging.\n"
"      --no-adapter-trim                Don't peform Illumina adapter trimming, which removes reads with AGATCGGAAGAGC present.\n"
"      --assemble-by-tag                Separate the assemblies and read-to-contig mapping by the given read tag. Useful for 10X Genomics data (e.g. --aseembly-by-tag BX).\n"
"  Assembly params\n"
"      --write-asqg                     Output an ASQG graph file for each 5000bp window. Default: false\n"
"  -e, --error-rate                     Fractional difference two reads can have to overlap. See SGA param. 0 is fast, but requires exact. Default: 0.05\n"
"  -c, --chunk-size                     Amount of genome to read in at once. High numbers have fewer I/O rounds, but more memory. Default 1000000 (1M). Suggested 50000000 (50M) or 'chr' for exomes\n"
"\n";

void runSnowman(int argc, char** argv) {

  parseRunOptions(argc, argv);
  
  if (opt::refilter) {
    makeVCFs();
    return;
  }

  std::cerr << 
    "-----------------------------------------------------------------" << std::endl << 
    "--- Running Snowman somatic indel and rearrangement detection ---" << std::endl <<
    "-----------------------------------------------------------------" << std::endl;
  std::cerr << 
    "***************************** PARAMS ****************************" << std::endl << 
    "    DBSNP Database file: " << opt::dbsnp << std::endl << 
    "    Max cov to assemble: " << opt::max_cov << std::endl << 
    "    ErrorRate: " << (opt::assemb::error_rate < 0.001f ? "EXACT (0)" : std::to_string(opt::assemb::error_rate)) << std::endl << 
    "    Remove clipped reads with adapters? " << (opt::adapter_trim ? "TRUE" : "FALSE") << std::endl;
  if (opt::assemb::writeASQG)
    std::cerr << "    Writing ASQG files. Suggest running R/snow-asqg2pdf.R -i <my.asqg> -o graph.pdf" << std::endl;
  if (opt::disc_cluster_only)
    std::cerr << "    ######## ONLY DISCORDANT READ CLUSTERING. NO ASSEMBLY ##############" << std::endl;

  std::cerr <<
    "*****************************************************************" << std::endl;			  

  if (opt::disc_cluster_only) {
    static std::string rules = "global@nbases[0,0];!hardclip;!supplementary;!duplicate;!qcfail;%region@WG%discordant[0,800];mapq[1,1000]";
  }
  
#ifdef MICROBE
  // make the microbe BWAWrapper
  if (SnowTools::read_access_test(opt::microbegenome) && !opt::disc_cluster_only) {
    std::cerr << "...loading the microbial genome " << opt::microbegenome << std::endl;
    microbe_bwa = new SnowTools::BWAWrapper();
    microbe_bwa->retrieveIndex(opt::microbegenome);

    // open the microbe bam for writing  
    b_microbe_writer.SetWriteHeader(microbe_bwa->HeaderFromIndex());
    b_microbe_writer.OpenWriteBam(opt::analysis_id + ".microbe.bam"); // open and write header
  } else {
    microbe_bwa = 0;
  }
 
#endif

  findex = fai_load(opt::refgenome.c_str());  // load the reference
  if (opt::microbegenome.length())
    findex_viral = fai_load(opt::microbegenome.c_str());  // load the reference
  if (!opt::disc_cluster_only) {

    std::cerr << "...loading the human reference sequence" << std::endl;

#ifndef NO_LOAD_INDEX
    // make the BWAWrapper
    main_bwa = new SnowTools::BWAWrapper();
    main_bwa->retrieveIndex(opt::refgenome);
#endif
  }

  // open the tumor bam to get header info
  bwalker.OpenReadBam(opt::bam.begin()->first);

  // open the r2c writer
  if (opt::r2c) {
    bam_hdr_t * r2c_hdr = bam_hdr_dup(bwalker.header());
    r2c_writer.SetWriteHeader(r2c_hdr);
    r2c_writer.OpenWriteBam(opt::analysis_id + ".r2c.bam");
  }

  // open the extracted reads writer
  if (opt::write_extracted_reads) {
    bam_hdr_t * r2c_hdr = bam_hdr_dup(bwalker.header());
    er_writer.SetWriteHeader(r2c_hdr);
    er_writer.OpenWriteBam(opt::analysis_id + ".tumor.extracted.reads.bam");    
  }

#ifndef NO_LOAD_INDEX
  // open the contig bam for writing
  if (!opt::disc_cluster_only) {
    // open the all-contig bam for writing
    b_allwriter.SetWriteHeader(main_bwa->HeaderFromIndex());
    b_allwriter.OpenWriteBam(opt::analysis_id + ".contigs.all.bam"); // open and write header
  }
#endif

  bool chr_in_header = false; // is this header a chr1 header or 1 header
  // check if there is a "chr" marker in header
    for (int i = 0; i < bwalker.header()->n_targets; ++i) {
      if (bwalker.header()->target_name[i] && std::string(bwalker.header()->target_name[i]).find("chr") != std::string::npos) {
	chr_in_header = true;
	break;
      }
    }

  // open the blacklist
  if (opt::blacklist.length()) {

    std::cerr << "...reading blacklist from " << opt::blacklist << std::endl;
    blacklist.regionFileToGRV(opt::blacklist, 0, bwalker.header(), chr_in_header);
    blacklist.createTreeMap();
    std::cerr << "...read in " << blacklist.size() << " blacklist regions " << std::endl;
  }

  // open the DBSnpFilter
  if (opt::dbsnp.length()) {
    std::cerr << "...loading the DBsnp database" << std::endl;
    dbsnp_filter = new SnowTools::DBSnpFilter(opt::dbsnp);
    std::cerr << (*dbsnp_filter) << std::endl;
  }

  // set the MiniRules to be applied to each region
  mr = new SnowTools::MiniRulesCollection(opt::rules);
  if (opt::verbose > 1)
    std::cerr << *mr;

  // learn some parameters
  learnParameters(); 
  std::cerr << "...found read length of " << opt::readlen << ". Min Overlap is " << opt::assemb::minOverlap << std::endl;

  // parse the indel mask
  if (opt::indel_mask == opt::blacklist)
    indel_blacklist_mask = blacklist;
  else if (opt::indel_mask.length()) {
    std::cerr << "...loading the indel greylist mask" << std::endl;
    indel_blacklist_mask.regionFileToGRV(opt::indel_mask, 0, bwalker.header(), chr_in_header);
    indel_blacklist_mask.createTreeMap();
    std::cerr << "...read in " << blacklist.size() << " blacklist regions " << std::endl;
  }

  // parse the region file, count number of jobs
  int num_jobs = countJobs(file_regions, regions_torun); 
  if (num_jobs)
    std::cerr << "...running on " << num_jobs << " big-chunk regions with chunk size of " << SnowTools::AddCommas<int>(opt::chunk) << std::endl;
  else
    std::cerr << "Chunk was <= 0: READING IN WHOLE GENOME AT ONCE" << std::endl;
  if (opt::verbose > 1)
    for (auto& i : regions_torun)
      std::cerr << " REGION TO RUN " << i << std::endl;

  // override the number of threads if need
  num_jobs = (num_jobs == 0) ? 1 : num_jobs;
  opt::numThreads = std::min(num_jobs, opt::numThreads);

  // open the mutex
  if (pthread_mutex_init(&snow_lock, NULL) != 0) {
      printf("\n mutex init failed\n");
      return;
  }

  // open the files
  std::string n1 = opt::analysis_id + ".alignments.txt.gz";
  std::string n2 = opt::analysis_id + ".bps.txt.gz";
  //std::string ns = opt::analysis_id + ".bps.secondary.txt.gz";
  std::string n3 = opt::analysis_id + ".discordant.txt.gz";
  std::string n4 = opt::analysis_id + ".tumor.coverage.bedgraph.gz";
  std::string n5 = opt::analysis_id + ".normal.coverage.bedgraph.gz";
  std::string n6 = opt::analysis_id + ".cigarmap.txt.gz";
  std::string n7 = opt::analysis_id + ".r2c_corrected.fa";
  std::string n8 = opt::analysis_id + ".log";
  all_align.open(n1.c_str(), std::ios::out);
  os_allbps.open(n2.c_str(), std::ios::out);
  //os_allbps_secondary.open(ns.c_str(), std::ios::out);
  all_disc_stream.open(n3.c_str(), ios::out);
  os_cigmap.open(n6.c_str(), std::ios::out);
  if (opt::r2c) r2c_c.open(n7.c_str(), std::ios::out); 
  log_file.open(n8, std::ios::out);
  
  // write the headers
  os_allbps << SnowTools::BreakPoint::header() << endl;
  //os_allbps_secondary << SnowTools::BreakPoint::header() << endl;
  all_disc_stream << SnowTools::DiscordantCluster::header() << endl;

  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);

  // send the jobs to the queue
  std::cerr << "---- Starting detection pipeline --- on " << opt::numThreads << " thread" << std::endl;
  sendThreads(regions_torun);

  if (microbe_bwa)
    delete microbe_bwa;
  if (main_bwa)
    delete main_bwa;  

  // close the files
  all_align.close();
  os_allbps.close();
  //os_allbps_secondary.close();
  all_disc_stream.close();
  if (opt::r2c) r2c_c.close();
  os_cigmap.close();
  log_file.close();

  // make the VCF file
  makeVCFs();

  std::cerr << SnowTools::displayRuntime(start) << std::endl;
}

void makeVCFs() {

  if (!bwalker.header()) {
    // open the tumor bam to get header info
    bwalker.OpenReadBam(opt::bam.begin()->first);
  }
    

  // make the VCF file
  if (opt::verbose)
    std::cerr << "...loading the bps files for conversion to VCF" << std::endl;

  std::string file = opt::analysis_id + ".bps.txt.gz";
  if (!SnowTools::read_access_test(file))
    file = opt::analysis_id + ".bps.txt";

  // primary VCFs
  if (SnowTools::read_access_test(file)) {
    if (opt::verbose)
      std::cerr << "...making the primary VCFs (unfiltered and filtered) from file " << file << std::endl;
    VCFFile snowvcf(file, opt::refgenome.c_str(), opt::analysis_id, bwalker.header());
    std::string basename = opt::analysis_id + ".snowman.unfiltered.";
    snowvcf.include_nonpass = true;
    snowvcf.writeIndels(basename, opt::zip);
    snowvcf.writeSVs(basename, opt::zip);

    basename = opt::analysis_id + ".snowman.";
    snowvcf.include_nonpass = false;
    snowvcf.writeIndels(basename, opt::zip);
    snowvcf.writeSVs(basename, opt::zip);

  } else {
    std::cerr << "Failed to make VCF. Could not file bps file " << file << std::endl;
  }

}

// parse the command line options
void parseRunOptions(int argc, char** argv) {
  bool die = false;

  if (argc <= 2) 
    die = true;

  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case OPT_REFILTER : opt::refilter = true; break;
      case 'p': arg >> opt::numThreads; break;
      case 'm': arg >> opt::assemb::minOverlap; break;
      case 'a': arg >> opt::analysis_id; break;
      case 'B': arg >> opt::blacklist; break;
      case 'M': arg >> opt::indel_mask; break;
      case 'Y': arg >> opt::microbegenome; break;
	//case 'q': arg >> opt::pon; break;
      case 'z': opt::zip = false; break;
      case 'h': die = true; break;
      case 'x': opt::r2c = true; break;
      case 'c': 
	tmp = "";
	arg >> tmp;
	if (tmp.find("chr") != std::string::npos) {
	  opt::chunk = 250000000; break;
	} else {
	  opt::chunk = stoi(tmp); break;
	}
    case OPT_ASQG: opt::assemb::writeASQG = true; break;
    case OPT_ADAPTER_TRIM: opt::adapter_trim = false; break;
    case OPT_NO_ASSEMBLE_NORMAL: opt::no_assemble_normal = true; break;
    case OPT_MAX_COV: arg >> opt::max_cov;  break;
    case OPT_READ_TRACK: opt::no_reads = false; break;
	case 't': 
	  tmp = "";
	  arg >> tmp;
	  opt::bam[tmp] = "t"; 
	  break;
	case 'n': 
	  tmp = "";
	  arg >> tmp;
	  opt::bam[tmp] = "n"; 
	  break;
      case 'v': arg >> opt::verbose; break;
      case 'k': arg >> opt::regionFile; break;
      case 'e': arg >> opt::assemb::error_rate; break;
      case 'G': arg >> opt::refgenome; break;
	//case 'M': arg >> opt::microbegenome; break;
      case 'D': arg >> opt::dbsnp; break;
      case 'r': arg >> opt::rules; break;
      case OPT_DISCORDANT_ONLY: opt::disc_cluster_only = true; break;
      case OPT_WRITE_EXTRACTED_READS: opt::write_extracted_reads = true; break;
    case OPT_GEMCODE_AWARE: arg >> opt::gemcode; break;
      default: die= true; 
    }
  }

  // check that BAM files exist
  for (auto& b : opt::bam)
    if (!SnowTools::read_access_test(b.first)) {
      std::cerr << "Error: BAM file " << b.first << " is not readable / existant" << std::endl;
      exit(EXIT_FAILURE);
    }
      

  // check that we input something
  if (opt::bam.size() == 0 && !die && !opt::refilter) {
    std::cerr << "Must add a bam file " << std::endl;
    exit(EXIT_FAILURE);
  }

  if (opt::numThreads <= 0) {
    std::cerr << "run: invalid number of threads: " << opt::numThreads << std::endl;
    die = true;
  }

  if (die) 
    {
      std::cerr << "\n" << RUN_USAGE_MESSAGE;
      exit(EXIT_FAILURE);
    }
}

// just get a count of how many jobs to run. Useful for limiting threads. Also set the regions
int countJobs(SnowTools::GRC &file_regions, SnowTools::GRC &run_regions) {

  // open the region file if it exists
  bool rgfile = SnowTools::read_access_test(opt::regionFile);
  if (rgfile)
    file_regions.regionFileToGRV(opt::regionFile, 0);
  // parse as a samtools string eg 1:1,000,000-2,000,000
  else if (opt::regionFile.find(":") != std::string::npos && opt::regionFile.find("-") != std::string::npos)
    file_regions.add(SnowTools::GenomicRegion(opt::regionFile, bwalker.header()));
  // add all chromosomes
  else {
    for (int i = 0; i < bwalker.header()->n_targets; i++) {
      int region_id = bam_name2id(bwalker.header(), bwalker.header()->target_name[i]);
      
      if (opt::verbose > 1)
	std::cerr << "chr id from header " << region_id << " name " << bwalker.header()->target_name[i] << " len " << bwalker.header()->target_len[i] << std::endl;

      if (region_id < 23) // don't add outsdie of 1-X
	file_regions.add(SnowTools::GenomicRegion(region_id, 30000, bwalker.header()->target_len[i]));
    }
  }

  if (opt::verbose > 1)
    for (auto& i : file_regions)
      std::cerr << "file regions " << i << std::endl;

  // check if the mask was successful
  if (file_regions.size() == 0) {
    std::cerr << "ERROR: Cannot read region file: " << opt::regionFile << " or something wrong with tumor bam header" << std::endl;
    exit(EXIT_FAILURE);
  }

  // divide it up
  if (opt::chunk > 0) // if <= 0, whole genome at once
    for (auto& r : file_regions) {
      SnowTools::GRC test(opt::chunk, WINDOW_PAD, r);
      run_regions.concat(test);
    }
  return run_regions.size();
  
}

bool runBigChunk(const SnowTools::GenomicRegion& region)
{
  std::vector<SnowTools::AlignedContig> alc;
  BamReadVector all_contigs;
  BamReadVector all_microbial_contigs;

  // start some counters
  //int num_n_reads = 0, num_t_reads = 0;
  SnowTimer st;
  st.start();

  // setup for the BAM walkers
  std::unordered_map<std::string, SnowmanBamWalker> walkers;

  // loop through the input bams and get reads
  size_t tcount = 0;
  size_t ncount = 0;
  for (auto& b : opt::bam)
    {
      // opt::bam is <filename, type>
      walkers[b.first] = SnowmanBamWalker(b.first);
      walkers[b.first].adapter_trim = opt::adapter_trim;
      walkers[b.first].max_cov = opt::max_cov;
      walkers[b.first].disc_only = opt::disc_cluster_only;
      walkers[b.first].prefix = b.second;
      walkers[b.first].blacklist = blacklist;
      if (!region.isEmpty())
	walkers[b.first].setBamWalkerRegion(region);
      walkers[b.first].SetMiniRulesCollection(*mr);
	
      // read the reads
      walkers[b.first].readBam();
     
      if (b.second == "t") {
	tcount += walkers[b.first].reads.size();
      } else {
	ncount += walkers[b.first].reads.size();
      }
    }

  size_t num_t_reads = tcount;
  size_t num_n_reads = ncount;

  if (opt::verbose > 1)
    std::cerr << "...read in " << tcount << "/" << ncount << " Tumor/Normal reads from "  << std::endl;

  st.stop("r");

  // collect all of the cigar maps
  CigarMap cigmap_n, cigmap_t;
  for (auto& b : opt::bam) {
    if (b.second.at(0) == 'n')
      for (auto& c : walkers[b.first].cigmap)
	cigmap_n[c.first] += c.second;
    else if (b.second.at(0) == 't')
      for (auto& c : walkers[b.first].cigmap)
	cigmap_t[c.first] += c.second;
  }

#ifdef GET_MATES
  // collect the normal mate regions
  MateRegionVector normal_mate_regions;
  for (auto& b : opt::bam)
    if (b.second.at(0) == 'n')
      normal_mate_regions.concat(walkers[b.first].mate_regions);
  normal_mate_regions.createTreeMap();

  // get the mates from somatic 3+ mate regions
  // These will be regions to grab exta 
  size_t LOOKUP_LIM = opt::disc_cluster_only ? 2 : 3;
  MateRegionVector somatic_mate_regions;
  for (auto& b : opt::bam)
    if (b.second.at(0) == 't')
      for (auto& i : walkers[b.first].mate_regions){
	//std::cerr << "...checking somatic mate region " << i << std::endl;
	//std::cerr << "   icount " << i.count << " i.chr " << i.chr << " indel_blacklist_mask.size() == 0 || " << 
	//  " indel_blacklist_mask.findOverlapping(i) == 0) :: " <<  
	//  (indel_blacklist_mask.size() == 0 || indel_blacklist_mask.findOverlapping(i) == 0) << std::endl;
	if (i.count >= LOOKUP_LIM && !normal_mate_regions.findOverlapping(i) && 
	    i.chr < 24 &&
	    (blacklist.size() == 0 || blacklist.findOverlapping(i) == 0))
	  somatic_mate_regions.add(i); //concat(walkers[b.first].mate_regions);
      }
  
  // reduce it down
  somatic_mate_regions.mergeOverlappingIntervals();

  //if (opt::verbose > 1)
  for (auto& i : somatic_mate_regions)
    log_file << "   somatic mate region " << i << " somatic count " << i.count << " origin " << i.partner << std::endl;

  if (opt::verbose > 3)
    std::cerr << "...grabbing mate reads from " << somatic_mate_regions.size() << " regions spanning " << SnowTools::AddCommas(somatic_mate_regions.width()) << std::endl;

  // add these regions to the walker and get the reads
  tcount = 0; ncount = 0;
  if (somatic_mate_regions.size())
  for (auto& b : opt::bam)
    {
      int oreads = walkers[b.first].reads.size();
      walkers[b.first].setBamWalkerRegions(somatic_mate_regions.asGenomicRegionVector());
      walkers[b.first].get_coverage = false;
      walkers[b.first].disc_only = opt::disc_cluster_only;
      walkers[b.first].max_cov = opt::max_cov;
      walkers[b.first].get_mate_regions = false;
      walkers[b.first].readBam();
      if (b.second == "t") {
	tcount += (walkers[b.first].reads.size() - oreads);
      } else {
	ncount += (walkers[b.first].reads.size() - oreads);
      }
    }

  num_t_reads += tcount;
  num_n_reads += ncount;

  if (somatic_mate_regions.size() && opt::verbose > 1)
    std::cerr << "           " << tcount << "/" << ncount << " Tumor/Normal mate reads in " << somatic_mate_regions.size() << " regions spanning " << SnowTools::AddCommas<int>(somatic_mate_regions.width()) << " bases" << std::endl;

  st.stop("m");
#endif

  // put all of the reads together and dedupe
  std::set<std::string> dedup;
  BamReadVector bav_this;
  for (auto& b : walkers)
    for (auto& r : b.second.reads)
      if (!dedup.count(r.GetZTag("SR"))) {
	dedup.insert(r.GetZTag("SR"));
	bav_this.push_back(r);
      }


  // do the kmer filtering
  if (!opt::disc_cluster_only) {
    KmerFilter kmer;
    kmer.correctReads(bav_this);
    st.stop("k");
  }
  
  // do the discordant read clustering
  if (opt::verbose > 1)
    std::cerr << "...doing the discordant read clustering" << std::endl;
  SnowTools::DiscordantClusterMap dmap = SnowTools::DiscordantCluster::clusterReads(bav_this, region);
  
  // 
  if (opt::verbose > 3)
    for (auto& i : dmap)
      std::cerr << i.first << " " << i.second << std::endl;

  GRC grv_small = makeAssemblyRegions(region);

  if (opt::verbose > 1 && !opt::disc_cluster_only)
    std::cerr << "running the assemblies for region " << region <<  std::endl;

  if (!opt::disc_cluster_only) {
    if (opt::verbose > 1)
      std::cerr << "...doing the assemblies" << std::endl;
    for (auto& g : grv_small) 
      {
	// set the contig uid
      std::string name = "c_" + std::to_string(g.chr+1) + "_" + std::to_string(g.pos1) + "_" + std::to_string(g.pos2);

      // check that we don't have too many reads
      if (bav_this.size() > (size_t)(region.width() * 20) && region.width() > 20000) {
	log_file << bav_this.size() << "\t" << g << std::endl;
	continue;
      }

      // print message about assemblies
      if (opt::verbose > 1 && bav_this.size() > 1)
	std::cerr << "Doing assemblies on " << name << std::endl;
      else if (opt::verbose > 1)
	std::cerr << "Skipping assembly (< 2 reads) for " << name << std::endl;
      if (bav_this.size() < 2)
	continue;

      ContigVector all_contigs_this;

      // 10X test
      if (opt::gemcode.length()) {

	std::unordered_map<std::string, SnowmanAssemblerEngine> bx_engines;

	// sort into gemcodes
	std::unordered_map<std::string, SnowTools::BamReadVector> bx_reads;
	for (auto& rr : bav_this) {
	  std::string gcode; rr.GetZTag(opt::gemcode);
	  if (gcode.empty())
	    gcode = "Default";
	  bx_reads[gcode].push_back(rr);
	}

	// make assembler engines
	for (auto& brv : bx_reads) {
	  bx_engines[brv.first] = SnowmanAssemblerEngine(name + "BX_" + brv.first, opt::assemb::error_rate, opt::assemb::minOverlap, opt::readlen);	  
	  std::unordered_map<std::string, SnowmanAssemblerEngine>::iterator en = bx_engines.find(brv.first);
	  en->second.fillReadTable(brv.second);
	  en->second.performAssembly();
	  all_contigs_this.insert(all_contigs_this.end(), en->second.getContigs().begin(), en->second.getContigs().end());
	}
	
      // normal non-10X data
      } else {

	SnowmanAssemblerEngine engine(name, opt::assemb::error_rate, opt::assemb::minOverlap, opt::readlen);
	// do the assemblies
	if (opt::assemb::writeASQG)
	  engine.setToWriteASQG();
	engine.fillReadTable(bav_this);
	
	engine.performAssembly();
	
	all_contigs_this = engine.getContigs();
	
	if (opt::verbose > 1)
	  std::cerr << "Assembled " << engine.getContigs().size() << " contigs for " << name << std::endl; 
      }
      
      std::vector<SnowTools::AlignedContig> this_alc;
      
      // align the contigs to the genome
      if (opt::verbose > 1)
	std::cerr << "...aligning contigs to genome and reads to contigs" << std::endl;
      SnowTools::USeqVector usv;
      
#ifndef NO_LOAD_INDEX
      for (auto& i : all_contigs_this) {
	BamReadVector ct_alignments;
	main_bwa->alignSingleSequence(i.getSeq(), i.getID(), ct_alignments, 0.90);
	
	if (opt::verbose > 3)
	  for (auto& i : ct_alignments)
	    std::cerr << " aligned contig: " << i << std::endl;
	
#ifdef MICROBE
	
	BamReadVector ct_plus_microbe;

	if (microbe_bwa && !hasRepeat(i.getSeq())) {

	  // do the microbial alignment
	  BamReadVector microbial_alignments;
	  microbe_bwa->alignSingleSequence(i.getSeq(), i.getID(), microbial_alignments, 0.90);

	  // if the microbe alignment is large enough and doesn't overlap human...
	  for (auto& j : microbial_alignments) {
	    // keep only long microbe alignments with decent mapq
	    if (j.NumMatchBases() >= MICROBE_MATCH_MIN && j.MapQuality() >= 10) { 
	      if (overlapSize(j, ct_alignments) <= 20) { // keep only those where most do not overlap human
		assert(microbe_bwa->ChrIDToName(j.ChrID()).length());
		j.AddZTag("MC", microbe_bwa->ChrIDToName(j.ChrID()));
		all_microbial_contigs.push_back(j);
		ct_plus_microbe.push_back(j);
	      }
	    }
	  }
	}
	
#endif
	// add in the chrosome name tag for human alignments
	if (main_bwa)
	  for (auto& r : ct_alignments) {
	    assert(main_bwa->ChrIDToName(r.ChrID()).length());
	    r.AddZTag("MC", main_bwa->ChrIDToName(r.ChrID()));
	  }
	
	// remove human alignments that are not as good as microbe
	// that is, remove human alignments that intersect with microbe. 
	// We can do this because we already removed microbial alignments taht 
	// intersect too much with human. Thus, we are effectively removing human 
	// alignments that are contained within a microbial alignment

	BamReadVector human_alignments;
	for (auto& j : ct_alignments) {
	  // keep human alignments that have < 50% overlap with microbe and have >= 25 bases matched
	  if (overlapSize(j, ct_plus_microbe) < 0.5 * j.Length() && j.NumMatchBases() >= 25) { 
	    human_alignments.push_back(j);
	    all_contigs.push_back(j);
	  }
	}
	// add in the microbe alignments
	human_alignments.insert(human_alignments.end(), ct_plus_microbe.begin(), ct_plus_microbe.end());
	
	// make the aligned contig object
	if (!human_alignments.size())
	  continue;

	SnowTools::AlignedContig ac(human_alignments);

	// assign the local variable to each
	ac.checkLocal(g);

	if (ac.hasLocal()/* && ac.hasVariant()*/) {
	  this_alc.push_back(ac);
	  usv.push_back({i.getID(), i.getSeq()});	  
	}
      }
#else 
      for (auto& i : engine.getContigs()) 
	usv.push_back({i.getID(), i.getSeq()});
#endif

      if (!usv.size())
	continue;

      // Align the reads to the contigs with BWA-MEM
      SnowTools::BWAWrapper bw;
      bw.constructIndex(usv);

      if (opt::verbose > 3)
	std::cerr << "...aligning " << bav_this.size() << " reads to " << this_alc.size() << " contigs " << std::endl;
      alignReadsToContigs(bw, usv, bav_this, this_alc);

      // Get contig coverage, discordant matching to contigs, etc
      for (auto& a : this_alc) {
	a.assignSupportCoverage();
	a.splitCoverage();	
	// add discordant reads support to each of the breakpoints
	a.addDiscordantCluster(dmap);
	// add in the cigar matches
	a.checkAgainstCigarMatches(cigmap_n, cigmap_t);
	// check the indel breakpoints against the indel blacklist. 
	// simply set the blacklist flag for these breaks if they hit
	a.blacklist(indel_blacklist_mask);
	// add to the final structure
	alc.push_back(a);
      }
      
    } // stop region loop
  }

  st.stop("as");
  if (opt::verbose > 1)
    std::cerr << "...done assembling, post processing" << std::endl;

  // get the breakpoints
  std::vector<SnowTools::BreakPoint> bp_glob;
  
  if (opt::verbose > 1)
    std::cerr << "...getting the breakpoints" << std::endl;
  for (auto& i : alc) {
    //if (i.m_bamreads.size() && !i.m_skip) { //m_bamreads.size() is zero for contigs with ...
    std::vector<SnowTools::BreakPoint> allbreaks = i.getAllBreakPoints();
    //std::vector<SnowTools::BreakPoint> allbreaks_2 = i.getAllBreakPointsSecondary();
    bp_glob.insert(bp_glob.end(), allbreaks.begin(), allbreaks.end());
    //bp_glob_secondary.insert(bp_glob_secondary.end(), allbreaks_2.begin(), allbreaks_2.end());
    //}
  }

  if (opt::verbose > 1)
    std::cerr << "...repeat sequence filtering" << std::endl;
  
  // repeat sequence filter
  /*if (findex && false) {
    for (auto& i : bp_glob)
      i.repeatFilter(findex);
    for (auto& i : bp_glob_secondary)
      i.repeatFilter(findex);
      }*/
  
  if (dbsnp_filter && opt::dbsnp.length()) {
    if (opt::verbose > 1)
      std::cerr << "...DBSNP filtering" << std::endl;
    for (auto & i : bp_glob) {
      dbsnp_filter->queryBreakpoint(i);
    }
  }

  // add in the discordant clusters as breakpoints
  for (auto& i : dmap) {
    // DiscordantCluster not associated with assembly BP and has 2+ read support
    if (!i.second.hasAssociatedAssemblyContig() && (i.second.tcount + i.second.ncount) > 1) {
      SnowTools::BreakPoint tmpbp(i.second, main_bwa);
      assert(tmpbp.b1.gr < tmpbp.b2.gr);
      bp_glob.push_back(tmpbp);
    }
  }
  
  if (opt::verbose > 1)
    std::cerr << "...deduplicating breakpoints"<< std::endl;

  // de duplicate the breakpoints
  std::sort(bp_glob.begin(), bp_glob.end());
  SnowTools::BPVec tmprr = bp_glob; //debug
  //std::cerr << "..num breaks before " << bp_glob.size() << std::endl;
  bp_glob.erase( std::unique( bp_glob.begin(), bp_glob.end() ), bp_glob.end() );
  //std::cerr << "..num breaks after " << bp_glob.size() << std::endl;
  
  //debug
  //if (tmprr.size() != bp_glob.size()) {
  //  for (auto& o : bp_glob)
  //    std::cerr << " after " << o << std::endl;
  //  for (auto& o : tmprr)
  //    std::cerr << " before " << o << std::endl;
  //}
  //std::sort(bp_glob_secondary.begin(), bp_glob_secondary.end());
  //bp_glob_secondary.erase( std::unique( bp_glob_secondary.begin(), bp_glob_secondary.end() ), bp_glob_secondary.end() );
  
  // add the coverage data to breaks for allelic fraction computation
  SnowTools::STCoverage * t_cov = nullptr;
  SnowTools::STCoverage * n_cov = nullptr;
  for (auto& b : opt::bam) {
    if (b.second.at(0) == 't')
      t_cov = &walkers[b.first].cov;
    else 
      n_cov = &walkers[b.first].cov;
  }
  for (auto& i : bp_glob)
    i.addAllelicFraction(t_cov, n_cov);
 
  // score them and set somatic / germline
  for (auto& i : bp_glob)
    i.scoreBreakpoint();


  // add teh ref and alt tags
  for (auto& i : bp_glob)
    i.setRefAlt(findex, findex_viral);
 

  ////////////////////////////////////
  // MUTEX LOCKED
  ////////////////////////////////////
  // write to the global contig out
  pthread_mutex_lock(&snow_lock);  
  
  // dump the cigmaps
  for (auto& i : cigmap_n) 
    os_cigmap << i.first << "\t" << i.second << "\tN" << endl;
  for (auto& i : cigmap_t) 
    os_cigmap << i.first << "\t" << i.second << "\tT" << endl;
  
  // print the alignment plots
  size_t contig_counter = 0;
  for (auto& i : alc)
    if (i.hasVariant()) {
      all_align << i << std::endl;
      ++contig_counter;
    }
  
  // send the microbe to file
  for (auto& b : all_microbial_contigs)
    b_microbe_writer.WriteAlignment(b);

  // send the discordant to file
  for (auto& i : dmap)
    if (std::max(i.second.mapq1, i.second.mapq2) >= 5)
      all_disc_stream << i.second.toFileString(!opt::no_reads) << std::endl;

  // send breakpoints to file
  for (auto& i : bp_glob)
    if (i.hasMinimal() || true) {
      os_allbps << i.toFileString(opt::no_reads) << std::endl;
      if (i.confidence == "PASS" && i.nsplit == 0 && i.ncigar == 0 && i.dc.ncount == 0) {
	std::cout << "SOM  " << i.toPrintString() << std::endl;
      } else if (i.confidence == "PASS") {
	std::cout << "GER " << i.toPrintString() << std::endl;
      }
    }

  //for (auto& i : bp_glob_secondary)
  //  os_allbps_secondary << i.toFileString(opt::no_reads) << std::endl;
  
#ifndef NO_LOAD_INDEX
  // write ALL contigs
  if (!opt::disc_cluster_only)  
  for (auto& i : all_contigs)
    b_allwriter.WriteAlignment(i);
#endif
  
  // write extracted reads
  if (opt::write_extracted_reads) 
    for (auto& b : opt::bam)
      if (b.second == "t")
	for (auto& r : walkers[b.first].reads) 
	  er_writer.WriteAlignment(r);
  
  // write all the to-assembly reads
  if (!opt::disc_cluster_only && opt::r2c)
    for (auto& b : opt::bam) {
      if (b.second == "t")
	for (auto& r : walkers[b.first].reads) {
	  r2c_writer.WriteAlignment(r);
	  
	//write the corrected
	  std::string new_seq  = r.GetZTag("KC");
	  if (!new_seq.length()) {
	    //new_seq = r.QualityTrimmedSequence(4, dum);
	    new_seq = r.QualitySequence();
	  }
	  r2c_c << ">" << r.GetZTag("SR") << std::endl << new_seq << std::endl;
	}
    }
  
  st.stop("pp");
  // display the run time
  if (opt::verbose > 0 && (num_t_reads + num_n_reads) > 0 && !region.isEmpty()) {
    string print1 = SnowTools::AddCommas<int>(region.pos1);
    string print2 = SnowTools::AddCommas<int>(region.pos2);
    char buffer[180];
    sprintf (buffer, "Ran %2s:%11s-%11s | T: %5d N: %5d C: %5d | ", 
	     bwalker.header()->target_name[region.chr],
	     print1.c_str(),print2.c_str(),
	     (int)num_t_reads, (int)num_n_reads, 
	     (int)contig_counter);
    
    std::cerr << string(buffer) << st << " | ";
    std::cerr << SnowTools::displayRuntime(start);
    std::cerr << endl;
  } else if (opt::verbose > 0 && (num_t_reads + num_n_reads) > 0) {
    char buffer[180];
    sprintf (buffer, "Ran Whole Genome | T: %5d N: %5d C: %5d | ", 
	     (int)num_t_reads, (int)num_n_reads, 
	     (int)contig_counter);
    
    std::cerr << string(buffer) << st << " | ";
    std::cerr << SnowTools::displayRuntime(start);
    std::cerr << endl;

  }

  ////////////////////////////////////
  // MUTEX UNLOCKED
  ////////////////////////////////////
  pthread_mutex_unlock(&snow_lock);
  
  return true;
}

void sendThreads(SnowTools::GRC& regions_torun) {

  // Create the queue and consumer (worker) threads
  wqueue<SnowmanWorkItem*>  queue;
  std::vector<ConsumerThread<SnowmanWorkItem>*> threadqueue;
  for (int i = 0; i < opt::numThreads; i++) {
    ConsumerThread<SnowmanWorkItem>* threadr = new ConsumerThread<SnowmanWorkItem>(queue, opt::verbose > 0);
    threadr->start();
    threadqueue.push_back(threadr);
  }

  // send the jobs
  size_t count = 0;
  for (auto& i : regions_torun) {
    SnowmanWorkItem * item     = new SnowmanWorkItem(SnowTools::GenomicRegion(i.chr, i.pos1, i.pos2), ++count);
    queue.add(item);
  }
  if (regions_torun.size() == 0) { // whole genome 
    SnowmanWorkItem * item     = new SnowmanWorkItem(SnowTools::GenomicRegion(), ++count);
    queue.add(item);
  }
    
  
  // wait for the threads to finish
  for (int i = 0; i < opt::numThreads; i++) 
    threadqueue[i]->join();
  
}


void learnParameters() {

  SnowTools::BamRead r;
  size_t count = 0;
  
  bool rule;
  while (bwalker.GetNextRead(r, rule) && ++count < 1000) 
    opt::readlen = std::max(r.Length(), opt::readlen);
  
  if (opt::assemb::minOverlap < 0)
    opt::assemb::minOverlap = 0.4 * opt::readlen;
}

GRC makeAssemblyRegions(const SnowTools::GenomicRegion& region) {

  // set the regions to run
  GRC grv_small;
  if (region.isEmpty()) {  // whole genome, so divide up the whole thing
    for (size_t c = 0; c < 23; ++c)
      grv_small.concat(GRC(/*LITTLECHUNK*/opt::chunk, WINDOW_PAD, SnowTools::GenomicRegion(c, WINDOW_PAD + 1, SnowTools::CHR_LEN[c] - WINDOW_PAD - 1)));
  }
  else if (region.width() >= /*LITTLECHUNK*/opt::chunk) // divide into smaller chunks
    grv_small = GRC(/*LITTLECHUNK*/opt::chunk, WINDOW_PAD, region);
  else
    grv_small.add(region);

  //if (region.isEmpty() && file_regions.size()) { // if whole genome, restrict to just regions in regions file
  //  GRC test = grv_small.intersection(file_regions);
  //  for (auto& i : test)
  //    std::cerr << i << std::endl;
  //}
  return grv_small;

}

 void alignReadsToContigs(SnowTools::BWAWrapper& bw, const SnowTools::USeqVector& usv, BamReadVector& bav_this, std::vector<SnowTools::AlignedContig>& this_alc) {

   if (!usv.size())
     return;

  for (auto i : bav_this) {

    BamReadVector brv;
    std::string seqr = i.GetZTag("KC");
    if (!seqr.length())
      seqr = i.QualitySequence();

    assert(seqr.length());
    bw.alignSingleSequence(seqr, i.Qname(), brv, 0.80);

    if (brv.size() == 0) 
      continue;

    // make sure we have only one alignment per contig
    std::set<std::string> cc;
    
    // check which ones pass
    BamReadVector bpass;
    for (auto& r : brv) {

      // make sure alignment score is OK
      if (r.NumMatchBases() * 0.9 > r.GetIntTag("AS")/* && i.GetZTag("SR").at(0) == 't'*/)
	continue;
      
      bool length_pass = (r.PositionEnd() - r.Position()) >= (seqr.length() * 0.95);
      bool mapq_pass = r.MapQuality();
      int ins_bases = r.MaxInsertionBases();
      int del_bases = r.MaxDeletionBases();

      ins_bases = 0;
      del_bases = 0;

      if (length_pass && mapq_pass && ins_bases == 0 && del_bases == 0 && !cc.count(usv[r.ChrID()].name)) {
	bpass.push_back(r);
	cc.insert(usv[r.ChrID()].name);
      }
    }
    
    // annotate the original read
    for (auto& r : bpass) {
      if (r.ReverseFlag())
	i.SmartAddTag("RC","1");
      else 
	i.SmartAddTag("RC","0");

      i.SmartAddTag("SL", std::to_string(r.Position()));
      i.SmartAddTag("SE", std::to_string(r.PositionEnd()));
      i.SmartAddTag("TS", std::to_string(r.AlignmentPosition()));
      i.SmartAddTag("TE", std::to_string(r.AlignmentEndPosition()));
      i.SmartAddTag("SC", r.CigarString());
      i.SmartAddTag("CN", usv[r.ChrID()].name/*getContigName()*/);

      for (auto& a : this_alc) {
	if (a.getContigName() == usv[r.ChrID()].name) {
	  
	  a.m_bamreads.push_back(i);
	  
	  // add the coverage to the aligned contig
	  int cc = r.Position();
	  std::string srr = i.GetZTag("SR");
	  if (srr.length()) {
	    if (srr.at(0) == 't')
	      while (cc <= r.PositionEnd() && cc < (int)a.tum_cov.size())
		++a.tum_cov[cc++];
	    else
	      while (cc <= r.PositionEnd() && cc < (int)a.norm_cov.size())
		++a.norm_cov[cc++];
	  }
	}
      }	  

    } // end passing bwa-aligned read loop 
  } // end main read loop
}

bool hasRepeat(const std::string& seq) {
  
  if ((seq.find(POLYT) == std::string::npos) && 
      (seq.find(POLYA) == std::string::npos) && 
      (seq.find(POLYC) == std::string::npos) && 
      (seq.find(POLYG) == std::string::npos) && 
      (seq.find(POLYCG) == std::string::npos) && 
      (seq.find(POLYAT) == std::string::npos) && 
      (seq.find(POLYTC) == std::string::npos) && 
      (seq.find(POLYAG) == std::string::npos) && 
      (seq.find(POLYCA) == std::string::npos) && 
      (seq.find(POLYTG) == std::string::npos))
    return false;
  return true;
  
}

int overlapSize(const BamRead& query, const BamReadVector& subject) {

  // get the amount covered by subjet sequences
  typedef std::pair<int32_t, int32_t> AP;
  typedef std::pair<AP, std::string> AP_wseq;
  std::vector<AP_wseq> align_pos;
  std::string convention_seq = "";
  if (subject.size()) { 
    convention_seq = subject[0].Sequence(); // set the orientation convention
    for(auto& j : subject) {
      if (j.Sequence() == convention_seq)
	align_pos.push_back(AP_wseq(AP(j.AlignmentPosition(),j.AlignmentEndPosition()), j.Sequence()));
      else
	align_pos.push_back(AP_wseq(AP(j.AlignmentPositionReverse(),j.AlignmentEndPositionReverse()), j.Sequence()));
    }
  }
  
  bool same_orientation = query.Sequence() == convention_seq; 
  int al1 = same_orientation ? query.AlignmentPosition() : query.AlignmentPositionReverse();
  int al2 = same_orientation ? query.AlignmentEndPosition() : query.AlignmentEndPositionReverse();
  
  int max_overlap = 0;
  for (auto& k : align_pos) { // check each subject alignment to see if it overlaps query
    AP kk = k.first;
    max_overlap = std::max(max_overlap, std::max(0, std::min(kk.second, al2) - std::max(kk.first,al1)));
    //  out = true;
    //if (al1 >= kk.first && al2 <= kk.second)
    //  query_contained = true;
  }
  return max_overlap;
  
}
