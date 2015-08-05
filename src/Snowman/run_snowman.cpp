#include "run_snowman.h"

#include <getopt.h>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>

#include "vcf.h"
#include "bwa/bwa.h"

#include "SnowTools/SnowUtils.h"
#include "SnowTools/BWAWrapper.h"
#include "SnowTools/MiniRules.h"
#include "SnowmanBamWalker.h"
#include "SnowmanAssemblerEngine.h"
#include "SnowTools/AlignedContig.h"
#include "SnowTools/DiscordantCluster.h"
#include "SnowTools/DBSnpFilter.h"

#define LITTLECHUNK 6000 
#define WINDOW_PAD 300
#define MATE_READ_LIMIT 3000
#define GET_MATES 1
#define MICROBE 1

typedef std::unordered_map<std::string, std::string> BamMap;

static faidx_t * findex;

static SnowTools::BamWalker bwalker;
static SnowTools::BamWalker bwriter;
static SnowTools::BamWalker r2c_writer;
static SnowTools::BamWalker r2c_corrected_writer;
static SnowTools::BamWalker b_microbe_writer;
static SnowTools::BamWalker b_allwriter;
static SnowTools::BWAWrapper * microbe_bwa;
static SnowTools::BWAWrapper * main_bwa;
static SnowTools::MiniRulesCollection * mr;
static SnowTools::GRC blacklist;
static SnowTools::GRC indel_blacklist_mask;
static SnowTools::DBSnpFilter * dbsnp_filter;

// output files
static ogzstream * all_align;
static ogzstream * all_tum_cov;
static ogzstream * all_norm_cov;
static ogzstream * os_allbps;
static ogzstream * all_disc_stream;
static ogzstream * os_cigmap; 

static ofstream * mates_file;

static pthread_mutex_t snow_lock;
static struct timespec start;

namespace opt {

  namespace assemb {
    static unsigned minOverlap = 35;
    static float error_rate = 0.00; 
    static bool writeASQG = false;
  }

  //static int isize = 1000;
  static bool no_assemble_normal = false;
  static std::string indel_mask = ""; //"/xchip/gistic/Jeremiah/Projects/HengLiMask/um75-hs37d5.bed.gz";
  //static bool output_cov = false;
  static bool no_reads = true;
  static int32_t readlen;
  static bool no_r2c = false;
  static bool zip = false;
  static std::string pon = "";

  // parameters for filtering reads
  static std::string rules = "global@nbases[0,0];!hardclip;!supplementary;!duplicate;!qcfail;phred[4,100];%region@WG%discordant[0,1000];mapq[1,1000]%mapq[1,1000];clip[5,1000]%ins[1,1000];mapq[1,100]%del[1,1000];mapq[1,1000]";

  static int max_cov = 500;
  static int chunk = 1000000;

  // runtime parameters
  static int verbose = 1;
  static int numThreads = 1;

  // data
  static BamMap bam;
  static std::string refgenome = SnowTools::REFHG19;  
  static std::string microbegenome = "/xchip/gistic/Jeremiah/Projects/SnowmanFilters/viral.1.1.genomic.fna";  
  static std::string analysis_id = "no_id";

  //subsample
  float subsample = 1.0;

  static std::string regionFile = "";
  static std::string blacklist = ""; //"/xchip/gistic/Jeremiah/Projects/HengLiMask/um75-hs37d5.bed.gz";

  static std::string dbsnp = ""; // /xchip/gistic/Jeremiah/SnowmanFilters/dbsnp_138.b37_indel.vcf

  // filters on when / how to assemble
  static bool disc_cluster_only = false;

}

enum { 
  OPT_ASQG,
  OPT_DISC_CLUSTER_ONLY,
  OPT_READ_TRACK,
  OPT_NO_ASSEMBLE_NORMAL,
  OPT_MAX_COV
};

static const char* shortopts = "hzxt:n:p:v:r:G:r:e:g:k:c:a:q:m:b:M:D:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "tumor-bam",               required_argument, NULL, 't' },
  { "indel-mask",              required_argument, NULL, 'm' },
  { "panel-of-normals",        required_argument, NULL, 'q' },
  { "id-string",               required_argument, NULL, 'a' },
  { "normal-bam",              required_argument, NULL, 'n' },
  { "threads",                 required_argument, NULL, 'p' },
  { "chunk-size",              required_argument, NULL, 'c' },
  { "region-file",             required_argument, NULL, 'k' },
  { "rules",                   required_argument, NULL, 'r' },
  { "reference-genome",        required_argument, NULL, 'G' },
  { "microbial-genome",        required_argument, NULL, 'M' },
  { "dbsnp-vcf",               required_argument, NULL, 'D' },
  { "g-zip",                  no_argument, NULL, 'z' },
  { "read-tracking",           no_argument, NULL, OPT_READ_TRACK },
  { "no-assemble-normal",       no_argument, NULL, OPT_NO_ASSEMBLE_NORMAL },
  { "no-r2c-bam",              no_argument, NULL, 'x' },
  { "write-asqg",              no_argument, NULL, OPT_ASQG   },
  { "error-rate",              required_argument, NULL, 'e'},
  { "verbose",                 required_argument, NULL, 'v' },
  { "blacklist",                 required_argument, NULL, 'b' },
  { "max-coverage",                 required_argument, NULL, OPT_MAX_COV },
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
"  Required input\n"
"  -G, --reference-genome               Path to indexed reference genome to be used by BWA-MEM. Default is Broad hg19 (/seq/reference/...)\n"
"  -t, --tumor-bam                      Tumor BAM file\n"
"  Optional input\n"                       
"  -n, --normal-bam                     Normal BAM file\n"
"  -r, --rules                          VariantBam style rules string to determine which reads to do assembly on. See documentation for default.\n"
"  -m, --min-overlap                    Minimum read overlap, an SGA parameter. Default: 0.4* readlength\n"
"  -k, --region-file                    Set a region txt file. Format: one region per line, Ex: 1,10000000,11000000\n"
"  -q, --panel-of-normals               Panel of normals gzipped txt file generated from snowman pon\n"
"  -m, --indel-mask                     BED-file with blacklisted regions for indel calling. Default /xchip/gistic/Jeremiah/Projects/HengLiMask/um75-hs37d5.bed.gz\n"
"  -b, --blacklist                      BED-file with blacklisted regions to not extract any reads from. Default /xchip/gistic/Jeremiah/Projects/HengLiMask/um75-hs37d5.bed.gz\n"
"  -z, --g-zip                          Gzip and tabix the output VCF files. Default: off\n"
"      --disc-cluster-only              Only run the discordant read clustering module, skip assembly. Default: off\n"
"      --read-tracking                  Track supporting reads. Increases file sizes.\n"
"      --max-coverage                   Maximum weird read coverage to send to assembler (per BAM). Subsample reads to achieve max if overflow. Defualt 500\n"
"  -M, --microbial-genome               Path to indexed reference genome of microbial sequences to be used by BWA-MEM to filter reads.\n"
"  Assembly params\n"
"      --write-asqg                     Output an ASQG graph file for each 5000bp window. Default: false\n"
"  -e, --error-rate                     Fractional difference two reads can have to overlap. See SGA param. 0 is fast, but requires exact. Default: 0.05\n"
"  -c, --chunk-size                     Amount of genome to read in at once. High numbers have fewer I/O rounds, but more memory. Default 1000000 (1M). Suggested 50000000 (50M) or 'chr' for exomes\n"
"\n";

void runSnowman(int argc, char** argv) {

  parseRunOptions(argc, argv);

  std::cerr << 
    "-----------------------------------------------------------------" << std::endl << 
    "--- Running Snowman somatic indel and rearrangement detection ---" << std::endl <<
    "-----------------------------------------------------------------" << std::endl;
  std::cerr << 
    "***************************** PARAMS ****************************" << std::endl << 
    "    DBSNP Database file: " << opt::dbsnp << std::endl << 
    "    Max cov to assemble: " << opt::max_cov << std::endl << 
    "    ErrorRate: " << (opt::assemb::error_rate < 0.001f ? "EXACT (0)" : std::to_string(opt::assemb::error_rate)) << std::endl;
  if (opt::assemb::writeASQG)
    std::cerr << "    Writing ASQG files. Suggest running R/snow-asqg2pdf.R -i <my.asqg> -o graph.pdf" << std::endl;
  std::cerr <<
    "*****************************************************************" << std::endl;			  
  
#ifdef MICROBE
  // make the microbe BWAWrapper
  microbe_bwa = new SnowTools::BWAWrapper();
  microbe_bwa->retrieveIndex(opt::microbegenome);

  // open the microbe bam for writing
  //b_microbe_writer.SetWriteHeader(microbe_bwa->HeaderFromIndex());
  //b_microbe_writer.OpenWriteBam(opt::analysis_id + ".microbe.bam"); // open and write header
 
#endif

  std::cerr << "...loading the human reference sequence" << std::endl;
  findex = fai_load(opt::refgenome.c_str());  // load the reference

  // make the BWAWrapper
  main_bwa = new SnowTools::BWAWrapper();
  main_bwa->retrieveIndex(opt::refgenome);

  // open the tumor bam to get header info
  bwalker.OpenReadBam(opt::bam.begin()->first);

  // open the r2c writer
  bam_hdr_t * r2c_hdr = bam_hdr_dup(bwalker.header());
  r2c_writer.SetWriteHeader(r2c_hdr);
  r2c_writer.OpenWriteBam(opt::analysis_id + ".r2c.bam");

  // open the r2c corrected writer
  bam_hdr_t * r2c_corrected_hdr = bam_hdr_dup(bwalker.header());
  r2c_corrected_writer.SetWriteHeader(r2c_corrected_hdr);
  r2c_corrected_writer.OpenWriteBam(opt::analysis_id + ".r2c.corrected.bam");

  // open the contig bam for writing
  bwriter.SetWriteHeader(main_bwa->HeaderFromIndex());
  bwriter.OpenWriteBam(opt::analysis_id + ".contigs.bam"); // open and write header

  // open the all-contig bam for writing
  b_allwriter.SetWriteHeader(main_bwa->HeaderFromIndex());
  b_allwriter.OpenWriteBam(opt::analysis_id + ".contigs.all.bam"); // open and write header

  // open the blacklist
  if (opt::blacklist.length()) {
    std::cerr << "...reading blacklist from " << opt::blacklist << std::endl;
    blacklist.regionFileToGRV(opt::blacklist, 0, bwalker.header());
    std::cerr << "...read in " << blacklist.size() << " blacklist regions " << std::endl;
  }

  // open the DBSnpFilter
  if (opt::dbsnp.length()) {
    std::cerr << "...loading the DB Snp database" << std::endl;
    dbsnp_filter = new SnowTools::DBSnpFilter(opt::dbsnp);
    std::cerr << (*dbsnp_filter) << std::endl;
  }

  // set the MiniRules to be applied to each region
  mr = new SnowTools::MiniRulesCollection(opt::rules);
  if (opt::verbose > 1)
    std::cerr << *mr;

  // learn some parameters
  //learnParameters(); //debug
  opt::readlen = 101;
  std::cerr << "...found read length of " << opt::readlen << std::endl;

  // parse the indel mask
  if (opt::indel_mask.length()) {
    std::cerr << "...loading the indel blacklist mask" << std::endl;
    indel_blacklist_mask.regionFileToGRV(opt::indel_mask, 0, bwalker.header());
    indel_blacklist_mask.createTreeMap();
  }

  // parse the region file, count number of jobs
  SnowTools::GRC file_regions, regions_torun;
  int num_jobs = countJobs(file_regions, regions_torun); 
  std::cerr << "...running on " << num_jobs << " big-chunk regions with chunk size of " << SnowTools::AddCommas<int>(opt::chunk) << std::endl;
  if (opt::verbose > 1)
    for (auto& i : regions_torun)
      std::cerr << " REGION TO RUN " << i << std::endl;

  // override the number of threads if need
  opt::numThreads = std::min(num_jobs, opt::numThreads);

  // open the mutex
  if (pthread_mutex_init(&snow_lock, NULL) != 0) {
      printf("\n mutex init failed\n");
      return;
  }

  // open the files
  std::string n1 = opt::analysis_id + ".alignments.txt.gz";
  std::string n2 = opt::analysis_id + ".bps.txt.gz";
  std::string n3 = opt::analysis_id + ".discordant.txt.gz";
  std::string n4 = opt::analysis_id + ".tumor.coverage.bedgraph.gz";
  std::string n5 = opt::analysis_id + ".normal.coverage.bedgraph.gz";
  std::string mates_string = opt::analysis_id + ".matelookups.txt";
  std::string n6 = opt::analysis_id + ".cigarmap.txt.gz";
  all_align = new ogzstream(n1.c_str(), std::ios::out);
  os_allbps = new ogzstream(n2.c_str(), std::ios::out);
  all_disc_stream  = new ogzstream(n3.c_str(), ios::out);
  all_tum_cov = new ogzstream(n4.c_str(), std::ios::out);
  all_norm_cov = new ogzstream(n5.c_str(), std::ios::out);
  os_cigmap        = new ogzstream(n6.c_str(), std::ios::out);
  mates_file = new ofstream();
  mates_file->open(mates_string);

  // write the headers
  (*os_allbps) << SnowTools::BreakPoint::header() << endl;
  (*all_disc_stream) << SnowTools::DiscordantCluster::header() << endl;

  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);

  // send the jobs to the queue
  std::cerr << "---- Starting detection pipeline --- on " << opt::numThreads << " thread" << std::endl;
  sendThreads(regions_torun);

  // close the files
  all_align->close();
  os_allbps->close();
  all_disc_stream->close();
  all_tum_cov->close();
  all_norm_cov->close();
  mates_file->close(); 
  os_cigmap->close();
  delete all_align;
  delete os_allbps;
  delete all_disc_stream;
  delete all_tum_cov;
  delete all_norm_cov;
  delete mates_file;
  delete os_cigmap;
  

  // make the VCF file
  if (opt::verbose)
    std::cerr << "...making the VCF files" << endl;

  VCFFile snowvcf(opt::analysis_id + ".bps.txt.gz", opt::refgenome.c_str(), '\t', opt::analysis_id);

  // write the indel one
  std::string basename = opt::analysis_id + ".broad-snowman.DATECODE.";
  snowvcf.writeIndels(basename, opt::zip);
  snowvcf.writeSVs(basename, opt::zip);

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
      case 'p': arg >> opt::numThreads; break;
      case 'a': arg >> opt::analysis_id; break;
      case 'b': arg >> opt::blacklist; break;
      case 'm': arg >> opt::indel_mask; break;
      case 'q': arg >> opt::pon; break;
      case 'z': opt::zip = false; break;
      case 'h': die = true; break;
      case 'x': opt::no_r2c = true; break;
      case 'c': 
	tmp = "";
	arg >> tmp;
	if (tmp.find("chr") != std::string::npos) {
	  opt::chunk = 250000000; break;
	} else {
	  opt::chunk = stoi(tmp); break;
	}
    case OPT_ASQG: opt::assemb::writeASQG = true; break;
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
      case 'M': arg >> opt::microbegenome; break;
      case 'D': arg >> opt::dbsnp; break;
      case 'r': arg >> opt::rules; break;
      case OPT_DISC_CLUSTER_ONLY: opt::disc_cluster_only = true; break;
      default: die= true; 
    }
  }

  // check that we input something
  if (opt::bam.size() == 0) {
    std::cerr << "Must add a bam file " << std::endl;
    exit(EXIT_FAILURE);
  }

  // check file validity
  //if (opt::regionFile.length()) { // && !SnowTools::read_access_test(opt::regionFile)) {
  //  std::cerr << "Region file does not exist or is not readable: " << opt::regionFile << std::endl;
  //  exit(EXIT_FAILURE);
  //}
    
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
      if (region_id < 25) // don't add outsdie of 1-Y
	file_regions.add(SnowTools::GenomicRegion(region_id, 1, bwalker.header()->target_len[i]));
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
  for (auto& r : file_regions) {
    SnowTools::GRC test(opt::chunk, 1000, r);
    run_regions.concat(test);
  }
  return run_regions.size();
  
  /*
  unsigned jj = 0; 
  int startr, endr;;
  int kk = 0;
  
  //set amount to modulate
  int thispad = 1000;
  if (rgfile)
    thispad = 0;

  std::cerr << "file_regions.size() " << file_regions.size() << std::endl;
  
  // loop through each region
  bool stoploop = false;
  while (jj < file_regions.size()) {

    int fr_pos1 = file_regions.at(jj).pos1;
    int fr_pos2 = file_regions.at(jj).pos2;
    int fr_chr = file_regions.at(jj).chr;

    //if (jj % 1000 == 0)
    std::cerr << "jj " << jj << std::endl;

    // if regions are greater than chunk, breakup
    if ( (fr_pos2 - fr_pos1) > opt::chunk) {
      startr = std::max(1,fr_pos1-thispad);
      endr = startr + opt::chunk;

      do {
	SnowTools::GenomicRegion grr(fr_chr, startr, endr);
	run_regions.add(grr);
	
	if (endr == fr_pos2)
	  stoploop = true;
	
	std::cerr << "run_regions.size() "  << run_regions.size() << std::endl;
	kk++;
	endr   = std::min(fr_pos2, (kk+1)*opt::chunk + fr_pos1 + thispad);
	startr = std::min(fr_pos2,  kk*opt::chunk + fr_pos1);
	//cov_endr   = min(file_regions[jj].pos2, (kk+1)*opt::chunk + file_regions[jj].pos1);
	
      } while (!stoploop);
      
    } else { // region size is below chunk
      run_regions.add(file_regions.at(jj));
      //cov_regions.push_back(file_regions[jj]);
    }
    jj++;
    kk = 0;
    stoploop = false;
  } // end big while

  return run_regions.size();
  */
}

bool runBigChunk(const SnowTools::GenomicRegion& region)
{

  std::vector<SnowTools::AlignedContig> alc;
  BamReadVector all_contigs;

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
      walkers[b.first].max_cov = opt::max_cov;
      walkers[b.first].coverage_region = region;
      walkers[b.first].prefix = b.second;
      //walkers[b.first].addBlacklist(blacklist);
      walkers[b.first].setBamWalkerRegion(region);
      walkers[b.first].SetMiniRulesCollection(*mr);
	
      // read the reads
      walkers[b.first].readBam(microbe_bwa);
     
      if (b.second == "t") {
	tcount += walkers[b.first].reads.size();
      } else {
	ncount += walkers[b.first].reads.size();
      }
    }

  size_t num_t_reads = tcount;
  size_t num_n_reads = ncount;

  std::cerr << "...read in " << tcount << "/" << ncount << " Tumor/Normal reads from " << region << std::endl;

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
  MateRegionVector somatic_mate_regions;
  for (auto& b : opt::bam)
    if (b.second.at(0) == 't')
      for (auto& i : walkers[b.first].mate_regions)
	if (!normal_mate_regions.findOverlapping(i) && 
	    i.chr < 24 && 
	    (indel_blacklist_mask.size() == 0 || indel_blacklist_mask.findOverlapping(i) == 0))
	  somatic_mate_regions.add(i); //concat(walkers[b.first].mate_regions);

  //if (opt::verbose > 1)
  for (auto& i : somatic_mate_regions)
    (*mates_file) << "   somatic mate region " << i << " somatic count " << i.count << std::endl;

  // add these regions to the walker and get the reads
  tcount = 0; ncount = 0;
  if (somatic_mate_regions.size())
  for (auto& b : opt::bam)
    {
      int oreads = walkers[b.first].reads.size();
      walkers[b.first].setBamWalkerRegions(somatic_mate_regions.asGenomicRegionVector());
      walkers[b.first].get_coverage = false;
      walkers[b.first].max_cov = opt::max_cov;
      walkers[b.first].setReadKeepLimit(MATE_READ_LIMIT);
      walkers[b.first].get_mate_regions = false;
      walkers[b.first].readBam(microbe_bwa);
      if (b.second == "t") {
	tcount += (walkers[b.first].reads.size() - oreads);
      } else {
	ncount += (walkers[b.first].reads.size() - oreads);
      }
    }


  num_t_reads += tcount;
  num_n_reads += ncount;

  if (somatic_mate_regions.size())
    std::cerr << "           " << tcount << "/" << ncount << " Tumor/Normal mate reads in " << somatic_mate_regions.size() << " regions spanning " << somatic_mate_regions.width() << " bases" << std::endl;

  st.stop("m");
#endif
  
  // put all of the reads together and dedupe
  std::set<std::string> dedup;
  BamReadVector bav_join;
  for (auto& b : walkers)
    for (auto& r : b.second.reads)
      if (!dedup.count(r.GetZTag("SR"))) {
	dedup.insert(r.GetZTag("SR"));
	bav_join.push_back(r);
      }

  // do the discordant read clustering
  SnowTools::DiscordantClusterMap dmap = SnowTools::DiscordantCluster::clusterReads(bav_join, region);

  // 
  if (opt::verbose > 2)
    for (auto& i : dmap)
      std::cerr << i.first << " " << i.second << std::endl;

  // set the regions to run
  GRC grv_small;

  if (region.width() > LITTLECHUNK) // divide into smaller chunks
    grv_small = GRC(LITTLECHUNK, WINDOW_PAD, region);
  else
    grv_small.add(region);

  if (opt::verbose > 1)
    std::cerr << "running the assemblies for region " << region <<  std::endl;

  for (auto& g : grv_small) 
    {
      // set the contig uid
      std::string name = "c_" + std::to_string(g.chr+1) + "_" + std::to_string(g.pos1) + "_" + std::to_string(g.pos2);

      // get the local reference at this region and index
      std::string local_ref = getRefSequence(g, findex);
      SnowTools::BWAWrapper local_bwa;
      SnowTools::USeqVector v = { {name, local_ref} };
      local_bwa.constructIndex(v);

      // figure out the mate regions for this
      SnowTools::GRC this_regions;
      this_regions.add(g);
      for (auto& i : dmap) {
	SnowTools::GenomicRegion grm = i.second.GetMateRegionOfOverlap(g);
	if (!grm.isEmpty()) {
	  grm.pad(100);
	  this_regions.add(grm);
	}
      }
      this_regions.mergeOverlappingIntervals();
      this_regions.createTreeMap();
     
      // get the reads for this region
      BamReadVector bav_this;
      int num_total = 0;
      int num_kept = 0;
      for (auto& r : bav_join) {

	  SnowTools::GenomicRegion read(r.ChrID(), r.Position(), r.Position());
	  SnowTools::GenomicRegion mate(r.MateChrID(), r.MatePosition(), r.MatePosition());

	  if (this_regions.findOverlapping(read) || this_regions.findOverlapping(mate)) {
	    //if (g.getOverlap(read) || g.getOverlap(mate)) {

	    // align locally
	    if (false) {
	      BamReadVector tmp_loc;
	      int dum;
	      std::string seqr = r.GetZTag("KC");
	      if (!seqr.length())
		seqr = r.QualityTrimmedSequence(4, dum);
	      local_bwa.alignSingleSequence(seqr, r.Qname(), tmp_loc, false);
	      
	      //std::cout << std::endl << r << std::endl;
	      //if (tmp_loc.size()) 
	      //  std::cout << tmp_loc[0] << std::endl;
	      //else
	      //  std::cout << "NO LOCAL ALIGNMENT" << std::endl;
	      
	      ++num_total;
	      ++num_kept;
	      
	      // if it mismatches at all, add
	      if (tmp_loc.size() == 0)
		bav_this.push_back(r);
	      else if (tmp_loc[0].CigarSize() != 1) // remove full matches
		bav_this.push_back(r);
	      //else if (std::abs(r.InsertSize()) > 1000 || !r.MappedFlag())
	      //  bav_this.push_back(r);
	      else
		--num_kept;


	    } else {
	      bav_this.push_back(r);
	    }
	      //else
	      //  std::cout << "REMOVING " << r << std::endl;
	      
	  }
      }
      if (opt::verbose > 1)
	std::cerr << "Keeping " << num_kept << " of " << num_total << " after local realignment filter " << std::endl;
      // do the assemblies
      if (opt::verbose > 1)
	std::cerr << "Doing assemblies on " << name << std::endl;
      SnowmanAssemblerEngine engine(name, opt::assemb::error_rate, opt::assemb::minOverlap, opt::readlen);
      if (opt::assemb::writeASQG)
	engine.setToWriteASQG();
      engine.fillReadTable(bav_this);
      engine.performAssembly();
      if (opt::verbose > 1)
	std::cerr << "Assembled " << engine.getContigs().size() << " contigs for " << name << std::endl; 
      st.stop("as");

      //for (auto& i : engine.getContigs())//debug
      //	std::cerr << i.getSeq() << std::endl;

      // map to the reference
      for (auto& i : engine.getContigs()) {

	BamReadVector ct_alignments;
	main_bwa->alignSingleSequence(i.getSeq(), i.getID(), ct_alignments, false);
	all_contigs.insert(all_contigs.begin(), ct_alignments.begin(), ct_alignments.end());

	// make the aligned contig object
	SnowTools::AlignedContig ac(ct_alignments);
	
	// assign the local variable to each
	ac.checkLocal(g);

	st.stop("bw");

	// for contigs with min quality, realign reads, get breaks
	if (ac.getMaxMapq() >= 40 && ac.hasLocal()) {

	  // align reads with BWA
	  ac.alignReads(bav_this);
	  ac.splitCoverage();	
	  
	  // add discordant reads support to each of the breakpoints
	  ac.addDiscordantCluster(dmap);
	  
	  // add in the cigar matches
	  ac.checkAgainstCigarMatches(cigmap_n, cigmap_t);

	  // check the indel breakpoints against the indel blacklist. 
	  // simply set the blacklist flag for these breaks if they hit
	  ac.blacklist(indel_blacklist_mask);

	  // add to the final structure
	  alc.push_back(ac);

	  st.stop("sw");

	  // print it out
	  //std::cerr << ac << std::endl;
	}

      } // stop contig loop
    } // stop region loop

  // get the breakpoints
  std::vector<SnowTools::BreakPoint> bp_glob;

  for (auto& i : alc) {
    //if (i.m_bamreads.size() && !i.m_skip) { //m_bamreads.size() is zero for contigs with ...
    std::vector<SnowTools::BreakPoint> allbreaks = i.getAllBreakPoints();
    bp_glob.insert(bp_glob.end(), allbreaks.begin(), allbreaks.end());
    //}
  }

  // repeat sequence filter
  for (auto& i : bp_glob)
    i.repeatFilter(findex);

  // db snp filter
  for (auto & i : bp_glob) {
    dbsnp_filter->queryBreakpoint(i);
  }

  // add in the discordant clusters as breakpoints
  for (auto& i : dmap) {
    // DiscordantCluster not associated with assembly BP and has 2+ read support
    if (!i.second.hasAssociatedAssemblyContig() && (i.second.tcount + i.second.ncount) > 1) {
      SnowTools::BreakPoint tmpbp(i.second);
      bp_glob.push_back(tmpbp);
    }
  }

  // de duplicate the breakpoints
  std::sort(bp_glob.begin(), bp_glob.end());
  bp_glob.erase( std::unique( bp_glob.begin(), bp_glob.end() ), bp_glob.end() );


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

  ////////////////////////////////////
  // MUTEX LOCKED
  ////////////////////////////////////
  // write to the global contig out
  pthread_mutex_lock(&snow_lock);  

  // dump the cigmaps
  for (auto& i : cigmap_n) 
    (*os_cigmap) << i.first << "\t" << i.second << "\tN" << endl;
  for (auto& i : cigmap_t) 
    (*os_cigmap) << i.first << "\t" << i.second << "\tT" << endl;


  // dump the coverage track
  //for (auto& b : opt::bam) {
  //  if (b.second == "t")
  //    walkers[b.first].cov.ToBedgraph(all_tum_cov, walkers[b.first].header());
  //  else if (b.second == "n")
  //    walkers[b.first].cov.ToBedgraph(all_norm_cov, walkers[b.first].header());
  //}

  // print the alignment plots
  size_t contig_counter = 0;
  for (auto& i : alc)
    if (i.hasVariant()) {
      (*all_align) << i << std::endl;
      ++contig_counter;
    }

#ifdef MICROBE
  // send the microbe to file
  //for (auto& b : bb_microbe)
  // b_microbe_writer.WriteAlignment(b);
#endif

  // send the discordant to file
  for (auto& i : dmap)
    if (i.second.getMeanMapq(false) >= 10 && i.second.getMeanMapq(true) >= 10)
      (*all_disc_stream) << i.second.toFileString(!opt::no_reads) << std::endl;

  // send breakpoints to file
  for (auto& i : bp_glob)
    if (i.hasMinimal() || true) {
      (*os_allbps) << i.toFileString(opt::no_reads) << std::endl;
      if (i.confidence == "PASS" && i.nsplit == 0 && i.ncigar == 0 && i.dc.ncount == 0)
	std::cerr << "SOMATIC " << i.toPrintString() << std::endl;
      else if (i.confidence == "PASS")
	std::cerr << "GERMLINE " << i.toPrintString() << std::endl;
    }

  // write the contigs
  for (auto& i : alc)
    i.writeToBAM(bwriter);

  // write ALL contigs
  for (auto& i : all_contigs)
    b_allwriter.WriteAlignment(i);
  
  // write all the to-assembly reads
  for (auto& b : opt::bam) {
    if (b.second == "t")
      for (auto& r : walkers[b.first].reads) {
	r2c_writer.WriteAlignment(r);
	
	//write the corrected
	std::string new_seq  = r.GetZTag("KC");
	if (new_seq.length())
	  r.SetSequence(new_seq);
	r2c_corrected_writer.WriteAlignment(r);
      }
  }
  
  // display the run time
  if (opt::verbose > 0 && (num_t_reads + num_n_reads) > 0) {
    string print1 = SnowTools::AddCommas<int>(region.pos1);
    string print2 = SnowTools::AddCommas<int>(region.pos2);
    char buffer[140];
    sprintf (buffer, "Ran chr%2s:%11s-%11s | T: %5d N: %5d C: %5d| ", 
	     bwalker.header()->target_name[region.chr],
	     print1.c_str(),print2.c_str(),
	     (int)num_t_reads, (int)num_n_reads, 
	     (int)contig_counter);
    
    std::cerr << string(buffer) << st << " | ";
    std::cerr << SnowTools::displayRuntime(start);
    //int perc = SnowTools::percentCalc<int>(CHR_LEN[refID] + pos2, 3100000000);
    //cout << " % of WG done: " << perc << "%" << endl;
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
  SnowTools::GenomicRegion it;
  for (auto& i : regions_torun) {
    SnowmanWorkItem * item     = new SnowmanWorkItem(SnowTools::GenomicRegion(i.chr, i.pos1, i.pos2), ++count);
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
  while (bwalker.GetNextRead(r, rule) && ++count < 10000) 
      opt::readlen = std::max(r.Length(), opt::readlen);

}

