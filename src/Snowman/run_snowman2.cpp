 #include "run_snowman.h"

#include <getopt.h>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <map>
#include <vector>

#include "bwa/bwa.h"

//#include "boost/filesystem/path.hpp" 
//#include "boost/archive/text_oarchive.hpp"
//#include "boost/archive/text_iarchive.hpp"

#include "SnowTools/MiniRules2.h"
#include "SnowTools/BLATWrapper.h"

#include "vcf.h"
#include "PONFilter.h"
#include "DBSnpFilter.h"
#include "SnowmanUtils.h"
#include "LearnBamParams.h"

#define MIN_CONTIG_MATCH 35
#define MATE_LOOKUP_MIN 3
#define SECONDARY_CAP 25

#define LITTLECHUNK 25000 
#define WINDOW_PAD 500
#define MICROBE_MATCH_MIN 50
#define GET_MATES 1
#define MICROBE 1

static BamParams t_params;
static BamParams n_params;

static SnowTools::GRC bad_mate_regions;

static std::set<std::string> prefixes;

static SnowTools::PONFilter * pon_filter = nullptr;

// Declares a default Aligner
//static StripedSmithWaterman::Aligner sw_aligner;
// Declares a default filter
//static StripedSmithWaterman::Filter sw_filter;

struct bidx_delete {
  void operator()(void* x) { hts_idx_destroy((hts_idx_t*)x); }
};

typedef std::map<std::string, std::shared_ptr<hts_idx_t>> BamIndexMap;
typedef std::map<std::string, std::string> BamMap;

static BamIndexMap bindices;
static faidx_t * findex;
static faidx_t * findex_viral;
typedef std::unordered_map<std::string, size_t> SeqHash;
static SeqHash over_represented_sequences;

static SnowTools::BamWalker bwalker, r2c_writer, er_writer, b_microbe_writer, b_allwriter, blat_allwriter;
static SnowTools::BWAWrapper * microbe_bwa = nullptr;
static SnowTools::BWAWrapper * main_bwa = nullptr;
static SnowTools::BLATWrapper * main_blat = nullptr;
static SnowTools::MiniRulesCollection * mr;
static SnowTools::GRC blacklist, indel_blacklist_mask;
static SnowTools::DBSnpFilter * dbsnp_filter;

// output files
static ogzstream all_align, os_allbps, all_disc_stream, os_cigmap, r2c_c;
static std::ofstream log_file, bad_bed;

static SnowTools::GRC file_regions, regions_torun;

static pthread_mutex_t snow_lock;
static struct timespec start;

// replace function
std::string myreplace(std::string &s,
                      std::string toReplace,
                      std::string replaceWith)
{
  return(s.replace(s.find(toReplace), toReplace.length(), replaceWith));
}

namespace opt {

  namespace assemb {
    static int minOverlap = -1;
    static float error_rate = 0.025; 
    static bool writeASQG = false;
  }

  static int sample_number = 0;

  static std::string args = "snowman ";

  static bool exome = false;

  static bool run_blat = false;

  static bool kmer_filter = false;

  static std::string ooc = "/xchip/gistic/Jeremiah/blat/11.ooc";

  static int num_assembly_rounds = 3;
  static bool write_extracted_reads = false;

  static std::string pon;
  static bool no_assemble_normal = false;
  static std::string indel_mask = ""; //"/xchip/gistic/Jeremiah/Projects/HengLiMask/um75-hs37d5.bed.gz";

  static bool no_reads = true;
  static int32_t readlen;
  static bool r2c = false;
  static bool zip = false;
  //  static std::string pon = "";

  static int gap_open_penalty = 6;
  static int gap_extension_penalty = 1;
  static int mismatch_penalty = 4;
  static int zdrop = 100;
  static int bandwidth = 100;
  static float reseed_trigger = 1.5;

  // parameters for filtering reads
  //static std::string rules = "global@!hardclip;!duplicate;!qcfail;phred[4,100];length[LLL,1000]%region@WG%!isize[0,800]%ic%clip[10,1000]%ins[1,1000];mapq[0,100]%del[1,1000];mapq[1,1000]%mapped;!mate_mapped;mapq[1,1000]%mate_mapped;!mapped";
  //static std::string rules = "global@!hardclip;!duplicate;!qcfail;phred[4,100]%region@WG%discordant[0,800]%ic%clip[10,1000]%ins[1,1000];mapq[0,100]%del[1,1000];mapq[1,1000]%mapped;!mate_mapped;mapq[1,1000]%mate_mapped;!mapped";
  static std::string rules = "{\"global\" : {\"duplicate\" : false, \"qcfail\" : false}, \"\" : { \"rules\" : [{\"isize\" : [800,0]},{\"rr\" : true},{\"ff\" : true}, {\"rf\" : true}, {\"ic\" : true}, {\"clip\" : [5, 1000], \"phred\" : [4,100]}, {\"ins\" : [1,1000]}, {\"del\" : [1,1000]}, {\"mapped\": true , \"mate_mapped\" : false}, {\"mate_mapped\" : true, \"mapped\" : false}]}}";
  //static std::string rules = "global@!duplicate;!qcfail%region@WG%discordant[0,600]%ic%clip[5,1000];phred[4,100]%ins[1,1000]%del[1,1000]%mapped;!mate_mapped%mate_mapped;!mapped";
  static int num_to_sample = 2000000;
  static std::string motif_exclude;
  
  static int max_cov = 100;
  static int chunk = LITTLECHUNK; // 1000000;

  // runtime parameters
  static int verbose = 1;
  static int numThreads = 1;

  // data
  static BamMap bam;
  static std::string refgenome = SnowTools::REFHG19;  
  static std::string microbegenome; //"/xchip/gistic/Jeremiah/Projects/SnowmanFilters/viral.1.1.genomic_ns.fna";  
  static std::string analysis_id = "no_id";

  //subsample
  float subsample = 1.0;

  static std::string regionFile;
  static std::string blacklist; // = "/xchip/gistic/Jeremiah/Projects/HengLiMask/um75-hs37d5.bed.gz";

  static std::string dbsnp; // = "/xchip/gistic/Jeremiah/SnowmanFilters/dbsnp_138.b37_indel.vcf";

  // filters on when / how to assemble
  static bool disc_cluster_only = false;

  static bool refilter = false;
  static bool adapter_trim = true;

  static std::string gemcode;

  static std::string normal_bam;
  static std::string tumor_bam;

  static size_t mate_lookup_min = 3;

  static double lod = 8;
  static double lod_db = 7;
  static double lod_no_db = 2.5;
  static double lod_germ = 3;

  static bool interchrom_lookup = true;

}

enum { 
  OPT_ASQG,
  OPT_LOD,
  OPT_DLOD,
  OPT_NDLOD,
  OPT_LODGERM,
  OPT_DISC_CLUSTER_ONLY,
  OPT_READ_TRACK,
  OPT_NO_ASSEMBLE_NORMAL,
  OPT_MAX_COV,
  OPT_DISCORDANT_ONLY,
  OPT_REFILTER,
  OPT_WRITE_EXTRACTED_READS,
  OPT_ADAPTER_TRIM,
  OPT_GEMCODE_AWARE,
  OPT_NUM_ASSEMBLY_ROUNDS,
  OPT_NUM_TO_SAMPLE,
  OPT_EXOME,
  OPT_NOBLAT,
  OPT_KMER,
  OPT_GAP_OPEN,
  OPT_GAP_EXTENSION,
  OPT_MISMATCH,
  OPT_ZDROP,
  OPT_BANDWIDTH,
  OPT_RESEED_TRIGGER
};

static const char* shortopts = "hzxt:n:p:v:r:G:r:e:g:k:c:a:m:B:D:Y:S:P:L:O:I";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "tumor-bam",               required_argument, NULL, 't' },
  //  { "indel-mask",              required_argument, NULL, 'M' },
  { "panel-of-normals",        required_argument, NULL, 'P' },
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
  { "motif-filter",            required_argument, NULL, 'S' },
  { "mate-lookup-min",         required_argument, NULL, 'L' },
  { "blat-ooc",                required_argument, NULL, 'O' },
  { "g-zip",                 no_argument, NULL, 'z' },
  { "no-interchrom-lookup",                 no_argument, NULL, 'I' },
  { "read-tracking",         no_argument, NULL, OPT_READ_TRACK },
  { "gap-open-penalty",         required_argument, NULL, OPT_GAP_OPEN },
  { "gap-extension-penalty",         required_argument, NULL, OPT_GAP_EXTENSION },
  { "mismatch-penalty",         required_argument, NULL, OPT_MISMATCH },
  { "z-dropoff",         required_argument, NULL, OPT_ZDROP },
  { "reseed-trigger",         required_argument, NULL, OPT_RESEED_TRIGGER },
  { "bandwidth",         required_argument, NULL, OPT_BANDWIDTH },
  { "write-extracted-reads", no_argument, NULL, OPT_WRITE_EXTRACTED_READS },
  { "no-assemble-normal",    no_argument, NULL, OPT_NO_ASSEMBLE_NORMAL },
  { "exome",                 no_argument, NULL, OPT_EXOME },
  { "lod",                   required_argument, NULL, OPT_LOD },
  { "lod-dbsnp",             required_argument, NULL, OPT_DLOD },
  { "lod-no-dbsnp",          required_argument, NULL, OPT_NDLOD },
  { "lod-germ",              required_argument, NULL, OPT_LODGERM },
  { "discordant-only",       no_argument, NULL, OPT_DISCORDANT_ONLY },
  { "num-to-sample",         required_argument, NULL, OPT_NUM_TO_SAMPLE },
  { "refilter",              no_argument, NULL, OPT_REFILTER },
  { "r2c-bam",               no_argument, NULL, 'x' },
  { "write-asqg",            no_argument, NULL, OPT_ASQG   },
  { "kmer-correct",              no_argument, NULL, OPT_KMER   },
  { "error-rate",            required_argument, NULL, 'e'},
  { "verbose",               required_argument, NULL, 'v' },
  { "blacklist",             required_argument, NULL, 'B' },
  { "max-coverage",          required_argument, NULL, OPT_MAX_COV },
  { "no-adapter-trim",       no_argument, NULL, OPT_ADAPTER_TRIM },
  { "num-assembly-rounds",   required_argument, NULL, OPT_NUM_ASSEMBLY_ROUNDS },
  { "assemble-by-tag",       required_argument, NULL, OPT_GEMCODE_AWARE },
  { "with-blat",             no_argument, NULL, OPT_NOBLAT },
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
"      --refilter                       Recreate the VCF files from the bps.txt.gz file. Need to supply correct id-string.\n"
"  Required input\n"
"  -G, --reference-genome               Path to indexed reference genome to be used by BWA-MEM. Default is Broad hg19 (/seq/reference/...)\n"
"  -t, --tumor-bam                      Tumor BAM file\n"
"  Variant filtering and classification\n"
"      --lod                            LOD cutoff to classify somatic indel as real / artifact (tests AF=0 vs AF=MaxLikelihood(AF)) [8]\n"
"      --lod-dbsnp                      LOD cutoff to classify indel as somatic (at DBsnp site) [7]\n"
"      --lod-no-dbsnp                   LOD cutoff to classify indel as somatic (not at DBsnp site) [2.5]\n"
"      --lod-germ                       LOD cutoff for germline indels that variant is real (tests AF=0 vs AF=0.5)[3]\n"
"  Optional input\n"                       
"  -n, --normal-bam                     Normal BAM file\n"
"  -r, --rules                          VariantBam style rules string to determine which reads to do assembly on. See documentation for default.\n"
"  -m, --min-overlap                    Minimum read overlap, an SGA parameter. Default: 0.4* readlength\n"
"  -k, --region-file                    Set a region txt file. Format: one region per line, Ex: 1,10000000,11000000\n"
  //"  -P, --panel-of-normals               Panel of normals gzipped txt file generated from snowman pon\n"
  //"  -M, --indel-mask                     BED-file with graylisted regions for stricter indel calling.\n"
"  -D, --dbsnp-vcf                      DBsnp database (VCF) to compare indels against\n"
"  -B, --blacklist                      BED-file with blacklisted regions to not extract any reads from.\n"
"  -V, --microbial-genome               Path to indexed reference genome of microbial sequences to be used by BWA-MEM to filter reads.\n"
"  -z, --g-zip                          Gzip and tabix the output VCF files. Default: off\n"
"  -L, --min-lookup-min                 Minimum number of somatic reads required to attempt mate-region lookup [3]\n"
"      --no-interchrom-lookup           Don't do mate lookup for inter-chromosomal translocations. This increased run time at cost of power for translocations.\n"
"      --kmer-correct                   Perform k-mer based error correction on the reads. Useful for noisy data, need to balance with error rate (w/kmer can have lower error rate). [off]\n"
"      --r2c-bam                        Output a BAM of reads that aligned to a contig, and fasta of kmer corrected sequences\n"
"      --discordant-only                Only run the discordant read clustering module, skip assembly. Default: off\n"
"      --read-tracking                  Track supporting reads. Increases file sizes.\n"
"      --max-coverage                   Maximum weird read coverage to send to assembler (per BAM). Subsample reads to achieve max if overflow. [500]\n"
"      --write-extracted-reads          For the tumor BAM, write the extracted reads (those sent to assembly) to a BAM file. Good for debugging.\n"
"      --no-adapter-trim                Don't peform Illumina adapter trimming, which removes reads with AGATCGGAAGAGC present.\n"
"      --assemble-by-tag                Separate the assemblies and read-to-contig mapping by the given read tag. Useful for 10X Genomics data (e.g. --aseembly-by-tag BX).\n"
"      --no-assemble-normal             Don't kmer correct or assembly normal reads. Still maps normal reads to contigs. Increases speed for somatic-only queries.\n"
"      --num-assembly-rounds            Number of times to run the assembler per window. > 1 will bootstrap the assembly with error rate = 0.05. [1]\n"
"      --num-to-sample                  When learning about BAM, number of reads to sample. Default 1,000,000.\n"
"      --motif-filter                   Add a motif file to exclude reads that contain this motif (e.g. illumina adapters).\n"
"      --with-blat                      Run BLAT on assembled contigs that have clips/indels to reference. Single-thread only. In-development. Outputs *.blat.bam\n"
"  Assembly params\n"
"      --write-asqg                     Output an ASQG graph file for each assembly window.\n"
"  -e, --error-rate                     Fractional difference two reads can have to overlap. See SGA param. 0 is fast, but requires exact matches and error correcting. [0.02]\n"
"  -c, --chunk-size                     Size of a local assembly window (in bp). [25000]\n"
"  Alignment params\n"
"      --gap-open-penalty               Set the BWA-MEM gap open penalty for contig to genome alignments. BWA-MEM -O [6]\n"
"      --gap-extension-penalty          Set the BWA-MEM gap extension penalty for contig to genome alignments. BWA-MEM -E [1]\n"
"      --mismatch-penalty               Set the BWA-MEM mismatch penalty for contig to genome alignments. BWA-MEM -b [4]\n"
"      --bandwidth                      Set the BWA-MEM SW alignment bandwidth for contig to genome alignments. BWA-MEM -w [100]\n"
"      --z-dropoff                      Set the BWA-MEM SW alignment Z-dropoff for contig to genome alignments. BWA-MEM -d [100]\n"
"      --reseed-trigger                 Set the BWA-MEM reseed trigger for reseeding mems for contig to genome alignments. BWA-MEM -r [1.5]\n"
"\n";

void runSnowman(int argc, char** argv) {

  parseRunOptions(argc, argv);

  // test
  //run_test_assembly();

  SnowmanUtils::fopen(opt::analysis_id + ".log", log_file);

  if (opt::refilter) {
    // read in all of the breakpoints
    makeVCFs();
    return;
  }

  std::stringstream ss;
  ss << 
    "-----------------------------------------------------------------" << std::endl << 
    "--- Running Snowman somatic indel and rearrangement detection ---" << std::endl <<
    "-----------------------------------------------------------------" << std::endl;
  ss << 
    "***************************** PARAMS ****************************" << std::endl << 
    "    DBSNP Database file: " << opt::dbsnp << std::endl << 
    "    BLAT OOC file:       " << opt::ooc << std::endl << 
    "    Indel PON file:      " << opt::pon << std::endl << 
    "    Max cov to assemble: " << opt::max_cov << std::endl << 
    "    ErrorRate: " << (opt::assemb::error_rate < 0.001f ? "EXACT (0)" : std::to_string(opt::assemb::error_rate)) << std::endl << 
    "    Num assembly rounds: " << opt::num_assembly_rounds << std::endl << 
    "    Num reads to sample: " << opt::num_to_sample << std::endl << 
    "    Remove clipped reads with adapters? " << (opt::adapter_trim ? "TRUE" : "FALSE") << std::endl << 
    "    Minimum number of reads for mate lookup " << opt::mate_lookup_min << std::endl <<
    "    LOD cutoff (artifact vs real) :  " << opt::lod << std::endl << 
    "    LOD cutoff (somatic vs germline, at DBSNP):  " << opt::lod_db << std::endl << 
    "    LOD cutoff (somatic vs germlin, no DBSNP):  " << opt::lod_no_db << std::endl << 
    "    LOD cutoff (germline, AF>=0.5 vs AF=0):  " << opt::lod_germ << std::endl <<
    "    Gap open penalty: " << opt::gap_open_penalty << std::endl << 
    "    Gap extension penalty: " << opt::gap_extension_penalty << std::endl <<
    "    Mismatch penalty: " << opt::mismatch_penalty << std::endl <<
    "    Aligment bandwidth: " << opt::bandwidth << std::endl <<
    "    Z-dropoff: " << opt::zdrop << std::endl <<
    "    Reseed trigger: " << opt::reseed_trigger << std::endl;

  if (!opt::run_blat)
    ss << "    Running with BWA-MEM only (no BLAT)" << std::endl;
  if (opt::assemb::writeASQG)
    ss << "    Writing ASQG files. Suggest running R/snow-asqg2pdf.R -i <my.asqg> -o graph.pdf" << std::endl;
  if (opt::write_extracted_reads)
    ss << "    Writing extracted reads and fasta of kmer-corrected reads." << std::endl;
  if (opt::r2c)
    ss << "    Writing read-to-contig files." << std::endl;
  if (opt::disc_cluster_only)
    ss << "    ######## ONLY DISCORDANT READ CLUSTERING. NO ASSEMBLY ##############" << std::endl;
  if (!opt::interchrom_lookup)
    ss << "    ######## NOT LOOKING UP MATES FOR INTERCHROMOSOMAL #################" << std::endl;
  ss <<
    "*****************************************************************" << std::endl;	  

  SnowmanUtils::print(ss, log_file, opt::verbose > 0);

  if (opt::disc_cluster_only) 
    static std::string rules = "global@nbases[0,0];!hardclip;!supplementary;!duplicate;!qcfail;%region@WG%discordant[0,800];mapq[1,1000]";
  
  // open the microbe genome
  if (!opt::microbegenome.empty()) {
    ss << "...loading the microbe reference sequence" << std::endl;
    SnowmanUtils::print(ss, log_file, opt::verbose > 0);
    microbe_bwa = new SnowTools::BWAWrapper();
    findex_viral = SnowmanUtils::__open_index_and_writer(opt::microbegenome, microbe_bwa, opt::analysis_id + ".microbe.bam", b_microbe_writer, findex_viral);  
  }
  
  // open the tumor bam to get header info
  bwalker.OpenReadBam(opt::tumor_bam);
  if (opt::r2c) // open the r2c writer
    SnowmanUtils::__openWriterBam(bwalker, opt::analysis_id + ".r2c.bam", r2c_writer);
  if (opt::write_extracted_reads) // open the extracted reads writer
    SnowmanUtils::__openWriterBam(bwalker, opt::analysis_id + ".tumor.extracted.reads.bam", er_writer);    

  // open the blacklists
  SnowmanUtils::__open_bed(opt::blacklist, blacklist, bwalker.header());
  SnowmanUtils::__open_bed(opt::indel_mask, indel_blacklist_mask, bwalker.header());

  // open the DBSnpFilter
  if (opt::dbsnp.length()) {
    std::cerr << "...loading the DBsnp database" << std::endl;
    dbsnp_filter = new SnowTools::DBSnpFilter(opt::dbsnp);
    std::cerr << "..." << (*dbsnp_filter) << std::endl;
  }

  // open the PONFilter 
  if (!opt::pon.empty()) {
    std::cerr << "...loading the PON Filter" << std::endl;
    if (!SnowTools::read_access_test(opt::pon)) {
      std::cerr << "!!!!!!!!! ERROR !!!!!!!!!" << std::endl << "   Cannot read PON file: " << opt::pon
		<< "!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    } else {
      pon_filter = new SnowTools::PONFilter(opt::pon);
      std::cerr << (*pon_filter) << std::endl;
    }
  }
  
  // open the Bam Headeres
  std::cerr << "...loading BAM indices" << std::endl;
  for (auto& b : opt::bam) 
    bindices[b.first] = std::shared_ptr<hts_idx_t>(hts_idx_load(b.second.c_str(), HTS_FMT_BAI), bidx_delete());

  // set the prefixes
  for (auto& b : opt::bam)
    prefixes.insert(b.first);

  // learn some parameters
  LearnBamParams tparm(opt::tumor_bam);
  tparm.learnParams(t_params, opt::num_to_sample);
  ss << opt::tumor_bam << std::endl << t_params << std::endl;
  opt::readlen = t_params.readlen;
  if (opt::assemb::minOverlap < 0) {
    opt::assemb::minOverlap = 0.4 * t_params.readlen;
  }
  ss << "...found read length of " << opt::readlen << ". Min Overlap is " << opt::assemb::minOverlap << std::endl;
  ss << "...isize distribution (tumor):  " << t_params.mean_isize << " (" << t_params.sd_isize << ") Median: " << t_params.median_isize << std::endl;
  ss << "...isize distribution (normal): " << n_params.mean_isize << " (" << n_params.sd_isize << ") Median: " << n_params.median_isize << std::endl;
  if (opt::assemb::minOverlap < 30)
    opt::assemb::minOverlap = 30;
  if (!opt::normal_bam.empty()) {
    LearnBamParams nparm(opt::normal_bam);
    nparm.learnParams(n_params, opt::num_to_sample);
    ss << opt::normal_bam << std::endl << n_params << std::endl;
    if (t_params.readlen != n_params.readlen)
      ss << "!!!!! Tumor and normal readlengths differ. Tumor: " << t_params.readlen << " Normal: " << n_params.readlen << std::endl;
  }
  SnowmanUtils::print(ss, log_file, opt::verbose);

  // get the seed length for printing
  int seedLength, seedStride;
  SnowmanAssemblerEngine enginetest("test", opt::assemb::error_rate, opt::assemb::minOverlap, opt::readlen);
  enginetest.calculateSeedParameters(opt::readlen, opt::assemb::minOverlap, seedLength, seedStride);
  ss << "...calculated seed size for error rate of " << opt::assemb::error_rate << " and read length " << opt::readlen << " is " << seedLength << std::endl;
  std::cerr << "..seed size w/error rate " << opt::assemb::error_rate << " and read length " << opt::readlen << " is " << seedLength << std::endl;


  // parse the region file, count number of jobs
  int num_jobs = SnowmanUtils::countJobs(opt::regionFile, file_regions, regions_torun,
					 bwalker.header(), opt::chunk, WINDOW_PAD); 
  if (num_jobs)
    std::cerr << "...running on " << SnowTools::AddCommas(num_jobs) << " regions with chunk size of " << SnowTools::AddCommas<int>(opt::chunk) << std::endl;
  else
    std::cerr << "Chunk was <= 0: READING IN WHOLE GENOME AT ONCE" << std::endl;

  //int len_cutoff = 0.5 * opt::readlen;
  //opt::rules = myreplace(opt::rules, "LLL", std::to_string(len_cutoff));

  // set the MiniRules to be applied to each region
  ss << opt::rules << std::endl;
  mr = new SnowTools::MiniRulesCollection(opt::rules);
  ss << *mr << std::endl;
  log_file << ss.str();
  if (opt::verbose > 1)
    std::cerr << ss.str();
  ss.str(std::string());

  // override the number of threads if need
  num_jobs = (num_jobs == 0) ? 1 : num_jobs;
  opt::numThreads = std::min(num_jobs, opt::numThreads);

  // open the human reference
  ss << "...loading the human reference sequence for BWA" << std::endl;
  SnowmanUtils::print(ss, log_file, opt::verbose > 0);
  main_bwa = new SnowTools::BWAWrapper();
  main_bwa->setGapOpen(opt::gap_open_penalty);
  main_bwa->setGapExtension(opt::gap_extension_penalty);
  main_bwa->setMismatchPenalty(opt::mismatch_penalty);
  main_bwa->setZDropoff(opt::zdrop);
  main_bwa->setBandwidth(opt::bandwidth);
  main_bwa->setReseedTrigger(opt::reseed_trigger);

  findex = SnowmanUtils::__open_index_and_writer(opt::refgenome, main_bwa, opt::analysis_id + ".contigs.bam", b_allwriter, findex);  

  // open the BLAT reference
  if (opt::run_blat) {
    ss << "...loading the human reference sequence for BLAT" << std::endl;
    SnowmanUtils::print(ss, log_file, opt::verbose > 0);
    main_blat = new SnowTools::BLATWrapper();

    // add a header made from the index, so BLAT can convert between reference names and ID #s
    bam_hdr_t * bwa_header = main_bwa->HeaderFromIndex(); // creates new one in memory, need to destroy bwa_header

    main_blat->addHeader(bwa_header);
    main_blat->loadIndex(opt::refgenome, opt::ooc);

    // serialize it
    //std::ofstream brf("blatref.dat");
    //boost::archive::text_oarchive oa(brf);
    //oa << (*main_blat);

    // open the BLAT writer
    blat_allwriter.SetWriteHeader(bwa_header); // passes off to smart pointer, which will destory it
    blat_allwriter.OpenWriteBam(opt::analysis_id + ".blat.bam"); // open and write header
  }

  // open the mutex
  if (pthread_mutex_init(&snow_lock, NULL) != 0) {
      printf("\n mutex init failed\n");
      return;
  }

  // open the files
  SnowmanUtils::fopen(opt::analysis_id + ".alignments.txt.gz", all_align);
  SnowmanUtils::fopen(opt::analysis_id + ".bps.txt.gz", os_allbps);
  SnowmanUtils::fopen(opt::analysis_id + ".discordant.txt.gz", all_disc_stream);
  SnowmanUtils::fopen(opt::analysis_id + ".cigarmap.txt.gz", os_cigmap);
  if (opt::write_extracted_reads) 
    SnowmanUtils::fopen(opt::analysis_id + ".corrected.fa", r2c_c); 
  SnowmanUtils::fopen(opt::analysis_id + ".bad_regions.bed", bad_bed);
  
  // write the headers
  os_allbps << SnowTools::BreakPoint::header();
  for (auto& b : opt::bam) 
    os_allbps << "\t" << b.first << "_" << b.second; //boost::filesystem::path(b.second).filename().string();
  os_allbps << std::endl;
  all_disc_stream << SnowTools::DiscordantCluster::header() << std::endl;

  // put args into string for VCF later
  for (int i = 0; i < argc; ++i)
    opt::args += std::string(argv[i]) + " ";

  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);

  // send the jobs to the queue
  std::cerr << std::endl << "---- Starting detection pipeline --- on " << 
    opt::numThreads << " thread" << (opt::numThreads > 1 ? "s" : "") << std::endl;
  sendThreads(regions_torun);

  if (microbe_bwa)
    delete microbe_bwa;
  //if (main_bwa)
  //  delete main_bwa;  
  
  // close the files
  all_align.close();
  os_allbps.close();
  all_disc_stream.close();
  if (opt::r2c) 
    r2c_c.close();
  os_cigmap.close();
  log_file.close();
  bad_bed.close();

  // make the VCF file
  makeVCFs();

  std::cerr << SnowTools::displayRuntime(start) << std::endl;
}

void makeVCFs() {

  if (opt::bam.size() == 0) {
    std::cerr << "makeVCFs error: must supply a BAM via -t to get header from" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (!bwalker.header()) {    // open the tumor bam to get header info
    bwalker.OpenReadBam(opt::bam.begin()->second);
  }
    
  // make the VCF file
  if (opt::verbose)
    std::cerr << "...loading the bps files for conversion to VCF" << std::endl;

  std::string file = opt::analysis_id + ".bps.txt.gz";
  if (!SnowTools::read_access_test(file))
    file = opt::analysis_id + ".bps.txt";

  // make the header
  VCFHeader header;
  header.filedate = SnowmanUtils::fileDateString();
  header.source = opt::args;
  header.reference = opt::refgenome;

  bam_hdr_t * bwa_header = nullptr;
  if (main_bwa)
    bwa_header = main_bwa->HeaderFromIndex();
  if (main_bwa)
    delete main_bwa;  
  if (bwa_header)
    for (int i = 0; i < bwa_header->n_targets; ++i)
      header.addContigField(bwa_header->target_name[i], bwa_header->target_len[i]);

  if (bwa_header)
    bam_hdr_destroy(bwa_header);

  for (auto& b : opt::bam) {
    //boost::filesystem::path bpf(b.second);
    std::string fname = b.second; //bpf.filename();
    header.addSampleField(fname);
    header.colnames += "\t" + fname; 
  }


  // primary VCFs
  if (SnowTools::read_access_test(file)) {
    if (opt::verbose)
      std::cerr << "...making the primary VCFs (unfiltered and filtered) from file " << file << std::endl;
    VCFFile snowvcf(file, opt::analysis_id, bwalker.header(), header);
    std::string basename = opt::analysis_id + ".snowman.unfiltered.";

    std::cerr << "...writing unfiltered VCFs" << std::endl;
    snowvcf.include_nonpass = true;
    snowvcf.writeIndels(basename, opt::zip, opt::bam.size() == 1);
    snowvcf.writeSVs(basename, opt::zip, opt::bam.size() == 1);

    std::cerr << "...writing filtered VCFs" << std::endl;
    basename = opt::analysis_id + ".snowman.";
    snowvcf.include_nonpass = false;
    snowvcf.writeIndels(basename, opt::zip, opt::bam.size() == 1);
    snowvcf.writeSVs(basename, opt::zip, opt::bam.size() == 1);

  } else {
    std::cerr << "Failed to make VCF. Could not file bps file " << file << std::endl;
  }

}

// parse the command line options
void parseRunOptions(int argc, char** argv) {
  bool die = false;

  if (argc <= 2) 
    die = true;

  bool help = false;
  std::stringstream ss;

  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case OPT_REFILTER : opt::refilter = true; break;
    case OPT_NOBLAT : opt::run_blat = true; break;
      case 'p': arg >> opt::numThreads; break;
      case 'm': arg >> opt::assemb::minOverlap; break;
      case 'a': arg >> opt::analysis_id; break;
      case 'B': arg >> opt::blacklist; break;
      case 'L': arg >> opt::mate_lookup_min; break;
      case 'M': arg >> opt::indel_mask; break;
    case 'I': opt::interchrom_lookup = false; break;
      case 'Y': arg >> opt::microbegenome; break;
      case 'P': arg >> opt::pon; break;
      case 'O': arg >> opt::ooc; break;
      case 'z': opt::zip = false; break;
      case 'h': help = true; break;
    case OPT_GAP_OPEN : arg >> opt::gap_open_penalty; break;
    case OPT_MISMATCH : arg >> opt::mismatch_penalty; break;
    case OPT_GAP_EXTENSION : arg >> opt::gap_extension_penalty; break;
    case OPT_ZDROP : arg >> opt::zdrop; break;
    case OPT_BANDWIDTH : arg >> opt::bandwidth; break;
    case OPT_RESEED_TRIGGER : arg >> opt::reseed_trigger; break;
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
    case OPT_EXOME: opt::exome = true; break;
    case OPT_LOD: arg >> opt::lod; break;
    case OPT_DLOD: arg >> opt::lod_db; break;
    case OPT_NDLOD: arg >> opt::lod_no_db; break;
    case OPT_LODGERM: arg >> opt::lod_germ; break;
    case OPT_MAX_COV: arg >> opt::max_cov;  break;
    case OPT_NUM_TO_SAMPLE: arg >> opt::num_to_sample;  break;
    case OPT_READ_TRACK: opt::no_reads = false; break;
	case 't': 
	  tmp = SnowmanUtils::__bamOptParse(opt::bam, arg, opt::sample_number++, "t");
	  if (opt::tumor_bam.empty())
	    opt::tumor_bam = tmp;
	  break;
	case 'n': 
	  tmp = SnowmanUtils::__bamOptParse(opt::bam, arg, opt::sample_number++, "n");
	  if (opt::normal_bam.empty())
	    opt::normal_bam = tmp;
	  break;
      case 'v': arg >> opt::verbose; break;
      case 'k': arg >> opt::regionFile; break;
      case 'e': arg >> opt::assemb::error_rate; break;
      case 'G': arg >> opt::refgenome; break;
	//case 'M': arg >> opt::microbegenome; break;
      case 'D': arg >> opt::dbsnp; break;
    case 'S': arg >> opt::motif_exclude; break;
      case 'r': 
	arg >> opt::rules; 
	if (opt::rules == "all")
	  opt::rules = "";
	break;
      case OPT_DISCORDANT_ONLY: opt::disc_cluster_only = true; break;
      case OPT_WRITE_EXTRACTED_READS: opt::write_extracted_reads = true; break;
    case OPT_NUM_ASSEMBLY_ROUNDS: arg >> opt::num_assembly_rounds; break;
    case OPT_KMER: opt::kmer_filter = true; break;
    case OPT_GEMCODE_AWARE: arg >> opt::gemcode; break;
      default: die= true; 
    }
  }

  // check that BAM files exist
  for (auto& b : opt::bam)
    if (!SnowTools::read_access_test(b.second)) {
      std::cerr << "Error: BAM file " << b.second << " is not readable / existant" << std::endl;
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

  if (die || help) 
    {
      std::cerr << "\n" << RUN_USAGE_MESSAGE;
      if (die)
	exit(EXIT_FAILURE);
      else 
	exit(EXIT_SUCCESS);	
    }
}

bool runBigChunk(const SnowTools::GenomicRegion& region)
{
  std::vector<SnowTools::AlignedContig> alc;
  BamReadVector all_contigs, all_microbial_contigs;

  BamReadVector blat_alignments;

  // start some counters
  SnowmanUtils::SnowTimer st;
  st.start();

  // setup for the BAM walkers
  std::map<std::string, SnowmanBamWalker> walkers;

  // loop through the input bams and get reads
  int num_t_reads = 0, num_n_reads = 0, tcount = 0, ncount = 0;
  for (auto& b : opt::bam) { 
    walkers[b.first] = __make_walkers(b.first, b.second, region, tcount, ncount);
  }
  num_t_reads = tcount; num_n_reads = ncount;
  st.stop("r");

  // collect all of the cigar maps
  CigarMap cigmap_n, cigmap_t;
  std::unordered_map<std::string, CigarMap> cigmap;
  for (auto& b : opt::bam) {
    cigmap[b.first] = walkers[b.first].cigmap;
    if (b.first.at(0) == 'n')
      for (auto& c : walkers[b.first].cigmap)
	cigmap_n[c.first] += c.second;
    else if (b.first.at(0) == 't')
      for (auto& c : walkers[b.first].cigmap)
	cigmap_t[c.first] += c.second;
  }

  MateRegionVector all_somatic_mate_regions;
  all_somatic_mate_regions.add(MateRegion(region.chr, region.pos1, region.pos2)); // add the origional, don't want to double back

  SnowTools::GRC badd;  
#ifdef GET_MATES
  const int MAX_MATE_ROUNDS = 3;
  for (int jjj = 0; jjj <  MAX_MATE_ROUNDS; ++jjj) {
    MateRegionVector normal_mate_regions = __collect_normal_mate_regions(walkers);
    
    // get the mates from somatic 3+ mate regions
    MateRegionVector tmp_somatic_mate_regions = __collect_somatic_mate_regions(walkers, normal_mate_regions);
    
    // keep regions that haven't been visited before
    MateRegionVector somatic_mate_regions;
    for (auto& s : tmp_somatic_mate_regions) {
      bool valid = true;
      
      // check if its not bad from mate region
      for (auto& br : bad_mate_regions)
	if (br.getOverlap(s))
	  valid = false;

      // check if we are allowed to lookup interchromosomal
      if (!opt::interchrom_lookup && s.chr != region.chr)
	valid = false;

      // new region overlaps with one already seen
      if (valid)
	for (auto& ss : all_somatic_mate_regions) 
	  if (s.getOverlap(ss)) 
	    valid = false;
      
      if (valid && (s.count > opt::mate_lookup_min * 2 || jjj == 0)) { // be more strict about higher rounds
	somatic_mate_regions.add(s);
	all_somatic_mate_regions.add(s);
      }
    }
    
    
    for (auto& i : somatic_mate_regions)
      log_file << "   somatic mate region " << i << " somatic count " << i.count << " origin " << i.partner << " on round " << (jjj+1) << std::endl;
    if (opt::verbose > 3)
      std::cerr << "...grabbing mate reads from " << somatic_mate_regions.size() << " regions spanning " << SnowTools::AddCommas(somatic_mate_regions.width()) << std::endl;
    
    
    // add these regions to the walker and get the reads
    tcount = 0; ncount = 0;
    
    if (somatic_mate_regions.size()) {
      for (auto& b : opt::bam)
	{
	  int oreads = walkers[b.first].reads.size();
	  walkers[b.first].m_limit = 1000;
	  walkers[b.first].setBamWalkerRegions(somatic_mate_regions.asGenomicRegionVector(), bindices[b.first]);
	  walkers[b.first].get_coverage = false;
	  walkers[b.first].readlen = (b.first == "n") ? n_params.readlen : t_params.readlen;
	  walkers[b.first].max_cov = opt::max_cov;
	  walkers[b.first].get_mate_regions = jjj != MAX_MATE_ROUNDS;
	  badd.concat(walkers[b.first].readBam(&log_file, dbsnp_filter));
	  if (b.first == "t") {
	    tcount += (walkers[b.first].reads.size() - oreads);
	  } else {
	    ncount += (walkers[b.first].reads.size() - oreads);
	  }
	}
    }
    
    num_t_reads += tcount;
    num_n_reads += ncount;
    
    if (somatic_mate_regions.size() && opt::verbose > 1)
      std::cerr << "           " << tcount << "/" << ncount << " Tumor/Normal mate reads in " << somatic_mate_regions.size() << " regions spanning " << SnowTools::AddCommas<int>(somatic_mate_regions.width()) << " bases after round " << (jjj+1) << std::endl;
    
  }
  
  
  st.stop("m");
#endif
  
  // put all of the reads together and dedupe
  // also get rid of reads where normal weird reads are too high
  std::set<std::string> dedup;
  
  //SnowTools::GRC excluded_bad_regions = __get_exclude_on_badness(walkers, region);;
  SnowTools::GRC excluded_bad_regions = SnowTools::GRC();
  std::unordered_set<std::string> excluded_bad_reads;
  
  BamReadVector bav_this;
  for (auto& b : walkers)
    for (auto& r : b.second.reads) {
      bool bad_pass = true;
      if (excluded_bad_regions.size())
	bad_pass = !excluded_bad_regions.findOverlapping(r.asGenomicRegion());
      if (!bad_pass && opt::write_extracted_reads)
	excluded_bad_reads.insert(r.GetZTag("SR"));
      if (!dedup.count(r.GetZTag("SR")) && bad_pass) {
	dedup.insert(r.GetZTag("SR"));
	bav_this.push_back(r);
      }
    }
  
  // do the kmer filtering
  BamReadVector bav_this_tum, bav_this_norm;
  if (!opt::disc_cluster_only) {
    KmerFilter kmer;
    int kmer_corrected = 0;

    if (!opt::no_assemble_normal && opt::kmer_filter)
      kmer_corrected = kmer.correctReads(bav_this);
    else {
      for (auto& r : bav_this)
	if (r.GetZTag("SR").at(0) == 't')
	  bav_this_tum.push_back(r);
	else
	  bav_this_norm.push_back(r);
      if (opt::kmer_filter)
	kmer_corrected = kmer.correctReads(bav_this_tum);
      bav_this = bav_this_tum;
      bav_this.insert(bav_this.end(), bav_this_norm.begin(), bav_this_norm.end());
    }
    st.stop("k");
    if (opt::verbose > 2 && opt::kmer_filter)
      std::cerr << "...kmer corrected " << kmer_corrected << " reads of " << bav_this.size() << std::endl;

    // write them to fasta
    if (opt::write_extracted_reads) {
      for (auto& r : bav_this)
	r2c_c << ">" << r.GetZTag("SR") << std::endl << r.QualitySequence() << std::endl;
    }

  }

  // do the discordant read clustering
  if (opt::verbose > 1)
    std::cerr << "...doing the discordant read clustering" << std::endl;

  //for (auto& r : bav_this)
  //  if (r.Qname() == "D0ENMACXX111207:1:2207:17771:106943")
  //    std::cerr << " FOUND HER " << std::endl;

  SnowTools::DiscordantClusterMap dmap = SnowTools::DiscordantCluster::clusterReads(bav_this, region);
  
  // remove the hardclips
  BamReadVector bav_tmp;
  for (auto& i : bav_this)
    if (i.NumHardClip() == 0) 
      bav_tmp.push_back(i);
  bav_this = bav_tmp;

  //for (auto& i : bav_this)
  //  std::cout << i << std::endl;

  // 
  if (opt::verbose > 3)
    for (auto& i : dmap)
      std::cerr << i.first << " " << i.second << std::endl;

  GRC grv_small = makeAssemblyRegions(region);

  if (opt::verbose > 1)
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
	log_file << "TOO MANY READS IN REGION " << bav_this.size() << "\t" << g << std::endl;
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
	  en->second.fillReadTable(brv.second) ;
	  en->second.performAssembly(opt::num_assembly_rounds);
	  all_contigs_this.insert(all_contigs_this.end(), en->second.getContigs().begin(), en->second.getContigs().end());
	}
	
      // normal non-10X data
      } else {

	SnowmanAssemblerEngine engine(name, opt::assemb::error_rate, opt::assemb::minOverlap, opt::readlen);
	// do the assemblies
	if (opt::assemb::writeASQG)
	  engine.setToWriteASQG();
	if (!opt::no_assemble_normal)
	  engine.fillReadTable(bav_this);
	else
	  engine.fillReadTable(bav_this_tum);	

	engine.performAssembly(opt::num_assembly_rounds);

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

	if ((int)i.getSeq().length() < (t_params.readlen * 1.15))
	  continue;
	
	BamReadVector ct_alignments;
	bool hardclip = false;
	main_bwa->alignSingleSequence(i.getSeq(), i.getID(), ct_alignments, hardclip, 0.90, SECONDARY_CAP);	

	if (opt::verbose > 3)
	  for (auto& i : ct_alignments)
	    std::cerr << " aligned contig: " << i << std::endl;
	
#ifdef MICROBE
	
	BamReadVector ct_plus_microbe;

	if (microbe_bwa && !SnowmanUtils::hasRepeat(i.getSeq())) {
	  
	  // do the microbial alignment
	  BamReadVector microbial_alignments;
	  bool hardclip = false;
	  microbe_bwa->alignSingleSequence(i.getSeq(), i.getID(), microbial_alignments, hardclip, 0.90, SECONDARY_CAP);

	  // if the microbe alignment is large enough and doesn't overlap human...
	  for (auto& j : microbial_alignments) {
	    // keep only long microbe alignments with decent mapq
	    if (j.NumMatchBases() >= MICROBE_MATCH_MIN && j.MapQuality() >= 10) { 
	      if (SnowmanUtils::overlapSize(j, ct_alignments) <= 20) { // keep only those where most do not overlap human
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

	// align with BLAT
	if (main_blat) {
	  std::set<std::string> sset;
	  for (auto& r : ct_alignments) {
	    std::string q = r.Qname();
	    if (!sset.count(q) && (r.NumClip() >= 20 || r.MaxDeletionBases() || r.MaxInsertionBases()))
	      main_blat->querySequence(q + "_BLAT", r.Sequence(), blat_alignments);
	    sset.insert(q);
	  }
	}
	
	// same for BLAT alignments
	if (main_blat)
	  for (auto& r : blat_alignments) {
	    //assert(main_bwa->ChrIDToName(r.ChrID()).length());
	    if (r.ChrID() >= 0) {
	      r.AddZTag("MC", main_bwa->ChrIDToName(r.ChrID()));
	    } else {
	      std::cerr << " FAILED TO PARSE CHR FOR BLAT " << std::endl;
	    }
	  }

	// remove human alignments that are not as good as microbe
	// that is, remove human alignments that intersect with microbe. 
	// We can do this because we already removed microbial alignments taht 
	// intersect too much with human. Thus, we are effectively removing human 
	// alignments that are contained within a microbial alignment

	BamReadVector human_alignments;
	for (auto& j : ct_alignments) {
	  // keep human alignments that have < 50% overlap with microbe and have >= 25 bases matched
	  if (SnowmanUtils::overlapSize(j, ct_plus_microbe) < 0.5 * j.Length() && j.NumMatchBases() >= MIN_CONTIG_MATCH) { 
	    human_alignments.push_back(j);
	    all_contigs.push_back(j);
	  }
	}

	// add the BLAT hits
	//for (auto& j : blat_alignments) {
	//   human_alignments.push_back(j);
	//   all_contigs.push_back(j);
	//}

	if (!human_alignments.size())
	  continue;

	// add in the microbe alignments
	human_alignments.insert(human_alignments.end(), ct_plus_microbe.begin(), ct_plus_microbe.end());

	// make the aligned contig object
	if (!human_alignments.size())
	  continue;

        SnowTools::AlignedContig ac(human_alignments, prefixes);

	// assign the local variable to each
	ac.checkLocal(g);

	//if ((ac.hasLocal() || ac.  && ac.hasVariant()) {
	this_alc.push_back(ac);
	usv.push_back({i.getID(), i.getSeq()});	  
        //}
      }
#else 
      for (auto& i : engine.getContigs()) 
	usv.push_back({i.getID(), i.getSeq()});
#endif

      if (!usv.size())
	continue;

      // Align the reads to the contigs with BWA-MEM
      SnowTools::BWAWrapper bw;
      //for (auto& i : usv) 
      //	std::cerr << i.name << " " << i.seq << std::endl;
      bw.constructIndex(usv);

      if (opt::verbose > 3)
	std::cerr << "...aligning " << bav_this.size() << " reads to " << this_alc.size() << " contigs " << std::endl;
      alignReadsToContigs(bw, usv, bav_this, this_alc);

      // Get contig coverage, discordant matching to contigs, etc
      for (auto& a : this_alc) {
	a.assignSupportCoverage();
	a.splitCoverage();	
	// now that we have all the break support, check that the complex breaks are OK
	a.refilterComplex(); 
	// add discordant reads support to each of the breakpoints
	a.addDiscordantCluster(dmap);
	// add in the cigar matches
	//a.checkAgainstCigarMatches(cigmap_n, cigmap_t, cigpos_n);
	a.checkAgainstCigarMatches(cigmap);
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
  for (auto& i : bp_glob)
    i.repeatFilter();
  
  // do the PON filter
  if (pon_filter) {
    for (auto & i : bp_glob) 
      i.checkPon(pon_filter);
  }

  if (dbsnp_filter && opt::dbsnp.length()) {
    if (opt::verbose > 1)
      std::cerr << "...DBSNP filtering" << std::endl;
    for (auto & i : bp_glob) {
      dbsnp_filter->queryBreakpoint(i);
    }
  }

  // filter against blacklist
  for (auto& i : bp_glob) {
    i.checkBlacklist(blacklist);
  }

  // add in the discordant clusters as breakpoints
  for (auto& i : dmap) {
    // DiscordantCluster not associated with assembly BP and has 2+ read support
    if (!i.second.hasAssociatedAssemblyContig() && (i.second.tcount + i.second.ncount) > 1) {
      SnowTools::BreakPoint tmpbp(i.second, main_bwa, dmap);
      //assert(tmpbp.b1.gr < tmpbp.b2.gr);
      bp_glob.push_back(tmpbp);
    }
  }
  
  if (opt::verbose > 1)
    std::cerr << "...deduplicating breakpoints"<< std::endl;

  // de duplicate the breakpoints
  std::sort(bp_glob.begin(), bp_glob.end());
  bp_glob.erase( std::unique( bp_glob.begin(), bp_glob.end() ), bp_glob.end() );
  
  if (opt::verbose > 2)
    std::cerr << "...integrating coverage counts" << std::endl;

  // add the coverage data to breaks for allelic fraction computation
  std::unordered_map<std::string, SnowTools::STCoverage*> covs;
  std::unordered_map<std::string, SnowTools::STCoverage*> clip_covs;
  for (auto& i : opt::bam) {
    covs[i.first] = &walkers[i.first].cov;
    clip_covs[i.first] = &walkers[i.first].clip_cov;
  }

  if (opt::verbose > 2)
    std::cerr << "...adding allelic fractions" << std::endl;

  for (auto& i : bp_glob)
    i.addCovs(covs, clip_covs);

  // score them and set somatic / germline
  if (opt::verbose > 2)
    std::cerr << "...scoring breakpoints" << std::endl;

  for (auto& i : bp_glob)
    i.scoreBreakpoint(opt::lod, opt::lod_db, opt::lod_no_db, opt::lod_germ);

  // label somatic breakpoints that intersect directly with normal as NOT somatic
  if (opt::verbose > 2)
    std::cerr << "...removing somatic overlap with normal" << std::endl;
  std::unordered_set<std::string> norm_hash;
  for (auto& i : bp_glob) // hash the normals
    if (!i.somatic_score && i.confidence == "PASS" && i.evidence == "INDEL") {
      norm_hash.insert(i.b1.hash());
      norm_hash.insert(i.b2.hash());
      norm_hash.insert(i.b1.hash(1));
      norm_hash.insert(i.b1.hash(-1));
      //norm_hash.insert(i.b1.hash(2));
      //norm_hash.insert(i.b1.hash(-2));
    }
  for (auto& i : bp_glob)  // find somatic that intersect with norm. Set somatic = 0;
    if (i.somatic_score && (norm_hash.count(i.b1.hash()) || norm_hash.count(i.b2.hash()))) {
      i.somatic_score = 0;
    }

  // remove somatic calls if they have a germline normal SV in them or indels with 
  // 2+germline normal in same contig
  if (opt::verbose > 2)
    std::cerr << "...removing somatic calls that have germline SV" << std::endl;
  std::unordered_set<std::string> bp_hash;
  for (auto& i : bp_glob) // hash the normals
    if (!i.somatic_score && i.evidence != "INDEL")
      bp_hash.insert(i.cname);
  for (auto& i : bp_glob)  // find somatic that intersect with norm. Set somatic = 0;
    if (i.somatic_score && i.num_align > 1 && bp_hash.count(i.cname))
      i.somatic_score = 0;

  // add the ref and alt tags
  if (opt::verbose > 2)
    std::cerr << "...setting ref and alt" << std::endl;
  for (auto& i : bp_glob)
    i.setRefAlt(findex, findex_viral);

  ////////////////////////////////////
  // MUTEX LOCKED
  ////////////////////////////////////
  // write to the global contig out
  pthread_mutex_lock(&snow_lock);  

  // update the bad mate regions
  bad_mate_regions.concat(badd);
    
  if (opt::verbose > 2)
    std::cerr << "...dumping files" << std::endl;

  // print out the bad reagions
  for (auto& i : excluded_bad_regions)
    bad_bed << bwalker.header()->target_name[i.chr] << "\t" << i.pos1 << "\t" << i.pos2 << std::endl;
  
  // dump the cigmaps
  for (auto& i : cigmap_n) 
    os_cigmap << i.first << "\t" << i.second << "\tN" << std::endl;
  for (auto& i : cigmap_t) 
    os_cigmap << i.first << "\t" << i.second << "\tT" << std::endl;

  if (opt::verbose > 2)
    std::cerr << "...writing alignments plot" << std::endl;
  
  // print the alignment plots
  size_t contig_counter = 0;
  for (auto& i : alc) {
    if (i.hasVariant()) {
      all_align << i << std::endl;
      ++contig_counter;
    }
  }

  if (opt::verbose > 2)
    std::cerr << "...dumping discordant" << std::endl;
  
  // send the microbe to file
  for (auto& b : all_microbial_contigs)
    b_microbe_writer.writeAlignment(b);

  // send the discordant to file
  for (auto& i : dmap)
    if (std::max(i.second.mapq1, i.second.mapq2) >= 5)
      all_disc_stream << i.second.toFileString(!opt::no_reads) << std::endl;

  // send breakpoints to file
  for (auto& i : bp_glob)
    if (i.hasMinimal() && i.confidence != "NOLOCAL") {
      os_allbps << i.toFileString(opt::no_reads) << std::endl;
      if (i.confidence == "PASS" && i.n.split == 0 && i.n.cigar == 0 && i.dc.ncount == 0 && opt::verbose > 1) {
	std::cout << "SOM  " << i << std::endl;
      } else if (i.confidence == "PASS" && opt::verbose > 1) {
	std::cout << "GER " << i << std::endl;
      }
    }

  // write ALL contigs
  if (opt::verbose > 2)
    std::cerr << "...writing contigs" << std::endl;

  if (!opt::disc_cluster_only) { 
    for (auto& i : all_contigs)
      b_allwriter.writeAlignment(i);
    for (auto& i : blat_alignments)
      blat_allwriter.writeAlignment(i);

  }
  
  // write extracted reads
  if (opt::write_extracted_reads) {
    std::cerr << walkers[opt::tumor_bam].reads.size() << std::endl;
    for (auto& r : walkers[opt::tumor_bam].reads) 
      if (!excluded_bad_reads.count(r.GetZTag("SR")))
	er_writer.writeAlignment(r);
  }
  
  // write all the to-assembly reads
  if (!opt::disc_cluster_only && opt::r2c) {
    for (auto& r : walkers[opt::tumor_bam].reads) {
      r2c_writer.writeAlignment(r);
      
      //write the corrected
      //std::string new_seq  = r.GetZTag("KC");
      //if (!new_seq.length()) {
      //new_seq = r.QualitySequence();
      //}
      //r2c_c << ">" << r.GetZTag("SR") << std::endl << new_seq << std::endl;
    }
  }

  
  st.stop("pp");

  // display the run time
  std::string rt = SnowmanUtils::runTimeString(num_t_reads, num_n_reads, contig_counter, region, bwalker.header(), st, start);
  log_file << rt;
  if (opt::verbose > 1)
    std::cerr << rt;

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

  return grv_small;

}

void alignReadsToContigs(SnowTools::BWAWrapper& bw, const SnowTools::USeqVector& usv, BamReadVector& bav_this, std::vector<SnowTools::AlignedContig>& this_alc) {
  
  if (!usv.size())
    return;
  
  for (auto i : bav_this) {
    
    BamReadVector brv;

    std::string seqr = i.QualitySequence();
    
    bool hardclip = false;
    assert(seqr.length());
    bw.alignSingleSequence(seqr, i.Qname(), brv, hardclip, 0.60, 10000);

    if (brv.size() == 0) 
      continue;
    
    // make sure we have only one alignment per contig
    std::set<std::string> cc;

    // check which ones pass
    BamReadVector bpass;
    for (auto& r : brv) {

      /*if (i.Qname() == "H75YKCCXX150615:7:2209:12997:14863") {
      	std::cerr << " rbrv " << brv.size() << " NM " << ((double)r.NumMatchBases() * 0.7) << " AS " << r.GetIntTag("AS") << std::endl;
      	std::cerr << " seqr length " << ((double)seqr.length() * 0.75) << " pos " << (r.PositionEnd() - r.Position()) << std::endl;
	std::cerr << " length pass " << ((r.PositionEnd() - r.Position()) >= ((double)seqr.length() * 0.75)) << "ccpass  " << (!cc.count(usv[r.ChrID()].name)) << std::endl;
	std::cerr << " should it pass score " << ((double)r.NumMatchBases() * 0.7 > r.GetIntTag("AS")) << std::endl;
	}*/
	
      // make sure alignment score is OK
      if ((double)r.NumMatchBases() * 0.5 > r.GetIntTag("AS")/* && i.GetZTag("SR").at(0) == 't'*/)
      	continue;
      
      bool length_pass = (r.PositionEnd() - r.Position()) >= ((double)seqr.length() * 0.75);

      if (length_pass && !cc.count(usv[r.ChrID()].name)) {

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
	    while (cc <= r.PositionEnd() && cc < (int)a.tum_cov.size()) 
	      ++a.cov[srr.substr(0,4)][cc++];

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

SnowmanBamWalker __make_walkers(const std::string& p, const std::string& b, const SnowTools::GenomicRegion& region, 
				int& tcount, int& ncount) {
  
  int readlen = p.at(0) == 't' ? t_params.readlen : n_params.readlen;
  SnowmanBamWalker walk(b, main_bwa, readlen, p, blacklist);
  walk.adapter_trim = opt::adapter_trim;
  walk.max_cov = opt::max_cov;

  if (!region.isEmpty()) 
    walk.setBamWalkerRegion(region, bindices[p]);
  walk.SetMiniRulesCollection(*mr);
  walk.m_limit = 100000; //8000;
  walk.readBam(&log_file, nullptr); 

  if (p.at(0) == 't') {
    tcount += walk.reads.size();
  } else {
    ncount += walk.reads.size();
  }

  if (opt::verbose > 1)
    std::cerr << "...read in " << walk.reads.size() << " reads from " << p << std::endl;
  
  return walk;

}

MateRegionVector __collect_normal_mate_regions(std::map<std::string, SnowmanBamWalker>& walkers) {
  
  MateRegionVector normal_mate_regions;
  for (auto& b : opt::bam)
    if (b.first.at(0) == 'n')
      normal_mate_regions.concat(walkers[b.first].mate_regions);
  normal_mate_regions.createTreeMap();
  
  return normal_mate_regions;

}

MateRegionVector __collect_somatic_mate_regions(std::map<std::string, SnowmanBamWalker>& walkers, MateRegionVector& bl) {


  MateRegionVector somatic_mate_regions;
  for (auto& b : opt::bam)
    if (b.first.at(0) == 't')
      for (auto& i : walkers[b.first].mate_regions){
	if (i.count >= opt::mate_lookup_min && !bl.findOverlapping(i) && 
	    i.chr < 24 
	    && (blacklist.size() == 0 || blacklist.findOverlapping(i) == 0))
	  somatic_mate_regions.add(i); 
      }
  
  // reduce it down
  somatic_mate_regions.mergeOverlappingIntervals();

  return somatic_mate_regions;

}

/*SnowTools::GRC __get_exclude_on_badness(std::map<std::string, SnowmanBamWalker>& walkers, const SnowTools::GenomicRegion& region) {
  
  SnowTools::STCoverage * weird_normal = &(walkers[opt::normal_bam].bad_cov);
  SnowTools::STCoverage * cov_normal   = &(walkers[opt::normal_bam].cov);

  SnowTools::GRC excluded_bad_regions;
  std::unordered_set<std::string> excluded_bad_reads;
  
  std::vector<int> excluded_on_badness(region.pos2 - region.pos1, 0);
  if (!opt::normal_bam.empty() && !opt::exome) 
    for (int i = region.pos1; i < region.pos2; ++i) {

      int vpos = i - region.pos1;
      double w_cov1 = weird_normal->getCoverageAtPosition(region.chr, i);
      double n_cov1 = cov_normal->getCoverageAtPosition(region.chr, i);
      
      // filter out bad positions by looking at bad reads in normal
      if (std::max(n_cov1,w_cov1) >= 10)
	if (( w_cov1 >= n_params.mean_cov * n_params.frac_bad * 4) || (n_cov1 >= 6 * n_params.mean_cov)) {
	  for (int j = -5; j <= 5; ++j)
	    if (vpos + j >= 0 && vpos + j < (int)excluded_on_badness.size())
	      excluded_on_badness[vpos+j] = w_cov1;      
	}
      
    }

    // print out the bad reagions
  if (excluded_on_badness.size()) {
    int start = -1, end = -1; 
    for (size_t i = 0; i < excluded_on_badness.size(); ++i) {
      if (excluded_on_badness[i] == 0)
	continue;
      else if (start < 0) 
	start = i;
      else if (i - end >= 2 ) {
      // && end - start >= 50) {
	excluded_bad_regions.add(SnowTools::GenomicRegion(region.chr, start+region.pos1, end+region.pos1));
	//bad_bed << bwalker.header()->target_name[region.chr] << "\t" << (start+region.pos1) << "\t" << (end+region.pos1) << std::endl;
	start = i;
      }
      end = i;
    }
    if (start > 0)
      excluded_bad_regions.add(SnowTools::GenomicRegion(region.chr, start+region.pos1, end+region.pos1));    
  }
  excluded_bad_regions.createTreeMap();

  return excluded_bad_regions;

  
}*/

/*
void run_test_assembly() {

  BamReadVector bav_test;
  std::string inputfile = "/xchip/gistic/Jeremiah/Data/HCC1143_mem.bam";
  std::string dum2 = "dum";
  SnowmanBamWalker walktest(inputfile, main_bwa, 101, dum2, blacklist);
  SnowTools::GenomicRegion grtest("1:1,000,000-1,000,500", walktest.header());
  walktest.setBamWalkerRegion(grtest);
  walktest.readBam(nullptr, nullptr);
  for (auto& r : walktest.reads)
    bav_test.push_back(r);
  SnowmanAssemblerEngine enginetest("test", opt::assemb::error_rate, 35, 101);
  enginetest.fillReadTable(bav_test);

  int seedLength, seedStride;
  enginetest.calculateSeedParameters(101, 35, seedLength, seedStride);
  std::cerr << "seed length " << seedLength << std::endl;;

  BamReadVector bav_original = bav_test;
  for (int i = 0; i < 100; ++i) {
    KmerFilter kmer;
    //int kmer_corrected = kmer.correctReads(bav_test);
    if (i % 1 == 0)
      std::cerr << "Assembling " << bav_test.size() << " read on " << SnowTools::AddCommas(i) << std::endl;
    enginetest.performAssembly(1);
    enginetest.clearContigs();
    //bav_test = bav_original;
  }
  
  exit(0);

  }*/
