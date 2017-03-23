#include "run_svaba.h"

#include <getopt.h>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <map>
#include <vector>
#include <cassert>

#include "SeqLib/ReadFilter.h"
#include "KmerFilter.h"
#include "vcf.h"
#include "DBSnpFilter.h"
#include "svabaUtils.h"
#include "LearnBamParams.h"
#include "SeqLib/BFC.h"

#define THREAD_READ_LIMIT  20000
#define THREAD_CONTIG_LIMIT 250

// minimum number of reads to support even reporting dscrd cluster 
// (if not assocaited with assembly contig)
#define MIN_DSCRD_READS_DSCRD_ONLY 3 

// useful replace function
std::string myreplace(std::string &s,
                      std::string toReplace,
                      std::string replaceWith)
{
  return(s.replace(s.find(toReplace), toReplace.length(), replaceWith));
}

// output files
static ogzstream all_align, os_allbps, os_discordant, os_corrected;
static std::ofstream log_file, bad_bed;
static std::stringstream ss; // initalize a string stream once

#define WRITELOG(msg, toerr, tolog)   \
  { if (tolog) log_file  << (msg) << std::endl;	\
    if (toerr) std::cerr << (msg) << std::endl; };

#define ERROR_EXIT(msg) { std::cerr << (msg) << std::endl; exit(EXIT_FAILURE); }

#define MIN_CONTIG_MATCH 35
#define MATE_LOOKUP_MIN 3
#define SECONDARY_CAP 10
#define MAX_MATE_ROUNDS 1
#define MATE_REGION_LOOKUP_LIMIT 400
#define MAX_NUM_MATE_WINDOWS 50000000

#define GERMLINE_CNV_PAD 10
#define WINDOW_PAD 500
#define MICROBE_MATCH_MIN 50
#define GET_MATES 1
#define MICROBE 1
#define LARGE_INTRA_LOOKUP_LIMIT 50000
#define SECONDARY_FRAC 0.90

// if a local alignment has < MIN_CLIP_FOR_LOCAL clips
// then it has a good local (and is not an SV candidate contig)
#define MIN_CLIP_FOR_LOCAL 40
// if local alignment to assembly has > MAX_NM_FOR_LOCAL
// NM, then dont' consider it a strong local match
#define MAX_NM_FOR_LOCAL 10 

static SeqLib::RefGenome * ref_genome, * ref_genome_viral;
static std::unordered_map<std::string, BamParamsMap> params_map; // key is bam id (t000), value is map with read group as key
static SeqLib::BamHeader bwa_header, viral_header;

static std::set<std::string> prefixes;

static int min_dscrd_size_for_variant = 0; // set a min size for what we can call with discordant reads only. 
// something like max(mean + 3*sd) for all read groups

static std::unordered_map<std::string, int> min_isize_for_disc;

static SeqLib::BamHeader b_header; // header for main bam
static SeqLib::BamReader b_reader; // reader for the main bam
static SeqLib::BamWriter er_writer, b_microbe_writer, b_contig_writer;
static SeqLib::BWAWrapper * microbe_bwa = nullptr;
static SeqLib::BWAWrapper * main_bwa = nullptr;
static SeqLib::Filter::ReadFilterCollection * mr;
static SeqLib::GRC blacklist, germline_svs, simple_seq;
static DBSnpFilter * dbsnp_filter;
static SeqLib::GRC file_regions, regions_torun;

// mutex and time
static pthread_mutex_t snow_lock;
static struct timespec start;

// learned value 
static int max_mapq_possible;
static std::string args = "svaba "; // hold string of what the input args were
static int32_t readlen = 0;

namespace opt {

  // SGA options
  namespace sga {
    static int minOverlap = 0;
    static float error_rate = 0; 
    static bool writeASQG = false;
    static int num_assembly_rounds = 3;
  }

  // error correction options
  static std::string ec_correct_type = "f";
  static double ec_subsample = 0.50;

  // run in single end mode?
  bool single_end = false;

  // output options
  static bool write_extracted_reads = false;
  static bool write_corrected_reads = false;
  static bool zip = false;
  static bool read_tracking = false; // turn on output of qnames
  static bool all_contigs = false;   // output all contigs
  static bool no_unfiltered = false; // don't output unfiltered variants

  // discordant clustering params
  static double sd_disc_cutoff = 3.92;
  static bool disc_cluster_only = false;

  // BWA MEM params
  namespace bwa {
    static int gap_open_penalty = 32;
    static int gap_extension_penalty = 1;
    static int mismatch_penalty = 18;
    static int sequence_match_score = 2;
    static int zdrop = 100;
    static int bandwidth = 1000;
    static float reseed_trigger = 1.5;
    static int clip3_pen = 5;
    static int clip5_pen = 5;
  }

  // parameters for filtering / getting reads
  static std::string rules = "{\"global\" : {\"duplicate\" : false, \"qcfail\" : false}, \"\" : { \"rules\" : [FRRULES,{\"rr\" : true},{\"ff\" : true}, {\"rf\" : true}, {\"ic\" : true}, {\"clip\" : 5, \"length\" : READLENLIM}, {\"ins\" : true}, {\"del\" : true}, {\"mapped\": true , \"mate_mapped\" : false}, {\"mate_mapped\" : true, \"mapped\" : false}]}}";  
  static int max_cov = 100;
  static size_t mate_lookup_min = 3;
  static bool interchrom_lookup = true;
  static int32_t max_reads_per_assembly = -1; // set default of 10000 in parseRunOptions

  // additional optional params
  static int chunk = 25000;
  static std::string regionFile;  // region to run on
  static std::string analysis_id = "no_id";
  static int num_to_sample = 2000000;  // num to learn from (eg isize distribution)

  // runtime parameters
  static int verbose = 0;
  static int numThreads = 1;
  static bool hp = false; // should run in highly-parallel mode? (no file dump til end)

  // data
  static BamMap bam;
  static std::string refgenome = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
  static std::string microbegenome; // = "/xchip/gistic/Jeremiah/Projects/SnowmanFilters/viral.1.1.genomic_ns.fna";
  static std::string simple_file; //  file of simple repeats as a filter
  static std::string blacklist; // = "/xchip/gistic/Jeremiah/Projects/HengLiMask/um75-hs37d5.bed.gz";
  static std::string germline_sv_file;
  static std::string dbsnp; // = "/xchip/gistic/Jeremiah/SnowmanFilters/dbsnp_138.b37_indel.vcf";
  static std::string main_bam = "-"; // the main bam

  // optimize defaults for single sample mode
  static bool germline = false; 
  
  // indel probability cutoffs
  static double lod = 8; // LOD that variant is not ref
  static double lod_db = 6; // same, but at DB snp site (want lower bc we have prior)
  static double lod_somatic = 6; // LOD that normal is REF
  static double lod_somatic_db = 10; // same, but at DBSNP (want higher bc we have prior that its germline)
  static double scale_error = 1; // how much to emphasize erorrs. 1 is standard. 0 is assume no errors

}

enum { 
  OPT_ASQG,
  OPT_HP,
  OPT_READLEN,
  OPT_LOD,
  OPT_LOD_DB,
  OPT_LOD_SOMATIC,
  OPT_LOD_SOMATIC_DB,
  OPT_DISC_CLUSTER_ONLY,
  OPT_READ_TRACK,
  OPT_EC_SUBSAMPLE,
  OPT_DISCORDANT_ONLY,
  OPT_WRITE_EXTRACTED_READS,
  OPT_NUM_ASSEMBLY_ROUNDS,
  OPT_NUM_TO_SAMPLE,
  OPT_GAP_OPEN,
  OPT_MATCH_SCORE,
  OPT_GAP_EXTENSION,
  OPT_MISMATCH,
  OPT_ZDROP,
  OPT_BANDWIDTH,
  OPT_RESEED_TRIGGER,
  OPT_CLIP5,
  OPT_CLIP3,
  OPT_GERMLINE,
  OPT_SCALE_ERRORS,
  OPT_NO_UNFILTERED
};

static const char* shortopts = "hzIAt:n:p:v:r:G:e:k:c:a:m:B:D:Y:S:L:s:V:R:K:E:C:x:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "tumor-bam",               required_argument, NULL, 't' },
  { "germline",                no_argument, NULL, OPT_GERMLINE },  
  { "all-contigs",             no_argument, NULL, 'A' },  
  { "panel-of-normals",        required_argument, NULL, 'P' },
  { "id-string",               required_argument, NULL, 'a' },
  { "hp",                      no_argument, NULL, OPT_HP },
  { "normal-bam",              required_argument, NULL, 'n' },
  { "threads",                 required_argument, NULL, 'p' },
  { "no-unfiltered",           no_argument, NULL, OPT_NO_UNFILTERED },
  { "chunk-size",              required_argument, NULL, 'c' },
  { "region-file",             required_argument, NULL, 'k' },
  { "rules",                   required_argument, NULL, 'r' },
  { "reference-genome",        required_argument, NULL, 'G' },
  { "microbial-genome",        required_argument, NULL, 'Y' },
  { "min-overlap",             required_argument, NULL, 'm' },
  { "dbsnp-vcf",               required_argument, NULL, 'D' },
  { "disc-sd-cutoff",          required_argument, NULL, 's' },
  { "mate-lookup-min",         required_argument, NULL, 'L' },
  { "germline-sv-database",    required_argument, NULL, 'V' },
  { "simple-seq-database",     required_argument, NULL, 'R' },
  { "g-zip",                   no_argument, NULL, 'z' },
  { "no-interchrom-lookup",    no_argument, NULL, 'I' },
  { "read-tracking",           no_argument, NULL, OPT_READ_TRACK },
  { "gap-open-penalty",        required_argument, NULL, OPT_GAP_OPEN },
  { "readlen",                 required_argument, NULL, OPT_READLEN },
  { "ec-subsample",            required_argument, NULL, 'E' },
  { "bwa-match-score",         required_argument, NULL, OPT_MATCH_SCORE },
  { "gap-extension-penalty",   required_argument, NULL, OPT_GAP_EXTENSION },
  { "mismatch-penalty",        required_argument, NULL, OPT_MISMATCH },
  { "z-dropoff",               required_argument, NULL, OPT_ZDROP },
  { "reseed-trigger",          required_argument, NULL, OPT_RESEED_TRIGGER },
  { "penalty-clip-3",          required_argument, NULL, OPT_CLIP3 },
  { "penalty-clip-5",          required_argument, NULL, OPT_CLIP5 },
  { "bandwidth",               required_argument, NULL, OPT_BANDWIDTH },
  { "write-extracted-reads",   no_argument, NULL, OPT_WRITE_EXTRACTED_READS },
  { "lod",                     required_argument, NULL, OPT_LOD },
  { "lod-dbsnp",               required_argument, NULL, OPT_LOD_DB },
  { "lod-somatic",             required_argument, NULL, OPT_LOD_SOMATIC },
  { "lod-somatic-dbsnp",       required_argument, NULL, OPT_LOD_SOMATIC_DB },
  { "scale-errors",            required_argument, NULL, OPT_SCALE_ERRORS },
  { "discordant-only",         no_argument, NULL, OPT_DISCORDANT_ONLY },
  { "num-to-sample",           required_argument, NULL, OPT_NUM_TO_SAMPLE },
  { "write-asqg",              no_argument, NULL, OPT_ASQG   },
  { "ec-correct-type",         required_argument, NULL, 'K'},
  { "error-rate",              required_argument, NULL, 'e'},
  { "verbose",                 required_argument, NULL, 'v' },
  { "blacklist",               required_argument, NULL, 'B' },
  { "max-coverage",            required_argument, NULL, 'C' },
  { "max-reads",               required_argument, NULL, 'x' },
  { "num-assembly-rounds",     required_argument, NULL, OPT_NUM_ASSEMBLY_ROUNDS },
  { NULL, 0, NULL, 0 }
};

static const char *RUN_USAGE_MESSAGE =
"Usage: svaba run -t <BAM/SAM/CRAM> -G <reference> -a myid [OPTIONS]\n\n"
"  Description: SV and indel detection using rolling SGA assembly and BWA-MEM realignment\n"
"\n"
"  General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 0 \n"
"  -h, --help                           Display this help and exit\n"
"  -p, --threads                        Use NUM threads to run svaba. Default: 1\n"
"  -a, --id-string                      String specifying the analysis ID to be used as part of ID common.\n"
"  Main input\n"
"  -G, --reference-genome               Path to indexed reference genome to be used by BWA-MEM.\n"
"  -t, --case-bam                       Case BAM/CRAM/SAM file (eg tumor). Can input multiple.\n"
"  -n, --control-bam                    (optional) Control BAM/CRAM/SAM file (eg normal). Can input multiple.\n"
"  -k, --region                         Run on targeted intervals. Accepts BED file or Samtools-style string\n"
"      --germline                       Sets recommended settings for case-only analysis (eg germline). (-I, -L5, assembles NM >= 3 reads)\n"
"  Variant filtering and classification\n"
"      --lod                            LOD cutoff to classify indel as non-REF (tests AF=0 vs AF=MaxLikelihood(AF)) [8]\n"
"      --lod-dbsnp                      LOD cutoff to classify indel as non-REF (tests AF=0 vs AF=MaxLikelihood(AF)) at DBSnp indel site [5]\n"
"      --lod-somatic                    LOD cutoff to classify indel as somatic (tests AF=0 in normal vs AF=ML(0.5)) [2.5]\n"
"      --lod-somatic-dbsnp              LOD cutoff to classify indel as somatic (tests AF=0 in normal vs AF=ML(0.5)) at DBSnp indel site [4]\n"
"      --scale-errors                   Scale the priors that a site is artifact at given repeat count. 0 means assume low (const) error rate [1]\n"
"  Additional options\n"                       
"  -L, --mate-lookup-min                Minimum number of somatic reads required to attempt mate-region lookup [3]\n"
"  -s, --disc-sd-cutoff                 Number of standard deviations of calculated insert-size distribution to consider discordant. [3.92]\n"
"  -c, --chunk-size                     Size of a local assembly window (in bp). Set 0 for whole-BAM in one assembly. [25000]\n"
"  -x, --max-reads                      Max total read count to read in from assembly region. Set 0 to turn off. [10000]\n"
"  -C, --max-coverage                   Max read coverage to send to assembler (per BAM). Subsample reads if exceeded. [500]\n"
"      --no-interchrom-lookup           Skip mate lookup for inter-chr candidate events. Reduces power for translocations but less I/O.\n"
"      --discordant-only                Only run the discordant read clustering module, skip assembly. \n"
"      --num-assembly-rounds            Run assembler multiple times. > 1 will bootstrap the assembly. [2]\n"
"      --num-to-sample                  When learning about inputs, number of reads to sample. [2,000,000]\n"
"      --hp                             Highly parallel. Don't write output until completely done. More memory, but avoids all thread-locks.\n"
"  Output options\n"
"  -z, --g-zip                          Gzip and tabix the output VCF files. [off]\n"
"  -A, --all-contigs                    Output all contigs that were assembled, regardless of mapping or length. [off]\n"
"      --read-tracking                  Track supporting reads by qname. Increases file sizes. [off]\n"
"      --write-extracted-reads          For the case BAM, write reads sent to assembly to a BAM file. [off]\n"
"  Optional external database\n"
"  -D, --dbsnp-vcf                      DBsnp database (VCF) to compare indels against\n"
"  -B, --blacklist                      BED-file with blacklisted regions to not extract any reads from.\n"
"  -Y, --microbial-genome               Path to indexed reference genome of microbial sequences to be used by BWA-MEM to filter reads.\n"
"  -V, --germline-sv-database           BED file containing sites of known germline SVs. Used as additional filter for somatic SV detection\n"
"  -R, --simple-seq-database            BED file containing sites of simple DNA that can confuse the contig re-alignment.\n"
"  Assembly and EC params\n"
"  -m, --min-overlap                    Minimum read overlap, an SGA parameter. Default: 0.4* readlength\n"
"  -e, --error-rate                     Fractional difference two reads can have to overlap. See SGA. 0 is fast, but requires error correcting. [0]\n"
"  -K, --ec-correct-type                (f) Fermi-kit BFC correction, (s) Kmer-correction from SGA, (0) no correction (then suggest non-zero -e) [f]\n"
"  -E, --ec-subsample                   Learn from fraction of non-weird reads during error-correction. Lower number = faster compute [0.5]\n"
"      --write-asqg                     Output an ASQG graph file for each assembly window.\n"
"  BWA-MEM alignment params\n"
"      --bwa-match-score                Set the BWA-MEM match score. BWA-MEM -A [2]\n"
"      --gap-open-penalty               Set the BWA-MEM gap open penalty for contig to genome alignments. BWA-MEM -O [32]\n"
"      --gap-extension-penalty          Set the BWA-MEM gap extension penalty for contig to genome alignments. BWA-MEM -E [1]\n"
"      --mismatch-penalty               Set the BWA-MEM mismatch penalty for contig to genome alignments. BWA-MEM -b [18]\n"
"      --bandwidth                      Set the BWA-MEM SW alignment bandwidth for contig to genome alignments. BWA-MEM -w [1000]\n"
"      --z-dropoff                      Set the BWA-MEM SW alignment Z-dropoff for contig to genome alignments. BWA-MEM -d [100]\n"
"      --reseed-trigger                 Set the BWA-MEM reseed trigger for reseeding mems for contig to genome alignments. BWA-MEM -r [1.5]\n"
"      --penalty-clip-3                 Set the BWA-MEM penalty for 3' clipping. [5]\n"
"      --penalty-clip-5                 Set the BWA-MEM penalty for 5' clipping. [5]\n"
"\n";

void runsvaba(int argc, char** argv) {

  parseRunOptions(argc, argv);

  // open the output streams
  svabaUtils::fopen(opt::analysis_id + ".log", log_file);
  //  svabaUtils::fopen(opt::analysis_id + ".bad_mate_regions.bed", bad_bed);

  // will check later if reads have different max mapq or readlen
  bool diff_read_len = false;
  bool diff_mapq = false;

  // set the germline parameters 
  if (opt::germline) {
    if (!opt::rules.empty())
      opt::rules = "{\"global\" : {\"duplicate\" : false, \"qcfail\" : false}, \"\" : { \"rules\" : [FRRULES,{\"rr\" : true},{\"ff\" : true}, {\"rf\" : true}, {\"ic\" : true}, {\"clip\" : 5, \"length\" : READLENLIM}, {\"ins\" : true}, {\"del\" : true}, {\"mapped\": true , \"mate_mapped\" : false}, {\"mate_mapped\" : true, \"mapped\" : false}, {\"nm\" : [3,0]}]}}";
    opt::interchrom_lookup = false;
    opt::mate_lookup_min = 5;
  } else {
    if (opt::sd_disc_cutoff==3.96)
      opt::sd_disc_cutoff=6;
  }

  // set the rules to skip read learning if doing stdin
  if (opt::main_bam == "-" && !opt::rules.empty()) 
    opt::rules = "{\"global\" : {\"duplicate\" : false, \"qcfail\" : false}, \"\" : { \"rules\" : [{\"isize\" : 800}, {\"rr\" : true},{\"ff\" : true}, {\"rf\" : true}, {\"ic\" : true}, {\"clip\" : 5, \"length\" : 30}, {\"ins\" : true}, {\"del\" : true}, {\"mapped\": true , \"mate_mapped\" : false}, {\"mate_mapped\" : true, \"mapped\" : false}, {\"nm\" : [3,0]}]}}";
  
  std::cerr << 
    "-----------------------------------------------------------" << std::endl << 
    "--- Running svaba SV and indel detection on " << SeqLib::AddCommas(opt::numThreads) << 
    " threads --" <<(opt::numThreads >= 10 ? "" : "-") << std::endl <<
    "---    (inspect *.log for real-time progress updates)   ---" << std::endl << 
    "-----------------------------------------------------------" << std::endl;
  
  ss << 
    "***************************** PARAMS ****************************" << std::endl << 
    "    DBSNP Database file: " << opt::dbsnp << std::endl << 
    "    Max cov to assemble: " << opt::max_cov << std::endl <<
    "    Error correction mode: " << opt::ec_correct_type << std::endl << 
    "    Subsample-rate for correction learning: " + std::to_string(opt::ec_subsample) << std::endl;
    ss << 
      "    ErrorRate: " << (opt::sga::error_rate < 0.001f ? "EXACT (0)" : std::to_string(opt::sga::error_rate)) << std::endl << 
      "    Num assembly rounds: " << opt::sga::num_assembly_rounds << std::endl;
  ss << 
    "    Num reads to sample: " << opt::num_to_sample << std::endl << 
    "    Discordant read extract SD cutoff:  " << opt::sd_disc_cutoff << std::endl << 
    "    Discordant cluster std-dev cutoff:  " << opt::sd_disc_cutoff << std::endl << 
    "    Minimum number of reads for mate lookup " << opt::mate_lookup_min << std::endl <<
    "    LOD cutoff (non-REF):            " << opt::lod << std::endl << 
    "    LOD cutoff (non-REF, at DBSNP):  " << opt::lod_db << std::endl << 
    "    LOD somatic cutoff:              " << opt::lod_somatic << std::endl << 
    "    LOD somatic cutoff (at DBSNP):   " << opt::lod_somatic_db << std::endl <<
    "    BWA-MEM params:" << std::endl <<
    "      Gap open penalty: " << opt::bwa::gap_open_penalty << std::endl << 
    "      Gap extension penalty: " << opt::bwa::gap_extension_penalty << std::endl <<
    "      Mismatch penalty: " << opt::bwa::mismatch_penalty << std::endl <<
    "      Aligment bandwidth: " << opt::bwa::bandwidth << std::endl <<
    "      Z-dropoff: " << opt::bwa::zdrop << std::endl <<
    "      Clip 3 penalty: " << opt::bwa::clip3_pen << std::endl <<
    "      Clip 5 penalty: " << opt::bwa::clip5_pen << std::endl <<
    "      Reseed trigger: " << opt::bwa::reseed_trigger << std::endl <<
    "      Sequence match score: " << opt::bwa::sequence_match_score << std::endl;

  if (opt::sga::writeASQG)
    ss << "    Writing ASQG files. Suggest running R/snow-asqg2pdf.R -i <my.asqg> -o graph.pdf" << std::endl;
  if (opt::write_extracted_reads)
    ss << "    Writing fasta of error-corrected reads." << std::endl;
  if (opt::disc_cluster_only)
    ss << "    ######## ONLY DISCORDANT READ CLUSTERING. NO ASSEMBLY ##############" << std::endl;
  if (!opt::interchrom_lookup)
    ss << "    ######## NOT LOOKING UP MATES FOR INTERCHROMOSOMAL #################" << std::endl;
  ss <<
    "*****************************************************************" << std::endl;	  
  WRITELOG(ss.str(), opt::verbose >= 1, true);
  ss.str(std::string());
  
  // make one anyways, we check if its empty later
  ref_genome_viral = new SeqLib::RefGenome;
  microbe_bwa = nullptr;
  
  // open the microbe genome
  
  if (!opt::microbegenome.empty()) {
    WRITELOG("...loading the microbe reference sequence", opt::verbose > 0, true)
    microbe_bwa = new SeqLib::BWAWrapper();
    svabaUtils::__open_index_and_writer(opt::microbegenome, microbe_bwa, opt::analysis_id + ".microbe.bam", b_microbe_writer, ref_genome_viral, viral_header);  
  }

  // open the main bam to get header info
  if (!b_reader.Open(opt::main_bam)) {
    if (opt::main_bam == "-")
      std::cerr << "ERROR: Cannot read from stdin" << std::endl;
    else
      std::cerr << "ERROR: Cannot open main bam file: " << opt::main_bam << std::endl;
    exit(EXIT_FAILURE);
  }

  // then open the main header
  b_header = b_reader.Header();
  if (b_header.isEmpty()) {
    std::cerr << "ERROR: empty header in main bam file" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  // open some writer bams
  if (opt::write_extracted_reads) // open the extracted reads writer
    svabaUtils::__openWriterBam(b_header, opt::analysis_id + ".extracted.reads.bam", er_writer);    

  // open the blacklists
  svabaUtils::__open_bed(opt::blacklist, blacklist, b_header);
  if (blacklist.size())
    ss << "...loaded " << blacklist.size() << " blacklist regions from " << opt::blacklist << std::endl;

  // open the germline sv database
  svabaUtils::__open_bed(opt::germline_sv_file, germline_svs, b_header);
  if (germline_svs.size())
    ss << "...loaded " << germline_svs.size() << " germline SVs from " << opt::germline_sv_file << std::endl;

  // open the simple seq database
  svabaUtils::__open_bed(opt::simple_file, simple_seq, b_header);
  if (simple_seq.size())
    ss << "...loaded " << simple_seq.size() << " simple sequence regions from " << opt::simple_file << std::endl;

  // open the DBSnpFilter
  if (opt::dbsnp.length()) {
    WRITELOG("...loading the DBsnp database", opt::verbose > 0, true)
      dbsnp_filter = new DBSnpFilter(opt::dbsnp, b_header);
    WRITELOG("...loaded DBsnp database", opt::verbose > 0, true)
  }

  // needed for aligned contig
  for (auto& b : opt::bam)
    prefixes.insert(b.first);

  if (opt::main_bam == "-")
    opt::single_end = true;

  // parse the region file, count number of jobs
  int num_jobs = svabaUtils::countJobs(opt::regionFile, file_regions, regions_torun,
					 b_header, opt::chunk, WINDOW_PAD); 

  // no learning for stdin or single-end mode
  if (opt::single_end) 
    goto afterlearn;

  // learn bam
  min_dscrd_size_for_variant = 0; // set a min size for what we can call with discordant reads only. 
  for (auto& b : opt::bam) {
    LearnBamParams parm(b.second);
    params_map[b.first] = BamParamsMap();
    parm.learnParams(params_map[b.first], opt::num_to_sample);
    for (auto& i : params_map[b.first]) {
      readlen = std::max(readlen, i.second.readlen);
      max_mapq_possible = std::max(max_mapq_possible, i.second.max_mapq);
      min_dscrd_size_for_variant = std::max(min_dscrd_size_for_variant, (int)std::floor(i.second.mean_isize + i.second.sd_isize * opt::sd_disc_cutoff)); 
    }

    ss << "BAM PARAMS FOR: " << b.first << "--" << b.second << std::endl;
    for (auto& i : params_map[b.first])
      ss << i.second << std::endl;
    ss << " min_dscrd_size_for_variant " << min_dscrd_size_for_variant << std::endl;
  }

  // check if differing read lengths or max mapq
  for (auto& a : params_map) {
    for (auto& i : a.second) {
      if (i.second.readlen != readlen)
	diff_read_len = true;
      if (i.second.max_mapq != max_mapq_possible)
	diff_mapq = true;
    } 
  }
  for (auto& a : params_map) {
    for (auto& i : a.second) {
      if (diff_read_len)
	std::cerr << "!!!! WARNING. Multiple readlengths mixed: " << i.first << "--" << i.second.readlen << " max readlen " << readlen << std::endl;
      if (diff_mapq)
	std::cerr << "!!!! WARNING. Multiple max mapq mixed: " << i.first << "--" << i.second.max_mapq << " max possible " << max_mapq_possible << std::endl;
    }
  }
  ss << "...min discordant-only variant size " << min_dscrd_size_for_variant << std::endl;

  // set the min overlap
  if (!opt::sga::minOverlap) 
    opt::sga::minOverlap = (0.6 * readlen) < 30 ? 30 : 0.6 * readlen;

  ss << "...found read length of " << readlen << ". Min Overlap is " << opt::sga::minOverlap << std::endl;
  ss << "...max read MAPQ detected: " << max_mapq_possible << std::endl;
  WRITELOG(ss.str(), opt::verbose > 1, true);
  ss.str(std::string());

afterlearn: 

  // if didn't learn anything, then make sure we set an overlap
  if (!readlen && !opt::sga::minOverlap) {
    std::cerr << "ERROR: Didn't learn from reads (stdin assembly?). Need to explicitly set readlen with --readlen" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  // get the seed length for printing
  int seedLength, seedStride;
  svabaAssemblerEngine enginetest("test", opt::sga::error_rate, opt::sga::minOverlap, readlen);
  enginetest.calculateSeedParameters(readlen, opt::sga::minOverlap, seedLength, seedStride);
  WRITELOG("...calculated seed size for error rate of " + std::to_string(opt::sga::error_rate) + " and read length " +
	   std::to_string(readlen) + " is " + std::to_string(seedLength), opt::verbose, true);

  // open the human reference
  WRITELOG("...loading the human reference sequence for BWA", opt::verbose, true);
  main_bwa = new SeqLib::BWAWrapper();
  main_bwa->SetAScore(opt::bwa::sequence_match_score);
  main_bwa->SetGapOpen(opt::bwa::gap_open_penalty);
  main_bwa->SetGapExtension(opt::bwa::gap_extension_penalty);
  main_bwa->SetMismatchPenalty(opt::bwa::mismatch_penalty);
  main_bwa->SetZDropoff(opt::bwa::zdrop);
  main_bwa->SetBandwidth(opt::bwa::bandwidth);
  main_bwa->SetReseedTrigger(opt::bwa::reseed_trigger);
  main_bwa->Set3primeClippingPenalty(opt::bwa::clip3_pen);
  main_bwa->Set5primeClippingPenalty(opt::bwa::clip5_pen);

  // open the reference for reading seqeuence
  ref_genome = new SeqLib::RefGenome;

  svabaUtils::__open_index_and_writer(opt::refgenome, main_bwa, opt::analysis_id + ".contigs.bam", b_contig_writer, ref_genome, bwa_header);
  if (ref_genome->IsEmpty()) {
    std::cerr << "ERROR: Unable to open index file: " << opt::refgenome << std::endl;
    exit(EXIT_FAILURE);
   }
  
  if (num_jobs) {
    WRITELOG("...running on " + SeqLib::AddCommas(num_jobs) + " chunks", opt::verbose, true);
  } else {
    WRITELOG("Chunk was <= 0: READING IN WHOLE GENOME AT ONCE", opt::verbose, true);
  }

  // loop through and construct the readgroup rules
  std::stringstream ss_rules;
  std::unordered_set<std::string> rg_seen;

  // set the insert size bounds for each read group
  for (auto& a : params_map) {
    for (auto& i : a.second) {
      if (!rg_seen.count(i.second.read_group)) {
	int mi = std::floor(i.second.mean_isize + i.second.sd_isize * opt::sd_disc_cutoff);
	ss_rules << "{\"isize\" : [ " << mi << ",0], \"rg\" : \"" << i.second.read_group << "\"},";
	rg_seen.insert(i.second.read_group);
	min_isize_for_disc.insert(std::pair<std::string, int>(i.second.read_group, mi));
      }
    } 
  }

  // format the rules JSON from the above string
  if (opt::rules.find("FRRULES") != std::string::npos) {
    std::string string_rules = ss_rules.str();
    if (!string_rules.empty()) // cut last comma
      string_rules = string_rules.substr(0, string_rules.length() - 1);
    else
      string_rules = "[0,0]";
    opt::rules = myreplace(opt::rules, "FRRULES", string_rules);
  }

  // set min length for clips
  if (opt::rules.find("READLENLIM") != std::string::npos) {
    if (readlen == 0)
      readlen = 30; // set some small default, in case we didn't learn the bam
    opt::rules = myreplace(opt::rules, "READLENLIM", std::to_string((int) (readlen * 0.4)));
  }

  // set the ReadFilterCollection to be applied to each region
  WRITELOG(opt::rules, opt::verbose > 1, true);
  mr = new SeqLib::Filter::ReadFilterCollection(opt::rules, bwa_header);
  WRITELOG(*mr, opt::verbose > 1, true);

  // override the number of threads if need
  num_jobs = (num_jobs == 0) ? 1 : num_jobs;
  opt::numThreads = std::min(num_jobs, opt::numThreads);

  // open the mutex
  if (pthread_mutex_init(&snow_lock, NULL) != 0) {
    std::cerr << "\n mutex init failed\n";
    exit(EXIT_FAILURE);
  }

  // open the text files
  svabaUtils::fopen(opt::analysis_id + ".alignments.txt.gz", all_align);
  svabaUtils::fopen(opt::analysis_id + ".bps.txt.gz", os_allbps);
  svabaUtils::fopen(opt::analysis_id + ".discordant.txt.gz", os_discordant);
  if (opt::write_extracted_reads) 
    svabaUtils::fopen(opt::analysis_id + ".corrected.fa.gz", os_corrected); 
  
  // write the headers to the text files
  os_allbps << BreakPoint::header();
  for (auto& b : opt::bam) 
    os_allbps << "\t" << b.first << "_" << b.second;
  os_allbps << std::endl;
  os_discordant << DiscordantCluster::header() << std::endl;

  // put args into string for VCF later
  for (int i = 0; i < argc; ++i)
    args += std::string(argv[i]) + " ";

  // start the timer
#ifndef __APPLE__
  clock_gettime(CLOCK_MONOTONIC, &start);
#endif

  // send the jobs to the queue
  WRITELOG("--- Loaded non-read data. Starting detection pipeline", true, true);
  sendThreads(regions_torun);

  if (microbe_bwa)
    delete microbe_bwa;

  // dump the bad bed regions
  /*SeqLib::GRC bad_mate_regions;
  for (const auto& i : )
    bad_mate_regions.Concat(i.second);
  bad_mate_regions.MergeOverlappingIntervals();
  bad_mate_regions.CoordinateSort();
  for (auto& i : bad_mate_regions)
    try { // throws out of range if there is mismatch between main header and other headers. just be graceful here
      bad_bed << i.ChrName(b_header) << "\t" << i.pos1 << "\t" << i.pos2 << "\t" << i.strand << std::endl;
    } catch (...) {
    }
  bad_bed.close();
  */

  // close the files
  all_align.close();
  os_allbps.close();
  os_discordant.close();
  if (opt::write_corrected_reads) 
    os_corrected.close();
  log_file.close();

  // more clean up 
  if (ref_genome)
    delete ref_genome;
  if (ref_genome_viral)
    delete ref_genome_viral;
  
  // make the VCF file
  makeVCFs();
  
#ifndef __APPLE__
  //  std::cerr << SeqLib::displayRuntime(start) << std::endl;
#endif
}

void makeVCFs() {

  if (opt::bam.size() == 0) {
    std::cerr << "makeVCFs error: must supply a BAM via -t to get header from" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (!b_reader.IsOpen())     // open the tumor bam to get header info
    b_reader.Open(opt::bam.begin()->second);
    
  // make the VCF file
  WRITELOG("...loading the bps files for conversion to VCF", opt::verbose, true);

  std::string file = opt::analysis_id + ".bps.txt.gz";
  //if !SeqLib::read_access_test(file))
  //  file = opt::analysis_id + ".bps.txt";

  // make the header
  VCFHeader header;
  header.filedate = svabaUtils::fileDateString();
  header.source = args;
  header.reference = opt::refgenome;

  if (main_bwa)
    delete main_bwa;  
  if (!bwa_header.isEmpty())
    for (int i = 0; i < bwa_header.NumSequences(); ++i)
      header.addContigField(bwa_header.IDtoName(i),bwa_header.GetSequenceLength(i));

  for (auto& b : opt::bam) {
    std::string fname = b.second; //bpf.filename();
    header.addSampleField(fname);
    header.colnames += "\t" + fname; 
  }


  // check if it has a matched control. If so, output "somatic / germline" vcfs
  bool case_control_run = false;
  for (auto& b : opt::bam)
    if (b.first.at(0) == 'n')
      case_control_run = true;

  // primary VCFs
  if (SeqLib::read_access_test(file)) {
    WRITELOG("...making the primary VCFs (unfiltered and filtered) from file " + file, opt::verbose, true);
    VCFFile snowvcf(file, opt::analysis_id, b_header, header, !opt::no_unfiltered);

    if (!opt::no_unfiltered) {
      std::string basename = opt::analysis_id + ".svaba.unfiltered.";
      snowvcf.include_nonpass = true;
      WRITELOG("...writing unfiltered VCFs", opt::verbose, true);
      snowvcf.writeIndels(basename, opt::zip, !case_control_run);
      snowvcf.writeSVs(basename, opt::zip,    !case_control_run);
    }

    WRITELOG("...writing filtered VCFs", opt::verbose, true);
    std::string basename = opt::analysis_id + ".svaba.";
    snowvcf.include_nonpass = false;
    snowvcf.writeIndels(basename, opt::zip, !case_control_run);
    snowvcf.writeSVs(basename, opt::zip,    !case_control_run);

  } else {
    WRITELOG("ERROR: Failed to make VCF. Could not file bps file " + file, true, true);
  }

}

// parse the command line options
void parseRunOptions(int argc, char** argv) {
  bool die = false;

  if (argc <= 2) 
    die = true;

  bool help = false;
  std::stringstream ss;

  int sample_number = 0;

  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'E' : arg >> opt::ec_subsample; break;
    case 'p': arg >> opt::numThreads; break;
    case 'm': arg >> opt::sga::minOverlap; break;
    case 'a': arg >> opt::analysis_id; break;
    case 'B': arg >> opt::blacklist; break;
    case 'V': arg >> opt::germline_sv_file; break;
    case 'R': arg >> opt::simple_file; break;
    case 'L': arg >> opt::mate_lookup_min; break;
    case 'I': opt::interchrom_lookup = false; break;
    case 'Y': arg >> opt::microbegenome; break;
    case 'z': opt::zip = true; break;
    case 'h': help = true; break;
    case OPT_GAP_OPEN : arg >> opt::bwa::gap_open_penalty; break;
    case 'x' : arg >> opt::max_reads_per_assembly; break;
    case 'A' : opt::all_contigs = true; break;
    case OPT_MATCH_SCORE : arg >> opt::bwa::sequence_match_score; break;
    case OPT_READLEN : arg >> readlen; break;
    case OPT_MISMATCH : arg >> opt::bwa::mismatch_penalty; break;
    case OPT_GERMLINE : opt::germline = true; break;
    case OPT_GAP_EXTENSION : arg >> opt::bwa::gap_extension_penalty; break;
    case OPT_ZDROP : arg >> opt::bwa::zdrop; break;
    case OPT_BANDWIDTH : arg >> opt::bwa::bandwidth; break;
    case OPT_RESEED_TRIGGER : arg >> opt::bwa::reseed_trigger; break;
    case OPT_CLIP3 : arg >> opt::bwa::clip3_pen; break;
    case OPT_CLIP5 : arg >> opt::bwa::clip5_pen; break;
      case 'c': 
	tmp = "";
	arg >> tmp;
	if (tmp.find("chr") != std::string::npos) {
	  opt::chunk = 250000000; break;
	} else {
	  opt::chunk = stoi(tmp); break;
	}
    case OPT_ASQG: opt::sga::writeASQG = true; break;
    case OPT_LOD: arg >> opt::lod; break;
    case OPT_NO_UNFILTERED: opt::no_unfiltered = true; break;
    case OPT_LOD_DB: arg >> opt::lod_db; break;
    case OPT_LOD_SOMATIC: arg >> opt::lod_somatic; break;
    case OPT_LOD_SOMATIC_DB: arg >> opt::lod_somatic_db; break;
    case OPT_HP: opt::hp = true; break;
    case OPT_SCALE_ERRORS: arg >> opt::scale_error; break;
    case 'C': arg >> opt::max_cov;  break;
    case OPT_NUM_TO_SAMPLE: arg >> opt::num_to_sample;  break;
    case OPT_READ_TRACK: opt::read_tracking = true; break;
	case 't': 
	  tmp = svabaUtils::__bamOptParse(opt::bam, arg, sample_number++, "t");
	  if (opt::main_bam == "-")
	    opt::main_bam = tmp;
	  break;
	case 'n': 
	  tmp = svabaUtils::__bamOptParse(opt::bam, arg, sample_number++, "n");
	  break;
      case 'v': arg >> opt::verbose; break;
      case 'k': arg >> opt::regionFile; break;
      case 'e': arg >> opt::sga::error_rate; break;
      case 'G': arg >> opt::refgenome; break;
      case 'D': arg >> opt::dbsnp; break;
    case 's': arg >> opt::sd_disc_cutoff; break;
      case 'r': 
	arg >> opt::rules; 
	if (opt::rules == "all")
	  opt::rules = "";
       	break;
      case OPT_DISCORDANT_ONLY: opt::disc_cluster_only = true; break;
      case OPT_WRITE_EXTRACTED_READS: opt::write_extracted_reads = true; break;
    case OPT_NUM_ASSEMBLY_ROUNDS: arg >> opt::sga::num_assembly_rounds; break;
    case 'K': arg >> opt::ec_correct_type; break;
      default: die= true; 
    }
  }

  if (opt::chunk <= 0 || opt::main_bam == "-")
    opt::max_reads_per_assembly = INT_MAX;
  else if (opt::max_reads_per_assembly < 0) 
    opt::max_reads_per_assembly = 10000; //set a default

      

  if (!(opt::ec_correct_type == "s" || opt::ec_correct_type == "f" || opt::ec_correct_type == "0")) {
    WRITELOG("ERROR: Error correction type must be one of s, f, or 0", true, true);
    exit(EXIT_FAILURE);
  }

  // check that we input something
  if (opt::bam.size() == 0 && !die) {
    WRITELOG("Must add a bam file with -t flag. stdin with -t -", true, true);
    exit(EXIT_FAILURE);
  }

  if (opt::numThreads <= 0) {
    WRITELOG("Invalid number of threads from -p flag: " + SeqLib::AddCommas(opt::numThreads), true, true);
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

bool runWorkUnit(const SeqLib::GenomicRegion& region, svabaWorkUnit& wu, long unsigned int thread_id) {
  
  WRITELOG("Running region " + region.ToString() + " on thread " + std::to_string(thread_id), opt::verbose > 1, true);

  for (auto& w : wu.walkers)
    set_walker_params(w.second);

  // create a new BFC read error corrector for this
  SeqPointer<SeqLib::BFC> bfc;
  if (opt::ec_correct_type == "f") {
    bfc = SeqPointer<SeqLib::BFC>(new SeqLib::BFC());
    for (auto& w : wu.walkers)
      w.second.bfc = bfc;
  }
  
  // setup structures to store the final data for this region
  std::vector<AlignedContig> alc;
  SeqLib::BamRecordVector all_contigs, all_microbial_contigs;

  // start a timer
  svabaUtils::svabaTimer st;
  st.start();

  // setup for the BAM walkers
  CountPair read_counts = {0,0};

  // read in alignments from the main region
  for (auto& w : wu.walkers) {

    // set the region to jump to
    if (!region.IsEmpty()) {
      w.second.SetRegion(region);
    } else { // whole BAM analysis. If region file set, then set regions
      if (file_regions.size()) {
	w.second.SetMultipleRegions(file_regions);
      } else { // no regions, literally take every read in bam
	// default is walk all regions, so leave empty
      }
    }

    // do the reading, and store the bad mate regions
    wu.badd.Concat(w.second.readBam(&log_file)); 
    wu.badd.MergeOverlappingIntervals();
    wu.badd.CreateTreeMap();

    
    // adjust the counts
    if (w.first.at(0) == 't') {
      read_counts.first += w.second.reads.size();
    } else {
      read_counts.second += w.second.reads.size();
    }

  }

  // collect all of the cigar strings in a hash
  std::unordered_map<std::string, SeqLib::CigarMap> cigmap;
  for (const auto& w : wu.walkers) 
    cigmap[w.first] = w.second.cigmap;

  // setup read collectors
  std::vector<char*> all_seqs;
  //SeqLib::BamRecordVector bav_this;
  svabaReadVector bav_this;

  // collect and clear reads from main round
  std::unordered_set<std::string> dedupe;
  collect_and_clear_reads(wu.walkers, bav_this, all_seqs, dedupe);

  // adjust counts and timer
  st.stop("r");

  // get the mate reads, if this is local assembly and has insert-size distro
  if (!region.IsEmpty() && !opt::single_end && min_dscrd_size_for_variant) {
    run_mate_collection_loop(region, wu.walkers, wu.badd);
    // collect the reads together from the mate walkers
    collect_and_clear_reads(wu.walkers, bav_this, all_seqs, dedupe);
    st.stop("m");
  }
  

  // do the discordant read clustering
  DiscordantClusterMap dmap, dmap_tmp;
  
  // if couldn't get insert size stats, skip discordant clustering
  if (!min_dscrd_size_for_variant || opt::single_end)
    goto afterdiscclustering;

  WRITELOG("...discordant read clustering", opt::verbose > 1, false);
  dmap = DiscordantCluster::clusterReads(bav_this, region, max_mapq_possible, &min_isize_for_disc);

  // tag FR clusters that are below min_dscrd_size_for_variant AND low support
  for (auto& d : dmap) {
    bool below_size = 	d.second.m_reg1.strand == '+' && d.second.m_reg2.strand == '-' && 
      (d.second.m_reg2.pos1 - d.second.m_reg1.pos2) < min_dscrd_size_for_variant && 
      d.second.m_reg1.chr == d.second.m_reg2.chr;

    // low support and low size, completely ditch it
    if (below_size && (d.second.tcount + d.second.ncount) < 4)
      continue;
    else
      dmap_tmp.insert(std::pair<std::string, DiscordantCluster>(d.first, d.second));
  }
  dmap = dmap_tmp;

  // print out results
  if (opt::verbose > 3)
    for (auto& i : dmap) 
      WRITELOG(i.first + " " + i.second.toFileString(false), true, false);

 afterdiscclustering:

  // skip all the assembly stuff?
  if (opt::disc_cluster_only)
    goto afterassembly;
  
  /////////////////////
  //// ASSEMBLY 
  /////////////////////

  // remove the hardclips, don't assemble them
  remove_hardclips(bav_this);

  if (opt::disc_cluster_only)
    goto afterassembly;

  // check that we don't have too many reads
  if (bav_this.size() > (size_t)(region.Width() * 20) && region.Width() > 20000) {
    std::stringstream ssss;
    WRITELOG("TOO MANY READS IN REGION " + SeqLib::AddCommas(bav_this.size()) + "\t" + region.ToString(), opt::verbose, false);
    goto afterassembly;
  }

  // print message about assemblies
  if (bav_this.size() > 1) {
    WRITELOG("Doing assemblies on " + region.ToString(), opt::verbose > 1, false);
  } else if (bav_this.size() < 3) { 
    WRITELOG("Skipping assembly (<= 2 reads) on " + region.ToString(), opt::verbose > 1, false);
    goto afterassembly;
  }

  // do the kmer correction, in place
  if (opt::ec_correct_type == "s") {
    correct_reads(all_seqs, bav_this);
  } else if (opt::ec_correct_type == "f" && bav_this.size() >= 8) {
    
    assert(bfc);
    int learn_reads_count = bfc->NumSequences();
    
    bfc->Train();
    bfc->clear();  // clear memory and reads. Keeps training data

    st.stop("t");

    // reload with the reads to be corrected
    for (auto& s : bav_this)
      bfc->AddSequence(s.Seq().c_str(), "", "");

    // error correct 
    bfc->ErrorCorrect();

    // retrieve the sequences
    std::string s, name_dum;
    for (auto& r : bav_this) {
      assert(bfc->GetSequence(s, name_dum));
      r.SetSeq(s);
    }
    
    double kcov = bfc->GetKCov();
    int kmer    = bfc->GetKMer();

    WRITELOG("...BFC attempted correct " + std::to_string(bav_this.size()) + " train: " + 
	     std::to_string(learn_reads_count) + " kcov: " + std::to_string(kcov) + 
	     " kmer: " + std::to_string(kmer), opt::verbose > 1, true);      

  }

  st.stop("k");
  
  // do the assembly, contig realignment, contig local realignment, and read realignment
  // modifes bav_this, alc, all_contigs and all_microbial_contigs
  run_assembly(region, bav_this, alc, all_contigs, all_microbial_contigs, dmap, cigmap, wu.ref_genome);

afterassembly:

  // clear it out, not needed anymore
  if (bfc)
    bfc->clear();
  
  st.stop("as");
  WRITELOG("...done assembling, post processing", opt::verbose > 1, false);

  // get the breakpoints
  std::vector<BreakPoint> bp_glob;
  
  for (auto& i : alc) {
    std::vector<BreakPoint> allbreaks = i.getAllBreakPoints();
    bp_glob.insert(bp_glob.end(), allbreaks.begin(), allbreaks.end());
  }

  if (dbsnp_filter && opt::dbsnp.length()) {
    WRITELOG("...DBSNP filtering", opt::verbose > 1, false);
    for (auto & i : bp_glob) 
      dbsnp_filter->queryBreakpoint(i);
  }

  // filter against blacklist
  for (auto& i : bp_glob) 
    i.checkBlacklist(blacklist);

  // add in the discordant clusters as breakpoints
  for (auto& i : dmap) {
    // dont send DSCRD if FR and below size
    bool below_size = 	i.second.m_reg1.strand == '+' && i.second.m_reg2.strand == '-' && 
      (i.second.m_reg2.pos1 - i.second.m_reg1.pos2) < min_dscrd_size_for_variant && 
      i.second.m_reg1.chr == i.second.m_reg2.chr;
    // DiscordantCluster not associated with assembly BP and has 2+ read support
    if (!i.second.hasAssociatedAssemblyContig() && 
	(i.second.tcount + i.second.ncount) >= MIN_DSCRD_READS_DSCRD_ONLY && i.second.valid() && !below_size) {
      BreakPoint tmpbp(i.second, main_bwa, dmap, region);
      bp_glob.push_back(tmpbp);
    }
  }
  
  // de duplicate the breakpoints
  std::sort(bp_glob.begin(), bp_glob.end());
  bp_glob.erase( std::unique( bp_glob.begin(), bp_glob.end() ), bp_glob.end() );

  // add the coverage data to breaks for allelic fraction computation
  std::unordered_map<std::string, STCoverage*> covs;
  for (auto& i : opt::bam) 
    covs[i.first] = &wu.walkers[i.first].cov;

  for (auto& i : bp_glob)
    i.addCovs(covs);

  for (auto& i : bp_glob) {
    i.readlen = readlen; // set the readlength
    i.scoreBreakpoint(opt::lod, opt::lod_db, opt::lod_somatic, opt::lod_somatic_db, opt::scale_error, min_dscrd_size_for_variant);
  }

  // label somatic breakpoints that intersect directly with normal as NOT somatic
  std::unordered_set<std::string> norm_hash;
  for (auto& i : bp_glob) // hash the normals
    if (!i.somatic_score && i.confidence == "PASS" && i.evidence == "INDEL") {
      norm_hash.insert(i.b1.hash());
      norm_hash.insert(i.b2.hash());
      norm_hash.insert(i.b1.hash(1));
      norm_hash.insert(i.b1.hash(-1));
    }

  // find somatic that intersect with norm. Set somatic = 0;
  for (auto& i : bp_glob)  
    if (i.somatic_score && i.evidence == "INDEL" && (norm_hash.count(i.b1.hash()) || norm_hash.count(i.b2.hash()))) {
      i.somatic_score = -3;
    }

  // remove indels at repeats that have multiple variants
  std::unordered_map<std::string, size_t> ccc;
  for (auto& i : bp_glob) {
    if (i.evidence == "INDEL" && i.repeat_seq.length() > 6) {
      ++ccc[i.b1.hash()];
    }
  }
  for (auto& i : bp_glob) {
    if (ccc[i.b1.hash()] > 1)
      i.confidence = "REPVAR";
  }

  // remove somatic calls if they have a germline normal SV in them or indels with 
  // 2+germline normal in same contig
  std::unordered_set<std::string> bp_hash;
  for (auto& i : bp_glob) { // hash the normals
    if (!i.somatic_score && i.evidence != "INDEL" && i.confidence == "PASS") {
      bp_hash.insert(i.cname);
    }
  }
  for (auto& i : bp_glob)  // find somatic that intersect with norm. Set somatic = 0;
    if (i.somatic_score && i.num_align > 1 && bp_hash.count(i.cname)) {
      i.somatic_score = -2;
    }

  // remove somatic SVs that overlap with germline svs
  if (germline_svs.size()) {
    for (auto& i : bp_glob) {
      if (i.somatic_score && i.b1.gr.chr == i.b2.gr.chr && i.evidence != "INDEL") {
	SeqLib::GenomicRegion gr1 = i.b1.gr;
	SeqLib::GenomicRegion gr2 = i.b2.gr;
	gr1.Pad(GERMLINE_CNV_PAD);
	gr2.Pad(GERMLINE_CNV_PAD);
	if (germline_svs.OverlapSameInterval(gr1, gr2)) {
	  i.somatic_score = -1;
	}
      }
    }
      
  }

  // add the ref and alt tags
  // WHY IS THIS NOT THREAD SAFE?
  for (auto& i : bp_glob)
    i.setRefAlt(wu.ref_genome, wu.vir_genome);

  // transfer local versions to thread store
  for (const auto& a : alc)
    if (a.hasVariant()) {
      wu.m_alc.push_back(a);
      wu.m_bamreads_count += a.NumBamReads();
    }
  for (const auto& d : dmap)
    wu.m_disc_reads += d.second.reads.size();
  
  wu.m_contigs.insert(wu.m_contigs.end(), all_contigs.begin(), all_contigs.end());
  wu.m_vir_contigs.insert(wu.m_vir_contigs.end(), all_microbial_contigs.begin(), all_microbial_contigs.end());
  wu.m_disc.insert(dmap.begin(), dmap.end());
  for (const auto& a : alc)
    wu.m_bamreads_count += a.NumBamReads();
  for (auto& i : bp_glob) 
    if ( i.hasMinimal() && (i.confidence != "NOLOCAL" || i.complex_local ) ) 
      wu.m_bps.push_back(i);
  
  // dump if getting to much memory
  if (wu.MemoryLimit(THREAD_READ_LIMIT, THREAD_CONTIG_LIMIT) && !opt::hp) {
    WRITELOG("writing contigs etc on thread " + std::to_string(thread_id) + " with limit hit of " + std::to_string(wu.m_bamreads_count), opt::verbose > 1, true);
    pthread_mutex_lock(&snow_lock);    
    WriteFilesOut(wu); 
    pthread_mutex_unlock(&snow_lock);
  }
  
  // write extracted reads
  if (opt::write_extracted_reads) {
    pthread_mutex_lock(&snow_lock);    
    for (auto& r : bav_this)
      er_writer.WriteRecord(r);
    pthread_mutex_unlock(&snow_lock);
  }
  
  // write the raw error corrected reads to a fasta
  if (opt::write_corrected_reads) {
    pthread_mutex_lock(&snow_lock);    
    for (auto& r : bav_this) {
      std::string seq;
      r.GetZTag("KC", seq);
      if (seq.empty())
	seq = r.QualitySequence();
      //os_corrected << ">" << SRTAG(r) << std::endl << seq << std::endl;
      os_corrected << ">" << r.SR() << std::endl << seq << std::endl;
    }
    pthread_mutex_unlock(&snow_lock);
  }

  st.stop("pp");
  
  // display the run time
  WRITELOG(svabaUtils::runTimeString(read_counts.first, read_counts.second, alc.size(), region, b_header, st, start), opt::verbose > 1, true);

  // clear out the reads and reset the walkers
  for (auto& w : wu.walkers) {
    w.second.clear(); 
    w.second.m_limit = opt::max_reads_per_assembly;
  }

  return true;
}

void sendThreads(SeqLib::GRC& regions_torun) {

  // Create the queue and consumer (worker) threads
  wqueue<svabaWorkItem*>  queue;
  std::vector<ConsumerThread<svabaWorkItem>*> threadqueue;
  for (int i = 0; i < opt::numThreads; i++) {
    ConsumerThread<svabaWorkItem>* threadr = new ConsumerThread<svabaWorkItem>(queue, opt::verbose > 0,
										   opt::refgenome, opt::microbegenome,
										   opt::bam);
    threadr->start();
    threadqueue.push_back(threadr);
  }

  // send the jobs
  size_t count = 0;
  for (auto& i : regions_torun) {
    svabaWorkItem * item     = new svabaWorkItem(SeqLib::GenomicRegion(i.chr, i.pos1, i.pos2), ++count);
    queue.add(item);
  }
  if (!regions_torun.size()) { // whole genome 
    svabaWorkItem * item     = new svabaWorkItem(SeqLib::GenomicRegion(), ++count);
    queue.add(item);
  }
  
  // wait for the threads to finish
  for (int i = 0; i < opt::numThreads; ++i) 
    threadqueue[i]->join();

  // write and free remaining items stored in the thread
  pthread_mutex_lock(&snow_lock);
  for (int i = 0; i < opt::numThreads; ++i) 
    WriteFilesOut(threadqueue[i]->wu); 
  pthread_mutex_unlock(&snow_lock);

}

SeqLib::GRC makeAssemblyRegions(const SeqLib::GenomicRegion& region) {

  // set the regions to run
  SeqLib::GRC grv_small;
  if (region.IsEmpty()) {  // whole genome, so divide up the whole thing
    for (size_t c = 0; c < 23; ++c)
      grv_small.Concat(SeqLib::GRC(opt::chunk, WINDOW_PAD, SeqLib::GenomicRegion(c, WINDOW_PAD + 1, b_header.GetSequenceLength(c) - WINDOW_PAD - 1)));
    //grv_small.Concat(SeqLib::GRC(opt::chunk, WINDOW_PAD, SeqLib::GenomicRegion(c, WINDOW_PAD + 1, SeqLib::CHR_LEN[c] - WINDOW_PAD - 1)));
  }
  else if (region.Width() >= opt::chunk) // divide into smaller chunks
    grv_small = SeqLib::GRC(opt::chunk, WINDOW_PAD, region);
  else
    grv_small.add(region);

  return grv_small;

}

void alignReadsToContigs(SeqLib::BWAWrapper& bw, const SeqLib::UnalignedSequenceVector& usv, 
			 svabaReadVector& bav_this, std::vector<AlignedContig>& this_alc, const SeqLib::RefGenome *  rg) {
  
  if (!usv.size())
    return;

  // get the reference info
  SeqLib::GRC g;
  for (auto& a : this_alc)
    for (auto& i : a.getAsGenomicRegionVector()) {
      i.Pad(100);
      g.add(i);
    }
  g.MergeOverlappingIntervals();

  // get the reference sequence
  std::vector<std::string> ref_alleles;
  for (auto& i : g)
    if (i.chr < 24) //1-Y
      try {
	std::string tmpref = rg->QueryRegion(i.ChrName(bwa_header), i.pos1, i.pos2);
	ref_alleles.push_back(tmpref); 
      } catch (...) {
	std::cerr << "Caught exception for ref_allele on " << i << std::endl;
      }
  // make the reference allele BWAWrapper
  SeqLib::BWAWrapper bw_ref;
  SeqLib::UnalignedSequenceVector usv_ref;
  int aa = 0;
  for (auto& i : ref_alleles) {
    if (!i.empty())
      usv_ref.push_back({std::to_string(aa++), i, std::string()}); // name, seq, qual
  }
  if (!usv_ref.size())
    bw_ref.ConstructIndex(usv_ref);
  
  // set up custom alignment parameters, mean
  bw_ref.SetGapOpen(16); // default 6
  bw.SetGapOpen(16); // default 6
  bw_ref.SetMismatchPenalty(9); // default 2
  bw.SetMismatchPenalty(9); // default 4

  for (auto i : bav_this) {
    
    SeqLib::BamRecordVector brv, brv_ref;

    // try the corrected seq first
    //std::string seqr = i.GetZTag("KC");
    //  if (seqr.empty())
    //	seqr = i.QualitySequence();
    std::string seqr = i.Seq();
    
    bool hardclip = false;
    assert(seqr.length());
    bw.AlignSequence(seqr, i.Qname(), brv, hardclip, 0.60, 10000);

    if (brv.size() == 0) 
      continue;

    // get the maximum non-reference alignment score
    int max_as = 0;
    for (auto& r : brv) {
      int thisas = 0;
      r.GetIntTag("AS", thisas);
      max_as = std::max(max_as, thisas);
    }

    // align to the reference alleles
    if (!bw_ref.IsEmpty())
      bw_ref.AlignSequence(seqr, i.Qname(), brv_ref, hardclip, 0.60, 10);

    // get the maximum reference alignment score
    int max_as_r = 0;
    for (auto& r : brv_ref) {
      int thisas= 0;
      r.GetIntTag("AS", thisas);
      max_as_r = std::max(max_as_r, thisas);
    }
    
    // reject if better alignment to reference
    if (max_as_r > max_as) {
      //std::cerr << " Alignment Rejected for " << max_as_r << ">" << max_as << "  " << i << std::endl;
      //std::cerr << "                        " << max_as_r << ">" << max_as << "  " << brv_ref[0] << std::endl;
      continue;
    }

    // convert to a svabaReadVector
    svabaReadVector brv_svaba;
    for (auto& r : brv)
      brv_svaba.push_back(svabaRead(r, i.Prefix()));
    brv.clear();

    // make sure we have only one alignment per contig
    std::set<std::string> cc;

    // check which ones pass
    SeqLib::BamRecordVector bpass;
    for (auto& r : brv_svaba) {

      // make sure alignment score is OK
      int thisas = 0;
      r.GetIntTag("AS", thisas);
      if ((double)r.NumMatchBases() * 0.5 > thisas/* && i.GetZTag("SR").at(0) == 't'*/)
      	continue;
      
      bool length_pass = (r.PositionEnd() - r.Position()) >= ((double)seqr.length() * 0.75);

      if (length_pass && !cc.count(usv[r.ChrID()].Name)) {
	bpass.push_back(r);
	cc.insert(usv[r.ChrID()].Name);
      }
    }

    // annotate the original read
    for (auto& r : bpass) {

      r2c this_r2c; // alignment of this read to this contig
      if (r.ReverseFlag())
	this_r2c.rc = true;
      //i.SmartAddTag("RC","1");
      //else 
      //i.SmartAddTag("RC","0");

      this_r2c.AddAlignment(r);
      i.AddR2C(usv[r.ChrID()].Name, this_r2c);
      
      //i.SmartAddTag("SL", std::to_string(r.Position()));
      //i.SmartAddTag("SE", std::to_string(r.PositionEnd()));
      //i.SmartAddTag("TS", std::to_string(r.AlignmentPosition()));
      //i.SmartAddTag("TE", std::to_string(r.AlignmentEndPosition()));
      //i.SmartAddTag("SC", r.CigarString());
      //i.SmartAddTag("CN", usv[r.ChrID()].Name);

      //i.AddZTag("SR", i.SR().substr(0,4));
      //i.AddZTag("GV", i.Seq());

      // add the read to the right contig (loop to check for right contig)
      for (auto& a : this_alc) {
	if (a.getContigName() != usv[r.ChrID()].Name)
	  continue;
	a.AddAlignedRead(i);
      }
      
    } // end passing bwa-aligned read loop 
  } // end main read loop
}

void set_walker_params(svabaBamWalker& walk) {

  walk.main_bwa = main_bwa; // set the pointer
  walk.blacklist = blacklist;
  walk.do_kmer_filtering = (opt::ec_correct_type == "s" || opt::ec_correct_type == "f");
  walk.simple_seq = &simple_seq;
  walk.kmer_subsample = opt::ec_subsample;
  walk.max_cov = opt::max_cov;
  walk.m_mr = mr;  // set the read filter pointer
  walk.m_limit = opt::max_reads_per_assembly;

}

CountPair run_mate_collection_loop(const SeqLib::GenomicRegion& region, WalkerMap& wmap, SeqLib::GRC& badd) {

  SeqLib::GRC this_bad_mate_regions; // store the newly found bad mate regions
  
  CountPair counts = {0,0};

  MateRegionVector all_somatic_mate_regions;
  all_somatic_mate_regions.add(MateRegion(region.chr, region.pos1, region.pos2)); // add the origional, don't want to double back
  
  for (int jjj = 0; jjj <  MAX_MATE_ROUNDS; ++jjj) {
    
    MateRegionVector normal_mate_regions;
    if (wmap.size() > 1) // have at least one control bam
      normal_mate_regions = __collect_normal_mate_regions(wmap);
    
    // get the mates from somatic 3+ mate regions
    // that don't overlap with normal mate region
    MateRegionVector tmp_somatic_mate_regions = __collect_somatic_mate_regions(wmap, normal_mate_regions);
    
    // no more regions to check
    if (!tmp_somatic_mate_regions.size())
      break;
    
    // keep regions that haven't been visited before
    MateRegionVector somatic_mate_regions;
    for (auto& s : tmp_somatic_mate_regions) {
      
      // check if its not bad from mate region
      if (badd.size())
	if (badd.CountOverlaps(s))
	  continue;

      // check if we are allowed to lookup interchromosomal
      if (!opt::interchrom_lookup && (s.chr != region.chr || std::abs(s.pos1 - region.pos1) < LARGE_INTRA_LOOKUP_LIMIT))
	continue;

      // new region overlaps with one already seen from another round
      for (auto& ss : all_somatic_mate_regions) 
	if (s.GetOverlap(ss)) 
	  continue;
      
      if (s.count > opt::mate_lookup_min * 2 || (jjj == 0)) { // be more strict about higher rounds and inter-chr
	somatic_mate_regions.add(s);
	all_somatic_mate_regions.add(s);
      }
      
      // don't add too many regions
      if (all_somatic_mate_regions.size() > MAX_NUM_MATE_WINDOWS) {
	all_somatic_mate_regions.clear(); // its a bad region. Don't even look up any
	break;
      }

    }
    
    // print out to log
    for (auto& i : somatic_mate_regions) 
      WRITELOG("...mate region " + i.ToString() + " case read count that triggered lookup: " + 
	       std::to_string(i.count) + " on mate-lookup round " + std::to_string(jjj+1), opt::verbose > 1, true);
    
    // collect the reads for this round
    std::pair<int,int> mate_read_counts = collect_mate_reads(wmap, somatic_mate_regions, jjj, this_bad_mate_regions);

    // update the counts
    counts.first += mate_read_counts.first;
    counts.second += mate_read_counts.second;

    WRITELOG("\t<case found, control found>: <" + SeqLib::AddCommas(counts.first) +
	     "," + SeqLib::AddCommas(counts.second) + "> on round " + std::to_string(jjj+1) , opt::verbose > 2, true);
    
    if (counts.first + counts.second == 0) // didn't get anything on first round, so quit
      break; 

  } // mate collection round loop

  // update this threads tally of bad mate regions
  badd.Concat(this_bad_mate_regions);
  badd.MergeOverlappingIntervals();
  badd.CreateTreeMap();
  WRITELOG("\tTotal of " + SeqLib::AddCommas(badd.size()) + " bad mate regions for this thread", opt::verbose > 1, true);

  return counts;
}


MateRegionVector __collect_normal_mate_regions(WalkerMap& walkers) {
  
  MateRegionVector normal_mate_regions;
  for (auto& b : opt::bam)
    if (b.first.at(0) == 'n')
      normal_mate_regions.Concat(walkers[b.first].mate_regions);

  normal_mate_regions.MergeOverlappingIntervals();
  normal_mate_regions.CreateTreeMap();
  
  return normal_mate_regions;

}

MateRegionVector __collect_somatic_mate_regions(WalkerMap& walkers, MateRegionVector& bl) {


  MateRegionVector somatic_mate_regions;
  for (auto& b : opt::bam)
    if (b.first.at(0) == 't')
      for (auto& i : walkers[b.first].mate_regions) {
	if (i.count >= opt::mate_lookup_min && !bl.CountOverlaps(i)
	    && (!blacklist.size() || !blacklist.CountOverlaps(i)))
	  somatic_mate_regions.add(i); 
      }
  
  // reduce it down
  somatic_mate_regions.MergeOverlappingIntervals();

  return somatic_mate_regions;

}

void correct_reads(std::vector<char*>& learn_seqs, svabaReadVector& brv) {

  if (!learn_seqs.size())
    return;

  if (opt::ec_correct_type == "s") {
    
    KmerFilter kmer;
    int kmer_corrected = 0;
 
    // make the index for learning correction
    kmer.makeIndex(learn_seqs);

    // free the training sequences
    // this memory was alloced w/strdup in collect_and_clear_reads
    for (size_t i = 0; i < learn_seqs.size(); ++i)
      if (learn_seqs[i]) 
    	free(learn_seqs[i]);
    
    // do the correction
    kmer_corrected = kmer.correctReads(brv); 
    //kmer_corrected = kmer.correctReads(brv); 

    WRITELOG("...SGA kmer corrected " + std::to_string(kmer_corrected) + " reads of " + std::to_string(brv.size()), opt::verbose > 1, true);  
  } 
}

void remove_hardclips(svabaReadVector& brv) {
  svabaReadVector bav_tmp;
  for (auto& i : brv)
    if (i.NumHardClip() == 0) 
      bav_tmp.push_back(i);
  brv = bav_tmp;
}

void run_assembly(const SeqLib::GenomicRegion& region, svabaReadVector& bav_this, std::vector<AlignedContig>& master_alc, 
		  SeqLib::BamRecordVector& master_contigs, SeqLib::BamRecordVector& master_microbial_contigs, DiscordantClusterMap& dmap,
		  std::unordered_map<std::string, SeqLib::CigarMap>& cigmap, SeqLib::RefGenome* refg) {

  // get the local region
  std::string lregion;
  if (!region.IsEmpty()) {
    try {
      lregion = refg->QueryRegion(bwa_header.IDtoName(region.chr), region.pos1, region.pos2);
    } catch (...) {
      WRITELOG(" Caught exception for lregion with reg " + region.ToString(), true, true);
      lregion = "";
    }
  }

  // make a BWA wrapper from the locally retrieved sequence
  SeqLib::UnalignedSequenceVector local_usv = {{"local", lregion, std::string()}};
  SeqLib::BWAWrapper local_bwa;
  if (local_usv[0].Seq.length() > 200) // have to have pulled some ref sequence
    local_bwa.ConstructIndex(local_usv);

  std::stringstream region_string;
  region_string << region;
  WRITELOG("...running assemblies for region " + region_string.str(), opt::verbose > 1, false);
  
  // set the contig uid
  std::string name = "c_" + std::to_string(region.chr+1) + "_" + std::to_string(region.pos1) + "_" + std::to_string(region.pos2);
  
  // where to store contigs
  SeqLib::UnalignedSequenceVector all_contigs_this;
  
  // setup the engine
  svabaAssemblerEngine engine(name, opt::sga::error_rate, opt::sga::minOverlap, readlen);
  if (opt::sga::writeASQG)
    engine.setToWriteASQG();
  engine.fillReadTable(bav_this);
  
  // do the actual assembly
  engine.performAssembly(opt::sga::num_assembly_rounds);
  
  // retrieve contigs
  all_contigs_this = engine.getContigs();
  WRITELOG("...assembled " + std::to_string(all_contigs_this.size()) + " contigs for " + name, opt::verbose > 1, true);

  // store the aligned contig struct
  std::vector<AlignedContig> this_alc;
      
  // align the contigs to the genome
  WRITELOG("...aliging contigs to genome", opt::verbose > 1, false);

  SeqLib::UnalignedSequenceVector usv;
  
  for (auto& i : all_contigs_this) {
    
    // if too short, skip
    if ((int)i.Seq.length() < (readlen * 1.15) && !opt::all_contigs)
      continue;
    
    bool hardclip = false;	

    //// LOCAL REALIGNMENT
    // align to the local region
    SeqLib::BamRecordVector local_ct_alignments;
    if (!local_bwa.IsEmpty())
      local_bwa.AlignSequence(i.Seq, i.Name, local_ct_alignments, hardclip, SECONDARY_FRAC, SECONDARY_CAP);
    
    // check if it has a non-local alignment
    bool valid_sv = true;
    for (auto& aa : local_ct_alignments) {
      if (aa.NumClip() < MIN_CLIP_FOR_LOCAL) // || aa.GetIntTag("NM") < MAX_NM_FOR_LOCAL)
	valid_sv = false; // has a non-clipped local alignment. can't be SV. Indel only
    }
    ////////////
    
    // do the main realignment
    SeqLib::BamRecordVector ct_alignments;
    main_bwa->AlignSequence(i.Seq, i.Name, ct_alignments, hardclip, SECONDARY_FRAC, SECONDARY_CAP);	

    if (opt::verbose > 3)
      for (auto& i : ct_alignments)
	std::cerr << " aligned contig: " << i << std::endl;
    
    // do the microbe realigenment
    SeqLib::BamRecordVector ct_plus_microbe;

    if (microbe_bwa && !svabaUtils::hasRepeat(i.Seq)) {
      
      // do the microbial alignment
      SeqLib::BamRecordVector microbial_alignments;
      bool hardclip = false;
      microbe_bwa->AlignSequence(i.Seq, i.Name, microbial_alignments, hardclip, SECONDARY_FRAC, SECONDARY_CAP);
      
      // if the microbe alignment is large enough and doesn't overlap human...
      for (auto& j : microbial_alignments) {
	// keep only long microbe alignments with decent mapq
	if (j.NumMatchBases() >= MICROBE_MATCH_MIN && j.MapQuality() >= 10) { 
	  if (svabaUtils::overlapSize(j, ct_alignments) <= 20) { // keep only those where most do not overlap human
	    assert(microbe_bwa->ChrIDToName(j.ChrID()).length());
	    j.AddZTag("MC", microbe_bwa->ChrIDToName(j.ChrID()));
	    master_microbial_contigs.push_back(j);
	    ct_plus_microbe.push_back(j);
	  }
	}
      }
    }
    
    // add in the chrosome name tag for human alignments
    if (main_bwa)
      for (auto& r : ct_alignments) {
	assert(main_bwa->ChrIDToName(r.ChrID()).length());
	r.AddZTag("MC", main_bwa->ChrIDToName(r.ChrID()));
	if (!valid_sv)
	  r.AddIntTag("LA", 1); // flag as having a valid local alignment. Can't be SV
      }
    
    // remove human alignments that are not as good as microbe
    // that is, remove human alignments that intersect with microbe. 
    // We can do this because we already removed microbial alignments taht 
    // intersect too much with human. Thus, we are effectively removing human 
    // alignments that are contained within a microbial alignment
    
    SeqLib::BamRecordVector human_alignments;
    for (auto& j : ct_alignments) {
      // keep human alignments that have < 50% overlap with microbe and have >= 25 bases matched
      if (svabaUtils::overlapSize(j, ct_plus_microbe) < 0.5 * j.Length() && j.NumMatchBases() >= MIN_CONTIG_MATCH) { 
	human_alignments.push_back(j);
        master_contigs.push_back(j);
      }
    }
    
    if (!human_alignments.size())
      continue;
    
    // add in the microbe alignments
    human_alignments.insert(human_alignments.end(), ct_plus_microbe.begin(), ct_plus_microbe.end());
    
    // make the aligned contig object
    if (!human_alignments.size())
      continue;

    // check simple sequence overlaps
    if (simple_seq.size())
      for (auto& k : human_alignments) {
	SeqLib::GRC ovl = simple_seq.FindOverlaps(k.AsGenomicRegion(), true);
	int msize = 0;
	for (auto& j: ovl) {
	  int nsize = j.Width() - k.MaxDeletionBases() - 1;
	  if (nsize > msize && nsize > 0)
	    msize = nsize;
	}
	k.AddIntTag("SZ", msize);
      }

    // make the aligned contigs
    AlignedContig ac(human_alignments, prefixes);
    
    // assign the local variable to each
    ac.checkLocal(region);
    
    this_alc.push_back(ac);
    usv.push_back({i.Name, i.Seq, std::string()});	  
  } // end loop through contigs

  assert(this_alc.size() == usv.size());

  // didnt get any contigs that made it all the way through
  if (!this_alc.size())
    return;
  
  // Align the reads to the contigs with BWA-MEM
  SeqLib::BWAWrapper bw;
  bw.ConstructIndex(usv);
  
  if (opt::verbose > 3)
    std::cerr << "...aligning " << bav_this.size() << " reads to " << this_alc.size() << " contigs " << std::endl;
  alignReadsToContigs(bw, usv, bav_this, this_alc, refg);
  
  // Get contig coverage, discordant matching to contigs, etc
  for (auto& a : this_alc) {
    
    // repeat sequence filter
    a.assessRepeats();
    
    a.splitCoverage();	
    // now that we have all the break support, check that the complex breaks are OK
    a.refilterComplex(); 
    // add discordant reads support to each of the breakpoints
    a.addDiscordantCluster(dmap);
    // add in the cigar matches
    a.checkAgainstCigarMatches(cigmap);
    // add to the final structure
    master_alc.push_back(a);
    
  }


}

CountPair collect_mate_reads(WalkerMap& walkers, const MateRegionVector& mrv, int round, SeqLib::GRC& this_bad_mate_regions) {

  CountPair counts = {0,0};  

  if (!mrv.size())
    return counts;
  
  for (auto& w : walkers) {

    int oreads = w.second.reads.size();
    w.second.m_limit = MATE_REGION_LOOKUP_LIMIT;

    // convert MateRegionVector to GRC
    SeqLib::GRC gg;
    for (auto& s : mrv) 
      gg.add(SeqLib::GenomicRegion(s.chr, s.pos1, s.pos2, s.strand));

    assert(w.second.SetMultipleRegions(gg));
    w.second.get_coverage = false;
    w.second.get_mate_regions = (round != MAX_MATE_ROUNDS);

    // clear out the current store of mate regions, since we 
    // already added these to the to-do pile
    w.second.mate_regions.clear();

    this_bad_mate_regions.Concat(w.second.readBam(&log_file)); 
    
    // update the counts
    if (w.first.at(0) == 't') 
      counts.first += (w.second.reads.size() - oreads);
    else
      counts.second += (w.second.reads.size() - oreads);

  }
  
  return counts;
}

//void collect_and_clear_reads(WalkerMap& walkers, SeqLib::BamRecordVector& brv, std::vector<char*>& learn_seqs, std::unordered_set<std::string>& dedupe) {
void collect_and_clear_reads(WalkerMap& walkers, svabaReadVector& brv, std::vector<char*>& learn_seqs, std::unordered_set<std::string>& dedupe) {

  // concatenate together all the reads from the different walkers
  for (auto& w : walkers) {
    for (auto& r : w.second.reads) {
      std::string sr = r.SR();
      if (!dedupe.count(sr)) {
	brv.push_back(r); 
        dedupe.insert(sr);
      }
    }

    // concat together all of the learning sequences
    if (opt::ec_correct_type != "s")
      assert(!w.second.all_seqs.size());
    for (auto& r : w.second.all_seqs) {
      learn_seqs.push_back(strdup(r));
      free(r); // free what was alloced in 
    }

    w.second.all_seqs.clear();
    w.second.reads.clear();
  }
}

void WriteFilesOut(svabaWorkUnit& wu) {

  // print the alignment plots
  for (const auto& i : wu.m_alc) 
    if (i.hasVariant()) 
      all_align << i << std::endl;

  // send the microbe to file
  for (const auto& b : wu.m_vir_contigs)
    b_microbe_writer.WriteRecord(b);
  
  // send the discordant to file
  for (auto& i : wu.m_disc)
    if (i.second.valid()) //std::max(i.second.mapq1, i.second.mapq2) >= 5)
      os_discordant << i.second.toFileString(opt::read_tracking) << std::endl;
  
  // write ALL contigs
  if (opt::verbose > 2)
    std::cerr << "...writing contigs" << std::endl;
  
  // write the contigs to a BAM
  if (!opt::disc_cluster_only) { 
    for (auto& i : wu.m_contigs) {
      i.RemoveTag("MC");
      b_contig_writer.WriteRecord(i);
    }
  }

  // send breakpoints to file
  for (auto& i : wu.m_bps) {
    if ( i.hasMinimal() && (i.confidence != "NOLOCAL" || i.complex_local))
      os_allbps << i.toFileString(!opt::read_tracking) << std::endl;
  }

  // clear them out
  wu.clear();

}

