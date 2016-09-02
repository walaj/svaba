#include "run_snowman.h"

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
#include "SnowmanUtils.h"
#include "LearnBamParams.h"

#include "SeqLib/BFC.h"

// output files
static ogzstream all_align, os_allbps, all_disc_stream, os_cigmap, r2c_c;
static std::ofstream log_file, bad_bed;
static std::stringstream ss;

#define PRETTYFORM(num) SeqLib::AddCommas(num)
#define MATETOGRC(m, g)				\
for (auto& i : (m))\
g.add(SeqLib::GenomicRegion(i.chr, i.pos1, i.pos2, i.strand));

#define WRITELOG(msg, toerr, tolog)   \
  { if (tolog) log_file  << (msg) << std::endl;	\
    if (toerr) std::cerr << (msg) << std::endl; };

#define READ_LOOKUP_MAIN_REGION_LIMIT 10000
#define MIN_CONTIG_MATCH 35
#define MATE_LOOKUP_MIN 3
#define SECONDARY_CAP 200
#define MAX_MATE_ROUNDS 3
#define MATE_REGION_LOOKUP_LIMIT 400

#define GERMLINE_CNV_PAD 10
#define LITTLECHUNK 25000 
#define WINDOW_PAD 500
#define MICROBE_MATCH_MIN 50
#define GET_MATES 1
#define MICROBE 1
#define LARGE_INTRA_LOOKUP_LIMIT 1000000
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

static BamIndexMap bindices;

static int min_dscrd_size_for_variant = 0; // set a min size for what we can call with discordant reads only. 
// something like max(mean + 3*sd) for all read groups

static std::unordered_map<std::string, int> min_isize_for_disc;

static SeqLib::BamReader b_reader;
static SeqLib::BamWriter r2c_writer, er_writer, b_microbe_writer, b_allwriter;
static SeqLib::BWAWrapper * microbe_bwa = nullptr;
static SeqLib::BWAWrapper * main_bwa = nullptr;
static SeqLib::Filter::ReadFilterCollection * mr;
static SeqLib::GRC blacklist, indel_blacklist_mask, germline_svs, simple_seq, bad_mate_regions;
static DBSnpFilter * dbsnp_filter;


static SeqLib::GRC file_regions, regions_torun;

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
    static float error_rate = 0; 
    static bool writeASQG = false;
  }

  static int sample_number = 0;

  static std::string args = "snowman ";

  static std::string kmer_correction = "f";

  static std::string simple_file;

  static int num_assembly_rounds = 2;
  static bool write_extracted_reads = false;

  static bool no_assemble_normal = false;

  static bool no_reads = true;
  static int32_t readlen;
  static bool r2c = false;
  static bool zip = false;
  static int max_mapq_possible;

  static double sd_disc_cutoff = 4;

  static int gap_open_penalty = 32; //16; //6;
  static int gap_extension_penalty = 1;
  static int mismatch_penalty = 18; //9; //4;
  static int sequence_match_score = 2;
  static int zdrop = 100;
  static int bandwidth = 1000;//100
  static float reseed_trigger = 1.5;
  static int clip3_pen = 5;
  static int clip5_pen = 5;

  // parameters for filtering reads
  static std::string rules = "{\"global\" : {\"duplicate\" : false, \"qcfail\" : false}, \"\" : { \"rules\" : [FRRULES,{\"rr\" : true},{\"ff\" : true}, {\"rf\" : true}, {\"ic\" : true}, {\"clip\" : 5, \"length\" : READLENLIM}, {\"ins\" : true}, {\"del\" : true}, {\"mapped\": true , \"mate_mapped\" : false}, {\"mate_mapped\" : true, \"mapped\" : false}]}}";  
  static int num_to_sample = 2000000;
  
  static int max_cov = 100;
  static int chunk = LITTLECHUNK; // 1000000;

  // runtime parameters
  static int verbose = 1;
  static int numThreads = 1;

  // data
  static BamMap bam;
  static std::string refgenome = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
  static std::string microbegenome; //"/xchip/gistic/Jeremiah/Projects/SnowmanFilters/viral.1.1.genomic_ns.fna";  
  static std::string analysis_id = "no_id";

  //subsample
  float subsample = 1.0;

  static std::string regionFile;
  static std::string blacklist; // = "/xchip/gistic/Jeremiah/Projects/HengLiMask/um75-hs37d5.bed.gz";
  static std::string germline_sv_file;

  static std::string dbsnp; // = "/xchip/gistic/Jeremiah/SnowmanFilters/dbsnp_138.b37_indel.vcf";

  // filters on when / how to assemble
  static bool disc_cluster_only = false;

  static bool refilter = false;

  static std::string gemcode;

  static std::string normal_bam;
  static std::string tumor_bam;

  static double kmer_subsample = 0.25;

  static bool germline = false; // optimize defaults for single sample mode
  
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
  OPT_KMER_SUBSAMPLE,
  OPT_NO_ASSEMBLE_NORMAL,
  OPT_MAX_COV,
  OPT_DISCORDANT_ONLY,
  OPT_REFILTER,
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
  OPT_GERMLINE
};

static const char* shortopts = "hzxt:n:p:v:r:G:r:e:g:k:c:a:m:B:D:Y:S:P:L:O:Is:V:R:K:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "tumor-bam",               required_argument, NULL, 't' },
  { "germline",               no_argument, NULL, OPT_GERMLINE },  
  //  { "indel-mask",          required_argument, NULL, 'M' },
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
  { "disc-sd-cutoff",          required_argument, NULL, 's' },
  { "mate-lookup-min",         required_argument, NULL, 'L' },
  { "germline-sv-database",    required_argument, NULL, 'V' },
  { "simple-seq-database",     required_argument, NULL, 'R' },
  { "g-zip",                   no_argument, NULL, 'z' },
  { "no-interchrom-lookup",    no_argument, NULL, 'I' },
  { "read-tracking",           no_argument, NULL, OPT_READ_TRACK },
  { "gap-open-penalty",        required_argument, NULL, OPT_GAP_OPEN },
  { "ec-subsample",          required_argument, NULL, OPT_KMER_SUBSAMPLE },
  { "bwa-match-score",         required_argument, NULL, OPT_MATCH_SCORE },
  { "gap-extension-penalty",   required_argument, NULL, OPT_GAP_EXTENSION },
  { "mismatch-penalty",        required_argument, NULL, OPT_MISMATCH },
  { "z-dropoff",               required_argument, NULL, OPT_ZDROP },
  { "reseed-trigger",          required_argument, NULL, OPT_RESEED_TRIGGER },
  { "penalty-clip-3",          required_argument, NULL, OPT_CLIP3 },
  { "penalty-clip-5",          required_argument, NULL, OPT_CLIP5 },
  { "bandwidth",               required_argument, NULL, OPT_BANDWIDTH },
  { "write-extracted-reads", no_argument, NULL, OPT_WRITE_EXTRACTED_READS },
  { "no-assemble-normal",    no_argument, NULL, OPT_NO_ASSEMBLE_NORMAL },
  { "lod",                   required_argument, NULL, OPT_LOD },
  { "lod-dbsnp",             required_argument, NULL, OPT_DLOD },
  { "lod-no-dbsnp",          required_argument, NULL, OPT_NDLOD },
  { "lod-germ",              required_argument, NULL, OPT_LODGERM },
  { "discordant-only",       no_argument, NULL, OPT_DISCORDANT_ONLY },
  { "num-to-sample",         required_argument, NULL, OPT_NUM_TO_SAMPLE },
  { "refilter",              no_argument, NULL, OPT_REFILTER },
  { "r2c-bam",               no_argument, NULL, 'x' },
  { "write-asqg",            no_argument, NULL, OPT_ASQG   },
  { "ec-correct-type",     required_argument, NULL, 'K'},
  { "error-rate",            required_argument, NULL, 'e'},
  { "verbose",               required_argument, NULL, 'v' },
  { "blacklist",             required_argument, NULL, 'B' },
  { "max-coverage",          required_argument, NULL, OPT_MAX_COV },
  { "num-assembly-rounds",   required_argument, NULL, OPT_NUM_ASSEMBLY_ROUNDS },
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
  //"  -M, --indel-mask                     BED-file with graylisted regions for stricter indel calling.\n"
"  -D, --dbsnp-vcf                      DBsnp database (VCF) to compare indels against\n"
"  -B, --blacklist                      BED-file with blacklisted regions to not extract any reads from.\n"
"  -Y, --microbial-genome               Path to indexed reference genome of microbial sequences to be used by BWA-MEM to filter reads.\n"
"  -V, --germline-sv-database           BED file containing sites of known germline SVs. Used as additional filter for somatic SV detection\n"
"  -R, --simple-seq-database            BED file containing sites of simple DNA that can confuse the contig re-alignment.\n"
"  -z, --g-zip                          Gzip and tabix the output VCF files. Default: off\n"
"  -L, --min-lookup-min                 Minimum number of somatic reads required to attempt mate-region lookup [3]\n"
"  -s, --disc-sd-cutoff                 Number of standard deviations of calculated insert-size distribution to consider discordant. [1.96]\n"
"      --no-interchrom-lookup           Don't do mate lookup for inter-chromosomal translocations. This increased run time at cost of power for translocations.\n"
"      --r2c-bam                        Output a BAM of reads that aligned to a contig, and fasta of error corrected sequences\n"
"      --discordant-only                Only run the discordant read clustering module, skip assembly. Default: off\n"
"      --read-tracking                  Track supporting reads. Increases file sizes.\n"
"      --max-coverage                   Maximum weird read coverage to send to assembler (per BAM). Subsample reads to achieve max if overflow. [500]\n"
"      --write-extracted-reads          For the tumor BAM, write the extracted reads (those sent to assembly) to a BAM file. Good for debugging.\n"
"      --assemble-by-tag                Separate the assemblies and read-to-contig mapping by the given read tag. Useful for 10X Genomics data (e.g. --aseembly-by-tag BX).\n"
"      --no-assemble-normal             Don't error correct or assembly normal reads. Still maps normal reads to contigs. Increases speed for somatic-only queries.\n"
"      --num-assembly-rounds            Number of times to run the assembler per window. > 1 will bootstrap the assembly with error rate = 0.05. [1]\n"
"      --num-to-sample                  When learning about BAM, number of reads to sample. Default 1,000,000.\n"
"      --motif-filter                   Add a motif file to exclude reads that contain this motif (e.g. illumina adapters).\n"
"  Assembly params\n"
"      --write-asqg                     Output an ASQG graph file for each assembly window.\n"
"  -e, --error-rate                     Fractional difference two reads can have to overlap. See SGA param. 0 is fast, but requires exact matches and error correcting. [0]\n"
"  -c, --chunk-size                     Size of a local assembly window (in bp). [25000]\n"
"      --ec-correct-type                (f) Fermi-kit BFC correction (default), (s) Kmer-correction implemented in SGA, (0) no correction. Need to balance with error rate (w/EC can have lower error rate -e)\n"
"      --ec-subsample                   Learn from only a fraction of the reads during error-correction. Reduces memory and increases speed, especially for high-coverage samples [0.24]\n"
"  Alignment params\n"
"      --bwa-match-score                Set the BWA-MEM match score. BWA-MEM -A [1]\n"
"      --gap-open-penalty               Set the BWA-MEM gap open penalty for contig to genome alignments. BWA-MEM -O [6]\n"
"      --gap-extension-penalty          Set the BWA-MEM gap extension penalty for contig to genome alignments. BWA-MEM -E [1]\n"
"      --mismatch-penalty               Set the BWA-MEM mismatch penalty for contig to genome alignments. BWA-MEM -b [4]\n"
"      --bandwidth                      Set the BWA-MEM SW alignment bandwidth for contig to genome alignments. BWA-MEM -w [100]\n"
"      --z-dropoff                      Set the BWA-MEM SW alignment Z-dropoff for contig to genome alignments. BWA-MEM -d [100]\n"
"      --reseed-trigger                 Set the BWA-MEM reseed trigger for reseeding mems for contig to genome alignments. BWA-MEM -r [1.5]\n"
"      --penalty-clip-3\n"
"      --penalty-clip-5\n"
"\n";

void runSnowman(int argc, char** argv) {

  parseRunOptions(argc, argv);

  // test
  //run_test_assembly();

  SnowmanUtils::fopen(opt::analysis_id + ".log", log_file);
  SnowmanUtils::fopen(opt::analysis_id + ".bad_mate_regions.bed", bad_bed);

  if (opt::refilter) {
    // read in all of the breakpoints
    makeVCFs();
    return;
  }

  // set the germline parameters 
  if (opt::germline) {
    std::string rules = "{\"global\" : {\"duplicate\" : false, \"qcfail\" : false}, \"\" : { \"rules\" : [FRRULES,{\"rr\" : true},{\"ff\" : true}, {\"rf\" : true}, {\"ic\" : true}, {\"clip\" : 5, \"length\" : READLENLIM}, {\"ins\" : true}, {\"del\" : true}, {\"mapped\": true , \"mate_mapped\" : false}, {\"mate_mapped\" : true, \"mapped\" : false}, {\"nm\" : [3,0]}]}}";
    opt::interchrom_lookup = false;
  }
  
  ss << 
    "-----------------------------------------------------------------" << std::endl << 
    "--- Running Snowman somatic indel and rearrangement detection ---" << std::endl <<
    "-----------------------------------------------------------------" << std::endl;
  ss << 
    "***************************** PARAMS ****************************" << std::endl << 
    "    DBSNP Database file: " << opt::dbsnp << std::endl << 
    "    Max cov to assemble: " << opt::max_cov << std::endl <<
    "    Error correction mode: " << (opt::kmer_correction == "s" ? "SGA" : "BFC" << std::endl << 
    "    Subsample-rate for correction learning: " + std::to_string(opt::kmer_subsample) << std::endl;
    ss << 
      "    ErrorRate: " << (opt::assemb::error_rate < 0.001f ? "EXACT (0)" : std::to_string(opt::assemb::error_rate)) << std::endl << 
      "    Num assembly rounds: " << opt::num_assembly_rounds << std::endl;
  ss << 
    "    Num reads to sample: " << opt::num_to_sample << std::endl << 
    "    Discordant read extract SD cutoff:  " << opt::sd_disc_cutoff << std::endl << 
    "    Discordant cluster std-dev cutoff:  " << (opt::sd_disc_cutoff * 1.5) << std::endl << 
    "    Minimum number of reads for mate lookup " << opt::mate_lookup_min << std::endl <<
    "    LOD cutoff (artifact vs real) :  " << opt::lod << std::endl << 
    "    LOD cutoff (somatic vs germline, at DBSNP):  " << opt::lod_db << std::endl << 
    "    LOD cutoff (somatic vs germlin, no DBSNP):  " << opt::lod_no_db << std::endl << 
    "    LOD cutoff (germline, AF>=0.5 vs AF=0):  " << opt::lod_germ << std::endl <<
    "    BWA-MEM params:" << std::endl <<
    "      Gap open penalty: " << opt::gap_open_penalty << std::endl << 
    "      Gap extension penalty: " << opt::gap_extension_penalty << std::endl <<
    "      Mismatch penalty: " << opt::mismatch_penalty << std::endl <<
    "      Aligment bandwidth: " << opt::bandwidth << std::endl <<
    "      Z-dropoff: " << opt::zdrop << std::endl <<
    "      Clip 3 penalty: " << opt::clip3_pen << std::endl <<
    "      Clip 5 penalty: " << opt::clip5_pen << std::endl <<
    "      Reseed trigger: " << opt::reseed_trigger << std::endl <<
    "      Sequence match score: " << opt::sequence_match_score << std::endl;

  if (opt::assemb::writeASQG)
    ss << "    Writing ASQG files. Suggest running R/snow-asqg2pdf.R -i <my.asqg> -o graph.pdf" << std::endl;
  if (opt::write_extracted_reads)
    ss << "    Writing extracted reads and fasta of error-corrected reads." << std::endl;
  if (opt::r2c)
    ss << "    Writing read-to-contig files." << std::endl;
  if (opt::disc_cluster_only)
    ss << "    ######## ONLY DISCORDANT READ CLUSTERING. NO ASSEMBLY ##############" << std::endl;
  if (!opt::interchrom_lookup)
    ss << "    ######## NOT LOOKING UP MATES FOR INTERCHROMOSOMAL #################" << std::endl;
  ss <<
    "*****************************************************************" << std::endl;	  

  SnowmanUtils::print(ss, log_file, opt::verbose > 0);

  // make one anyways, we check if its empty later
  ref_genome_viral = new SeqLib::RefGenome;
  microbe_bwa = nullptr;
  
  // open the microbe genome
  
  if (!opt::microbegenome.empty()) {
    WRITELOG("...loading the microbe reference sequence", opt::verbose > 0, true)
    microbe_bwa = new SeqLib::BWAWrapper();
    SnowmanUtils::__open_index_and_writer(opt::microbegenome, microbe_bwa, opt::analysis_id + ".microbe.bam", b_microbe_writer, ref_genome_viral, viral_header);  
  }

  // open the tumor bam to get header info
  b_reader.Open(opt::tumor_bam);

  // open some writer bams
  if (opt::r2c) // open the r2c writer
    SnowmanUtils::__openWriterBam(b_reader.Header(), opt::analysis_id + ".r2c.bam", r2c_writer);
  if (opt::write_extracted_reads) // open the extracted reads writer
    SnowmanUtils::__openWriterBam(b_reader.Header(), opt::analysis_id + ".extracted.reads.bam", er_writer);    

  // open the blacklists
  SnowmanUtils::__open_bed(opt::blacklist, blacklist, b_reader.Header());
  if (blacklist.size())
    ss << "...loaded " << blacklist.size() << " blacklist regions from " << opt::blacklist << std::endl;

  // open the germline sv database
  SnowmanUtils::__open_bed(opt::germline_sv_file, germline_svs, b_reader.Header());
  if (germline_svs.size())
    ss << "...loaded " << germline_svs.size() << " germline SVs from " << opt::germline_sv_file << std::endl;

  // open the simple seq database
  SnowmanUtils::__open_bed(opt::simple_file, simple_seq, b_reader.Header());
  if (simple_seq.size())
    ss << "...loaded " << simple_seq.size() << " simple sequence regions from " << opt::simple_file << std::endl;

  // open the DBSnpFilter
  if (opt::dbsnp.length()) {
    WRITELOG("...loading the DBsnp database", opt::verbose > 0, true)
      dbsnp_filter = new DBSnpFilter(opt::dbsnp, b_reader.Header());
    WRITELOG("...loaded DBsnp database", opt::verbose > 0, true)
  }

  // open the Bam Headeres
  WRITELOG("...loading BAM indices", opt::verbose > 1, true);
  for (auto& b : opt::bam) 
    bindices[b.first] = SeqLib::SharedIndex(hts_idx_load(b.second.c_str(), HTS_FMT_BAI), idx_delete());

  // set the prefixes
  for (auto& b : opt::bam)
    prefixes.insert(b.first);

  // learn the bam
  min_dscrd_size_for_variant = 0; // set a min size for what we can call with discordant reads only. 
  for (auto& b : opt::bam) {
    LearnBamParams parm(b.second);
    params_map[b.first] = BamParamsMap();
    parm.learnParams(params_map[b.first], opt::num_to_sample);
    for (auto& i : params_map[b.first]) {
      opt::readlen = std::max(opt::readlen, i.second.readlen);
      opt::max_mapq_possible = std::max(opt::max_mapq_possible, i.second.max_mapq);
      min_dscrd_size_for_variant = std::max(min_dscrd_size_for_variant, (int)std::floor(i.second.mean_isize + i.second.sd_isize * opt::sd_disc_cutoff * 1.5)); 
    }

    ss << "BAM PARAMS FOR: " << b.first << "--" << b.second << std::endl;
    for (auto& i : params_map[b.first])
      ss << i.second << std::endl;
    ss << " min_dscrd_size_for_variant " << min_dscrd_size_for_variant << std::endl;
  }

  // check if differing read lengths or max mapq
  for (auto& a : params_map) {
    for (auto& i : a.second) {
      if (i.second.readlen != opt::readlen)
	std::cerr << "!!!! WARNING. Multiple readlengths mixed: " << i.first << "--" << i.second.readlen << std::endl;
      if (i.second.max_mapq != opt::max_mapq_possible)
	std::cerr << "!!!! WARNING. Multiple max mapq mixed: " << i.first << "--" << i.second.max_mapq << std::endl;
    }
  }

  ss << "...min discordant-only variant size " << min_dscrd_size_for_variant << std::endl;

  // set the min overlap
  if (opt::assemb::minOverlap < 0) 
    opt::assemb::minOverlap = (0.6 * opt::readlen) < 30 ? 30 : 0.5 * opt::readlen;

  ss << "...found read length of " << opt::readlen << ". Min Overlap is " << opt::assemb::minOverlap << std::endl;
  ss << "...max read MAPQ detected: " << opt::max_mapq_possible << std::endl;
  SnowmanUtils::print(ss, log_file, opt::verbose);
  
  // get the seed length for printing
  int seedLength, seedStride;
  SnowmanAssemblerEngine enginetest("test", opt::assemb::error_rate, opt::assemb::minOverlap, opt::readlen);
  enginetest.calculateSeedParameters(opt::readlen, opt::assemb::minOverlap, seedLength, seedStride);
  WRITELOG("...calculated seed size for error rate of " + std::to_string(opt::assemb::error_rate) + " and read length " +
	   std::to_string(opt::readlen) + " is " + std::to_string(seedLength), opt::verbose, true);

  // open the human reference
  WRITELOG("...loading the human reference sequence for BWA", opt::verbose, true);
  SnowmanUtils::print(ss, log_file, opt::verbose > 0);
  main_bwa = new SeqLib::BWAWrapper();
  main_bwa->SetAScore(opt::sequence_match_score);
  main_bwa->SetGapOpen(opt::gap_open_penalty);
  main_bwa->SetGapExtension(opt::gap_extension_penalty);
  main_bwa->SetMismatchPenalty(opt::mismatch_penalty);
  main_bwa->SetZDropoff(opt::zdrop);
  main_bwa->SetBandwidth(opt::bandwidth);
  main_bwa->SetReseedTrigger(opt::reseed_trigger);
  main_bwa->Set3primeClippingPenalty(opt::clip3_pen);
  main_bwa->Set5primeClippingPenalty(opt::clip5_pen);

  // open the reference for reading seqeuence
  ref_genome = new SeqLib::RefGenome;
  SnowmanUtils::__open_index_and_writer(opt::refgenome, main_bwa, opt::analysis_id + ".contigs.bam", b_allwriter, ref_genome, bwa_header);
  assert(!ref_genome->IsEmpty());
  
  // parse the region file, count number of jobs
  int num_jobs = SnowmanUtils::countJobs(opt::regionFile, file_regions, regions_torun,
					 bwa_header, opt::chunk, WINDOW_PAD); 
  if (num_jobs) {
    WRITELOG("...running on " + PRETTYFORM(num_jobs) + " chunks", opt::verbose, true);
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
    opt::rules = myreplace(opt::rules, "FRRULES", string_rules);
  }

  // set min length for clips
  if (opt::rules.find("READLENLIM") != std::string::npos) {
    opt::rules = myreplace(opt::rules, "READLENLIM", std::to_string((int) (opt::readlen * 0.4)));
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
      printf("\n mutex init failed\n");
      return;
  }

  // open the text files
  SnowmanUtils::fopen(opt::analysis_id + ".alignments.txt.gz", all_align);
  SnowmanUtils::fopen(opt::analysis_id + ".bps.txt.gz", os_allbps);
  SnowmanUtils::fopen(opt::analysis_id + ".discordant.txt.gz", all_disc_stream);
  if (opt::write_extracted_reads) 
    SnowmanUtils::fopen(opt::analysis_id + ".corrected.fa.gz", r2c_c); 
  
  // write the headers to the text files
  os_allbps << BreakPoint::header();
  for (auto& b : opt::bam) 
    os_allbps << "\t" << b.first << "_" << b.second;
  os_allbps << std::endl;
  all_disc_stream << DiscordantCluster::header() << std::endl;

  // put args into string for VCF later
  for (int i = 0; i < argc; ++i)
    opt::args += std::string(argv[i]) + " ";

  // start the timer
#ifndef __APPLE__
  clock_gettime(CLOCK_MONOTONIC, &start);
#endif

  // send the jobs to the queue
  WRITELOG("--- Starting detection pipeline --- on " + PRETTYFORM(opt::numThreads) + " thread" + (opt::numThreads > 1 ? "s" : ""), opt::verbose, true);
  sendThreads(regions_torun);

  // clean up everything
  if (microbe_bwa)
    delete microbe_bwa;

  // dump the bad bed regions
  bad_mate_regions.CoordinateSort();
  for (auto& i : bad_mate_regions)
    try { // throws out of range if there is mismatch between main header and other headers. just be graceful here
      bad_bed << i.ChrName(b_reader.Header()) << "\t" << i.pos1 << "\t" << i.pos2 << "\t" << i.strand << std::endl;
    } catch (...) {
    }
  bad_bed.close();

  // close the files
  all_align.close();
  os_allbps.close();
  all_disc_stream.close();
  if (opt::r2c) r2c_c.close();
  log_file.close();

  // make the VCF file
  makeVCFs();

  // more clean up 
  if (ref_genome)
    delete ref_genome;
  if (ref_genome_viral)
    delete ref_genome_viral;
  
  
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
  //if (!SeqLib::read_access_test(file))
  //  file = opt::analysis_id + ".bps.txt";

  // make the header
  VCFHeader header;
  header.filedate = SnowmanUtils::fileDateString();
  header.source = opt::args;
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

  // primary VCFs
  if (SeqLib::read_access_test(file)) {
    WRITELOG("...making the primary VCFs (unfiltered and filtered) from file " + file, opt::verbose, true);
    VCFFile snowvcf(file, opt::analysis_id, b_reader.Header(), header);
    std::string basename = opt::analysis_id + ".snowman.unfiltered.";

    WRITELOG("...writing unfiltered VCFs", opt::verbose, true);
    snowvcf.include_nonpass = true;
    snowvcf.writeIndels(basename, opt::zip, opt::bam.size() == 1);
    snowvcf.writeSVs(basename, opt::zip, opt::bam.size() == 1);

    WRITELOG("...writing filtered VCFs", opt::verbose, true);
    basename = opt::analysis_id + ".snowman.";
    snowvcf.include_nonpass = false;
    snowvcf.writeIndels(basename, opt::zip, opt::bam.size() == 1);
    snowvcf.writeSVs(basename, opt::zip, opt::bam.size() == 1);

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

  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case OPT_REFILTER : opt::refilter = true; break;
    case OPT_KMER_SUBSAMPLE : arg >> opt::kmer_subsample; break;
      case 'p': arg >> opt::numThreads; break;
      case 'm': arg >> opt::assemb::minOverlap; break;
      case 'a': arg >> opt::analysis_id; break;
      case 'B': arg >> opt::blacklist; break;
    case 'V': arg >> opt::germline_sv_file; break;
    case 'R': arg >> opt::simple_file; break;
      case 'L': arg >> opt::mate_lookup_min; break;
    case 'I': opt::interchrom_lookup = false; break;
      case 'Y': arg >> opt::microbegenome; break;
      case 'z': opt::zip = false; break;
      case 'h': help = true; break;
    case OPT_GAP_OPEN : arg >> opt::gap_open_penalty; break;
    case OPT_MATCH_SCORE : arg >> opt::sequence_match_score; break;
    case OPT_MISMATCH : arg >> opt::mismatch_penalty; break;
    case OPT_GERMLINE : opt::germline = true; break;
    case OPT_GAP_EXTENSION : arg >> opt::gap_extension_penalty; break;
    case OPT_ZDROP : arg >> opt::zdrop; break;
    case OPT_BANDWIDTH : arg >> opt::bandwidth; break;
    case OPT_RESEED_TRIGGER : arg >> opt::reseed_trigger; break;
    case OPT_CLIP3 : arg >> opt::clip3_pen; break;
    case OPT_CLIP5 : arg >> opt::clip5_pen; break;
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
    case OPT_NO_ASSEMBLE_NORMAL: opt::no_assemble_normal = true; break;
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
      case 'D': arg >> opt::dbsnp; break;
    case 's': arg >> opt::sd_disc_cutoff; break;
      case 'r': 
	arg >> opt::rules; 
	if (opt::rules == "all")
	  opt::rules = "";
       	break;
      case OPT_DISCORDANT_ONLY: opt::disc_cluster_only = true; break;
      case OPT_WRITE_EXTRACTED_READS: opt::write_extracted_reads = true; break;
    case OPT_NUM_ASSEMBLY_ROUNDS: arg >> opt::num_assembly_rounds; break;
    case 'K': arg >> opt::kmer_correction; break;
      default: die= true; 
    }
  }

  if (!(opt::kmer_correction == "s" || opt::kmer_correction == "f" || opt::kmer_correction = "0")) {
    WRITELOG("ERROR: Error correction type must be one of s, f, or 0", true, true);
    exit(EXIT_FAILURE);
  }

  // check that BAM files exist
  for (auto& b : opt::bam)
    if (!SeqLib::read_access_test(b.second)) {
      WRITELOG("ERROR: BAM file" + b.second + " is not readable / existant", true, true);
      exit(EXIT_FAILURE);
    }

  // check that we input something
  if (opt::bam.size() == 0 && !die && !opt::refilter) {
    WRITELOG("Must add a bam file with -t flag", true, true);
    exit(EXIT_FAILURE);
  }

  if (opt::numThreads <= 0) {
    WRITELOG("Invalid number of threads from -p flag: " + PRETTYFORM(opt::numThreads), true, true);
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

bool runBigChunk(const SeqLib::GenomicRegion& region)
{
  // setup structures to store the final data for this region
  std::vector<AlignedContig> alc;
  SeqLib::BamRecordVector all_contigs, all_microbial_contigs;
  SeqLib::GRC badd; // bad mate regions for this

  // start a timer
  SnowmanUtils::SnowTimer st;
  st.start();

  SeqLib::BFC * bfc = new SeqLib::BFC();

  // setup for the BAM walkers
  WalkerMap walkers;
  CountPair read_counts = {0,0};

  // loop through the input BAMs and get reads
  for (auto& b : opt::bam)
    walkers[b.first] = make_walkers(b.first, b.second, region, read_counts, badd, bfc);

  // collect all of the cigar strings in a hash
  std::unordered_map<std::string, SeqLib::CigarMap> cigmap;
  for (auto& w : walkers) 
    cigmap[w.first] = w.second.cigmap;

  // setup read collectors
  std::vector<char*> all_seqs;
  SeqLib::BamRecordVector bav_this;

  // collect and clear reads from main round
  std::unordered_set<std::string> dedupe;
  collect_and_clear_reads(walkers, bav_this, all_seqs, dedupe);

  // adjust counts and timer
  st.stop("r");

  // get the mate reads
  CountPair mate_counts = run_mate_collection_loop(region, walkers, badd);
  WRITELOG("...looked up <case, control> mate reads: <" + SeqLib::AddCommas(mate_counts.first) + "," + 
	   SeqLib::AddCommas(mate_counts.second) + ">", opt::verbose > 1, true);
  // collect the reads together from the mate walkers
  collect_and_clear_reads(walkers, bav_this, all_seqs, dedupe);

  st.stop("m");
  
  // do the discordant read clustering
  WRITELOG("...discordant read clustering", opt::verbose > 1, false);
  DiscordantClusterMap dmap = DiscordantCluster::clusterReads(bav_this, region, opt::max_mapq_possible, &min_isize_for_disc);

  // tag FR clusters that are below min_dscrd_size_for_variant AND low support
  DiscordantClusterMap dmap_tmp;
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
    WRITELOG("TOO MANY READS IN REGION " + PRETTYFORM(bav_this.size()) + "\t" + region.ToString(), opt::verbose, false);
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
  if (opt::kmer_correction == "s") {
    correct_reads(all_seqs, bav_this);
  } else if (opt::kmer_correction == "f") {
    
    int learn_reads_count = bfc->NumSequences();
    bfc->Train();
    bfc->clear();  // clear memory and reads

    // add the reads for correction
    bfc->ErrorCorrectToTag(bav_this, "KC");

    double kcov = bfc->GetKCov();
    int kmer    = bfc->GetKMer();

    WRITELOG("   BFC attempted correct " + std::to_string(bav_this.size()) + " train: " + 
	     std::to_string(learn_reads_count) + " kcov: " + std::to_string(kcov) + 
	     " kmer: " + std::to_string(kmer), opt::verbose > 1, true);      
    
    delete bfc;
    bfc = nullptr;
  }
  st.stop("k");
  
  // do the assembly, contig realignment, contig local realignment, and read realignment
  // modifes bav_this, alc, all_contigs and all_microbial_contigs
  run_assembly(region, bav_this, alc, all_contigs, all_microbial_contigs, dmap, cigmap);

afterassembly:
  
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
    if (!i.second.hasAssociatedAssemblyContig() && (i.second.tcount + i.second.ncount) > 1 && i.second.valid() && !below_size) {
      BreakPoint tmpbp(i.second, main_bwa, dmap);
      bp_glob.push_back(tmpbp);
    }
  }
  
  // de duplicate the breakpoints
  std::sort(bp_glob.begin(), bp_glob.end());
  bp_glob.erase( std::unique( bp_glob.begin(), bp_glob.end() ), bp_glob.end() );
  
  // add the coverage data to breaks for allelic fraction computation
  std::unordered_map<std::string, STCoverage*> covs;
  for (auto& i : opt::bam) 
    covs[i.first] = &walkers[i.first].cov;

  for (auto& i : bp_glob)
    i.addCovs(covs);

  for (auto& i : bp_glob) 
    i.scoreBreakpoint(opt::lod, opt::lod_db, opt::lod_no_db, opt::lod_germ, min_dscrd_size_for_variant);

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
      if (i.somatic_score && i.b1.gr.chr == i.b2.gr.chr) {
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

  ////////////////////////////////////
  // MUTEX LOCKED
  ////////////////////////////////////
  // write to the global contig out
  pthread_mutex_lock(&snow_lock);  
  
  // add the ref and alt tags
  // WHY IS THIS NOT THREAD SAFE?
  for (auto& i : bp_glob)
    i.setRefAlt(ref_genome, ref_genome_viral);

  // update the bad mate regions
  badd.MergeOverlappingIntervals();
  bad_mate_regions.Concat(badd);
  for (auto& r : badd) 
    WRITELOG("...added bad mate region " + r.ToString(), opt::verbose > 1, true);
  bad_mate_regions.MergeOverlappingIntervals();
  bad_mate_regions.CreateTreeMap();
  
  // print the alignment plots
  size_t contig_counter = 0;
  for (auto& i : alc) {
    if (i.hasVariant()) {
      all_align << i << std::endl;
      ++contig_counter;
    }
  }
  
  // send the microbe to file
  for (auto& b : all_microbial_contigs)
    b_microbe_writer.WriteRecord(b);
  
  // send the discordant to file
  for (auto& i : dmap)
    if (i.second.valid()) //std::max(i.second.mapq1, i.second.mapq2) >= 5)
      all_disc_stream << i.second.toFileString(!opt::no_reads) << std::endl;
  
  // write ALL contigs
  if (opt::verbose > 2)
    std::cerr << "...writing contigs" << std::endl;
  
  // write the contigs to a BAM
  if (!opt::disc_cluster_only) { 
    for (auto& i : all_contigs) {
      i.RemoveTag("MC");
      b_allwriter.WriteRecord(i);
    }
  }
  
  // write extracted reads
  if (opt::write_extracted_reads) {
    for (auto& r : bav_this) //walkers[opt::tumor_bam].reads) 
      //if (!excluded_bad_reads.count(r.GetZTag("SR")))
      er_writer.WriteRecord(r);
  }
  
  // write all the extracted reads to output
  if (!opt::disc_cluster_only && opt::r2c) {
    for (auto& r : walkers[opt::tumor_bam].reads) {
      r2c_writer.WriteRecord(r);
    }
  }
  
  // write the raw error corrected reads to a fasta
  if (opt::write_extracted_reads) 
    for (auto& r : bav_this) {
      std::string seq = r.GetZTag("KC");
      if (seq.empty())
	seq = r.QualitySequence();
      r2c_c << ">" << r.GetZTag("SR") << std::endl << seq << std::endl;
    }
  
  // send breakpoints to file
  for (auto& i : bp_glob) 
    if ( i.hasMinimal() && i.confidence != "NOLOCAL" ) 
      os_allbps << i.toFileString(opt::no_reads) << std::endl;

  st.stop("pp");
  
  // display the run time
  WRITELOG(SnowmanUtils::runTimeString(read_counts.first, read_counts.second, contig_counter, region, b_reader.Header(), st, start), opt::verbose > 1, true);

  ////////////////////////////////////
  // MUTEX UNLOCKED
  ////////////////////////////////////
  pthread_mutex_unlock(&snow_lock);
  
  return true;
}

void sendThreads(SeqLib::GRC& regions_torun) {

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
    SnowmanWorkItem * item     = new SnowmanWorkItem(SeqLib::GenomicRegion(i.chr, i.pos1, i.pos2), ++count);
    queue.add(item);
  }
  if (regions_torun.size() == 0) { // whole genome 
    SnowmanWorkItem * item     = new SnowmanWorkItem(SeqLib::GenomicRegion(), ++count);
    queue.add(item);
  }
    
  
  // wait for the threads to finish
  for (int i = 0; i < opt::numThreads; i++) 
    threadqueue[i]->join();

  // free the items
  // warning: deleting object of polymorphic class type ‘ConsumerThread<SnowmanWorkItem>’
  // which has non-virtual destructor might cause undefined behaviour
  //for (size_t i = 0; i < threadqueue.size(); ++i)
  //  delete threadqueue[i];
}

SeqLib::GRC makeAssemblyRegions(const SeqLib::GenomicRegion& region) {

  // set the regions to run
  SeqLib::GRC grv_small;
  if (region.IsEmpty()) {  // whole genome, so divide up the whole thing
    for (size_t c = 0; c < 23; ++c)
      grv_small.Concat(SeqLib::GRC(opt::chunk, WINDOW_PAD, SeqLib::GenomicRegion(c, WINDOW_PAD + 1, b_reader.Header().GetSequenceLength(c) - WINDOW_PAD - 1)));
    //grv_small.Concat(SeqLib::GRC(opt::chunk, WINDOW_PAD, SeqLib::GenomicRegion(c, WINDOW_PAD + 1, SeqLib::CHR_LEN[c] - WINDOW_PAD - 1)));
  }
  else if (region.Width() >= opt::chunk) // divide into smaller chunks
    grv_small = SeqLib::GRC(opt::chunk, WINDOW_PAD, region);
  else
    grv_small.add(region);

  return grv_small;

}

void alignReadsToContigs(SeqLib::BWAWrapper& bw, const SeqLib::UnalignedSequenceVector& usv, SeqLib::BamRecordVector& bav_this, std::vector<AlignedContig>& this_alc, SeqLib::RefGenome *  rg) {
  
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
  pthread_mutex_lock(&snow_lock);  
  for (auto& i : g)
    if (i.chr < 24) //1-Y
      try {
	std::string tmpref = rg->QueryRegion(i.ChrName(bwa_header), i.pos1, i.pos2);
	ref_alleles.push_back(tmpref); 
      } catch (...) {
	std::cerr << "Caught exception for ref_allele on " << i << std::endl;
      }
  pthread_mutex_unlock(&snow_lock);
  
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

    std::string seqr = i.QualitySequence();
    
    bool hardclip = false;
    assert(seqr.length());
    bw.AlignSequence(seqr, i.Qname(), brv, hardclip, 0.60, 10000);

    if (brv.size() == 0) 
      continue;

    // get the maximum non-reference alignment score
    int max_as = 0;
    for (auto& r : brv)
      max_as = std::max(max_as, r.GetIntTag("AS"));

    // align to the reference alleles
    if (!bw_ref.IsEmpty())
      bw_ref.AlignSequence(seqr, i.Qname(), brv_ref, hardclip, 0.60, 10);

    // get the maximum reference alignment score
    int max_as_r = 0;
    for (auto& r : brv_ref)
      max_as_r = std::max(max_as_r, r.GetIntTag("AS"));
    
    // reject if better alignment to reference
    if (max_as_r > max_as) {
      //std::cerr << " Alignment Rejected for " << max_as_r << ">" << max_as << "  " << i << std::endl;
      //std::cerr << "                        " << max_as_r << ">" << max_as << "  " << brv_ref[0] << std::endl;
      continue;
    }
    
    // make sure we have only one alignment per contig
    std::set<std::string> cc;

    // check which ones pass
    SeqLib::BamRecordVector bpass;
    for (auto& r : brv) {

      // make sure alignment score is OK
      if ((double)r.NumMatchBases() * 0.5 > r.GetIntTag("AS")/* && i.GetZTag("SR").at(0) == 't'*/)
      	continue;
      
      bool length_pass = (r.PositionEnd() - r.Position()) >= ((double)seqr.length() * 0.75);

      if (length_pass && !cc.count(usv[r.ChrID()].Name)) {
	bpass.push_back(r);
	cc.insert(usv[r.ChrID()].Name);
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
      i.SmartAddTag("CN", usv[r.ChrID()].Name);

      for (auto& a : this_alc) {
	if (a.getContigName() == usv[r.ChrID()].Name) {
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

SnowmanBamWalker make_walkers(const std::string& p, const std::string& b, const SeqLib::GenomicRegion& region, CountPair& counts, SeqLib::GRC& badd, SeqLib::BFC * bfc) {


  // setup the walker
  SnowmanBamWalker walk;
  walk.Open(b);
  if (opt::kmer_correction == "f")
    walk.bfc = bfc;
  walk.main_bwa = main_bwa;
  walk.prefix = p;
  walk.blacklist = blacklist;
  walk.do_kmer_filtering = (opt::kmer_correction == "s" || opt::kmer_correction == "f");
  walk.simple_seq = &simple_seq;
  walk.kmer_subsample = opt::kmer_subsample;
  walk.max_cov = opt::max_cov;
  walk.m_mr = mr; 
  walk.m_limit = READ_LOOKUP_MAIN_REGION_LIMIT;

  // set the loaded the index and jump to region
  assert(walk.SetPreloadedIndex(bindices[p]));
  if (!region.IsEmpty()) 
    walk.SetRegion(region);
  badd.Concat(walk.readBam(&log_file)); 

  // adjust the counts
  if (p.at(0) == 't') {
    counts.first += walk.reads.size();
  } else {
    counts.second += walk.reads.size();
  }

  WRITELOG("...read in " + std::to_string(walk.reads.size()) + " reads from region " + p, opt::verbose > 1, true);
  
  return walk;

}

CountPair run_mate_collection_loop(const SeqLib::GenomicRegion& region, WalkerMap& wmap, SeqLib::GRC& badd) {
  
  CountPair counts = {0,0};

  MateRegionVector all_somatic_mate_regions;
  all_somatic_mate_regions.add(MateRegion(region.chr, region.pos1, region.pos2)); // add the origional, don't want to double back
  
  for (int jjj = 0; jjj <  MAX_MATE_ROUNDS; ++jjj) {
    
    MateRegionVector normal_mate_regions;
    if (wmap.size() == 1) // only a case BAM
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
      if (bad_mate_regions.size())
	if (bad_mate_regions.CountOverlaps(s))
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
    }
    
    // print out to log
    for (auto& i : somatic_mate_regions) 
      WRITELOG("...mate region " + i.ToString() + " case read count " + std::to_string(i.count) + " on round " + std::to_string(jjj+1), opt::verbose > 1, true);
    
    // collect the reads for this round
    std::pair<int,int> mate_read_counts = collect_mate_reads(wmap, somatic_mate_regions, jjj, badd);

    // update the counts
    counts.first += mate_read_counts.first;
    counts.second += mate_read_counts.second;
    
  } // mate collection round loop

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

void correct_reads(std::vector<char*>& learn_seqs, SeqLib::BamRecordVector brv) {

  if (!learn_seqs.size()) // || opt::kmer_correction != "s")
    return;

  if (opt::kmer_correction == "s") {
    
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

    WRITELOG("...SGA kmer corrected " + std::to_string(kmer_corrected) + " reads of " + std::to_string(brv.size()), opt::verbose > 1, true);  
  } //else if (opt::kmer_correction == "f") {
    
    //std::cerr << " LEARNING FROM " << learn_seqs.size() << " reads, to correct " << brv.size() << std::endl;
    //bfc.TrainCorrection(learn_seqs);
    //bfc.ErrorCorrectToTag(brv, "KC");
  //}
    
}

void remove_hardclips(SeqLib::BamRecordVector& brv) {
  SeqLib::BamRecordVector bav_tmp;
  for (auto& i : brv)
    if (i.NumHardClip() == 0) 
      bav_tmp.push_back(i);
  brv = bav_tmp;
}

void run_assembly(const SeqLib::GenomicRegion& region, SeqLib::BamRecordVector& bav_this, std::vector<AlignedContig>& master_alc, 
		  SeqLib::BamRecordVector& master_contigs, SeqLib::BamRecordVector& master_microbial_contigs, DiscordantClusterMap& dmap,
		  std::unordered_map<std::string, SeqLib::CigarMap>& cigmap) {

  // get the local region
  pthread_mutex_lock(&snow_lock);  
  std::string lregion;
  try {
    lregion = ref_genome->QueryRegion(bwa_header.IDtoName(region.chr), region.pos1, region.pos2);
  } catch (...) {
    WRITELOG(" Caught exception for lregion with reg " + region.ToString(), true, true);
    lregion = "";
  }
  pthread_mutex_unlock(&snow_lock);

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
  SnowmanAssemblerEngine engine(name, opt::assemb::error_rate, opt::assemb::minOverlap, opt::readlen);
  if (opt::assemb::writeASQG)
    engine.setToWriteASQG();
  if (!opt::no_assemble_normal)
    engine.fillReadTable(bav_this);

  // do the actual assembly
  engine.performAssembly(opt::num_assembly_rounds);
  
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
    if ((int)i.Seq.length() < (opt::readlen * 1.15))
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
      //std::cout << aa << std::endl;
      if (aa.NumClip() < MIN_CLIP_FOR_LOCAL && aa.GetIntTag("NM") < MAX_NM_FOR_LOCAL)
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

    if (microbe_bwa && !SnowmanUtils::hasRepeat(i.Seq)) {
      
      // do the microbial alignment
      SeqLib::BamRecordVector microbial_alignments;
      bool hardclip = false;
      microbe_bwa->AlignSequence(i.Seq, i.Name, microbial_alignments, hardclip, SECONDARY_FRAC, SECONDARY_CAP);
      
      // if the microbe alignment is large enough and doesn't overlap human...
      for (auto& j : microbial_alignments) {
	// keep only long microbe alignments with decent mapq
	if (j.NumMatchBases() >= MICROBE_MATCH_MIN && j.MapQuality() >= 10) { 
	  if (SnowmanUtils::overlapSize(j, ct_alignments) <= 20) { // keep only those where most do not overlap human
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
      if (SnowmanUtils::overlapSize(j, ct_plus_microbe) < 0.5 * j.Length() && j.NumMatchBases() >= MIN_CONTIG_MATCH) { 
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
	SeqLib::GRC ovl = simple_seq.FindOverlaps(k.asGenomicRegion(), true);
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
  alignReadsToContigs(bw, usv, bav_this, this_alc, ref_genome);
  
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
    // check the indel breakpoints against the indel blacklist. 
    // simply set the blacklist flag for these breaks if they hit
    a.blacklist(indel_blacklist_mask);

    // add to the final structure
    master_alc.push_back(a);
    
  }


}

CountPair collect_mate_reads(WalkerMap& walkers, const MateRegionVector& mrv, int round, SeqLib::GRC& badd) {

  CountPair counts = {0,0};  

  if (!mrv.size())
    return counts;
  
  for (auto& w : walkers) {
    int oreads = w.second.reads.size();
    w.second.m_limit = MATE_REGION_LOOKUP_LIMIT;
    SeqLib::GRC gg;
    for (auto& s : mrv)
      gg.add(SeqLib::GenomicRegion(s.chr, s.pos1, s.pos2, s.strand));
    assert(w.second.SetMultipleRegions(gg));
    w.second.get_coverage = false;
    w.second.get_mate_regions = (round != MAX_MATE_ROUNDS);
    badd.Concat(w.second.readBam(&log_file));
    
    // update the counts
    if (w.first == "t") 
      counts.first += (w.second.reads.size() - oreads);
    else
      counts.second += (w.second.reads.size() - oreads);
  }
  
  return counts;
}

void collect_and_clear_reads(WalkerMap& walkers, SeqLib::BamRecordVector& brv, std::vector<char*>& learn_seqs, std::unordered_set<std::string>& dedupe) {

  // concatenate together all the reads from the different walkers
  for (auto& w : walkers) {
    for (auto& r : w.second.reads) {
      std::string sr = r.GetZTag("SR");
      if (!dedupe.count(sr)) {
	brv.push_back(r);
	dedupe.insert(sr);
      }
    }
    
    // concat together all of the learning sequences
    for (auto& r : w.second.all_seqs) {
      learn_seqs.push_back(strdup(r));
      free(r); // free what was alloced in snowmanbamwalker
    }

    w.second.all_seqs.clear();
    w.second.reads.clear();
  }
}
