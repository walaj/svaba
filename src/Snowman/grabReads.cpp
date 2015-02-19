#include "SVBamReader.h"
#include "grabReads.h"
#include "unistd.h"
#include <iostream>
#include <unordered_map>
#include <time.h>
#include <stdlib.h>
#include "SnowUtils.h"
#include "GenomicRegion.h"
#include "api/BamWriter.h"
#include <pthread.h>
#include "zfstream.h"
#include "SnowmanAssemble.h"
#include "SnowmanOverlapAlgorithm.h"
#include "SnowmanASQG.h"
#include "Util.h"
#include "faidx.h"
#include "SeqanTools.h"
#include "BWAWrapper.h"

#include <signal.h>

#include "EncodedString.h"
#include "AlignedContig.h"

#define LEARN_SAMPLE 100000
#define INTER_CHR_COV 5
#define CLIP_COV_LO 3
#define CLIP_COV_HI 6
#define MAPQ_THRESH 8

#define LITTLECHUNK 5000 
#define WINDOW_PAD 300

#define SGAPATH "/home/unix/jwala/sga-master/src/SGA/sga"
// discordant clustering parameters
// DISC_READ_LIMIT is limit to nuimber of disc reads BEFORE looking up region
// DISC_LOOKUP_LIMIT is limit to how many TOTAL reads during one disc region lookup before SVBam fails
#define DISC_READ_LIMIT 2500
#define DISC_LOOKUP_LIMIT 3000
#define DISC_CLUSTER_BUFF 500
#define MIN_DISC_PER_CLUSTER 3
#define MAX_DISC_CLUSTERS 8
#define DISC_LOOKUP_PAD 400 

using namespace std;
using namespace BamTools;

const int MAX_CHARS_PER_LINE = 512;
const int MAX_TOKENS_PER_LINE = 20;
const char* const DELIMITER = ",";

static int read_counter = 0;
static int contig_counter = 0;
static int discordant_counter = 0;

static pthread_mutex_t snow_lock;
//static AlignedContigVec * cont_total;

//static EncodedBAVector * reads_all_encoded;
static BamWriter * r2c_writer;
//static BamWriter * bam_writer;
static BamWriter * disc_writer;

static BamReader * treader_for_convert;

static ofstream * all_align_stream;
static ofstream * os_allbps; 
static ofstream * contigs_sam;

//bwa
static bwaidx_t *idx;

namespace opt {

  static string rules;

  static unsigned verbose = 1;
  static unsigned numThreads = 1;
  static unsigned minOverlap = 35;
  static int cutoff = 120;
  static string tbam;
  static string nbam;
  static int isize = 1000;
  static unsigned sleepDelay = 0;
  static std::string regionFile = "x";
  static std::string outdir = "./";
  static unsigned mapq = 0;
  static unsigned chunk = 1e5;
  static int qualthresh = 4;
  static unsigned skip = 6000; // 6000 is 30x coverage of 20000bp with 100bp reads
  static unsigned numBubbleRounds = 3;
  static float divergence = 0.05;
  static float gap_divergence = 0.05;
  static float error_rate = 0.05;
  static bool skip_disc = false;
  static bool writeASQG = false;
  static bool skip_supp = true;
  static bool skip_r2 = false; //debug
  static bool throw_disc = false;
  static int  memory = 8192;
  static int  memory_goal = 8192;
  static bool learn_params = false;

  static int nmlim = 15;
  static bool assemble_only_disc = false;
  static bool normal_assemble = false;
  static bool disc_cluster_only = false;
  
  // assembly params
  static int maxEdges = 128;
  static int numTrimRounds = 0; //
  static int trimLengthThreshold = -1; // doesn't matter
  static bool bPerformTR = false; // transitivie edge reducetion
  static bool bValidate = false;
  static int resolveSmallRepeatLen = -1; 
  static int maxIndelLength = 20;
  static bool bExact = true;
  static string outVariantsFile = ""; // dummy

  // 
  static bool ignore_skip_cent = false;
  static bool debug = false;
  static size_t contig_write_limit = 50000;
  
  static int min_read_cov = 3;
  static int min_clip = 10;

  static bool contig_bam = false;
  static bool no_r2c = false;

  static string refgenome = REFHG19;
}

enum { 
  OPT_NOR2 = 1,
  OPT_NODISC,
  OPT_NOSUPP,
  OPT_ASQG,
  OPT_IGNORE_SKIP_CENT,
  OPT_THROW_DISC,
  OPT_MEMORY,
  OPT_MEMORY_GOAL,
  OPT_LEARN,
  OPT_NM,
  OPT_ASSEMBLE_ONLY_DISC,
  OPT_ASSEMBLE_NORMAL,
  OPT_DEBUG,
  OPT_CONTIG_LIMIT,
  OPT_MIN_READ_COV,
  OPT_MIN_CLIP,
  OPT_CONTIG_BAM,
  OPT_NO_R2CBAM,
  OPT_DISC_CLUSTER_ONLY,
};

static const char* shortopts = "t:n:p:m:l:v:i:s:r:ho:w:c:q:y:g:b:r:e:g:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "tumor-bam",               required_argument, NULL, 't' },
  { "normal-bam",              required_argument, NULL, 'n' },
  { "output-directory",        required_argument, NULL, 'o' },
  { "threads (\"processes\")", required_argument, NULL, 'p' },
  { "min-overlap",             required_argument, NULL, 'm' },
  { "length-cutoff",           required_argument, NULL, 'l' },
  { "insert-size",             required_argument, NULL, 'i' },
  { "sleep-delay",             required_argument, NULL, 'i' },
  { "map-quality",             required_argument, NULL, 'w'},
  { "skip-num",                required_argument, NULL, 'y'},
  { "qual-thresh",             required_argument, NULL, 'q'},
  { "region-file",             required_argument, NULL, 'r' },
  { "reference-genome",        required_argument, NULL, 'g' },
  //  { "chunk-size",              required_argument, NULL, 'c' },
  { "nm-limit",              required_argument, NULL, OPT_NM },
  { "error-overlap",           required_argument, NULL, 'e' },
  { "write-asqg",              no_argument, NULL, OPT_ASQG   },
  { "assemble-only-disc",              no_argument, NULL, OPT_ASSEMBLE_ONLY_DISC   },
  { "no-disc-lookup",          no_argument, NULL, OPT_NODISC },
  { "no-supplementary"     ,   no_argument, NULL, OPT_NOSUPP },
  { "no-R2"     ,              no_argument, NULL, OPT_NOR2   },
  { "ignore-skip-centromere"     ,              no_argument, NULL, OPT_IGNORE_SKIP_CENT   },
  { "throw-discordant"     ,   no_argument, NULL, OPT_THROW_DISC   },
  { "sort-memory"     ,        required_argument, NULL, OPT_MEMORY   },
  { "memory-goal"     ,        required_argument, NULL, OPT_MEMORY_GOAL   },
  { "learn"     ,              no_argument, NULL, OPT_LEARN   },
  { "assemble-normal"     ,    no_argument, NULL, OPT_ASSEMBLE_NORMAL   },
  { "debug"     ,    no_argument, NULL, OPT_DEBUG   },
  { "contig-write-limit"     ,    required_argument, NULL, OPT_CONTIG_LIMIT   },
  { "min-read-cov"     ,    required_argument, NULL, OPT_MIN_READ_COV   },
  { "min-clip"     ,    required_argument, NULL, OPT_MIN_CLIP   },
  { "contig-bam"     ,   no_argument, NULL, OPT_CONTIG_BAM   },
  { "no-r2c-bam"     ,   no_argument, NULL, OPT_NO_R2CBAM   },
  { "disc-cluster-only"     ,   no_argument, NULL, OPT_DISC_CLUSTER_ONLY   },
  { "verbose",                 required_argument, NULL, 'v' },
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
"      --sort-memory                    Memory (in Mb) used to run \"samtools sort\". *Not* a global mem limit. Default: 4096\n"
"  -g, --reference-genome               Path to indexed reference genome to be used by BWA-MEM. Default is Broad hg19 (/seq/reference/...)\n"
  //"      --learn                          Learn the best parameters to use (-m, -l, -c, -i)\n"
"  Required input\n"
"  -t, --tumor-bam                      Tumor BAM file\n"
"  Optional input\n"                       
"  -n, --normal-bam                     Normal BAM file\n"
"  -m, --min-overlap                    Minimum read overlap, an SGA parameter. Default: 0.4* readlength\n"
  //"  -c, --chunk-size                     Size of region to read-in at once. Lower numbers require less memory, but more I/O. Default: 50e6\n"
"  -l, --length-cutoff                  Minimum length of the contig (discared otherwwise). Default: 1.5*readlength\n"
"  -r, --region-file                    Set a region txt file. Format: one region per line, Ex: 1,10000000,11000000\n"
"  -o, --output-directory               Where to place all of the intermediate and final files. Defaul: ./\n"
"      --no-r2c-bam                     Don't output a read to contig bam. Default: true\n"
"      --debug                          Don't delete temporaray files\n"
"      --contig-write-limit             Number of contigs+reads to store in memory before writing. Default: 50000 (about 1Gb max memory per thread)\n"
"      --disc-cluster-only              Only run the discordant read clustering module, skip assembly. Default: off\n"
"  Read filters\n"
"  -i, --insert-size                    Minimum insert size to consider reads discordant. *Set to -1 for no discordant*. Default: 1000\n"
"  -w, --map-quality                    Minimum MAPQ a read must have to move on to assembly. Default: 0\n"
"  -q, --qual-thresh                    Clip all bases with optical quality less than qual-thresh. Default: 4\n"
"      --skip-num                       Do not assemble regions with > skip-num reads in 20 kb window. Default: 6000 (30x cov w/101bp reads)\n"
"      --no-disc-lookup                 Do not retrieve discordant sequences when running on non-SV BAM\n"
  //"      --no-supplementary               Do not include supplementary alignments into the assembly\n"
"      --no-R2                          Do not include R2 tags (pair-mates from SVBam) in the assembly\n"
"      --ignore-skip-centromere         Set this flag to assemble in centromeric/telomeric regions. Default is to skip.\n"
"      --throw-discordant               Throw away discordant reads when they don't cluster to 3 or more. Useful for bad libraries.\n"
"      --nm-limit                       Skip reads that have more than nm mismatches (NM tag). Default: 3\n"
"      --assemble-only-disc             Assemble only regions where there is discordant read clusters\n"
"      --assemble-normal                Do assembly on normal as well, not just read mapping\n"
"      --min-clip                       Minimum number of clipped reads to included in assembly\n"
"      --min-read-cov                   Minimum number of reads required to build a valid contig\n"
"  Assembly params\n"
"      --write-asqg                     Output an ASQG graph file for each 5000bp window. Default: false\n"
"  -e, --error-overlap                  What is the error tolerance on read overlaps. Default: 0.05\n"
"\n";

static struct timespec start;
//static faidx_t * findex;

//char *bwa_pg;

// handle a ctrl C
void my_handler(int s){

  printf("Caught signal %d. Closing BAMs and ofstreams\n",s);
  r2c_writer->Close();
  disc_writer->Close();

  all_align_stream->close();
  os_allbps->close(); 
  contigs_sam->close();

  exit(1); 

}


bool runTaiga(int argc, char** argv) {

  // start the interrupt handler
  struct sigaction sigIntHandler;
  sigIntHandler.sa_handler = my_handler;
  sigemptyset(&sigIntHandler.sa_mask);
  sigIntHandler.sa_flags = 0;

  sigaction(SIGINT, &sigIntHandler, NULL);

  if (argc != 1000) // 1000 is a dummy to make re-call of runTaiga easier for contig bam
    parseTaigaOptions(argc, argv);
 
  // learn the parameters
  if (opt::learn_params)
    learnParameters();

  // parse the region file, count number of jobs
  GenomicRegionVector file_regions, regions_torun;
  int num_jobs = countJobs(file_regions, regions_torun); 

  // override the number of threads if need
  opt::numThreads = min(num_jobs, static_cast<int>(opt::numThreads));

  if (opt::verbose > 0 && !opt::disc_cluster_only) {
    string disc_msg = opt::skip_disc ? "ON" : "OFF";
    string supp_msg = opt::skip_supp ? "ON" : "OFF";
    string r2_msg   = opt::skip_r2   ? "ON" : "OFF";
    string chuck_msg   = opt::throw_disc   ? "ON" : "OFF";
    string learn_msg   = opt::learn_params   ? "ON" : "OFF";
    string only_msg   = opt::assemble_only_disc   ? "ON" : "OFF";
    string normal_msg   = opt::normal_assemble   ? "ON" : "OFF";
    string cent_msg     = opt::ignore_skip_cent ? "ON" : "OFF";
    string contig_msg     = opt::contig_bam ? "ON" : "OFF";
    cout << "Num threads:      " << opt::numThreads << endl;
    cout << "Min overlap:      " << opt::minOverlap << endl;
    cout << "Length cutoff:    " << opt::cutoff     << endl;
    cout << "Output directory: " << opt::outdir     << endl;
    //    cout << "Chunk size:       " << opt::chunk      << endl;
    cout << "Overlap error:    " << opt::error_rate << endl;
    cout << "ContigBAM assembly?:   " << contig_msg << endl;
    cout << "Contig write limit:   " << opt::contig_write_limit << endl;

    cout << "Read filters:\n"; 
    cout << "   Min read cov on contig:    " << opt::min_read_cov << endl;
    cout << "   Min clip length:           " << opt::min_clip << endl;
    cout << "   Min MAPQ:                  " << opt::mapq  << endl;
    cout << "   Max NM:                    " << opt::nmlim  << endl;
    cout << "   Insert size:               " << opt::isize << endl;
    cout << "   Min quality score:         " << opt::qualthresh << endl;
    cout << "   Skip supplementary reads:  " << supp_msg << endl;
    cout << "   Assemble only disc clusts: " << only_msg << endl;
    cout << "   Parameter learning:        " << learn_msg << endl;
    cout << "   Assemble normal also:      " << normal_msg << endl; 
    cout << "   Ignore skip centromere:    " << cent_msg << endl; 

    if (opt::isize < 0)
      cout << "   RUNNING WITH NO DISCORDANT READ CONSIDERATION. R2 LOOKUP ON CLIPPED READS ONLY" << endl;
    else {
      cout << "   Skip discordant lookup:    " << disc_msg << endl;
      cout << "   Skip R2 tag lookup:        " << r2_msg << endl;
      cout << "   Throw non-clust discs:     " << chuck_msg << endl;
    }
  } else if (opt::verbose > 0) {
    cout << "Num threads:      " << opt::numThreads << endl;
    cout << "Output directory: " << opt::outdir     << endl;
    cout << "RUNNING DISCORDANT CLUSTERING ONLY" << endl;
  }
 
  //get the header from the tumor bam file
  SamHeader sam;
  string samheader = SVBamReader::getSamHeader(opt::tbam, sam);

  // get the reference data
  RefVector ref;  
  SVBamReader::getRefVector(opt::tbam, ref);

  r2c_writer = new BamWriter();
  //  bam_writer = new BamWriter();
  disc_writer =new BamWriter();

  BamReader br;
  br.Open(opt::tbam);
  string sam_header = br.GetHeaderText();

  
  //debug
  // TODO remove the SO:coordinate because will be wrong.

  // set the r2c writer
  if (!r2c_writer->Open("r2c.bam", sam_header, ref)) {
    cerr << "Could not open the r2c.bam for reading." << endl;
    exit(EXIT_FAILURE);
  }
  
  

  // set the bam writer
  //  if (!bam_writer->Open(opt::outdir + "/contigs.bam", "", br.GetReferenceData()/*samheader, ref*/)) {
  // cerr << "Could not open the contigs.bam for reading." << endl;
  //  exit(EXIT_FAILURE);
  // }


  //debug
  //cout << samheader << endl;
  //bam_writer->Close();
  //exit(1);

  // set the discordant writer
  if (!disc_writer->Open(opt::outdir + "/discordant.bam", samheader, ref)) {
    cerr << "Could not open the discordant.bam for writing." << endl;
   exit(EXIT_FAILURE);
  }

  // make the folder for ascii plots
  //std::string dir = "mkdir -p " + opt::outdir + "/alignments";
  //std::system(dir.c_str());
  
  // setup the ASCII plots
  //ofstream som_align_stream(opt::outdir + "/alignments/alignments.somatic.txt", ios::out);
  //ofstream ger_align_stream(opt::outdir + "/alignments/alignments.germline.txt", ios::out);
  all_align_stream = new ofstream(opt::outdir + "alignments.all.txt", ios::out);
  os_allbps = new ofstream(opt::outdir + "allbps.txt", ios::out);
  contigs_sam = new ofstream(opt::outdir + "contigs.sam", ios::out);
  

  (*contigs_sam) << sam_header;

  if (opt::verbose > 0)
    cout << "Loading the reference BWT index file for: "  << opt::refgenome << endl;
  idx = bwa_idx_load(opt::refgenome.c_str(), BWA_IDX_ALL);
  if (idx == NULL) {
    std::cerr << "Could not load the reference: " << opt::refgenome << endl;
    exit(EXIT_FAILURE);
  }

  // Create the queue and consumer (worker) threads
  wqueue<TaigaWorkItem*>  queue;
  vector<ConsumerThread<TaigaWorkItem>*> threadqueue;
  for (unsigned i = 0; i < opt::numThreads; i++) {
    ConsumerThread<TaigaWorkItem>* threadr = new ConsumerThread<TaigaWorkItem>(queue, opt::verbose > 0);
    threadr->start();
    threadqueue.push_back(threadr);
  }

  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);

  // print the jobs
  if (opt::verbose > 0)
    for (GenomicRegionVector::const_iterator it = file_regions.begin(); it != file_regions.end(); it++)
      cout << "Input Regions: " << it->toStringOffset() << endl;
  
  // print regions to run
  if (opt::verbose > 1)
    for (GenomicRegionVector::const_iterator it = regions_torun.begin(); it != regions_torun.end(); it++)
      cout << "Per Thread Regions: " << it->toStringOffset() << endl;

  // global contig vector and reads vector
  //cont_total = new AlignedContigVec();
  //reads_all_encoded  = new EncodedBAVector();
  //disc_map = new DMap();

  if (pthread_mutex_init(&snow_lock, NULL) != 0) {
      printf("\n mutex init failed\n");
      return false;
  }

  // initialize the reader for converting chromosomes
  treader_for_convert = new BamReader();
  treader_for_convert->Open(opt::tbam);

  // send the jobs
  unsigned count = 0;
  for (GenomicRegionVector::const_iterator it = regions_torun.begin(); it != regions_torun.end(); it++) {
    count++;
    AlignedContigVec  * cont_out = new AlignedContigVec();
    DMap * disc_out = new DMap();
    TaigaWorkItem * item     = new TaigaWorkItem(it->chr, it->pos1, it->pos2, count, cont_out, disc_out);
    queue.add(item);
  }

  // wait for the threads to finish
  for (unsigned i = 0; i < opt::numThreads; i++) 
    threadqueue[i]->join();

  // display the run time
  struct timespec finish;
  clock_gettime(CLOCK_MONOTONIC, &finish);
  double elapsed = (finish.tv_sec - start.tv_sec);
  int t = clock()/CLOCKS_PER_SEC;
  int min = (int)floor(elapsed / 60.0);
  int sec = (int)(elapsed-min*60);
  char buffer[100];
  sprintf (buffer, "CPU time: %dm%ds    Wall time: %dm%ds", 
            (int)floor( ((double)t) /60.0), t % 60, min, sec);
  printf ("%s\n",buffer);

  return true;

}

//TODO turn to const
int runAll(BamAlignmentVector &bavd, string name, AlignedContigVec * cont_out, DMap * disc_out)
{
  // parse name to get genomic region
  GenomicRegion gr(name);
  // check that this doesn't overlap a blacklist regions
  bool black = gr.blacklistOverlap() != 0 && bavd.size() > 2000;
  if (black) {
    if (opt::verbose) 
      cout << "Blacklist region hit at: " << gr.toStringOffset() << " with num reads: " << bavd.size() << endl;
    return 0;
  }

  // grab discordant read sequences, figure out how to deal with it
  if (!opt::contig_bam)
    handleDiscordant(bavd, name, gr,disc_out);
  if ( bavd.size() == 0 )
    return 0;

  // de-duplicate the reads by readname and position
  BamAlignmentVector bav, bav_name, bavd_tum;
  SVBamReader::deduplicateReads(bavd, bav_name);
  SVBamReader::deduplicateReadsPos(bav_name, bav);

  // remove the normals from the piece sent to assembler
  if (!opt::normal_assemble && !opt::contig_bam) {
    string tmp;
    for (BamAlignmentVector::const_iterator it = bav.begin(); it != bav.end(); it++) {
      it->GetTag("JW", tmp);
      if (tmp.at(0) == 't')
	bavd_tum.push_back(*it);
    }
  }

  // make the reads tables
  ReadTable pRT;
  if (opt::contig_bam || opt::normal_assemble)
    //BamAlignmentVectorToReadTable(bav, pRT);
      pRT = ReadTable(bav);
  else 
    //BamAlignmentVectorToReadTable(bavd_tum, pRT);    
      pRT = ReadTable(bavd_tum);

  ContigVector contigs;
  // do the assembly
  ContigVector contigs1;
  doAssembly(&pRT, name, contigs1, 0);

  ReadTable pRTc;
  if (!opt::contig_bam) {
    pRTc = ReadTable(contigs1);
    //ContigsToReadTable(contigs, pRT);
    doAssembly(&pRTc, name, contigs, 1);
  }

  //filter contigs with BWA-MEM
  BWAWrapper wrap;
  BWAReadVec bwarv;
  SamRecordVec samv;

  if (!opt::contig_bam) {
    for (auto it = contigs.begin(); it != contigs.end(); it++) 
      bwarv.push_back(BWARead(it->getID(), it->getSeq()));
    wrap.addSequences(bwarv, idx, samv);

  }

  // make aligned contigs
  for (auto it = samv.begin(); it != samv.end(); it++) {
    AlignedContig ac(it->record, treader_for_convert);
    if (ac.m_align.size() > 1 && ac.getMinMapq() >= 10 && ac.getMaxMapq() >= 40) { 
      ac.alignReadsToContigs(bav);
      //ac.sortReads();
      ac.splitCoverage();
      ac.getBreakPairs();
      if (ac.maxSplit() > 1)
	cont_out->push_back(ac);
    }
  }

  return 0;

}

bool grabReads(int refID, int pos1, int pos2, AlignedContigVec * cont_out, DMap * bav_disc) {
 
  //  SeqRecordVector tsrv, nsrv;
  BamAlignmentVector tbav, nbav;

  ////// TUMOR
  std::clock_t startr;
  startr = std::clock();

  SVBamReader t_reader(opt::tbam, "t", opt::isize, opt::mapq, opt::qualthresh, opt::minOverlap, opt::skip_supp, opt::skip_r2, 
		       opt::min_clip, opt::verbose);
  t_reader.disc_cluster_only = opt::disc_cluster_only;
  t_reader.setNMLimit(opt::nmlim);
  // find the bai
  if (!t_reader.findBamIndex())
    cerr << "Failed to open BAM index in Tumor" << endl;
  // set the BAM region
  if (!t_reader.setBamRegion(refID, pos1, pos2))
    cerr << "Failed to set BAM position in Tumor" << endl;
  // add the reads to the BamAlignmentVector
  if (!opt::contig_bam) {
    if (!t_reader.bamToBAVec(tbav)) {
      cerr << "Failed to get BAM reads for Tumor" << endl;
      exit(EXIT_FAILURE);
      return false;
    }
  } else {
    if (!t_reader.contigBamToBAVec(tbav)) {
      cerr << "Failed to get BAM reads for Tumor Contig BAM" << endl;
      exit(EXIT_FAILURE);
      return false;
    }
  }
  int num_t_reads = tbav.size();
  double svbam_time_tumor = (std::clock() - startr) / (double)(CLOCKS_PER_SEC / 1000);

  ////// NORMAL
  // open the reads
  int num_n_reads = 0;
  startr = std::clock();
  if (opt::nbam.length() > 0) {
    SVBamReader n_reader(opt::nbam, "n", opt::isize, opt::mapq, opt::qualthresh, opt::minOverlap, opt::skip_supp, opt::skip_r2, 
			 opt::min_clip, opt::verbose);
    n_reader.disc_cluster_only = opt::disc_cluster_only;
    n_reader.setNMLimit(opt::nmlim);
    // find the bai
    n_reader.findBamIndex();
    // set the BAM region
    if (!n_reader.setBamRegion(refID, pos1, pos2))
      cerr << "Failed to set BAM position in Normal" << endl;
    // add the reads to the SeqRecordVector
    n_reader.bamToBAVec(nbav);
    num_n_reads = nbav.size();
  }
  double svbam_time_normal = (std::clock() - startr) / (double)(CLOCKS_PER_SEC / 1000);
  
  // return if no reads found
  if (num_t_reads == 0)
    return true;

  // verbose
  if (opt::verbose > 1) {
    string print1 = SnowUtils::AddCommas<int>(pos1);
    string print2 = SnowUtils::AddCommas<int>(pos2);
    printf("Running- chr%2d:%11s-%11s Tum: %5d Norm: %5d\n",refID+1,print1.c_str(),print2.c_str(),num_t_reads,num_n_reads);
  }

  // Divide up the reads and send them to the assembler
  startr = std::clock();
  // TODO make safe in case normal has differnt refids than tumor
  chunkReadsForAssembly(refID, pos1, LITTLECHUNK, WINDOW_PAD, cont_out, &tbav, &nbav, bav_disc);
  double assembly_time = (std::clock() - startr) / (double)(CLOCKS_PER_SEC / 1000);

  // create and open the .tmp.fa file for contigs
  startr = std::clock();

  // push reads and sort them
  EncodedBAVector this_reads;
  for (AlignedContigVec::iterator it = cont_out->begin(); it != cont_out->end(); it++) {
    //    BamAlignmentVector tmpvec = it->getBamAlignments();
    for (BamAlignmentVector::iterator jt = it->m_bamreads.begin(); jt != it->m_bamreads.end(); jt++) {
      //jt->AddTag("CN", "Z", it->getID());
      EncodedBA eba;
      eba.a = (*jt);
      eba.a.QueryBases = "";
      //eba.a.Qualities = "";
      eba.a.RemoveTag("TS");
      eba.s = DNAEncodedString(jt->QueryBases);
      this_reads.push_back(eba);
    }
  }
  // sort by read name
  ReadMap reads_map; 
  combineR2C(this_reads, reads_map);

  // combine discordant reads with breaks
  combineContigsWithDiscordantClusters(bav_disc, cont_out);

  // get the breakpoints
  BPVec bp_glob;
  for (auto it = cont_out->begin(); it != cont_out->end(); it++) 
    bp_glob.push_back(it->getGlobalBreak());

  // add discordant reads
  addDiscordantPairsBreakpoints(bp_glob, bav_disc);

  //set as somatic or not
  for (auto it = bp_glob.begin(); it != bp_glob.end(); it++) {
    if (it->isGoodSomatic(30, 3, 0))
      it->isSomatic = true;
    else if (it->isGoodGermline(30, 3))
      it->isGermline = true;
  }
    
  ////////////////////////////////////
  // MUTEX LOCKED
  ////////////////////////////////////
  // write to the global contig out
  pthread_mutex_lock(&snow_lock);  

  // send all bps to file
  stringstream out_allbps; out_allbps << opt::outdir << "/allbps.txt";
  for (BPVec::const_iterator it = bp_glob.begin(); it != bp_glob.end(); it++) 
    (*os_allbps) << it->toString() << endl;

  // send all alignments to ASCII plot
  for (auto it = cont_out->begin(); it != cont_out->end(); it++) 
    (*all_align_stream) << it->printAlignments() << endl;

  // add these reads to final reads structure
  if (!opt::contig_bam && !opt::no_r2c)
    for (ReadMap::iterator it = reads_map.begin(); it != reads_map.end(); it++) {
      //write the BAM immediately
      //cout << "writing read "  << endl;
      it->second.a.QueryBases = it->second.s.toString();


      string cn;
      int32_t al;
      int32_t sw;
      it->second.a.GetTag("CN", cn);
      it->second.a.GetTag("AL", al);
      it->second.a.GetTag("SW", sw);
      stringstream cig;
      for (auto& cc : it->second.a.CigarData) {
	cig << cc.Length << cc.Type;
      }

      r2c_writer->SaveAlignment(it->second.a);
      read_counter++;
    }

  /*
  // write the contigs
  /*  for (auto it = cont_out->begin(); it != cont_out->end(); it++) {
    for (auto jt = it->m_align.begin(); jt != it->m_align.end(); jt++) {
      //      cout << jt->align.RefID + ":" << jt->align.Position << endl;
      if (!bam_writer->SaveAlignment(jt->align))
	cerr << "Bam alignment not saved" << endl;
    }
  }
  */

  //debug
  for (AlignedContigVec::const_iterator jj = cont_out->begin(); jj != cont_out->end(); jj++) {
    contig_counter++;
    for (auto& a : jj->m_align) {
      stringstream cig;
      for (auto& cc : a.align.CigarData) {
	cig << cc.Length << cc.Type;
      }

      assert(a.align.QueryBases.length() == a.align.Qualities.length());
      
      (*contigs_sam) << a.align.Name << "\t" << a.align.AlignmentFlag << "\t" << treader_for_convert->GetReferenceData()[a.align.RefID].RefName << 
	"\t" << a.align.Position << "\t" << a.align.MapQuality << "\t" << cig.str() << "\t" << 
	"*" << "\t" << "0" << "\t" << 0 << "\t" << a.align.QueryBases << "\t" << a.align.Qualities << endl;
    }
  }

  // write the discordant reads
  for (auto it = bav_disc->begin(); it != bav_disc->end(); it++) {
    for (auto jt = it->second.reads.begin(); jt != it->second.reads.end(); jt++) {
      discordant_counter++;
     if (!disc_writer->SaveAlignment(*jt))
  	cerr << "Discordant alignment not saved" << endl;
   }
  }

  pthread_mutex_unlock(&snow_lock);
  /////////////////////////////////////
  // MUTEX UNLOCKED
  /////////////////////////////////////
  
  // delete the contig 
  delete cont_out;
  delete bav_disc;

  // display the run time
  double write_time = (std::clock() - startr) / (double)(CLOCKS_PER_SEC / 1000);

  if (opt::verbose > 0) {
    struct timespec finish;
    clock_gettime(CLOCK_MONOTONIC, &finish);
    double elapsed = (finish.tv_sec - start.tv_sec);
    int t = clock()/CLOCKS_PER_SEC;
    int min = (int)floor(elapsed / 60.0);
    int sec = (int)(elapsed-min*60);
    char buffer[140];
    double total_time = svbam_time_normal + svbam_time_tumor + assembly_time + write_time;
    if (total_time == 0)
      total_time = 1;
    int tumor_time  = static_cast<int>(floor(svbam_time_tumor /total_time * 100.0));
    int normal_time = static_cast<int>(floor(svbam_time_normal/total_time * 100.0));
    int assem_time  = static_cast<int>(floor(assembly_time    /total_time * 100.0));
    //    int w_time      = static_cast<int>(floor(write_time       /total_time * 100.0));
    string print1 = SnowUtils::AddCommas<int>(pos1);
    string print2 = SnowUtils::AddCommas<int>(pos2);
    sprintf (buffer, "Ran chr%2d:%11s-%11s T: %5d(%2d%%) N: %5d(%2d%%) ASBL: %2d%% #Contigs: %6d #Reads: %6d #Disc: %6d-- CPU: %dm%ds Wall: %dm%ds", 
	     refID+1,print1.c_str(),print2.c_str(),num_t_reads,tumor_time,
	     num_n_reads, normal_time, 
	     assem_time,
	     contig_counter, read_counter, discordant_counter,
	     //static_cast<int>(reads_all_encoded->size()),
	     static_cast<int>(floor( static_cast<double>(t) /60.0)), t % 60,
	     min, sec);
    printf ("%s\n",buffer);
  }
 
 return true;
}

void parseTaigaOptions(int argc, char** argv) {
  bool die = false;

  if (argc < 2) 
    die = true;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
      case 'p': arg >> opt::numThreads; break;
      case 'h': die = true; break;
      case OPT_ASQG: opt::writeASQG = true; break;
      case 'm': arg >> opt::minOverlap; break;
      case 'l': arg >> opt::cutoff; break;
      case 't': arg >> opt::tbam; break;
      case 'n': arg >> opt::nbam; break;
      case 'v': arg >> opt::verbose; break;
      case 'i': arg >> opt::isize; break;
      case 's': arg >> opt::sleepDelay; break;
      case 'r': arg >> opt::regionFile; break;
      case 'o': arg >> opt::outdir; break;
      case 'w': arg >> opt::mapq; break;
	//      case 'c': arg >> opt::chunk; break;
      case 'q': arg >> opt::qualthresh; break;
      case 'y': arg >> opt::skip; break;
      case 'b': arg >> opt::numBubbleRounds; break;
      case 'x': arg >> opt::divergence; break;
      case 'e': arg >> opt::error_rate; break;
      case 'g': arg >> opt::refgenome; break;
      case OPT_MIN_READ_COV: arg >> opt::min_read_cov; break;
      case OPT_MIN_CLIP: arg >> opt::min_clip; break;
      case OPT_NODISC: opt::skip_disc = true; break;
      case OPT_CONTIG_BAM: opt::contig_bam = true; break;
      case OPT_ASSEMBLE_ONLY_DISC: opt::assemble_only_disc = true; break;
      case OPT_NOSUPP: opt::skip_supp = true; break;
      case OPT_DISC_CLUSTER_ONLY: opt::disc_cluster_only = true; break;
      case OPT_NOR2: opt::skip_r2 = true; break;
      case OPT_NO_R2CBAM: opt::no_r2c = true; break;
      case OPT_IGNORE_SKIP_CENT: opt::ignore_skip_cent = true; break;
      case OPT_THROW_DISC: opt::throw_disc = true; break;
      case OPT_MEMORY: arg >> opt::memory; break;
      case OPT_LEARN: opt::learn_params = true; break;
      case OPT_ASSEMBLE_NORMAL: opt::normal_assemble = true; break;
      case OPT_MEMORY_GOAL: arg >> opt::memory_goal; break;
      case OPT_NM: arg >> opt::nmlim; break;
      case OPT_DEBUG: opt::debug = true; break;
      case OPT_CONTIG_LIMIT: arg >> opt::contig_write_limit; break;
    }
  }

  // clean the outdir
  //opt::outdir = SnowUtils::getDirPath(opt::outdir);

  // check file validity
  if (opt::regionFile != "x") {
    //if (!existTest(opt::regionFile)) {
    if (false) {
      cerr << "Region file does not exist: " << opt::regionFile << endl;
      exit(EXIT_FAILURE);
    }
  }
    
  if (opt::numThreads <= 0) {
      cout << "run: invalid number of threads: " << opt::numThreads << "\n";
      die = true;
  }

  if (opt::cutoff < 0)
      opt::cutoff = opt::minOverlap * 1.5;

  if (opt::tbam.length() == 0) {
    cout << "run: must supply a tumor bam"<< "\n";
    die = true;
  }

  if (die) 
    {
      cout << "\n" << RUN_USAGE_MESSAGE;
      exit(1);
    }
}

// Divide up the reads and send them to the assembler
void chunkReadsForAssembly(const int refID, const int pos1, const int chunk, const int pad, 
			   AlignedContigVec * cont_out, BamAlignmentVector * tbav, BamAlignmentVector *nbav,
			   DMap * disc_out) {

  IterPairMap tmap, nmap;
  // return vectors of iterators that points to positions for chunks
  if (tbav->size() > 0)
    getChunkReads(tbav, refID, pos1, chunk, pad, tmap);
  if (nbav->size() > 0)
    getChunkReads(nbav, refID, pos1, chunk, pad, nmap);
  // grab all the normal and tumor keys
  StringVec keys;
  for (IterPairMap::const_iterator it = tmap.begin(); it != tmap.end(); it++)
    keys.push_back(it->first);
  for (IterPairMap::const_iterator it = nmap.begin(); it != nmap.end(); it++)
    keys.push_back(it->first);
  // find the unique chunks
  sort( keys.begin(), keys.end() );
  keys.erase( unique( keys.begin(), keys.end() ), keys.end());

  for (StringVec::const_iterator it = keys.begin(); it < keys.end(); it++) {
     BamAlignmentVector fvec;
     if (tmap.count( (*it) ) > 0)
       fvec.insert(fvec.end(), tmap[(*it)].start, tmap[(*it)].end); 
     if (nmap.count( (*it) ) > 0)
       fvec.insert(fvec.end(), nmap[(*it)].start, nmap[(*it)].end);

     bool fvec_pass = fvec.size() < opt::skip && fvec.size() >= 3;
     if ( fvec_pass || (opt::contig_bam && fvec.size() > 0)) 
       runAll(fvec, *it, cont_out, disc_out);
     else if (opt::verbose > 1)
       cout << "Too many/few reads at: " << *it << " with " << fvec.size() << " reads\n";
  }
  
}

// 
void getChunkReads(const BamAlignmentVector * bav, const unsigned refID, const unsigned pos1, const unsigned chunk, const unsigned pad, 
                   IterPairMap &mmap) {

  // subset the SeqRecordVector into small chunks for assembly
  BamAlignmentVector::const_iterator myk            = bav->begin();
  BamAlignmentVector::const_iterator myk_curr_start = bav->begin();
  BamAlignmentVector::const_iterator myk_next_start = bav->begin();

  int nextstart = pos1 + chunk - pad;
  int currend   = pos1 + chunk;       
  bool start_trigger = false; // only update read pointer for the first read that crosses nextstart
  while (myk < bav->end()) {

    std::string myk_pos_str;
    myk->GetTag("HP", myk_pos_str);
    int myk_pos = std::stoi(myk_pos_str.substr(3,myk_pos_str.length()));
    //cerr << myk_pos_str << " Unsigned: " << myk_pos << " MYK " << bav->end() - myk << endl;

    // update the start of the next segment
    if (myk_pos >= nextstart && !start_trigger) {
      start_trigger = true;
      myk_next_start = myk;
    }

    if (myk_pos > currend) {  // went pass the boundary, update next boundary

      // count the number of reads
      int num_reads = myk-1 - myk_curr_start; // number of reads in completed intervals

      // add to output if there are 2+ reads
      if (num_reads >= 2) {
	stringstream sstm;
	int cstart = max(nextstart - (int)chunk, (int)pos1);
	sstm << "c_" << refID + 1 << ":" << cstart << "-" << currend << "_";
        IterPair p(myk_curr_start, myk-1);
        mmap.insert(pair<string, IterPair>(sstm.str(), p));
      } 

      // setup the next interval  
      while (myk_pos > currend) { // keep updating interval until it hits current read
	currend   += chunk;
        nextstart = nextstart + chunk;
      }
     
      // update the counters to start them at the next interval
      myk_curr_start = myk_next_start; // start at the next interval
      myk = myk_curr_start;            // set the counter back to the start for this range
      start_trigger = false;

    } else {  // didn't make it past current end
      myk++;
    }
  }

  // send the last one, since it didn't make it past the while
  int num_reads = myk-1 - myk_curr_start; // number of reads in completed intervals  
  if (num_reads >= 2) {
    stringstream sstm;
    int cstart = max(nextstart - (int)chunk, (int)pos1);
    sstm << "c_" << refID + 1 << ":" << cstart << "-" << currend << "_";
    IterPair p(myk_curr_start, myk-1); //myk-1 is the same as bav.end()
    mmap.insert(pair<string, IterPair>(sstm.str(), p));
  } 

}

// just get a count of how many jobs to run. Useful for limiting threads
// also set the regions
int countJobs(GenomicRegionVector &file_regions, GenomicRegionVector &run_regions) {

  // open the region file if it exists
  bool rgfile = opt::regionFile.compare("x") != 0;
  if (rgfile) 
    file_regions = GenomicRegion::regionFileToGRV(opt::regionFile, 0);
  else if (!opt::ignore_skip_cent)
    file_regions = GenomicRegion::non_centromeres; // from GenomicRegion.cpp
  else
    file_regions = GenomicRegion::getWholeGenome();

  int threadchunk = opt::chunk;
  unsigned jj = 0; 
  int startr, endr;;
  int kk = 0;
  
  //set amount to modulate
  int thispad = 1000;
  if (rgfile)
    thispad = 0;

  // loop through each region
  bool stoploop = false;
  while (jj < file_regions.size()) {

    // if regions are greater than chunk, breakup
    if ( (file_regions[jj].pos2 - file_regions[jj].pos1) > threadchunk) {
      startr = file_regions[jj].pos1-thispad; 
      endr = startr + threadchunk;

      do {
        GenomicRegion grr(file_regions[jj].chr, startr, endr);
		run_regions.push_back(grr);

		if (endr == file_regions[jj].pos2)
	  		stoploop = true;

		kk++;
        endr   = min(file_regions[jj].pos2, (kk+1)*threadchunk + file_regions[jj].pos1 + thispad);
        startr = min(file_regions[jj].pos2,  kk*threadchunk + file_regions[jj].pos1);

      } while (!stoploop);

    } else { // region size is below chunk
		run_regions.push_back(file_regions[jj]);
    }
    jj++;
    kk = 0;
    stoploop = false;
  } // end big while

  return run_regions.size();

}

// grab the pairmate reads
void grabPairmateReads(BamAlignmentVector &bav, const GenomicRegion window, DMap * bav_disc) {

  // get count of discordant reads
  size_t disc_countr = 0;
  for (BamAlignmentVector::const_iterator it = bav.begin(); it != bav.end(); it++) {
    bool good_chr = it->MateRefID < 25 && it->RefID < 25 && it->IsMapped() && it->IsMateMapped();
    if ( (abs(it->InsertSize) > opt::isize || (it->MateRefID != it->RefID)) && good_chr) 
      disc_countr++;
  }

  // need at least three disc reads to proceed
  if (disc_countr < 2) 
    return;

  // sort the reads by discordant position
  sort( bav.begin(), bav.end(), ByMatePosition()); 

  // get the regions
  RMap rmap;
  GenomicRegionVector disc_regions;
  clusterReads(bav, disc_regions, rmap, '+', '+');
  clusterReads(bav, disc_regions, rmap, '-', '-');
  clusterReads(bav, disc_regions, rmap, '+', '-');
  clusterReads(bav, disc_regions, rmap, '-', '+');

  if (opt::verbose > 2)
    for (GenomicRegionVector::const_iterator it = disc_regions.begin(); it != disc_regions.end(); it++)
      cout << "Before trim: " << *it << " TCount: " << it->tcount << " NCount: " << it->ncount << " Cluster " << it->cluster << endl;

  // remove ones with too few reads
  GenomicRegionVector tmp_vec;
  for (GenomicRegionVector::const_iterator it = disc_regions.begin(); it != disc_regions.end(); it++) {
    //bool cent    = it->centromereOverlap() == 0;
    bool window_overlap = it->getOverlap(window) == 2;
    //bool black   = it->blacklistOverlap() == 0;
    double som_ratio = static_cast<double>(it->tcount) / (static_cast<double>(it->ncount) + 0.01);

    bool som_pass = som_ratio > 10;
    bool clust_pass = it->tcount >= MIN_DISC_PER_CLUSTER;
    bool look_lim_pass = it->tcount <= DISC_LOOKUP_LIMIT;
    bool chr_pass = it->chr < 24;
    if ( som_pass && clust_pass && look_lim_pass && chr_pass && !window_overlap)
      tmp_vec.push_back(*it);
  }
  disc_regions = tmp_vec;

  // debug
  if (opt::verbose > 2)
    for (GenomicRegionVector::const_iterator it = disc_regions.begin(); it != disc_regions.end(); it++)
      cout << "After trim: " << *it << " TCount: " << it->tcount << " NCount: " << it->ncount << " from window " << window << endl;

  // if there are too many regions, forget it
  if (disc_regions.size() >= MAX_DISC_CLUSTERS)
    return;

  // open the tumor BAM
  SVBamReader t_disc_reader(opt::tbam, "td", opt::isize, opt::mapq, opt::qualthresh, opt::minOverlap, 
			    opt::min_clip, opt::verbose);    
  t_disc_reader.disc_cluster_only = opt::disc_cluster_only;
  t_disc_reader.setReadLimit(DISC_READ_LIMIT);
  t_disc_reader.setNMLimit(opt::nmlim);
  t_disc_reader.setSkipR2(true);
  if (!t_disc_reader.findBamIndex())
    cerr << "Failed to open BAM index in Tumor" << endl;

  // open the normal BAM if provided
  bool normal_here = opt::nbam.length() > 0;
  SVBamReader n_disc_reader(opt::nbam, "nd", opt::isize, opt::mapq, opt::qualthresh, opt::minOverlap, 
			    opt::min_clip, opt::verbose);
  n_disc_reader.disc_cluster_only = opt::disc_cluster_only;
  n_disc_reader.setNMLimit(opt::nmlim);
  n_disc_reader.setSkipR2(true);
  if (normal_here) {
    if (!n_disc_reader.findBamIndex())
      cerr << "Failed to open BAM index in Normal" << endl;
    n_disc_reader.setReadLimit(DISC_READ_LIMIT);
   }

  BamAlignmentVector tmp_disc;
  for (GenomicRegionVector::iterator it = disc_regions.begin(); it != disc_regions.end(); it++) {

    it->pad(DISC_LOOKUP_PAD);
    BamAlignmentVector tbav_disc, nbav_disc;

    if (!t_disc_reader.setBamRegion(*it))
      cerr << "Failed to set BAM position in Tumor for Disc Read on position: " << it->toString() << endl;
    if (normal_here && !n_disc_reader.setBamRegion(*it))
      cerr << "Failed to set BAM position in Normal for Disc Read" << it->toString() << endl;

    // read the tumor
    if (!t_disc_reader.bamToBAVec(tbav_disc)) 
      cerr << "Tumor bam disc reader failed on region " << *it << " in window " << window << endl;

    // read the normal
    if (normal_here && tbav_disc.size() > 0) // only run if tumor found something
      if (!n_disc_reader.bamToBAVec(nbav_disc)) 
	cerr << "Normal bam disc reader failed on region " << *it << " in window " << window << endl;

    // concatenate discorant region reads with original reads
    bav.insert(bav.end(), tbav_disc.begin(), tbav_disc.end());
    if (normal_here)
      bav.insert(bav.end(), nbav_disc.begin(), nbav_disc.end());

  }

  // add to disc vector. bav contains all original and discordant reads. There may be dups, but OK 
  // since I dedupe for assembly with SVBamReader::deduplicateReads, and dedupe for discordant.bam
  // in cleanDisc()
  for (BamAlignmentVector::iterator kt = bav.begin(); kt != bav.end(); kt++) {

      unordered_map<string, string>::const_iterator ff = rmap.find(kt->Name);    

      if (ff != rmap.end() && !kt->HasTag("IR")) {

	BamAlignment ba = *kt;
        ba.AddTag("DC", "Z", ff->second);
	
	// reduce JW to make BAM smaller
	string jw;
	if (!ba.GetTag("JW", jw))
	  cerr << "Missing JW tag"  << endl;
	ba.EditTag("JW", "Z", jw.substr(0, 2));

	ba.RemoveTag("TS");
	ba.RemoveTag("HP");
	ba.RemoveTag("J2");
	//ba.RemoveTag("RP");
	ba.RemoveTag("PG");
	ba.RemoveTag("AM");

	DMap::iterator gg = bav_disc->find(ff->second);
	if (gg == bav_disc->end()) { // make a new one
	  DiscordantCluster dcc(ff->second);
	  dcc.reads.push_back(ba);
	  (*bav_disc)[ff->second] = DiscordantCluster(ff->second);
	} else {
	  gg->second.reads.push_back(ba);
	}
      }
  }

  return;

}

// call the assembler
void doAssembly(ReadTable *pRT, string name, ContigVector &contigs, int pass) {

  if (pRT->getCount() == 0)
    return;

  // forward
  SuffixArray* pSAf = new SuffixArray(pRT, 1, false); //1 is num threads. false is isReverse
  RLBWT *pBWT= new RLBWT(pSAf, pRT);

  // reverse
  pRT->reverseAll();
  SuffixArray * pSAr = new SuffixArray(pRT, 1, true); //1 is num threads. false is isReverse
  RLBWT *pRBWT = new RLBWT(pSAr, pRT);
  pRT->reverseAll();

  pSAf->writeIndex();
  pSAr->writeIndex();

  double errorRate;
  int min_overlap = opt::minOverlap;
  if (pass > 0) 
    min_overlap = 60;

  //errorRate = 0;
  //if (pass > 0)
  errorRate = (double)opt::error_rate;
  //if (pass > 0)
  //  errorRate = 0.05;
  int cutoff = opt::cutoff;
  int seedLength = min_overlap;
  int seedStride = seedLength;
  bool bIrreducibleOnly = true; // default

  SnowmanOverlapAlgorithm* pOverlapper = new SnowmanOverlapAlgorithm(pBWT, pRBWT, 
                                                       errorRate, seedLength, 
                                                       seedStride, bIrreducibleOnly);
  pOverlapper->setExactModeOverlap(false);
  pOverlapper->setExactModeIrreducible(false);

  stringstream hits_stream;
  stringstream asqg_stream;

  SnowmanASQG::HeaderRecord headerRecord;
  headerRecord.setOverlapTag(min_overlap);
  headerRecord.setErrorRateTag(errorRate);
  headerRecord.setInputFileTag("");
  headerRecord.setContainmentTag(true); // containments are always present
  headerRecord.setTransitiveTag(!bIrreducibleOnly);
  headerRecord.write(asqg_stream);    

  pRT->setZero();

  size_t workid = 0;
  SeqItem si;
  while (pRT->getRead(si)) {
    SeqRecord read;
    read.id = si.id;
    read.seq = si.seq;
    OverlapBlockList obl;
    OverlapResult rr = pOverlapper->overlapReadInexact(read, min_overlap, &obl);

    pOverlapper->writeOverlapBlocks(hits_stream, workid, rr.isSubstring, &obl);

    SnowmanASQG::VertexRecord record(read.id, read.seq.toString());
    record.setSubstringTag(rr.isSubstring);
    record.write(asqg_stream);

    workid++;
  }
  string line;
  bool bIsSelfCompare = true;
  ReadInfoTable* pQueryRIT = new ReadInfoTable(pRT);

  while(getline(hits_stream, line)) {
    size_t readIdx;
    size_t totalEntries;
    bool isSubstring; 
    OverlapVector ov;
    OverlapCommon::parseHitsString(line, pQueryRIT, pQueryRIT, pSAf, pSAr, bIsSelfCompare, readIdx, totalEntries, ov, isSubstring);
    for(OverlapVector::iterator iter = ov.begin(); iter != ov.end(); ++iter)
    {
       SnowmanASQG::EdgeRecord edgeRecord(*iter);
       edgeRecord.write(asqg_stream);
    }

  }

  // optionally output the graph structure
  if (opt::writeASQG) {

    // write ASQG to file for visualization
    stringstream asqgfile;
    asqgfile << opt::outdir << "/" << name << "pass_" << pass << ".asqg";
    ofstream ofile(asqgfile.str(), ios::out);
    ofile << asqg_stream.str();
    ofile.close();

    // write the hits stream file
    stringstream hitsfile;
    hitsfile << opt::outdir << "/" << name << "pass_" << pass << ".hits";
    ofstream ofile_hits(hitsfile.str(), ios::out);
    ofile_hits << hits_stream.str();
    ofile_hits.close();

   }

    // Get the number of strings in the BWT, this is used to pre-allocated the read table
    delete pOverlapper;
    delete pBWT; 
    delete pRBWT;
    delete pSAf;
    delete pSAr;

    //
    // PERFORM THE ASSMEBLY
    //
 
    assemble(asqg_stream, min_overlap, opt::maxEdges, opt::bExact, 
	     opt::trimLengthThreshold, opt::bPerformTR, opt::bValidate, opt::numTrimRounds, 
	     opt::resolveSmallRepeatLen, opt::numBubbleRounds, opt::gap_divergence, 
	     opt::divergence, opt::maxIndelLength, cutoff, name, contigs);

    delete pQueryRIT;
    asqg_stream.str("");
    hits_stream.str("");

    // print out some results
    if (opt::verbose > 2) {
      if (contigs.size() >= 1) {
	cout << "PASS: " << pass << " Contig Count: " << contigs.size() << " at " << name << endl;
	for (ContigVector::iterator i = contigs.begin(); i != contigs.end(); i++) 
	  cout << "   " << i->getID() << " " << i->getSeq().length() << " " << i->getSeq() << std::endl;
      }
      cerr << "PASS " << pass << " Num contigs: " << contigs.size() << endl;
    }


}

// combine reads across different 
void combineR2C(EncodedBAVector &read_in, ReadMap &read_out) {

  // try with a hash table
  for (EncodedBAVector::const_iterator it = read_in.begin(); it != read_in.end(); it++) {
    string jw;
    jw = it->a.Name + to_string(it->a.AlignmentFlag);
    ReadMap::iterator ff = read_out.find(jw);
    if (ff == read_out.end()) {
      read_out.insert(pair<string, EncodedBA>(jw, *it));
    } 
  }
    
  return;
  
}

// learn the parameters
void learnParameters() {

  if (opt::verbose > 0)
    cout << "...learning parameters" << endl;

  // read in the 100,000 reads from the tumor bam
  BamReader t_reader;
  if (!t_reader.Open(opt::tbam)) {
    cerr << "Failed to open BAM: " << opt::tbam << endl;
    exit(EXIT_FAILURE);
  }
  // read in the 100,000 reads from the normal bam
  BamReader n_reader;
  bool nexist = opt::nbam.length() > 0;
  if (nexist && !n_reader.Open(opt::nbam)) {
    cerr << "Failed to open BAM: " << opt::nbam << endl;
    exit(EXIT_FAILURE);
  }

  // get tumor stuff
  vector<double> tmapq;
  vector<double> tisize;
  double tinter = 0;
  GenomicRegion tgr;
  int treadlen = 0;
  double tclipcov = 0;
  _learn_params(t_reader, tmapq, tisize, tinter, tclipcov, tgr, treadlen);

  // get normal stuff
  vector<double> nmapq = {60, 60};
  vector<double> nisize;
  double ninter = 0;
  GenomicRegion ngr;
  int nreadlen = 0;
  double nclipcov = 0;
  cout << "Normal Exist: " << nexist << endl;
  if (nexist)
    _learn_params(n_reader, nmapq, nisize, ninter, nclipcov, ngr, nreadlen);

  // print tumor
  char tbuffer[100];
  sprintf(tbuffer, "   Tumor--  MAPQ: %2.0f(%2.1f); InsertSize: %4.0f(%4.1f); Interchrom coverage: %2.2f;", // Clip coverage: %2.2f; ReadLength: %3dbp", 
	  tmapq[0], tmapq[1], tisize[0], tisize[1], tinter);
  cout << tbuffer << " from 100,000 reads on " << tgr.toStringOffset() << endl;

  // print normal
  if (nexist) {
    char nbuffer[100];
    sprintf(nbuffer, "   Normal-- MAPQ: %2.0f(%2.1f); InsertSize: %4.0f(%4.1f); Interchrom coverage: %2.2f;", //Clip coverage: %2.2f; ReadLength: %3dbp", 
	    nmapq[0], nmapq[1], nisize[0], nisize[1], ninter);
    cout << nbuffer << " from 100,000 reads on " << ngr.toStringOffset() << endl;
  }
  
  // set the parameters
  if (tinter > INTER_CHR_COV || ninter > INTER_CHR_COV) {
    cout << "Interchromosomal pair coverage exceeds " << INTER_CHR_COV << " threshold. Bad library?\n---Setting ignore discordant (isize -> -1)" << endl;
    opt::isize = -1;
  }
  if (nclipcov > CLIP_COV_LO || tclipcov > CLIP_COV_LO) {
    cout << "Coverage of clipped reads is high. Adjusting MAPQ threshold to 30 to compensate" << endl;
    unsigned newmap = 30;
    opt::mapq = max(newmap, opt::mapq);
  } else if (nclipcov > CLIP_COV_HI || tclipcov > CLIP_COV_HI) {
    cout << "Coverage of clipped reads is way too high. Adjusting MAPQ threshold to 60 to compensate" << endl;
    opt::mapq = 60;
  }
  if (nmapq[0] < MAPQ_THRESH || tmapq[0] < MAPQ_THRESH) {
    cout << "Mean mapping quality is very low, setting MAPQ to max(10, opt::mapq)" << endl;
    unsigned newmap = 10;
    opt::mapq = max(newmap, opt::mapq);
  }

  opt::minOverlap = static_cast<unsigned>(floor(static_cast<double>(treadlen) * 0.35));
  if (opt::isize > 0)
    opt::cutoff     = static_cast<int>     (floor(static_cast<double>(treadlen) * 1.35));
  else 
    opt::cutoff     = static_cast<int>     (floor(static_cast<double>(treadlen) * 1.15));    
}

//    
void _learn_params(BamReader &reader, vector<double> &mapq_result, vector<double> &isize_result,
		   double &inter_cov, double &clip_cov, GenomicRegion &gr, int &readlen) {
  
  // dummies, can probably delete
  mapq_result.push_back(1.0);
  mapq_result.push_back(1.0);
  isize_result.push_back(1.0);
  isize_result.push_back(1.0);
  inter_cov = 0.1;
  gr.chr = 1; gr.pos1 = 0; gr.pos2 = 1;
  readlen = 1;

  BamAlignment a;
  size_t count = 0;
  vector<uint16_t> mapq;
  vector<int32_t>  isize;
  vector<bool>     inter;

  // learning the tumor properties
  reader.GetNextAlignmentCore(a);
  int32_t start = a.Position;
  int32_t end   = a.Position;
  int32_t refstart = a.RefID;
  int32_t refend   = a.RefID;
  int32_t read_length = 0;
  int clip_count = 0;
  int disc_count = 0;
  
  // 
  while(reader.GetNextAlignmentCore(a) && count < LEARN_SAMPLE) {
    count++;
    mapq.push_back(a.MapQuality);
    isize.push_back(abs(a.InsertSize));
    //inter.push_back(a.RefID != a.MateRefID && a.MapQuality > opt::mapq);
    if (a.RefID != a.MateRefID && a.MapQuality > opt::mapq)
      disc_count++;
    read_length = max(read_length, a.Length);
    end = a.Position;
    refend = a.RefID;

    // check if this is a clipped read
    for (vector<CigarOp>::const_iterator j = a.CigarData.begin(); j != a.CigarData.end(); j++) {
      if (j->Type == 'S' && j->Length >= 10) {
	clip_count++;
	break;
      }
    }
  }
  reader.Close();
  gr.chr = refstart; gr.pos1 = start; gr.pos2 = GenomicRegion::convertPos(refend, end);
   
  // get the values
  double tmapq_sd = calc_sd<uint16_t>(mapq);
  double tisize_sd = calc_sd<int32_t>(isize);
  double tmapq_mean = accumulate(mapq.begin(), mapq.end(), 0.0) / mapq.size();
  double tisize_mean = accumulate(isize.begin(), isize.end(), 0.0) / isize.size();
  double spread = static_cast<double>(GenomicRegion::convertPos(refend, end)) - static_cast<double>(GenomicRegion::convertPos(refstart, start));
  inter_cov = disc_count * static_cast<double>(read_length) / spread;
  clip_cov  = clip_count * static_cast<double>(read_length) / spread;
  
  mapq_result[0] = tmapq_mean;
  mapq_result[1] = tmapq_sd;
  isize_result[0] = tisize_mean;
  isize_result[1] = tisize_sd;

  readlen = static_cast<int>(read_length);

  return;
}

void cleanR2C() {

  BamReader reader;

  string input_bam = opt::outdir + "/r2c_tmp.bam";
  string tmp_clean_bam = opt::outdir + "/r2c_tmp_clean.bam";
  string final_r2c_bam = opt::outdir + "/r2c_clean";

  if (!reader.Open(input_bam)) {
    cerr << "Failed to open merged tmp r2c bam: " << input_bam << endl;
    return;
  }

  // setup the writer
  BamWriter writer;

  //get the header
  SamHeader sam;
  string samheader = SVBamReader::getSamHeader(opt::tbam, sam);

  // get the reference data
  RefVector ref;  
  SVBamReader::getRefVector(opt::tbam, ref);

  if (!writer.Open(tmp_clean_bam, samheader, ref))  
    cerr << "Error initializing the BAM for: " << tmp_clean_bam << endl;

  BamAlignment a;
  unordered_map<string, bool> rm;

  int count = 0;
  int added_read_count = 0;

  size_t div = 500000;

  if (opt::verbose > 0)
    cout << "...processing reads for cleaning" << endl;
    
  // loop through and write
  while (reader.GetNextAlignment(a)) {
    
    count++;
    if (count % div == 0) {
      SnowUtils::percentCalc<int>(added_read_count, count);
      int perc  = static_cast<int>(floor((float)added_read_count / (float)count * 100.0));
      char buffer[100];
      int disp_count = count / 100000;
      sprintf (buffer, " Checking read  %2de5, keeping %7d (%2d%%)", 
	       disp_count, added_read_count, static_cast<int32_t>(perc));
      cout << buffer << "  ";
      SnowUtils::displayRuntime(start);
      cout << endl;
    }
    
    string pstring = to_string(a.Position);
    string jw = a.Name.substr(a.Name.length()-8, a.Name.length()) + to_string(a.AlignmentFlag) + pstring;
    if (rm.find(jw) == rm.end()) {
      rm.insert(pair<string, bool>(jw, true));
      writer.SaveAlignment(a);      
      added_read_count++;
    }

  }

  writer.Close();
  SnowUtils::displayRuntime(start); 
  cout << endl;

  // do the read sort and index
  string cmd = "samtools sort -m " + to_string(opt::memory) + "M " + tmp_clean_bam + " " + final_r2c_bam;
  if (opt::verbose > 0)
    cout << cmd << endl;
  system(cmd.c_str());
  cmd = "samtools index " + final_r2c_bam + ".bam";
  if (opt::verbose > 0)
    cout << cmd << endl;
  system(cmd.c_str());

  // delete tmp files
  //cmd = "rm " + opt::outdir + "/r2c_clean_qsort.bam " + opt::outdir + "/r2c_qsort.bam";
  if (!opt::debug) {
    cmd = "rm " + input_bam + " " + tmp_clean_bam;
    if (opt::verbose > 0)
      cout << cmd << endl;
    system(cmd.c_str());
  }
  SnowUtils::displayRuntime(start); 
  cout << endl;

  //
  return;

}

// given a set of reads, output structure that has either all the mate sequence,
// none of the discordant reads, or empty, depending on strictness of discordant 
// assembly options
void handleDiscordant(BamAlignmentVector &bavd, string name, GenomicRegion gr, DMap * bav_disc) {

  // grab pairmate region readsdd
  if (!opt::skip_disc && opt::isize > 0) {

    // grab the pairmate sequences from the BAM
    unsigned orig_reads = bavd.size();
    grabPairmateReads(bavd, gr, bav_disc);
    unsigned post_reads = bavd.size();
    if (post_reads != orig_reads && opt::verbose > 1)
      cout << " -- Disc add " << orig_reads << " to " << post_reads << " at " << name << endl;
    
    // throw away discordant reads if no cluster found. Useful for bad libraries
    if (post_reads == orig_reads && opt::throw_disc) {
      
      // if only do asssembly if discordant found, clear and return
      if (opt::assemble_only_disc) {
	bavd.clear();
	return; 
      }

      // otherwise, remove discordant NON-CLIPPED reads
      BamAlignmentVector tmp;
      for (BamAlignmentVector::const_iterator it = bavd.begin(); it != bavd.end(); it++) {
	bool unmap_pass = !it->IsMapped() || !it->IsMateMapped();
	bool clip_pass  = SVBamReader::getClipCount(*it) >= 5;
	if (unmap_pass || clip_pass) 
	  tmp.push_back(*it);
      }
      bavd = tmp;
      if (opt::verbose > 1)
	cout << " -- Disc throw: From " << orig_reads << " to " << bavd.size() << endl;
    }


  }


}

void cleanDiscBam() {

  BamReader reader;

  string input_bam = opt::outdir + "/discordant_tmp.bam";
  string tmp_clean_bam = opt::outdir + "/discordant_tmp_clean.bam";
  string final_disc_bam = opt::outdir + "/discordant_clean";

  if (!reader.Open(input_bam)) {
    cerr << "Failed to open merged tmp discordant bam: " << input_bam << endl;
    return;
  }

  // setup the writer
  BamWriter writer;
  SamHeader sam("none");
  // set the reference for the BAM
  RefVector ref;  
  for (int i = 0; i < 25; i++) {
    RefData rf(CHR_NAME[i], CHR_LEN[i]);
    ref.push_back(rf);      
  }
  if (!writer.Open(tmp_clean_bam, sam, ref))  
    cerr << "Error initializing the BAM for: " << tmp_clean_bam << endl;

  BamAlignment a;
  unordered_map<string, bool> rm;
  
  size_t count = 0;
  int added_read_count = 0;

  size_t div = 500000;

  if (opt::verbose > 0)
    cout << "...processing discordant reads for cleaning" << endl;
  // loop through and write

  while (reader.GetNextAlignment(a)) {
    
    count++;
    if (count % div == 0) {
      int perc  = static_cast<int>(floor((float)added_read_count / (float)count * 100.0));
      char buffer[100];
      int disp_count = count / 100000;
      sprintf (buffer, " Checking discordant read  %2de5, keeping %7d (%2d%%)", 
	       disp_count, added_read_count, static_cast<int32_t>(perc));
      cout << buffer << "  ";
      SnowUtils::displayRuntime(start);
      cout << endl;
    }
    
    string pstring = to_string(a.Position);
    string jw = a.Name.substr(a.Name.length()-8, a.Name.length()) + to_string(a.AlignmentFlag) + pstring;
    if (rm.find(jw) == rm.end()) {
      rm.insert(pair<string, bool>(jw, true));
      writer.SaveAlignment(a);      
      added_read_count++;
    }

  }

  writer.Close();
  SnowUtils::displayRuntime(start); 
  cout << endl;

  // do the read sort and index
  string cmd = "samtools sort -m " + to_string(opt::memory) + "M " + tmp_clean_bam + " " + final_disc_bam;
  if (opt::verbose > 0)
    cout << cmd << endl;
  system(cmd.c_str());
  cmd = "samtools index " + final_disc_bam + ".bam";
  if (opt::verbose > 0)
    cout << cmd << endl;
  system(cmd.c_str());

  // delete tmp files
  //cmd = "rm " + opt::outdir + "/r2c_clean_qsort.bam " + opt::outdir + "/r2c_qsort.bam";
  if (!opt::debug) {
    cmd = "rm " + input_bam + " " + tmp_clean_bam;
    if (opt::verbose > 0)
      cout << cmd << endl;
    system(cmd.c_str());
  }
  SnowUtils::displayRuntime(start); 
  cout << endl;

  //
  return;

}

// perform clustering of discordant reads, with supplied orientations
void clusterReads(BamAlignmentVector &bav, GenomicRegionVector &grv, RMap &rmap, char anchor_strand, char partner_strand) {

  if (bav.size() < 3)
    return;

  bool ancrev = anchor_strand == '-';
  bool parrev = partner_strand == '-';

  int dref = -1, dpos1 = -1, dpos2 = -1;
  int dtcount = 0, dncount = 0; 

  BamAlignmentVector::const_iterator it = bav.begin();

  // find the initial read with correct orientation
  for (; it != bav.end(); it++) {
    bool isNotR2 = !it->HasTag("IR");
    bool mapped = it->IsMapped() && it->IsMateMapped();
    bool oriented =  it->IsReverseStrand() == ancrev && it->IsMateReverseStrand() == parrev;
    bool discordant = it->InsertSize >= opt::isize || (it->RefID != it->MateRefID);
    if (mapped && oriented && discordant && isNotR2) {
      dref  = it->MateRefID;
      dpos1 = it->MatePosition;
      dpos2 = it->MatePosition;
      if (SVBamReader::IsTumorRead(*it))
	dtcount = 1;
      else
	dncount = 1;
      break;
    }
  }

  // no clusters found
  if (dref < 0)
    return;

  // start the first cluster
  grv.push_back(GenomicRegion(dref, dpos1, dpos2)); 

  // keep track of the anchor breakpoint
  GenomicRegion anc(it->RefID, it->Position, it->Position);
  anc.strand = anchor_strand;
  GenomicRegion par(dref, it->MatePosition, it->MatePosition);
  par.strand = partner_strand; 

  stringstream clusttag;
  for (it = it + 1; it != bav.end(); it++) {

    bool isNotR2 = !it->HasTag("IR");
    bool mapped = it->IsMapped() && it->IsMateMapped();
    bool oriented =  it->IsReverseStrand() == ancrev && it->IsMateReverseStrand() == parrev;
    bool discordant = it->InsertSize >= opt::isize || (it->RefID != it->MateRefID);

    if (mapped && oriented && discordant && isNotR2) {
      int diff = it->MatePosition - dpos2;
      bool different_chr = dref != it->MateRefID;

      assert(diff >= 0 || different_chr); // make sure ordering of mate position is true

      if ( (diff > DISC_CLUSTER_BUFF) || different_chr) { // set new cluster

	dref = it->MateRefID;
	
	// finish off the old cluster
	finalizeCluster(grv, rmap, dpos2, anc, par, dtcount, dncount);
	
        // start the new one
	GenomicRegion gr(it->MateRefID, it->MatePosition); //leaves pos2 undefined
	gr.strand = anchor_strand;
	gr.mapq.push_back(it->MapQuality);
	grv.push_back(gr);
	dtcount = 0;
	dncount = 0;
	anc.chr = it->RefID;
	anc.pos1 = it->Position;
	anc.pos2 = it->Position;
	par.chr = it->MateRefID;
	par.pos1 = it->MatePosition;
	par.pos2 = it->MatePosition;
      }
   
      if (SVBamReader::IsTumorRead(*it))
	dtcount++;
      else
	dncount++;
      
      // add the read name to the GenomicRegion
      grv.back().rname.insert(pair<string, size_t>(it->Name, 0));
      grv.back().mapq.push_back(it->MapQuality);

      dpos2 = it->MatePosition;
      
      // update the anchor and partner breakpoint position
      //int anc_tip_pos = it->IsReverseStrand() ? it->Position : it->GetEndPosition();
      //int par_tip_pos = it->IsMateReverseStrand() ? it->MatePosition : it->GetEndPosition();
      
      anc.pos1 = min(anc.pos1, it->Position);
      par.pos1 = min(par.pos1, it->MatePosition);
      anc.pos2 = max(anc.pos2, it->Position);
      par.pos2 = max(par.pos2, it->MatePosition);

      // discordant clusters shouldn't grow too big
      int lim = 20000;
      if (abs(anc.width()) >= lim || abs(par.width()) >= lim)
	cerr << anc << " " << par << endl;
      //assert(abs(anc.pos1 - anc.pos2) < lim);
      //assert(abs(par.pos1 - par.pos2) < lim);
    }
  }
    
  // finish the last one
  finalizeCluster(grv, rmap, dpos2, anc, par, dtcount, dncount);

  return;
}

void finalizeCluster(GenomicRegionVector &grv, RMap &rmap, int pos, GenomicRegion anc, GenomicRegion par, int dtcount, int dncount) {

  if (grv.size() == 0)
    return;

  // if cluster has only 1-2 reads, remove it
  if ( (dtcount + dncount) < 3) {
    grv.pop_back();
    return;
  }

  // back sure the mapq is sufficient
  //double mean_mapq = accumulate(grv.back().mapq.begin(), grv.back().mapq.end(), 0.0) / grv.back().mapq.size();
  //if (mean_mapq < 5) {
    //cerr << "...discarding due to low mapq " << mean_mapq << endl;
  //  grv.pop_back();
  //  return;
  // }
    
  // regions that are too big are probably false positives, skip
  if (abs(grv.back().width()) > 3000) {
    grv.pop_back();
    return;
  }

  // finish the last one
  grv.back().pos2 = pos;
  grv.back().tcount = dtcount;					       
  grv.back().ncount = dncount;					       

  // update the string
  stringstream clusttag;
  if (anc < par)
    clusttag << anc << "_" << par;
  else
    clusttag << par << "_" << anc;     
  grv.back().cluster = clusttag.str();

  // write the read map
  string ctag = clusttag.str();
  GMap nams = grv.back().rname;
  for (GMap::iterator jt = nams.begin(); jt != nams.end(); jt++)
    rmap.insert(pair<string, string>(jt->first, ctag));
  grv.back().rname.clear();

}

//
void combineContigsWithDiscordantClusters(DMap * dmap, AlignedContigVec * cont_out) {

  int padr = 400; 
  size_t count = 0;
  
  for (auto it = cont_out->begin(); it != cont_out->end(); it++) {
    count++;

    // check the global break
    GenomicRegion bp1(it->m_farbreak.refID1, it->m_farbreak.pos1, it->m_farbreak.pos1);
    bp1.pad(padr);
    GenomicRegion bp2(it->m_farbreak.refID2, it->m_farbreak.pos2, it->m_farbreak.pos2);
    bp2.pad(padr);
 
    for (DMap::iterator kt = dmap->begin(); kt != dmap->end(); kt++) {
    
      bool bp1reg1 = bp1.getOverlap(kt->second.reg1) != 0;
      bool bp2reg2 = bp2.getOverlap(kt->second.reg2) != 0;
      
      //debug
      bool bp1reg2 = bp1.getOverlap(kt->second.reg2) != 0;
      bool bp2reg1 = bp2.getOverlap(kt->second.reg1) != 0;

      bool pass = bp1reg1 && bp2reg2;

      //debug
      pass = pass || (bp2reg1 && bp1reg2);

      if (pass) {
		it->addDiscordantCluster(kt->second); // add disc cluster to contig

		cout << "FOUND AN OVERLAP ON CONTIG " << it->getContigName() << endl;
		// check that we haven't already added a cluster
		kt->second.contig = it->getContigName(); // add contig to disc cluster
		if (it->m_farbreak.dc.reg1.pos1 == 0) {
		  it->m_farbreak.dc = kt->second; // add cluster to global breakpoints
		} else if (it->m_farbreak.dc.ncount < kt->second.ncount) { // choose one with normal support
		  it->m_farbreak.dc = kt->second;
		} else if (it->m_farbreak.dc.tcount < kt->second.tcount) { // choose one with more tumor support
		  it->m_farbreak.dc = kt->second;
		}
      }
      
    }
  }
}


// transfer discordant pairs not mapped to a split break 
void addDiscordantPairsBreakpoints(BPVec &bp, DMap * dmap) {

  if (opt::verbose > 0)
    cout << "...transfering discordant clusters to breakpoints structure" << endl;

  for (DMap::const_iterator it = dmap->begin(); it != dmap->end(); it++) {
    if (it->second.contig == "") { // this discordant cluster isn't already associated with a break
      cout << "ADDING DISC"  << endl;
      BreakPoint tmpbp(it->second);
      bp.push_back(tmpbp);
    }
  }
}
