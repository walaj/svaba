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

static pthread_mutex_t snow_lock;
static ContigVector * cont_total;
static int out_count = 0;
static BamAlignmentVector * reads_all;
static BamAlignmentVector * disc_vec;

namespace opt
{
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
  //"      --memory-goal                    Set the chunk-size to aim for this max-memory. Not guarenteed. Default: 4096\n"
  //"      --learn                          Learn the best parameters to use (-m, -l, -c, -i)\n"
"  Required input\n"
"  -t, --tumor-bam                      Tumor BAM file\n"
"  Optional input\n"                       
"  -n, --normal-bam                     Normal BAM file\n"
"  -m, --min-overlap                    Minimum read overlap, an SGA parameter. Default: 0.4* readlength\n"
  //"  -s, --sleep-delay                    Delay (in seconds) between starting processes to avoid BAM IO stampede. Default: 0\n"
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
  //"      --write-asqg                     Output an ASQG graph file for each 5000bp window. Default: false\n"
"  -e, --error-overlap                  What is the error tolerance on read overlaps. Default: 0.05\n"
"\n";

static struct timespec start;

bool runTaiga(int argc, char** argv) {

  if (argc != 1000) // 1000 is a dummy to make re-call of runTaiga easier for contig bam
    parseTaigaOptions(argc, argv);

  //create marker file to let know started, terminate if already made
  /*std::stringstream odir;
  odir << opt::outdir << "/started.txt";
  if (existTest(odir.str()) && !opt::contig_bam) {
    std::cerr << "Snowman already started. Terminating...\n";
    return true;
  }
  std::ofstream marker(odir.str(), ios::out);
  marker << "Snowman started";
  marker.close();
  */

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
  cont_total = new ContigVector();
  reads_all  = new BamAlignmentVector();
  disc_vec   = new BamAlignmentVector();

  if (pthread_mutex_init(&snow_lock, NULL) != 0) {
      printf("\n mutex init failed\n");
      return false;
  }

  // send the jobs
  unsigned count = 0;
  for (GenomicRegionVector::const_iterator it = regions_torun.begin(); it != regions_torun.end(); it++) {
    count++;
    ContigVector  * cont_out = new ContigVector();
    TaigaWorkItem * item     = new TaigaWorkItem(it->chr, it->pos1, it->pos2, count, cont_out);
    queue.add(item);
  }

  // wait for the threads to finish
  for (unsigned i = 0; i < opt::numThreads; i++) 
    threadqueue[i]->join();

  // dump final contigs and reads
  clearMemWriteData();

  delete reads_all;
  delete cont_total;

  num_jobs = out_count;

  if (opt::verbose > 0)  
    cout << endl << "All threads done, merging temp files" << endl;

  if (num_jobs == 0) {
    cerr << "No regions run. Something went wrong..." << endl;
    exit(EXIT_FAILURE);
  }

  // merge the contig fasta files
  stringstream ss;
  if (opt::contig_bam) 
    ss << "cat " << opt::outdir << "/*contig_tmp.fa > " << opt::outdir << "/all_contigs_bootstrap.fa";
  else
    ss << "cat " << opt::outdir << "/*tmp.fa > " << opt::outdir << "/all_contigs.fa";    
  if (!opt::disc_cluster_only) {
    if (opt::verbose > 1)
      cout << ss.str() << endl;
    system(ss.str().c_str());
  }

  // remove the tmp fasta files
  stringstream rms;
  if (opt::contig_bam)
    rms << " rm -f " << opt::outdir << "/*contig_tmp.fa";
  else
    rms << " rm -f " << opt::outdir << "/*tmp.fa";    
  if (!opt::debug && !opt::disc_cluster_only)
    system(rms.str().c_str());

  // dedupe the contigs
  stringstream cmd_dedup;
  if (opt::contig_bam)
    cmd_dedup << "cd " << opt::outdir << "; " << SGAPATH << " index -t " << to_string(opt::numThreads) << " all_contigs_bootstrap.fa; " << 
      SGAPATH << " rmdup -t " << to_string(opt::numThreads) << " -e 0.01 all_contigs_bootstrap.fa; mv all_contigs_bootstrap.rmdup.fa all_contigs_bootstrap.fa; rm *bwt *sai *dups.fa";
  //else
  //  cmd_dedup << "cd " << opt::outdir << "; " << SGAPATH << " index -t " << to_string(opt::numThreads) << " all_contigs.fa; " << 
  //    SGAPATH << " rmdup -t " << to_string(opt::numThreads) << " -e 0.01 all_contigs.fa; mv all_contigs.rmdup.fa all_contigs.fa; rm *bwt *sai all_contigs.rmdup.dups.fa;";
  if (!opt::disc_cluster_only) {
    system(cmd_dedup.str().c_str());
    cout << cmd_dedup.str() << endl;
  }

  // merge the r2c BAM files
  sleep(2);
  stringstream bam_merge, disc_merge;
  if (!opt::disc_cluster_only && opt::verbose > 0)
    cout << "Number of tmp bam files to merge: " << num_jobs << endl;

  string rmstring, drmstring;
  if (opt::debug) {
    rmstring = "";
    drmstring = "";
  } else if (!opt::contig_bam && !opt::no_r2c) {
    rmstring = " rm " + opt::outdir + "/*reads2contig.bam;" ;
    drmstring= " rm " + opt::outdir + "/*dreads.bam;" ;
  } else if (!opt::contig_bam) {
    drmstring= " rm " + opt::outdir + "/*dreads.bam;" ;
  }

  if (num_jobs == 1) {
    bam_merge << "cd " << opt::outdir << "; mv r0_reads2contig.bam r2c_tmp.bam;";
    disc_merge<< "cd " << opt::outdir << "; mv r0_dreads.bam discordant_tmp.bam;";
  }
  else {
    bam_merge << "time samtools merge -f " << opt::outdir << "/r2c_tmp.bam " << opt::outdir << "/*reads2contig.bam;" << rmstring; 
    disc_merge<< "time samtools merge -f " << opt::outdir << "/discordant_tmp.bam " << opt::outdir << "/*dreads.bam;" << drmstring; 
  }

  // deal with discordant
  if (opt::verbose > 0 && !opt::contig_bam)
    cout << disc_merge.str() << endl;
  if (!opt::contig_bam)
    system(disc_merge.str().c_str());

  // deal with r2c
  bool run_r2c = !opt::contig_bam && !opt::no_r2c && !opt::disc_cluster_only;
  if (opt::verbose > 0 && run_r2c)
    cout << bam_merge.str() << endl;
  if (run_r2c)
    system(bam_merge.str().c_str());
  
  // peform the alignment of the contigs to the genome
  if (!opt::disc_cluster_only)
    runBWA();

  // clean up the r2c bam
  if (run_r2c) 
    cleanR2C();

  //clean up the discordant bam
  if (!opt::contig_bam) 
    cleanDiscBam();

  // re-run snowman on contigs
  if (!opt::contig_bam && !opt::disc_cluster_only) {
    opt::contig_bam = true;
    opt::tbam = opt::outdir + "/all_bwa.bam";
    opt::nbam = "";
    opt::regionFile = "x";
    opt::verbose = 0;
    opt::contig_write_limit = 500000;
    runTaiga(1000, argv);
  }
  
  // set the ending marker
  stringstream odir2;
  odir2 << opt::outdir << "/ended.txt";
  ofstream marker2(odir2.str(), ios::out);
  marker2 << "Snowman ended";
  marker2.close();

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
int runAll(BamAlignmentVector &bavd, string name, ContigVector * cont_out)
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
    handleDiscordant(bavd, name, gr);
  if ( bavd.size() == 0)
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
    BamAlignmentVectorToReadTable(bav, pRT);
  //    pRT = ReadTable(bav);
  else 
    BamAlignmentVectorToReadTable(bavd_tum, pRT);    
  //    pRT = ReadTable(bavd_tum);
  //pRT = opt::normal_assemble ? bav : bavd_tum;


  ContigVector contigs;
  // do the assembly
  ContigVector contigs1;
  doAssembly(&pRT, name, contigs1, 0);

  ReadTable pRTc;
  if (!opt::contig_bam) {
    //pRTc = ReadTable(contigs1);
    ContigsToReadTable(contigs, pRT);
    doAssembly(&pRTc, name, contigs, 1);
  }

  // do the matching of reads to contigs. Keep only contigs with good support
  if (!opt::contig_bam) {
    matchReads2Contigs(&contigs, bav, cont_out);
  } else {
    for (ContigVector::const_iterator it = contigs1.begin(); it != contigs1.end(); it++)
      cont_out->push_back(*it);
  }
  return 0;
}

bool grabReads(int refID, int pos1, int pos2, ContigVector * cont_out) {
 
  //  SeqRecordVector tsrv, nsrv;
  BamAlignmentVector tbav, nbav;

  ////// TUMOR
  // open the reads0
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
  unsigned num_t_reads = tbav.size();
  double svbam_time_tumor = (std::clock() - startr) / (double)(CLOCKS_PER_SEC / 1000);

  ////// NORMAL
  // open the reads
  unsigned num_n_reads = 0;
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
  chunkReadsForAssembly(refID, pos1, LITTLECHUNK, WINDOW_PAD, cont_out, &tbav, &nbav);
  double assembly_time = (std::clock() - startr) / (double)(CLOCKS_PER_SEC / 1000);

  // create and open the .tmp.fa file for contigs
  startr = std::clock();

  // push reads and sort them
  BamAlignmentVector this_reads, reads_reduced;
  for (ContigVector::iterator it = cont_out->begin(); it != cont_out->end(); it++) {
    BamAlignmentVector tmpvec = it->getBamAlignments();
    for (BamAlignmentVector::iterator jt = tmpvec.begin(); jt != tmpvec.end(); jt++) {
      //jt->AddTag("CN", "Z", it->getID());
      this_reads.push_back(*jt);
    }
  }

  // sort by read name
  ReadMap reads_map; 
  combineR2C(this_reads, reads_map);

  ////////////////////////////////////
  // MUTEX LOCKED
  ////////////////////////////////////
  // write to the global contig out
  pthread_mutex_lock(&snow_lock);  

  // add these reads to final reads structure
  if (!opt::contig_bam && !opt::no_r2c)
    for (ReadMap::const_iterator it = reads_map.begin(); it != reads_map.end(); it++)
      reads_all->push_back(it->second);

  // add contigs to final structure
  for (ContigVector::iterator it = cont_out->begin(); it != cont_out->end(); it++) {
    it->clearReads();
    cont_total->push_back(*it);
  }

  if (cont_total->size() > opt::contig_write_limit) 
    clearMemWriteData();

  pthread_mutex_unlock(&snow_lock);
  /////////////////////////////////////
  // MUTEX UNLOCKED
  /////////////////////////////////////
  
  // delete the contig 
  delete cont_out;

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
    int w_time      = static_cast<int>(floor(write_time       /total_time * 100.0));
    string print1 = SnowUtils::AddCommas<int>(pos1);
    string print2 = SnowUtils::AddCommas<int>(pos2);
    sprintf (buffer, "Ran chr%2d:%11s-%11s T: %5d N: %5d -- Tumor: %2d%% Normal: %2d%% Assembly: %2d%% Write: %2d%% -- CPU: %dm%ds Wall: %dm%ds", 
	     refID+1,print1.c_str(),print2.c_str(),num_t_reads,num_n_reads, 
	     tumor_time, normal_time, assem_time, w_time, 
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
			   ContigVector * cont_out, BamAlignmentVector * tbav, BamAlignmentVector *nbav) {

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
       runAll(fvec, *it, cont_out);
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


// Parse the region file specified by the -r flag. Must be comma delimited
bool parseRegionFile(GenomicRegionVector &gr) {

  ifstream fin;
  fin.open(opt::regionFile); // open a file
  if (!fin.good()) 
    return false; // exit if file not found

  // establish valid chr to run on
  unordered_map<string, int> VALID_CHR;
  for (unsigned i = 0; i < CHR_NAME_NUM.size(); i++)
    VALID_CHR.insert(pair<string, int>(CHR_NAME_NUM[i], 0));

  /*
  // open the optionally zipped region file in BED format
  gzifstream stream(opt::regionFile.c_str());
  if (!stream.is_open()) {
    std::cerr << "ERROR: Cannot open " + opt::regionFile << std::endl;
    exit(EXIT_FAILURE);
  }

  // parse the BED file
  BEDRow row;
  while (stream >> row) {
    cout << "BEDROW " << row.gr << endl;
  }
  */

  // read each line of the file
  while (!fin.eof()) {

    // read an entire line into memory
    char buf[MAX_CHARS_PER_LINE];
    fin.getline(buf, MAX_CHARS_PER_LINE);
  
    // parse the line into blank-delimited tokens
    int n = 0; // a for-loop index
  
    // array to store memory addresses of the tokens in buf
    const char* token[MAX_TOKENS_PER_LINE] = {}; // initialize to 0
  
    // parse the line
    token[0] = strtok(buf, DELIMITER); // first token
    if (token[0]) // zero if line is blank
      for (n = 1; n < MAX_TOKENS_PER_LINE; n++) {
        token[n] = strtok(0, DELIMITER); // subsequent tokens
  	if (!token[n]) 
	  break; // no more tokens
       }
  	
    // loop through the lines
    for (int i = 0; i < n; i += 3) {
	string mchr  = token[i];
  	string mpos1 = token[i+1];
  	string mpos2 = token[i+2];

  	// deal with x and y
  	bool isx = mchr.compare("X")==0;
  	bool isy = mchr.compare("Y")==0;
	mchr = isx ? "23" : mchr;
	mchr = isy ? "24" : mchr;

  	string error = "";
  	bool valid_chr = VALID_CHR.find(mchr) != VALID_CHR.end();
  	if (!valid_chr)
  	  error = "Invalid chromosome";

  	// create the genomic region
  	GenomicRegion this_gr;
  	try {
  	  this_gr.chr = atoi(mchr.c_str())-1;
  	  this_gr.pos1 = min(atoi(mpos1.c_str()), CHR_LEN[this_gr.chr]);
  	  this_gr.pos2 = min(atoi(mpos2.c_str()), CHR_LEN[this_gr.chr]);
  	} catch (...) {
  	  cerr << "Caught error with parsing region file for region " << mchr << ":" << mpos1 << "-" << mpos2 << endl;
  	}

  	// deteremine if this region overlaps with a centromere
  	int overlap_result = 0;
	if (!opt::ignore_skip_cent) {
	  overlap_result = this_gr.centromereOverlap();
	  if (overlap_result != 0) 
	    error = (overlap_result == 1) ? "Partial overlap with centromere" : "Full overlap with centromere";
	}
  	    
  	// keep only those regions that are valid
  	if (valid_chr && overlap_result != 2) 
  	  gr.push_back(this_gr);

  	// print out if verbose
	if (opt::verbose > 0 && error.length() > 0) 
  	  cout << "   " << this_gr.toStringOffset() << " " << error << endl;

    } // end i+=3 for
  } // end while

  return true;
  
}

// fill the worker queue with Taiga jobs tiling the whole genome
/*int sendWholeGenomeJobs(wqueue<TaigaWorkItem*>  &mqueue) {

  int threadchunk = opt::chunk;
  int jj = 1; 
  int endr = threadchunk;
  int startr = 1;
  int refr = 0;
  int num_jobs = 0;
  
  while (refr <= 23) { 
    while (endr <= CHR_LEN[refr] && startr < CHR_LEN[refr]) {
      jj++;
      num_jobs++;
      ContigVector * cont_out = new ContigVector();
      TaigaWorkItem * item = new TaigaWorkItem(opt::tbam, opt::nbam, refr, startr, endr, jj, cont_out);
      mqueue.add(item);
      endr = min(CHR_LEN[refr], (jj+1)*threadchunk);
      startr = min(CHR_LEN[refr], 1 + jj*threadchunk);
      sleep(opt::sleepDelay);
    }
    refr++;
    startr = 1;
    endr = min(threadchunk, CHR_LEN[refr]);
    jj = 1;
  }

  return num_jobs;
  }*/

// just get a count of how many jobs to run. Useful for limiting threads
// also set the regions
int countJobs(GenomicRegionVector &file_regions, GenomicRegionVector &run_regions) {

  // open the region file if it exists
  bool rgfile = opt::regionFile.compare("x") != 0;
  if (rgfile) 
    parseRegionFile(file_regions);
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

// contigs hold all the contigs for this window
// bav holds all the reads for this window
// cont_out holds all final contigs across windows
void matchReads2Contigs(ContigVector * contigs, BamAlignmentVector &bav, ContigVector  * cont_out) {

  if (opt::no_r2c) {
    for (ContigVector::iterator i = contigs->begin(); i != contigs->end(); i++) 
      cont_out->push_back(*i);
    return;
  }
  
  // MATCHING BY FIND
  int buff = 12;
  int pad = 10;

  for (ContigVector::iterator i = contigs->begin(); i != contigs->end(); i++) {
    StringMap name_map;
    for (BamAlignmentVector::iterator j = bav.begin(); j != bav.end(); j++) {

      string QB;
      
      if (!j->GetTag("TS", QB))
	QB = j->QueryBases;
      int seqlen = QB.length();
      size_t posa = i->getSeq().find(QB.substr(pad, buff)); // try first part of read
      size_t posb = i->getSeq().find(QB.substr(max(seqlen-buff-pad,0),buff)); // try second part of read
      bool hit1 = posa != string::npos;
      bool hit2 = posb != string::npos;
      string read_name;
      j->GetTag("JW", read_name);

      // PROCEED IF ALIGNS TO FORWARD
      if (hit1 || hit2) {
	int tpos = posa - pad;
	i->addRead(*j, std::max(tpos, 0), true);
	name_map.insert(pair<string, unsigned>(read_name, 0));
      }

      //OTHERWISE TRY REVERSE
      else {
	string rstring = QB;
	SnowUtils::rcomplement(rstring); 
	posa = i->getSeq().find(rstring.substr(pad,buff)); 
	posb = i->getSeq().find(rstring.substr(max(seqlen-buff-pad,0),buff)); 
	hit1 = posa != string::npos;
	hit2 = posb != string::npos;

	if (hit1 || hit2) {
	  int tpos = posa - pad;
	  //j->EditTag("TS", "Z", rstring); // edit the tag to be reverseComplemented
	  i->addRead(*j, std::max(tpos, 0), true);
	  name_map.insert(pair<string, unsigned>(read_name, 0));
	}
      }
      ////////////////////////////
    } // end read loop

    // add the alignment if there are 3+ tumor reads that made contig
    if (i->getContigTumorReadCount() >= opt::min_read_cov) { 

      // add pairmates if not already included
      for  (BamAlignmentVector::const_iterator it = bav.begin(); it != bav.end(); it++) {
	string tmp1, tmp2; 
	it->GetTag("JW", tmp1);
	it->GetTag("J2", tmp2);
	
	bool readInContig1 = name_map.find(tmp1) != name_map.end();
	bool readInContig2 = name_map.find(tmp2) != name_map.end();
	if (!readInContig1 && readInContig2) { //JW not in but J2 is, proceed
	  i->addRead(*it, -1, false);
	}
      }
      
      // clean the tags
      i->clearFinalTags();
      
      cont_out->push_back(*i); 
    }
    
#ifdef DEBUG_GRABREADS
    // debug
    if (i->getID().compare("contig_1:15001-20001_25") == 0 || true) {
      std::cout << i->getID() << std::endl;
      i->printBamAlignments();
     }
#endif

    // print verbose message
    if (opt::verbose > 2) 
	cerr << "Contig Size: " << i->getSeq().length() << " Read Support: " << i->getReadCount() << endl;

  } // end contig loop
}

// grab the pairmate reads
void grabPairmateReads(BamAlignmentVector &bav, const GenomicRegion window) {

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

  //debug 
  /*
  for (RMap::const_iterator it = rmap.begin(); it != rmap.end(); it++)
    cerr << "rmap: " << it->first << " " << it->second << endl;
  for (BamAlignmentVector::const_iterator it = bav.begin(); it != bav.end(); it++)
    if (it->Name == "D0ENMACXX111207:3:2303:3581:23759")
      cout << "BAV look: " << it->RefID+1 << ":" << it->Position << "--" << it->MateRefID+1 << ":" << it->MatePosition << " window: " << window << endl;
  for (BamAlignmentVector::const_iterator it = bav_disc_keep_fwd.begin(); it != bav_disc_keep_fwd.end(); it++)
    if (it->Name == "D0ENMACXX111207:3:2303:3581:23759")
      cout << "BAV look FWD: " << it->RefID+1 << ":" << it->Position << "--" << it->MateRefID+1 << ":" << it->MatePosition << " window: " << window << endl;
  for (BamAlignmentVector::const_iterator it = bav_disc_keep_rev.begin(); it != bav_disc_keep_rev.end(); it++)
    if (it->Name == "D0ENMACXX111207:3:2303:3581:23759")
      cout << "BAV REV look: " << it->RefID+1 << ":" << it->Position << "--" << it->MateRefID+1 << ":" << it->MatePosition << " window: " << window << endl;
  */

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

	tmp_disc.push_back(ba);
      }

  }

  // add discordant clusters to the global
  //MUTEX LOCK
  pthread_mutex_lock(&snow_lock);  
  for (BamAlignmentVector::const_iterator it = tmp_disc.begin(); it != tmp_disc.end(); it++) {
    disc_vec->push_back(*it);
  }
  pthread_mutex_unlock(&snow_lock);
  //MUTEX UNLOCKED

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
void combineR2C(BamAlignmentVector &read_in, ReadMap &read_out) {

  // try with a hash table
  for (BamAlignmentVector::const_iterator it = read_in.begin(); it != read_in.end(); it++) {
    string jw;
    jw = it->Name + to_string(it->AlignmentFlag);
    //it->GetTag("JW", jw);
    ReadMap::iterator ff = read_out.find(jw);
    if (ff == read_out.end()) {
      read_out.insert(pair<string, BamAlignment>(jw, *it));
    } else {
      /*string currcn, curral;
      ff->second.GetTag("CN", currcn);
      ff->second.GetTag("AL", curral);
      string newcn, newal;
      it->GetTag("CN", newcn);
      it->GetTag("AL", newal);
      currcn = currcn + "x" + newcn;
      curral = curral + "x" + newal;
      ff->second.EditTag("CN", "Z", currcn);
      ff->second.EditTag("AL", "Z", curral);*/
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


void runBWA() {

  string allbwa, infasta;
  if (opt::contig_bam)
    allbwa = opt::outdir + "/all_bwa_bootstrap.bam";
  else
    allbwa = opt::outdir + "/all_bwa.bam";    
  
  if (opt::contig_bam)
    infasta = opt::outdir + "/all_contigs_bootstrap.fa";
  else
    infasta = opt::outdir + "/all_contigs.fa";
  /*if (existTest(allbwa)) {
    if (opt::verbose > 0) {
      cout << "Already found all_bwa.bam, skipping BWA: " << allbwa << endl;
    }
    return;
  }*/
  
  string cmd = "bwa mem -t " + to_string(opt::numThreads) + " " + opt::refgenome + " " + infasta +
    " > " + opt::outdir + "/all_bwa_tmp.sam";
  if (opt::verbose > 0)
    cout << cmd << endl;
  system(cmd.c_str());

  // make BAM file from this
  cmd = "samtools view " + opt::outdir + "/all_bwa_tmp.sam -Sb -h > " + opt::outdir + "/all_bwa_tmp.bam";
  system(cmd.c_str());
  if (opt::verbose > 0)
    cout << cmd << endl;

  // sort the bam
  string allbwa_noext = allbwa;
  cmd = "samtools sort -m " + to_string(opt::memory) + "M " + opt::outdir + "/all_bwa_tmp.bam " + allbwa_noext.erase(allbwa_noext.length() - 4, 4);
  system(cmd.c_str());
  if (opt::verbose > 0)
    cout << cmd << endl;

  // index the bam
  cmd = "samtools index " + allbwa;
  system(cmd.c_str());
  if (opt::verbose > 0)
    cout << cmd << endl;

  // remove the SAM file and tmp bam
  cmd =  "rm " + opt::outdir + "/all_bwa_tmp.sam " + opt::outdir + "/all_bwa_tmp.bam " + infasta; 
  if (!opt::debug) {
    if (opt::verbose > 0)
      cout << cmd << endl;
    system(cmd.c_str());
  }

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
  string samheader = SVBamReader::getSamHeader(opt::tbam);

  // get the reference data
  RefVector ref;  
  SVBamReader::getRefVector(opt::tbam, ref);

  if (!writer.Open(tmp_clean_bam, samheader, ref))  
    cerr << "Error initializing the BAM for: " << tmp_clean_bam << endl;

  BamAlignment a;
  unordered_map<string, bool> rm;
  
  //BamAlignment curr_align, a;
  //reader.GetNextAlignment(curr_align);
  
  //string curr_cn, curr_al;
  //a.GetTag("CN", curr_cn);
  //a.GetTag("AL", curr_al);

  //string this_cn;
  //string this_al;

  size_t count = 0;
  int added_read_count = 0;

  size_t div = 500000;

  if (opt::verbose > 0)
    cout << "...processing reads for cleaning" << endl;
  // loop through and write

  while (reader.GetNextAlignment(a)) {
    
    count++;
    if (count % div == 0) {
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

    /*   
    if (a.AlignmentFlag == curr_align.AlignmentFlag && a.Name.compare(curr_align.Name)==0) {
      a.GetTag("CN", this_cn);
      a.GetTag("AL", this_al);
      
      curr_cn = curr_cn + "x" + this_cn;
      curr_al = curr_al + "x" + this_al;
    } else {
      
      // update the old one
      curr_align.EditTag("CN", "Z", curr_cn);
      curr_align.EditTag("AL", "Z", curr_al);

      // clean up
      if (opt::clear_tags) {
	curr_align.RemoveTag("J2");
	curr_align.RemoveTag("HP");
	curr_align.RemoveTag("Q2");
	curr_align.RemoveTag("NM");
	curr_align.RemoveTag("R2");
	curr_align.RemoveTag("AS"); 
	curr_align.RemoveTag("SA");
	curr_align.RemoveTag("MD");
	curr_align.RemoveTag("XM");
	curr_align.RemoveTag("XS");
	curr_align.RemoveTag("OQ");
	curr_align.RemoveTag("RP");
      }

      writer.SaveAlignment(curr_align);

      // make the next one
      curr_align = a;
      a.GetTag("CN", curr_cn);
      a.GetTag("AL", curr_al);

      added_read_count++;
    }
    */
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

void writeContigFasta(ContigVector *ct) {

  ofstream ostream;
  string ofile;
  if (opt::contig_bam)
    ofile = opt::outdir + "/" + "contigs_" + to_string(out_count) + ".contig_tmp.fa";
  else 
    ofile = opt::outdir + "/" + "contigs_" + to_string(out_count) + ".tmp.fa";
  ostream.open(ofile);
  
  // write the contigs to the .tmp.fa file
  for(ContigVector::const_iterator it = ct->begin(); it != ct->end(); ++it) {
    ostream << ">" << it->getID() << "\n";
    ostream << it->getSeq() << "\n";
  }
  ostream.close();
}

// grab the reads from the contig struct and write to BAM file
void writeReadsBam(BamAlignmentVector * reads) {

  // do a final combine of R2C reads before writing
  ReadMap readm;
  combineR2C(*reads, readm);

  // set the reference for the BAM
  RefVector ref;  
  for (int i = 0; i < 25; i++) {
    RefData rf(CHR_NAME[i], CHR_LEN[i]);
    ref.push_back(rf);      
  }
  
  // open the BAM for writing
  SamHeader sam("none");
  string outbam  = opt::outdir + "/" + "r" + to_string(out_count) + "_reads2contig.bam";
  if (opt::verbose > 1)
    cout << "Writing the read2contig BAM: " << outbam << endl;
  BamWriter writer;
  if (!writer.Open(outbam, sam, ref))  
    cerr << "Error initializing the BAM for: " << outbam << endl;
  
  // write the alignments to the BAM
  for (ReadMap::const_iterator it = readm.begin(); it != readm.end(); it++) 
    if (it->second.RefID < 25 && it->second.MateRefID < 25) { // dont write reads that arent on standard chr, for now
      writer.SaveAlignment(it->second);
    }
  
  writer.Close();
}

// given a set of reads, output structure that has either all the mate sequence,
// none of the discordant reads, or empty, depending on strictness of discordant 
// assembly options
void handleDiscordant(BamAlignmentVector &bavd, string name, GenomicRegion gr) {

  // grab pairmate region readsdd
  if (!opt::skip_disc && opt::isize > 0) {

    // grab the pairmate sequences from the BAM
    unsigned orig_reads = bavd.size();
    grabPairmateReads(bavd, gr);
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

// limit of contigs in memory hit, so clear data
void clearMemWriteData() {

  if (opt::verbose > 0) 
    cout << "Writing tmp fasta + bam, contig count: " << cont_total->size() << endl;
  
  // writing the contigs to a fasta
  if (cont_total->size() > 0)
    writeContigFasta(cont_total);
  
  // write the reads to a BAM
  if (reads_all->size() > 0 && !opt::contig_bam && !opt::no_r2c)
    writeReadsBam(reads_all);
  
  // write the disc_reads to a BAM
  if (disc_vec->size() > 0 && !opt::contig_bam)
    writeDiscBam(disc_vec);

  out_count++;
  cont_total->clear();
  reads_all->clear();
  disc_vec->clear();
  
}


/*void addDiscCluster(BamAlignment a1, BamAlignment a2, size_t cluster) {

  DiscordantPair dp(a1, a2, cluster);
  pthread_mutex_lock(&snow_lock);  
  disc_vec->push_back(dp);
  pthread_mutex_unlock(&snow_lock);

}*/


void writeDiscBam(BamAlignmentVector * disc) {

  // do a final combine of R2C reads before writing
  //ReadMap readm;
  //combineR2C(*reads, readm);

  // set the reference for the BAM
  RefVector ref;  
  for (int i = 0; i < 25; i++) {
    RefData rf(CHR_NAME[i], CHR_LEN[i]);
    ref.push_back(rf);      
  }
  
  // open the BAM for writing
  SamHeader sam("none");
  string outbam  = opt::outdir + "/" + "r" + to_string(out_count) + "_dreads.bam";
  if (opt::verbose > 1)
    cout << "Writing the discordant BAM: " << outbam << endl;
  BamWriter writer;
  if (!writer.Open(outbam, sam, ref))  
    cerr << "Error initializing the BAM for: " << outbam << endl;
  
  // write the alignments to the BAM
  for (BamAlignmentVector::const_iterator it = disc->begin(); it != disc->end(); it++) 
    if (it->RefID < 25 && it->MateRefID < 25) { // dont write reads that arent on standard chr, for now
      writer.SaveAlignment(*it);
    }
  
  writer.Close();
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

void finalizeCluster(GenomicRegionVector &grv, RMap &rmap, int pos, GenomicRegion anc, GenomicRegion par,
		     int dtcount, int dncount) {

  if (grv.size() == 0)
    return;

  // if cluster has only 1-2 reads, remove it
  if ( (dtcount + dncount) < 3) {
    grv.pop_back();
    return;
  }

  // back sure the mapq is sufficient
  double mean_mapq = accumulate(grv.back().mapq.begin(), grv.back().mapq.end(), 0.0) / grv.back().mapq.size();
  if (mean_mapq < 5) {
    //cerr << "...discarding due to low mapq " << mean_mapq << endl;
    grv.pop_back();
    return;
  }
    
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

/*
void SGAassemble(stringstream &asqg_stream, int minOverlap, int cutoff, string prefix, ContigVector &contigs) {

  Bigraph* pGraph = SGUtil::loadASQG(asqg_stream, minOverlap, true, opt::maxEdges);

  // this is the SGA default
  pGraph->setExactMode(true);

  // remove containments from graph
  SGContainRemoveVisitor containVisit;
  while(pGraph->hasContainment())
    pGraph->visit(containVisit);

  // Compact together unbranched chains of vertices
  pGraph->simplify();

  // do the bubble popping
  if(opt::numBubbleRounds > 0) {
    double maxBubbleDivergence = 0.05f;
    double maxBubbleGapDivergence = 0.01f;
    int maxIndelLength = 20;    

    SGSmoothingVisitor smoothingVisit(opt::outVariantsFile, maxBubbleGapDivergence, maxBubbleDivergence, maxIndelLength);
    int numSmooth = opt::numBubbleRounds;
    while(numSmooth-- > 0)
      pGraph->visit(smoothingVisit);
    pGraph->simplify();
  }
  
  pGraph->renameVertices(prefix);

  // extract the contigs
  SGVisitorContig av;
  pGraph->visit(av);

  // push onto final structure
  ContigVector tmp = av.m_ct;
  for (ContigVector::const_iterator it = tmp.begin(); it != tmp.end(); it++) {
    if (it->getLength() >= cutoff) 
      contigs.push_back(*it);
  }
  
  delete pGraph;

}

// load the ASQG from a stringstream
Bigraph* loadASQG(stringstream& pReader, const unsigned int minOverlap, bool allowContainments = false, size_t maxEdges = -1) {
  
  // Initialize graph
  Bigraph* pGraph = new Bigraph;

  int stage = 0;
  int line = 0;
  string recordLine;
  while(getline(pReader, recordLine))
    {
      ASQG::RecordType rt = ASQG::getRecordType(recordLine);
      switch(rt)
        {
	case ASQG::RT_HEADER:
	  {
	    if(stage != 0)
	      {
		std::cerr << "Error: Unexpected header record found at line " << line << "\n";
		exit(EXIT_FAILURE);
	      }
	    
	    ASQG::HeaderRecord headerRecord(recordLine);
	    const SQG::IntTag& overlapTag = headerRecord.getOverlapTag();
	    if(overlapTag.isInitialized())
	      pGraph->setMinOverlap(overlapTag.get());
	    else
	      pGraph->setMinOverlap(0);
	    
	    const SQG::FloatTag& errorRateTag = headerRecord.getErrorRateTag();
	    if(errorRateTag.isInitialized())
	      pGraph->setErrorRate(errorRateTag.get());
	    
	    const SQG::IntTag& containmentTag = headerRecord.getContainmentTag();
	    if(containmentTag.isInitialized())
	      pGraph->setContainmentFlag(containmentTag.get());
	    else
	      pGraph->setContainmentFlag(true); // conservatively assume containments are present
	    
	    const SQG::IntTag& transitiveTag = headerRecord.getTransitiveTag();
	    if(!transitiveTag.isInitialized())
	      {
		std::cerr << "Warning: ASQG does not have transitive tag\n";
		pGraph->setTransitiveFlag(true);
	      }
	    else
	      {
		pGraph->setTransitiveFlag(transitiveTag.get());
	      }
	    
	    break;
	  }
	case ASQG::RT_VERTEX:
	  {
	    // progress the stage if we are done the header
	    if(stage == 0)
	      stage = 1;
	    
	    if(stage != 1)
	      {
		std::cerr << "Error: Unexpected vertex record found at line " << line << "\n";
		exit(EXIT_FAILURE);
	      }
	    
	    ASQG::VertexRecord vertexRecord(recordLine);
	    const SQG::IntTag& ssTag = vertexRecord.getSubstringTag();
	    
	    //                Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(vertexRecord.getID(), vertexRecord.getSeq());
	    Vertex* pVertex = new Vertex(vertexRecord.getID(), vertexRecord.getSeq());
	    if(ssTag.isInitialized() && ssTag.get() == 1)
	      {
		// Vertex is a substring of some other vertex, mark it as contained
		pVertex->setContained(true);
		pGraph->setContainmentFlag(true);
	      }
	    pGraph->addVertex(pVertex);
	    break;
	  }
	case ASQG::RT_EDGE:
	  {
	    if(stage == 1)
	      stage = 2;
	    
	    if(stage != 2)
	      {
		std::cerr << "Error: Unexpected edge record found at line " << line << "\n";
		exit(EXIT_FAILURE);
	      }
	    
	    ASQG::EdgeRecord edgeRecord(recordLine);
	    const Overlap& ovr = edgeRecord.getOverlap();
	    
	    // Add the edge to the graph
	    if(ovr.match.getMinOverlapLength() >= (int)minOverlap)
	      SGAlgorithms::createEdgesFromOverlap(pGraph, ovr, allowContainments, maxEdges);
	    break;
	  }
        }
      ++line;
    }
  
  // Completely delete the edges for all nodes that were marked as super-repetitive in the graph
  SGSuperRepeatVisitor superRepeatVisitor;
  pGraph->visit(superRepeatVisitor);
  
  // Remove any duplicate edges
  SGDuplicateVisitor dupVisit;
  pGraph->visit(dupVisit);
  
  return pGraph;
  
}
*/

void ContigsToReadTable(const ContigVector &contigs, ReadTable &pRT) {
  
  //idx = 0;
  for (ContigVector::const_iterator it = contigs.begin(); it != contigs.end(); it++) {
    SeqItem si;
    si.seq = it->getSeq();
    si.id = it->getID();
    pRT.addRead(si);
  }
}

void BamAlignmentVectorToReadTable(const BamAlignmentVector &bav, ReadTable &pRT) {
  //idx = 0;
  //  m_pIndex = NULL; // not built by default
  for(BamAlignmentVector::const_iterator i = bav.begin(); i != bav.end(); i++) { 
    SeqItem si;
    i->GetTag("JW", si.id);
    //si.id = i->Name;
    //si.seq = i->QueryBases;
    std::string seqr;
    
    if (!i->GetTag("TS", seqr))
      seqr = i->QueryBases;
    si.seq = seqr;
    pRT.addRead(si); //m_table.push_back(si);
  }
}
