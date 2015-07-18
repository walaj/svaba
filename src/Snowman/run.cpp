#include "run.h"

#include <iostream>
#include <unordered_map>
#include <stdlib.h>
#include <cstdio>
#include <getopt.h>
#include <signal.h>
#include <memory>

#include "SnowTools/SnowUtils.h"
#include "SnowTools/SnowToolsCommon.h"
#include "SnowTools/gzstream.h"
#include "SnowTools/MiniRules.h"
#include "SnowmanAssemblerEngine.h"

#include "faidx.h"
#include "vcf.h"
#include "SnowToolsCoverage.h"
#include "SnowmanBamWalker.h"

#define LITTLECHUNK 3000 
#define WINDOW_PAD 300

#define DISC_PAD 400
#define MIN_PER_CLUST 2

#define MATES

using SnowTools::BamWalker;

#define CLOCK_COUNTER 1
#ifdef CLOCK_COUNTER
static ofstream * clock_counter;
#endif

// BamMap will store filename + type (e.g. /home/mynormalbam.bam, "n1")
typedef unordered_map<string, string> BamMap;

static int read_counter = 0;
static int contig_counter = 0;
static int discordant_counter = 0;
static pthread_mutex_t snow_lock;

static size_t interupt_counter = 0;

static unique_ptr<PON> pon;

// mask
SnowTools::GenomicRegionCollection<SnowTools::GenomicRegion> grv_mask;

// output writers
static ogzstream * all_align_stream;
static ogzstream * os_allbps; 
static ogzstream * os_cigmap; 
static ogzstream * contigs_all; 
static ofstream * contigs_sam;
static ogzstream * all_disc_stream;
static ogzstream *os_coverage;

// bwa index
static bwaidx_t * idx;

// a basic BamWalker to retreive BAM info
BamWalker bwalker;

// detected at least one contig
static bool hashits = false;

namespace opt {

  namespace assemb {
    static unsigned minOverlap = 35;
    //static unsigned numBubbleRounds = 3;
    //static float divergence = 0.05;
    //static float gap_divergence = 0.00;
    static float error_rate = 0.05; 
    static bool writeASQG = false;
    //static int maxEdges = 128;
    //static int numTrimRounds = 0; //
    //static int trimLengthThreshold = -1; // doesn't matter
    //static bool bPerformTR = false; // transitivie edge reducetion
    //static bool bValidate = false;
    //static int resolveSmallRepeatLen = -1; 
    //static int maxIndelLength = 20;
    //static bool bExact = true;
    //static string outVariantsFile = ""; // dummy
  }

  static int isize = 1000;

  static bool no_assemble_normal = false;

  static string indel_mask = ""; //"/xchip/gistic/Jeremiah/Projects/HengLiMask/um75-hs37d5.bed.gz";

  static bool output_cov = false;

  static bool no_reads = false;
  
  static int32_t readlen;

  static bool no_r2c = false;

  static bool zip = true;

  static string pon = "";

  // parameters for filtering reads
  //static string rules = "global@nbases[0,0];!nm[7,1000]!hardclip;!supplementary;!duplicate;!qcfail;phred[4,100];length[50,300]%region@WG%discordant[0,800];mapq[1,100]%mapq[1,1000];clip[5,100]%ins[1,1000]%del[1,1000]";
  static string rules = "global@nbases[0,0];!hardclip;!supplementary;!duplicate;!qcfail;phred[4,100];%region@WG%discordant[0,1000];mapq[1,1000]%mapq[1,1000];clip[5,1000]%ins[1,1000];mapq[1,100]%del[1,1000];mapq[1,1000]";

  static int chunk = 1000000;

  // runtime parameters
  static unsigned verbose = 1;
  static unsigned numThreads = 1;

  // data
  static BamMap bam;
  static string refgenome = SnowTools::REFHG19;  
  static string analysis_id = "no_id";

  //subsample
  float subsample = 1.0;

  // regions to run
  static std::string regionFile = "";

  // filters on when / how to assemble
  //static bool normal_assemble = false;
  static bool disc_cluster_only = false;

}

enum { 
  OPT_ASQG,
  OPT_DISC_CLUSTER_ONLY,
  OPT_NO_READS,
  OPT_NO_ASSEMBLE_NORMAL
};

static const char* shortopts = "hzxt:n:p:v:r:g:r:e:g:k:c:a:q:m:b:";
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
  { "reference-genome",        required_argument, NULL, 'g' },
  { "no-zip",                  no_argument, NULL, 'z' },
  { "no-read-tracking",        no_argument, NULL, OPT_NO_READS },
  { "no-assemble-normal",       no_argument, NULL, OPT_NO_ASSEMBLE_NORMAL },
  { "no-r2c-bam",              no_argument, NULL, 'x' },
  { "write-asqg",              no_argument, NULL, OPT_ASQG   },
  { "error-rate",              required_argument, NULL, 'e'},
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
"  -g, --reference-genome               Path to indexed reference genome to be used by BWA-MEM. Default is Broad hg19 (/seq/reference/...)\n"
"  -a, --id-string                      String specifying the analysis ID to be used as part of ID common.\n"
"  Required input\n"
"  -t, --tumor-bam                      Tumor BAM file\n"
"  Optional input\n"                       
"  -n, --normal-bam                     Normal BAM file\n"
"  -r, --rules                          VariantBam style rules string to determine which reads to do assembly on. See documentation for default.\n"
"  -m, --min-overlap                    Minimum read overlap, an SGA parameter. Default: 0.4* readlength\n"
"  -k, --region-file                    Set a region txt file. Format: one region per line, Ex: 1,10000000,11000000\n"
"  -q, --panel-of-normals               Panel of normals gzipped txt file generated from snowman pon\n"
"  -m, --indel-mask                     BED-file with blacklisted regions for indel calling. Default /xchip/gistic/Jeremiah/Projects/HengLiMask/um75-hs37d5.bed.gz\n"
"      --disc-cluster-only              Only run the discordant read clustering module, skip assembly. Default: off\n"
"      --no-read-tracking               Don't track supporting reads. Reduces file size.\n"
"  Assembly params\n"
"      --write-asqg                     Output an ASQG graph file for each 5000bp window. Default: false\n"
"  -e, --error-rate                     Fractional difference two reads can have to overlap. See SGA param. 0 is fast, but requires exact. Default: 0.05\n"
"  -c, --chunk-size                     Amount of genome to read in at once. High numbers have fewer I/O rounds, but more memory. Default 1000000 (1M). Suggested 50000000 (50M) or 'chr' for exomes\n"
"\n";

static struct timespec start;

static SnowTools::MiniRulesCollection * mr;

/*
// handle a ctrl C
void my_handler(int s){

  pthread_mutex_lock(&snow_lock);  

  interupt_counter++;
  if (interupt_counter > 2)
    exit(EXIT_FAILURE);

  // release the reference index 
  if (idx) {
    //bwa_idx_destroy(idx);
    //idx = NULL;
  }

  if (s > 0)
    printf("\nCaught signal %d. Closing BAMs and ofstreams\n",s);

    contigs_all->close();
  
  if (hashits) {

    all_align_stream->close();
    all_disc_stream->close();
    if (opt::output_cov)
      os_coverage->close();
    os_allbps->close(); 
    os_cigmap->close(); 
    contigs_sam->close();

#ifdef CLOCK_COUNTER
    clock_counter->close();
    delete clock_counter;
#endif

    delete all_align_stream;
    delete all_disc_stream;
    delete os_coverage;
    delete os_allbps;
    delete os_cigmap;
    delete contigs_sam;
    delete contigs_all;

    // make the VCF file
    if (opt::verbose)
      cout << "...making the VCF files" << endl;

    VCFFile snowvcf(opt::analysis_id + ".bps.txt.gz", opt::refgenome.c_str(), '\t', opt::analysis_id);

    //ogzstream out;
    //out.open("vars.vcf.gz", ios::out);
    //out << snowvcf;

    // write the indel one
    string basename = opt::analysis_id + ".broad-snowman.DATECODE.";
    snowvcf.writeIndels(basename, opt::zip);
    snowvcf.writeSVs(basename, opt::zip);


  } else {
    cout << "%%%%%%%%%%%%%%%%%%%%" << endl; 
    cout << "NO VARIANTS DETECTED" << endl;
    cout << "%%%%%%%%%%%%%%%%%%%%" << endl; 
  }

  if (s == 0) 
    cout << "******************************" << endl 
	 << "Snowman completed successfully" << endl
         << "******************************" << endl;
  else 
    cout << "Snowman stopped due to signal interrupt, but successfully wrote output files"<< endl;
  SnowTools::displayRuntime(start);
  cout << endl;
  if (s > 0)
    exit(EXIT_FAILURE); 
  else 
    exit(EXIT_SUCCESS);

}
*/

void sendThreads(GRC& regions_torun) {

  /*
  // Create the queue and consumer (worker) threads
  wqueue<SnowmanWorkItem*>  queue;
  vector<ConsumerThread<SnowmanWorkItem>*> threadqueue;
  for (unsigned i = 0; i < opt::numThreads; i++) {
    ConsumerThread<SnowmanWorkItem>* threadr = new ConsumerThread<SnowmanWorkItem>(queue, opt::verbose > 0);
    threadr->start();
    threadqueue.push_back(threadr);
  }

  // send the jobs
  size_t count = 0;
  SnowTools::GenomicRegion it;
  for (auto& i : regions_torun) {
    SnowmanWorkItem * item     = new SnowmanWorkItem(SnowTools::GenomicRegion(i.chr, i.pos1, i.pos2), ++count, idx);
    queue.add(item);
  }

  // wait for the threads to finish
  for (size_t i = 0; i < opt::numThreads; i++) 
    threadqueue[i]->join();
  
    if (idx) {
      bwa_idx_destroy(idx);
      idx = NULL;
    }
  */
}


// main function to kick-off snowman
bool runSnowman(int argc, char** argv) {

  /*
#ifdef SIG_CATCH
  // start the interrupt handle
  struct sigaction sigIntHandler;
  sigIntHandler.sa_handler = my_handler;
  sigemptyset(&sigIntHandler.sa_mask);
  sigIntHandler.sa_flags = 0;
  sigaction(SIGINT, &sigIntHandler, NULL);
#endif

  /*
  parseRunOptions(argc, argv);

  // open the tumor bam to get header info
  bwalker.OpenReadBam(opt::bam.begin()->first);

  // open all of the output files
  initializeFiles();

  exit(1);
  // learn some parameters
  learnParameters();

  // parse the region file, count number of jobs
  GRC file_regions, regions_torun;
  int num_jobs = countJobs(file_regions, regions_torun); 

  // override the number of threads if need
  opt::numThreads = min(num_jobs, static_cast<int>(opt::numThreads));

  if (opt::verbose > 0) {
    //cout << "Num threads:         " << opt::numThreads << endl;
    //cout << "Min overlap:         " << opt::assemb::minOverlap << endl;
    //cout << "Assembly error rate: " << opt::assemb::error_rate << endl;
    cout << "Indel mask BED:      " << opt::indel_mask << endl;
    cout << "Read length:         " << opt::readlen << endl;
    if (opt::no_assemble_normal)
      cout << "~~~~ Normal assembly turned off ~~~~" << endl;
    cout << "BGzipping and Tabixing output?: " << (opt::zip ? "ON" : "OFF") << endl;
    cout << "Output r2c.bam?:     " << (opt::no_r2c ? "OFF" : "ON") << endl;
    if (opt::no_reads)
      cout << "Skipping read tracking for bps and vcf files" << endl; 
  }

  // import the panel of normals if it exists
  if (SnowTools::read_access_test(opt::pon)) {
    if (opt::verbose) 
      cout << "...importing panel of normals -- " << opt::pon << endl;
    BreakPoint::readPON(opt::pon, pon);
  } else if (!SnowTools::read_access_test(opt::pon) && opt::pon.length()) {
    cerr << "!!!!!!!!!!" << endl << "!!!!!!!!!!" << endl << "!!!!!!!!!!" << endl;
    cerr << "Panel of normals file does not exist or is not readable -- " << opt::pon << endl;
    cerr << "!!!!!!!!!!" << endl << "!!!!!!!!!!" << endl << "!!!!!!!!!!" << endl;
    //exit(EXIT_FAILURE);
  }

  if (opt::disc_cluster_only) 
    cout << "RUNNING DISCORDANT CLUSTERING ONLY" << endl;

  if (opt::verbose > 0) {
    for (auto i : opt::bam)
      cout << i.first << " -- " << i.second << endl;
  }

  // read the indel mask
  if (!SnowTools::read_access_test(opt::indel_mask) && opt::indel_mask.length()) {
    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cerr << "indel mask " << opt::indel_mask << " does not exist / is not readable. Skipping indel masking."  << endl;
    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    opt::indel_mask = "";
  }

  if (opt::indel_mask.length()) 
    grv_mask.regionFileToGRV(opt::indel_mask, 0);

  // set the MiniRules
  mr = new SnowTools::MiniRulesCollection(opt::rules);
  if (opt::verbose > 1)
    cout << *mr;

  // open a mutex
  if (pthread_mutex_init(&snow_lock, NULL) != 0) {
      printf("\n mutex init failed\n");
      return false;
  }

  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);

  // print the jobs
  if (opt::verbose > 0) {
    if (SnowTools::read_access_test(opt::regionFile)) {
      for (auto& i : file_regions) 
	cout << "Input Regions: " << i << endl;
    }
    else {
      cout << "***Running across whole genome***" << endl;
    }
  }  

#ifndef NO_LOAD_IDX  
  // load the index reference genome onto the heap
  if (opt::verbose > 0)
    cout << "Loading the reference BWT index file for: "  << opt::refgenome << endl;

  //idx = unique_ptr<bwaidx_t>(bwa_idx_load(opt::refgenome.c_str(), BWA_IDX_ALL));
  idx = bwa_idx_load(opt::refgenome.c_str(), BWA_IDX_ALL);
  if (idx == NULL) {
    std::cerr << "Could not load the reference: " << opt::refgenome << endl;
    exit(EXIT_FAILURE);
  }
#endif  

  if (opt::verbose)
    cout << "Starting detection pipeline." << endl;

  sendThreads(regions_torun);
  
  // close the BAMs
  my_handler(0);

  if (opt::verbose > 0) {
    SnowTools::displayRuntime(start);
    cout << endl;
  }

  */
  return true;

}


bool grabReads(const SnowTools::GenomicRegion &region, bwaidx_t* idx) {

  /*
  // place to store the contigs that are made and their read alignments
  AlignedContigVec alc;

  // place to store the reads
  ReadMap this_r2c;

  // count the number of reads
  int num_n_reads = 0, num_t_reads = 0;
  SnowTimer st;

  std::unordered_map<std::string, SnowmanBamWalker> walkers;

  st.start();

  // loop through the input bams and get reads
  for (auto& b : opt::bam)
    {
      std::cerr << "debug starting walker " << std::endl;
      // opt::bam is <filename, type>
      walkers[b.first] = SnowmanBamWalker(b.first);
      walkers[b.first].setBamWalkerRegion(region);
      walkers[b.first].SetMiniRulesCollection(*mr);

      // read the reads
      walkers[b.first].readBam();
    }

  st.stop("r");

  // make master normal cigarmap
  CigarMap cigmap_n, cigmap_t;
  for (auto& b : opt::bam) {
    if (b.second.at(0) == 'n') {
      for (auto& c : walkers[b.first].cigmap)
	cigmap_n[c.first] += c.second;
    } else if (b.second.at(0) == 't') {
      for (auto& c : walkers[b.first].cigmap)
	cigmap_t[c.first] += c.second;
    }
  }

  // get the mates from somatic 3+ mate regions
  MateRegionVector somatic_mate_regions;
  for (auto& b : opt::bam)
    if (b.second.at(0) == 't')
      somatic_mate_regions.concat(walkers[b.first].mate_regions);
  
  // add these regions to the walker and get the reads
  for (auto& b : opt::bam)
    {
      walkers[b.first].setBamWalkerRegions(somatic_mate_regions.asGenomicRegionVector());
      walkers[b.first].get_coverage = false;
      walkers[b.first].get_mate_regions = false;
      walkers[b.first].readBam();
    }
  st.stop("m");
			
  // do the discordant read clustering. Put all reads in one vector, and dedupe doubles
  unordered_map<string,bool> dd;

  ReadVec reads_fdisc;
  for (auto& b : opt::bam) {
    for (auto& r : walkers[b.first].reads) {
      string tmp;
      r_get_Z_tag(r, "SR", tmp);
      if (dd.count(tmp) == 0) 
	reads_fdisc.push_back(r);
      dd[tmp] = true;
    }
  }

  // cluster the discordant reads
  DMap this_disc = clusterDiscordantReads(reads_fdisc, region);
  st.stop("cl");

  if (opt::verbose) {
    string print1 = SnowTools::AddCommas<int>(region.pos1);
    string print2 = SnowTools::AddCommas<int>(region.pos2);
    char buffer[140];
    sprintf (buffer, "...Got reads from chr%2s:%11s-%11s | T: %5d N: %5d", 
	     //treader_for_convert->GetReferenceData()[region.chr].RefName.c_str(),
	     "1",
	     print1.c_str(),print2.c_str(),
	     num_t_reads, num_n_reads);
    std::cout << string(buffer) << std::endl;
  }
  
  // set the regions to run
  GRC grv_small(LITTLECHUNK, WINDOW_PAD, region);

  // loop through the regions and do the assembly
  for (auto& g : grv_small) {
    
    ReadVec bav_join;
    for (auto& b : walkers)
      {
	for (auto& r : b.second.reads)
	  {
	    SnowTools::GenomicRegion read(r_id(r), r_pos(r), r_pos(r));
	    SnowTools::GenomicRegion mate(r_mid(r), r_mpos(r), r_mpos(r));
	    if (g.getOverlap(read) || g.getOverlap(mate))
	      bav_join.push_back(r);
	  }
      }

    // set the contig uid
    std::string name = "c_" + std::to_string(g.chr+1) + "_" + std::to_string(g.pos1) + "_" + std::to_string(g.pos2);

    // do the assemblies
    SnowmanAssemblerEngine engine(name, opt::assemb::error_rate, opt::assemb::minOverlap, opt::readlen);
    engine.fillReadTable(bav_join);
    engine.performAssembly();
    st.stop("as");
      
    // peform alignment of contigs to reference with BWA-MEM
    BWAWrapper wrap;
    BWAReadVec bwarv;
    SamRecordVec samv;
    for (auto& i : engine.getContigs()) 
      bwarv.push_back(BWARead(i.getID(), i.getSeq()));
    wrap.addSequences(bwarv, idx, samv);
    
    // make aligned contigs from SAM records
    size_t new_alc_start = alc.size();
    for (auto& r : samv) 
      alc.push_back(AlignedContig(r.record, /*treader_for_convert,*/ g, cigmap_n, cigmap_t));
    
      // dedupe the contigs
      for (size_t i = new_alc_start; i < alc.size(); i++) 
	for (size_t j = new_alc_start; j < alc.size(); j++) 
	  if (alc[i].hasVariant() && alc[j].hasVariant() && !alc[i].m_skip && !alc[j].m_skip && i != j)
	    if (alc[i].isWorse(alc[j])) 
	      alc[i].m_skip = true;
      

      // mark to remove contigs where the variants are all in the indel mask
      if (opt::indel_mask.length()) 
	for (size_t i = new_alc_start; i < alc.size(); i++)
	  alc[i].blacklist(grv_mask);

      st.stop("bw");     

      for (size_t i = new_alc_start; i < alc.size(); i++) {
	if (!alc[i].m_skip && alc[i].hasVariant() && alc[i].getMinMapq() >= 10 && alc[i].getMaxMapq() >= 40) { 
	  
	  alc[i].alignReadsToContigs(bav_join);
	  alc[i].splitCoverage();
	}
      }
      st.stop("sw");
      
  } // end the assembly regions loop

  // combine discordant reads with breaks
  combineContigsWithDiscordantClusters(this_disc, alc);
  st.stop("cl");

  // get the breakpoints
  //BPVec bp_glob;
  vector<BreakPoint> bp_glob;
  for (auto& i : alc) {
    if (i.m_bamreads.size() && !i.m_skip) { //m_bamreads.size() is zero for contigs with ...
      vector<BreakPoint> allbreaks = i.getAllBreakPoints();
      bp_glob.insert(bp_glob.end(), allbreaks.begin(), allbreaks.end());
    }
  }

  // add discordant reads
  addDiscordantPairsBreakpoints(bp_glob, this_disc);
    
  // store the reads in a map of read2contig alignment
  // TODO fancy combine
  for (auto& it : alc) 
    if (!it.m_skip)
      for (auto& r : it.m_bamreads) 
	this_r2c[to_string(r_flag(r)) + r_qname(r)] = r;

  // add in the reads from the discordant
  for (auto& i : this_disc) { 
    if (i.second.reads_mapq >= 10 && i.second.mates_mapq >= 10) {
      for (auto& r : i.second.reads)
	this_r2c[to_string(r_flag(r.second)) + r_qname(r.second)] = r.second;
      for (auto& r : i.second.mates)
	this_r2c[to_string(r_flag(r.second)) + r_qname(r.second)] = r.second;
    }
  }



  if (!hashits)
  for (auto& i : alc)
    if (i.m_bamreads.size() && !i.m_skip) {
      hashits = true;
      break;
    }

  // check the panel of normals
  for (auto& i : bp_glob)
    if (i.hasMinimal())
      i.checkPon(pon);

  // get the allelic fractions
  SnowToolsCoverage * tcov = NULL;
  SnowToolsCoverage * ncov = NULL;
  for (auto& b : opt::bam) {
    if (b.second.at(0) == 't')
      tcov = &walkers[b.first].cov;
    if (b.second.at(0) == 'n')
      ncov = &walkers[b.first].cov;
  }
  for (auto& i : bp_glob)
    if (i.hasMinimal() && i.isindel)
      i.addAllelicFraction(tcov, ncov);

  ////////////////////////////////////
  // MUTEX LOCKED
  ////////////////////////////////////
  // write to the global contig out
  pthread_mutex_lock(&snow_lock);  

  read_counter += this_r2c.size();
  discordant_counter += this_disc.size();
  
  // write the r2c
  //if (!opt::no_r2c)
  //  writeR2C(this_r2c);

  // send all bps to file
  for (auto& i : bp_glob) {
    if (i.hasMinimal()) {
      (*os_allbps) << i.toFileString(opt::no_reads) << endl;
      if (i.confidence == "PASS" && i.ncigar == 0 && i.nsplit == 0 && i.dc.ncount == 0)
	cout << i.toPrintString() << endl;
    }
  }

  // send the cigmaps to file
  for (auto& i : cigmap_n) 
    (*os_cigmap) << i.first << "\t" << i.second << "\tN" << endl;
  for (auto& i : cigmap_t) 
    (*os_cigmap) << i.first << "\t" << i.second << "\tT" << endl;


  // send all alignments to ASCII plot
  for (auto& i : alc) 
    if (i.m_bamreads.size() && !i.m_skip && i.hasVariant())
      (*all_align_stream) << i << endl;

#ifdef MATES
  // send all discordant clusters to txt
  for (auto &i : this_disc)
    if (i.second.reads_mapq >= 10 && i.second.mates_mapq >= 10)
      (*all_disc_stream) << i.second.toFileString(!opt::no_reads) << endl;
#endif

  // send all the variant contigs to a sam file
  for (auto& i : alc) {
    if (!i.m_skip) {
      if (i.m_bamreads.size()) {
	contig_counter++;
	(*contigs_sam) << i.samrecord;
      }
      (*contigs_all) << i.samrecord;
    }
  }

  // write the coverage file
  if (opt::output_cov)
    for (auto& b : opt::bam) 
      if (b.second.at(0) == 't')
	(*os_coverage) << walkers[b.first].cov;
  st.stop("wr");

  // display the run time
  if (opt::verbose > 0 && (num_t_reads + num_n_reads) > 0) {
    string print1 = SnowTools::AddCommas<int>(region.pos1);
    string print2 = SnowTools::AddCommas<int>(region.pos2);
    char buffer[140];
    sprintf (buffer, "Ran chr%2s:%11s-%11s | T: %5d N: %5d C: %5d R: %5d D: %5d | ", 
	     /*treader_for_convert->GetReferenceData()[region.chr].RefName.c_str()*/"dum",print1.c_str(),print2.c_str(),
	     num_t_reads, num_n_reads, 
	     contig_counter, read_counter, discordant_counter);

    cout << string(buffer) << st << " | ";
    SnowTools::displayRuntime(start);
    //int perc = SnowTools::percentCalc<int>(CHR_LEN[refID] + pos2, 3100000000);
    //cout << " % of WG done: " << perc << "%" << endl;
    cout << endl;
  }


  pthread_mutex_unlock(&snow_lock);
  /////////////////////////////////////
  // MUTEX UNLOCKED
  /////////////////////////////////////
*/ 
 return true;
}

// parse the command line options
void parseRunOptions(int argc, char** argv) {
  bool die = false;

  if (argc <= 2) 
    die = true;

  string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
      case 'p': arg >> opt::numThreads; break;
      case 'a': arg >> opt::analysis_id; break;
      case 'm': arg >> opt::indel_mask; break;
      case 'q': arg >> opt::pon; break;
      case 'z': opt::zip = false; break;
      case 'h': die = true; break;
      case 'x': opt::no_r2c = true; break;
      case 'c': 
	tmp = "";
	arg >> tmp;
	if (tmp.find("chr") != string::npos) {
	  opt::chunk = 250000000; break;
	} else {
	  opt::chunk = stoi(tmp); break;
	}
    case OPT_ASQG: opt::assemb::writeASQG = true; break;
    case OPT_NO_ASSEMBLE_NORMAL: opt::no_assemble_normal = true; break;
    case OPT_NO_READS: opt::no_reads = true; break;
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
      case 'g': arg >> opt::refgenome; break;
      case 'r': arg >> opt::rules; break;
      case OPT_DISC_CLUSTER_ONLY: opt::disc_cluster_only = true; break;
      default: die= true; 
    }
  }

  // check that we input something
  if (opt::bam.size() == 0) {
    cerr << "Must add a bam file " << endl;
    exit(EXIT_FAILURE);
  }

  // check file validity
  if (opt::regionFile.length() && !SnowTools::read_access_test(opt::regionFile)) {
    cerr << "Region file does not exist or is not readable: " << opt::regionFile << endl;
    exit(EXIT_FAILURE);
  }
    
  if (opt::numThreads <= 0) {
      cout << "run: invalid number of threads: " << opt::numThreads << "\n";
      die = true;
  }

  //  if (opt::tbam.length() == 0) {
  // cout << "run: must supply a tumor bam"<< "\n";
  //  die = true;
  //}

  if (die) 
    {
      cout << "\n" << RUN_USAGE_MESSAGE;
      exit(EXIT_FAILURE);
    }
}

// just get a count of how many jobs to run. Useful for limiting threads. Also set the regions
int countJobs(GRC &file_regions, GRC &run_regions) {
  /*
  // open the region file if it exists
  bool rgfile = SnowTools::read_access_test(opt::regionFile);
  if (rgfile)
    file_regions.regionFileToGRV(opt::regionFile, 0);
  else {
    //#ifdef HAVE_HTSLIB
    size_t dumcount = 0;
    for (int i = 0; i < bwalker.header()->n_targets; i++) {
      if (++dumcount <= 25)
    	file_regions.add(SnowTools::GenomicRegion(bam_name2id(bwalker.header(), bwalker.header()->target_name[i]), 1, bwalker.header()->target_len[i]));
    }
  }
  
  if (file_regions.size() == 0) {
    cerr << "ERROR: Cannot read region file: " << opt::regionFile << " or something wrong with tumor bam header" << endl;
    exit(EXIT_FAILURE);
  }

  //int threadchunk = opt::chunk;
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

    int fr_pos1 = file_regions.at(jj).pos1;
    int fr_pos2 = file_regions.at(jj).pos2;
    int fr_chr = file_regions.at(jj).chr;

    // if regions are greater than chunk, breakup
    if ( (fr_pos2 - fr_pos1) > opt::chunk) {
      startr = max(1,fr_pos1-thispad);
      endr = startr + opt::chunk;

      do {
	SnowTools::GenomicRegion grr(fr_chr, startr, endr);
	run_regions.add(grr);
	
	if (endr == fr_pos2)
	  stoploop = true;
	
	kk++;
	endr   = min(fr_pos2, (kk+1)*opt::chunk + fr_pos1 + thispad);
	startr = min(fr_pos2,  kk*opt::chunk + fr_pos1);
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


/**
 * Merge breaks found by assembly with breaks from discordant clusters
 *
 * Currently, for each AlignedContig, loops through all of the discordant
 * clusters and checks if it overlaps the event. If so, add the cluster
 * to the AlignedContig, and mark the DiscordantCluster as being associated
 * with the AlignedContig.
 *
 * @param dm Map of discordant clusters 
 * @param contigs Vector of AlignedContigs, which store the breakpoints information from assembly
 */
void combineContigsWithDiscordantClusters(DMap &dm, AlignedContigVec &contigs) {
  /*
  int padr = 400; 
  size_t count = 0;
  
  for (auto& i : contigs) {
    if (i.m_bamreads.size() && i.m_align.size() > 1 && !i.m_skip) {
      count++;
      
      // check the global break
      SnowTools::GenomicRegion bp1 = i.m_global_bp.gr1; 
      bp1.pad(padr);
      SnowTools::GenomicRegion bp2 = i.m_global_bp.gr2; 
      bp2.pad(padr);
      
      for (auto &kt : dm) { 
	
	bool bp1reg1 = bp1.getOverlap(kt.second.reg1) != 0;
	bool bp2reg2 = bp2.getOverlap(kt.second.reg2) != 0;
	
	//debug
	bool bp1reg2 = bp1.getOverlap(kt.second.reg2) != 0;
	bool bp2reg1 = bp2.getOverlap(kt.second.reg1) != 0;
	
	// TODO have discordant reads support non global breakpoints
	
	//debug
	//cout << bp1 << " " << bp2 << " bp1reg1: " << bp1reg1 << " bp1reg2: " << bp1reg2 << 
	//	" bp2reg1: " << bp2reg1 << " bp2reg2: " << bp2reg2 << " kt.second.reg1 " << kt.second.reg1 << " kt.second.reg2 " <<  kt.second.reg2 << endl;
	
	bool pass = bp1reg1 && bp2reg2;
	
	//debug
	pass = pass || (bp2reg1 && bp1reg2);
	
	if (pass) {
	  i.addDiscordantCluster(kt.second); // add disc cluster to contig
	  
	  // check that we haven't already added a cluster
	  kt.second.contig = i.getContigName(); // add contig to disc cluster
	  if (i.m_global_bp.dc.reg1.pos1 == 0) {
	    i.m_global_bp.dc = kt.second; // add cluster to global breakpoints} else if (i.m_global_bp.dc.ncount < kt.second.ncount) { // choose one with normal support
	    i.m_global_bp.dc = kt.second;
	  } else if (i.m_global_bp.dc.tcount < kt.second.tcount) { // choose one with more tumor support
	    i.m_global_bp.dc = kt.second;
	  }
	}
	
      }
    }
  }
  */
}

// transfer discordant pairs not mapped to a split break 
  /*
void addDiscordantPairsBreakpoints(vector<BreakPoint> &bp, DMap &dmap) {

  if (opt::verbose > 1)
    cout << "...transfering discordant clusters to breakpoints structure" << endl;
  for (auto& i : dmap) { 
    if (i.second.contig == "" && (i.second.tcount+i.second.ncount) > 1) { // this discordant cluster isn't already associated with a break
      BreakPoint tmpbp(i.second);
      bp.push_back(tmpbp);
    }
  }

}
  */

/**
 * Write the read support BAM file.
 *
 * Sorts the reads that are supplied by coordinate and adds them to the open 
 * BAM file. This method does not ensure that reads are not duplicated
 * within the BAM, as it is writing to an already open bam. See cleanR2CBam.
 *
 * @param  r2c Map of the reads to write
 */
//void writeR2C(ReadMap &r2c) {

  /*
  if (r2c.size() == 0)
    return; 

  // sort by position
  ReadVec r2c_vec;
  for (auto& r : r2c)
    r2c_vec.push_back(r.second);
  sort(r2c_vec.begin(), r2c_vec.end(), ByReadPosition());

  // add the reads
  //#ifdef HAVE_BAMTOOLS
  //for (auto& r : r2c_vec) {
  // r2c_writer->SaveAlignment(*r);
  //}
  //#endif

    //#ifdef HAVE_HTSLIB
  //for (auto& r : r2c_vec) {
  //  sam_write1(r2c_fop, hdr, r.get());
  //}
  //#endif
  */
//}

// cluster the discordant reads
  /*
DMap clusterDiscordantReads(ReadVec &bav, const SnowTools::GenomicRegion &interval) {

  // remove any reads that are not present twice or have sufficient isize
  unordered_map<string, int> tmp_map;
  for (auto& i : bav)
    if (tmp_map.count(r_qname(i)) ==0)
      tmp_map[r_qname(i)] = 1;
   else
      tmp_map[r_qname(i)]++;

  ReadVec bav_dd;
  for (auto&i : bav) {
    bool disc_r = (abs(r_isize(i)) >= opt::isize) || (r_mid(i) != r_id(i));
    if (tmp_map[r_qname(i)] == 2 && disc_r)
      bav_dd.push_back(i);
  }

  // sort by position
  sort(bav_dd.begin(), bav_dd.end(), ByReadPosition());

  // clear the tmp map. Now we want to use it to store if we already clustered read
  tmp_map.clear();

  vector<ReadVec> fwd, rev, fwdfwd, revrev, fwdrev, revfwd;
  ReadVec this_fwd, this_rev;

  pair<int, int> fwd_info, rev_info; // refid, pos
  fwd_info = {-1,-1};
  rev_info = {-1,-1};

  // cluster in the READ direction
  for (auto& i : bav_dd) {
    if (r_is_pmapped(i) && tmp_map.count(r_qname(i)) == 0) {

      tmp_map[r_qname(i)] = 0;
      // forward clustering
      if (!r_is_rstrand(i)) 
	_cluster(fwd, this_fwd, i, false);
      // reverse clustering 
      else 
      	_cluster(rev, this_rev, i, false);
    }
  }
  // finish the last clusters
  if (this_fwd.size() > 0)
    fwd.push_back(this_fwd);
  if (this_rev.size() > 0)
    rev.push_back(this_rev);

  // cluster in the MATE direction for FWD facing READ
  for (auto& v : fwd) {
    fwd_info = {-1,-1};
    rev_info = {-1,-1};
    this_fwd.clear(); this_rev.clear();
    // sort by mate position to prepare for second clustering
    sort(v.begin(), v.end(), ByMatePosition());
    for (auto& i : v) {
      // forward clustering
      if (!r_is_mrstrand(i)) 
	_cluster(fwdfwd, this_fwd, i, true);
      // reverse clustering 
      else 
	_cluster(fwdrev, this_rev, i, true);
    }
    // finish the last clusters
    if (this_fwd.size() > 0)
      fwdfwd.push_back(this_fwd);
    if (this_rev.size() > 0)
      fwdrev.push_back(this_rev);
  }

  // cluster in the MATE direction for REV facing READ
  for (auto& v : rev) {
    fwd_info = {-1,-1};
    rev_info = {-1,-1};
    this_fwd.clear(); this_rev.clear();
    // sort by mate position to prepare for second clustering
    sort(v.begin(), v.end(), ByMatePosition());
    for (auto& i : v) {
      if (!r_is_mrstrand(i) )
	_cluster(revfwd, this_fwd, i, true);
      // reverse clustering 
      else 
	_cluster(revrev, this_rev, i, true);
    }
    // finish the last clusters
    if (this_fwd.size() > 0)
      revfwd.push_back(this_fwd);
    if (this_rev.size() > 0)
      revrev.push_back(this_rev);
  }    

  DMap dd;
  _convertToDiscordantCluster(dd, fwdfwd, bav_dd);
  _convertToDiscordantCluster(dd, fwdrev, bav_dd);
  _convertToDiscordantCluster(dd, revfwd, bav_dd);
  _convertToDiscordantCluster(dd, revrev, bav_dd);

  // remove clusters that dont overlap with the window
  DMap dd_clean;
  for (auto& i : dd) {
    if (i.second.reg1.getOverlap(interval) > 0 || i.second.reg2.getOverlap(interval))
      dd_clean[i.first] = i.second;
  }

  // print the clusters
  if (opt::verbose > 3) {
    for (auto& i : dd_clean)
      cout << i.second << endl;
  }
  
  return dd_clean;

}
  */

// is this a read from a tumor
bool isTumorRead(const Read &a) {
  /*
  string tmp;
  r_get_SR(a, tmp);
  return tmp.at(0) == 't';

  */
}

/**
 * Cluster reads by alignment position 
 * 
 * Checks whether a read belongs to a cluster. If so, adds it. If not, ends
 * and stores cluster, adds a new one.
 *
 * @param cvec Stores the vector of clusters, which themselves are vectors of read pointers
 * @param clust The current cluster that is being added to
 * @param a Read to add to cluster
 * @param mate Flag to specify if we should cluster on mate position instead of read position
 * @return Description of the return value
 */
  /*
bool _cluster(vector<ReadVec> &cvec, ReadVec &clust, Read &a, bool mate) {

  // get the position of the previous read. If none, we're starting a new one so make a dummy
  pair<int,int> last_info;
  if (clust.size() == 0)
    last_info = {-1, -1};
  else if (mate)
    last_info = {r_mid(clust.back()), r_mpos(clust.back())};
  else
    last_info = {r_id(clust.back()), r_pos(clust.back())};

  // get the position of the current read
  pair<int,int> this_info;
  if (mate)
    this_info = {r_mid(a), r_mpos(a)};
  else 
    this_info = {r_id(a), r_pos(a)};

  // check if this read is close enough to the last
  if ( (this_info.first == last_info.first) && (this_info.second - last_info.second) <= DISC_PAD ) {
    clust.push_back(a);
    last_info = this_info;
    return true;
  // read does not belong to cluster. close this cluster and add to cvec
  } else {

    // if enough supporting reads, add as a cluster
    if (clust.size() >= MIN_PER_CLUST) 
      cvec.push_back(clust);

    // clear this cluster and start a new one
    clust.clear();
    clust.push_back(a);

    return false;
  }
}
  */

/**
 * Remove / combine duplicate reads in read support BAM.
 *
 * To save memory, the supporting read BAM is written on-the-fly.
 * This means that the final BAM has duplicates and is unsorted.
 * This reads in the r2c.bam file and de-duplicates by combining
 * the AL, CN, SW and DC tags with an 'x' separation. Writes the 
 * r2c_clean.bam file.
 */
/*
void cleanR2CBam() {

  if (opt::verbose)
    cout << "...cleaning R2C bam" << endl;

  // open the unsorted / uncleaned bam
  BamReader br;
  if (!br.Open("r2c.bam")) {
    cerr << "Cannot open r2c.bam" << endl;
    return;
  }

  // open a BAM for writing
  BamWriter bw;
  if (!bw.Open("r2c_clean.bam", br.GetHeaderText(), br.GetReferenceData())) {
    cerr << "Cannot open r2c_clean.bam" << endl;
    return;
  }

  unordered_map<string, Read> map;
  ReadVec vec;

  BamAlignment a;
  while (br.GetNextAlignment(a)) {
    string tmp;
    a.GetTag("SR",tmp);
    if (map.count(tmp) == 0) {
      Read b(new BamAlignment(a));
      map[tmp] = b;
    } else {
      
      string cn, dc, sw, al;
      a.GetTag("CN",cn);
      a.GetTag("DC",dc);
      a.GetTag("SW",sw);
      a.GetTag("AL",al);

      SnowTools::SmartAddTag(map[tmp],"CN",cn);
      SnowTools::SmartAddTag(map[tmp],"AL",al);
      SnowTools::SmartAddTag(map[tmp],"SW",sw);
      SnowTools::SmartAddTag(map[tmp],"DC",dc);

    }
  }

  // transfer to vector and sort
  for (auto& r : map)
    vec.push_back(r.second);
  sort(vec.begin(), vec.end(), ByReadPosition());

  // write it
  for (auto& r : vec)
    bw.SaveAlignment(*r);
  bw.Close();

  // index it
  br.Close();
  if (!br.Open("r2c_clean.bam")) {
    cerr << "Cannot open cleaned r2c_clean.bam" << endl;
    return;
  }
  br.CreateIndex();

  string cmd = "rm r2c.bam;";
  system(cmd.c_str());

}
*/
// convert a cluster of reads into a new DiscordantCluster object
/*
void _convertToDiscordantCluster(DMap &dd, vector<ReadVec> cvec, ReadVec &bav) {

  for (auto& v : cvec) {
    if (v.size() > 1) {
      DiscordantCluster d(v, bav); /// slow but works
      dd[d.id] = d;
    }
  }
}
*/

/**
 * Open the output files for writing
 */
void initializeFiles() {

  /*
  std::cerr << "...initializing files" << std::endl;
  
  // setup the ASCII plots
  string n1 = opt::analysis_id + ".alignments.txt.gz";
  string n2 = opt::analysis_id + ".discordant.txt.gz";
  string n3 = opt::analysis_id + ".bps.txt.gz";
  string n4 = opt::analysis_id + ".contigs.sam";
  string n5 = opt::analysis_id + ".contigs_all.sam.gz";
  string n6 = opt::analysis_id + ".cigarmap.txt.gz";
  string n7 = opt::analysis_id + ".coverage.bed.gz";
  all_align_stream = new ogzstream(n1.c_str(), ios::out);
  if (opt::output_cov)
    os_coverage      = new ogzstream(n7.c_str(), ios::out);
  all_disc_stream  = new ogzstream(n2.c_str(), ios::out);
  os_allbps        = new ogzstream(n3.c_str(), ios::out);
  os_cigmap        = new ogzstream(n6.c_str(), ios::out);
  contigs_sam      = new ofstream(n4.c_str(), ios::out);
  contigs_all      = new ogzstream(n5.c_str(), ios::out);
  
#ifdef CLOCK_COUNTER
  clock_counter = new ofstream("clock.counter.csv");
#endif

  // write the bp header
  (*os_allbps) << BreakPoint::header() << endl;
  
  // write the discordant reads header
  (*all_disc_stream) << DiscordantCluster::header() << endl;

  // write the tumor header to the contigs SAM file
  string r2c_name = opt::analysis_id + ".r2c.bam";

  // write the header for the contigs sam files
  //(*contigs_all) << string(hdr->text);

  // write the header for the contigs sam files
  //(*contigs_sam) << string(hdr->text);


  //#endif
  */

}

//
void learnParameters() {
  /*
  std::cerr << "...learning params" << std::endl;

  //while (treader_for_convert->GetNextAlignmentCore(a) && count++ < 10000) 
  //  opt::readlen = max(a.Length, opt::readlen);
  //treader_for_convert->Rewind();
  //opt::assemb::minOverlap = 0.30 * opt::readlen;
  //opt::assemb::minOverlap = (opt::assemb::minOverlap > 50) ? 50 : opt::assemb::minOverlap;
  //#endif 
  */
}

