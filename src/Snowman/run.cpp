#include "run.h"
//#include "VariantBamReader.h"
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
#include "Util.h" // TODO
#include "faidx.h"
#include "BWAWrapper.h"

#include <signal.h>
#include <memory>

#include "EncodedString.h"
#include "AlignedContig.h"

#define LITTLECHUNK 5000 
#define BIGCHUNK 1000000
#define WINDOW_PAD 300

#define DISC_PAD 400
#define MIN_PER_CLUST 2

// DISC_LOOKUP_LIMIT is limit to how many TOTAL reads during one disc region lookup before SVBam fails
#define DISC_READ_LIMIT 2500
#define DISC_LOOKUP_LIMIT 3000
#define DISC_CLUSTER_BUFF 500
#define MIN_DISC_PER_CLUSTER 3
#define MAX_DISC_CLUSTERS 8
#define DISC_LOOKUP_PAD 400 

#define MATES

using namespace std;
using namespace BamTools;

// BamMap will store filename + type (e.g. /home/mynormalbam.bam, "n1")
typedef unordered_map<string, string> BamMap;

static int read_counter = 0;
static int contig_counter = 0;
static pthread_mutex_t snow_lock;

static DMap * dmap;

// output writers
static BamWriter * r2c_writer;
//static BamWriter * disc_writer;

// store the reads that support contigs
static ReadMap * r2c_reads;

// reader for converting reference IDs
static BamReader * treader_for_convert;

// make a map to store the BAM readers for everything
//VariantBamReaderMap * vmap;

// writers for string data
static ofstream * all_align_stream;
static ofstream * os_allbps; 
static ofstream * contigs_sam;

// bwa index
static bwaidx_t *idx;

namespace opt {

  // assembly parameters
  namespace assemb {
    static unsigned minOverlap = 35;
    static unsigned numBubbleRounds = 3;
    static float divergence = 0.05;
    static float gap_divergence = 0.05;
    static float error_rate = 0.05;
    static bool writeASQG = false;
    static int maxEdges = 128;
    static int numTrimRounds = 0; //
    static int trimLengthThreshold = -1; // doesn't matter
    static bool bPerformTR = false; // transitivie edge reducetion
    static bool bValidate = false;
    static int resolveSmallRepeatLen = -1; 
    static int maxIndelLength = 20;
    static bool bExact = true;
    static string outVariantsFile = ""; // dummy
  }

  // parameters for filtering reads
  static string rules;
  static int isize = 800;

  // runtime parameters
  static unsigned verbose = 1;
  static unsigned numThreads = 1;

  // data
  static BamMap bam;
  static string refgenome = REFHG19;  

  // regions to run
  static std::string regionFile = "x";

  // filters on when / how to assemble
  //static bool normal_assemble = false;
  static bool disc_cluster_only = false;
  
  // output filter
  //static bool no_r2c = false;

}

enum { 
  OPT_ASQG,
  OPT_LEARN,
  OPT_ASSEMBLE_NORMAL,
  OPT_DISC_CLUSTER_ONLY,
};

static const char* shortopts = "ht:n:p:m:v:r:q:g:r:e:g:k:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "tumor-bam",               required_argument, NULL, 't' },
  { "normal-bam",              required_argument, NULL, 'n' },
  { "threads (\"processes\")", required_argument, NULL, 'p' },
  //{ "min-overlap",             required_argument, NULL, 'm' },
  { "region-file",             required_argument, NULL, 'k' },
  { "rules",                   required_argument, NULL, 'r' },
  { "reference-genome",        required_argument, NULL, 'g' },
  { "write-asqg",              no_argument, NULL, OPT_ASQG   },
  //{ "assemble-only-disc",      no_argument, NULL, OPT_ASSEMBLE_ONLY_DISC   },
  //{ "no-disc-lookup",          no_argument, NULL, OPT_NODISC },
  //{ "learn",                   no_argument, NULL, OPT_LEARN   },
  //{ "assemble-normal",         no_argument, NULL, OPT_ASSEMBLE_NORMAL   },
  //{ "disc-cluster-only",       no_argument, NULL, OPT_DISC_CLUSTER_ONLY   },
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
"  Required input\n"
"  -t, --tumor-bam                      Tumor BAM file\n"
"  Optional input\n"                       
"  -n, --normal-bam                     Normal BAM file\n"
"  -m, --min-overlap                    Minimum read overlap, an SGA parameter. Default: 0.4* readlength\n"
"  -k, --region-file                    Set a region txt file. Format: one region per line, Ex: 1,10000000,11000000\n"
"      --disc-cluster-only              Only run the discordant read clustering module, skip assembly. Default: off\n"
"  Assembly params\n"
"      --write-asqg                     Output an ASQG graph file for each 5000bp window. Default: false\n"
"\n";

static struct timespec start;

static MiniRulesCollection * mr;

// handle a ctrl C
void my_handler(int s){

  if (s > 0)
    printf("\nCaught signal %d. Closing BAMs and ofstreams\n",s);
  
  // 
  if (r2c_reads->size() > 0) {
    cout << "Writing the r2c bam" << endl;
    writeR2C(false);
  }

  r2c_reads->clear();
  r2c_writer->Close();

  sleep(1);
  cleanR2CBam();

  all_align_stream->close();
  os_allbps->close(); 
  contigs_sam->close();
  sleep(1);

  // convert SAM to BAM
  string cmd = "samtools view contigs.sam -h -Sb > contigs.tmp.bam; samtools sort contigs.tmp.bam contigs; rm contigs.tmp.bam; samtools index contigs.bam";
  cout << cmd << endl;
  system(cmd.c_str());

  if (s == 0) 
    cout << "******************************" << endl 
	 << "Snowman completed successfully" << endl
         << "******************************" << endl;
  else 
    cout << "Snowman stopped due to signal interrupt, but successfully wrote output files"<< endl;
  exit(1); 

}

// main function to kick-off snowman
bool runSnowman(int argc, char** argv) {

  // start the interrupt handle
  struct sigaction sigIntHandler;
  sigIntHandler.sa_handler = my_handler;
  sigemptyset(&sigIntHandler.sa_mask);
  sigIntHandler.sa_flags = 0;
  sigaction(SIGINT, &sigIntHandler, NULL);

  parseRunOptions(argc, argv);

  // initialize the reader for converting chromosomes
  treader_for_convert = new BamReader();
  
  if (!treader_for_convert->Open(opt::bam.begin()->first)) {
    cerr << "Cannot open " << opt::bam.begin()->first << endl;
    exit(EXIT_FAILURE);
  }

  // parse the region file, count number of jobs
  GenomicRegionVector file_regions, regions_torun;
  int num_jobs = countJobs(file_regions, regions_torun); 

  // override the number of threads if need
  opt::numThreads = min(num_jobs, static_cast<int>(opt::numThreads));

  if (opt::verbose > 0 && !opt::disc_cluster_only) {
    //string disc_msg = opt::skip_disc ? "ON" : "OFF";
    //string chuck_msg   = opt::throw_disc   ? "ON" : "OFF";
    //string learn_msg   = opt::learn_params   ? "ON" : "OFF";
    //string only_msg   = opt::assemble_only_disc   ? "ON" : "OFF";
    //string normal_msg   = opt::normal_assemble   ? "ON" : "OFF";
    //string cent_msg     = opt::ignore_skip_cent ? "ON" : "OFF";
    cout << "Num threads:      " << opt::numThreads << endl;
    //cout << "Min overlap:      " << opt::assemb::minOverlap << endl;
    //cout << "Overlap error:    " << opt::assemb::error_rate << endl;
    cout << "Read filters:" << endl; 
    //cout << "   Min quality score:         " << opt::qualthresh << endl;
    //cout << "   Assemble only disc clusts: " << only_msg << endl;
    //cout << "   Parameter learning:        " << learn_msg << endl;
    //cout << "   Assemble normal also:      " << normal_msg << endl; 
    //cout << "   Ignore skip centromere:    " << cent_msg << endl; 

  } else if (opt::verbose > 0) {
    cout << "Num threads:      " << opt::numThreads << endl;
    cout << "RUNNING DISCORDANT CLUSTERING ONLY" << endl;
  }

  if (opt::verbose > 0) {
    for (auto i : opt::bam)
      cout << i.first << " -- " << i.second << endl;
  }

  r2c_writer = new BamWriter();
  //disc_writer= new BamWriter();

  dmap = new DMap();

  //debug
  // TODO remove the SO:coordinate because will be wrong.
  // set the r2c writer
  if (!r2c_writer->Open("r2c.bam", treader_for_convert->GetHeaderText(), treader_for_convert->GetReferenceData())) {
    cerr << "Could not open the r2c.bam for reading." << endl;
    exit(EXIT_FAILURE);
  }

  // set the discordant writer
  //  if (!disc_writer->Open("discordant.bam", treader_for_convert->GetHeaderText(), treader_for_convert->GetReferenceData())) {
  // cerr << "Could not open the discordant.bam for writing." << endl;
  //exit(EXIT_FAILURE);
  //}
  
  // setup the ASCII plots
  all_align_stream = new ofstream("alignments.all.txt", ios::out);
  os_allbps        = new ofstream("allbps.txt", ios::out);
  contigs_sam      = new ofstream("contigs.sam", ios::out);

  // setup the r2c_reads
  r2c_reads = new ReadMap();

  // set the MiniRules
  mr = new MiniRulesCollection("global@nbases[0,0];!hardclip;!supplementary;!duplicate;!qcfail%region@WG%discordant[0,800]%phred[4,100];clip[5,100]");
  if (opt::verbose > 0)
    cout << *mr << endl;

  // write the tumor header to the contigs SAM file
  (*contigs_sam) << treader_for_convert->GetHeaderText();

  // load the index reference genome onto the heap
  if (opt::verbose > 0)
    cout << "Loading the reference BWT index file for: "  << opt::refgenome << endl;
  idx = bwa_idx_load(opt::refgenome.c_str(), BWA_IDX_ALL);
  if (idx == NULL) {
    std::cerr << "Could not load the reference: " << opt::refgenome << endl;
    exit(EXIT_FAILURE);
  }

  // Create the queue and consumer (worker) threads
  wqueue<SnowmanWorkItem*>  queue;
  vector<ConsumerThread<SnowmanWorkItem>*> threadqueue;
  for (unsigned i = 0; i < opt::numThreads; i++) {
    ConsumerThread<SnowmanWorkItem>* threadr = new ConsumerThread<SnowmanWorkItem>(queue, opt::verbose > 0);
    threadr->start();
    threadqueue.push_back(threadr);
  }

  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);

  // print the jobs
  if (opt::verbose > 0) {
    if (opt::regionFile != "x" ) {
      for (auto& i : file_regions) 
	cout << "Input Regions: " << i << endl;
    }
    else {
      cout << "Running across whole genome" << endl;
    }
  }  

  // open a mutex
  if (pthread_mutex_init(&snow_lock, NULL) != 0) {
      printf("\n mutex init failed\n");
      return false;
  }

  // send the jobs
  unsigned count = 0;
  for (GenomicRegionVector::const_iterator it = regions_torun.begin(); it != regions_torun.end(); it++) {
    count++;
    AlignedContigVec  * cont_out = new AlignedContigVec();
    SnowmanWorkItem * item     = new SnowmanWorkItem(it->chr, it->pos1, it->pos2, count, cont_out);
    queue.add(item);
  }

  // wait for the threads to finish
  for (unsigned i = 0; i < opt::numThreads; i++) 
    threadqueue[i]->join();

  // close the BAMs
  my_handler(0);

  // display the run time
  if (opt::verbose > 0) {
    SnowUtils::displayRuntime(start);
    cout << endl;
  }

  return true;

}

// chunk a large section of reads and read them all into memory. Then scatter to small assemblies with runAll
bool grabReads(int refID, int pos1, int pos2, AlignedContigVec * cont_out) {

  // count the number of reads
  int num_n_reads = 0, num_t_reads = 0;
  double svbam_time_tumor = 0, svbam_time_normal = 0, assemble_time = 0, bwa_time = 0;

  typedef unordered_map<string, unique_ptr<BamAndReads> > BARMap;
  unique_ptr<BARMap> bar(new BARMap);

  BamQC qc; // debug. This should be a heaped vector
  
  //vector<pair<size_t,int>> varbam_data; // num reads, time

  // TODO make grabReads take GR
  GenomicRegion gr(refID, pos1, pos2);
  gr.pos1 = max(1, gr.pos1);
    
  // Open the BAMs and get the reads. only add to first (anchor)
  for (auto& bam : opt::bam) {

    // open a new BAM file at this region
    unique_ptr<BamAndReads> b(new BamAndReads(gr, mr, opt::verbose, bam.first, bam.second));
    //zauto b = make_unique<BamAndReads>(gr, mr, opt::verbose, bam.first, bam.second);
    (*bar)[bam.first] = move(b);
    
    // grab the reads and assign them to assembly regions
    (*bar)[bam.first]->readBam();

    // grab the mate region reads
#ifdef MATES
    (*bar)[bam.first]->calculateMateRegions();
#endif

    // add the times
    assert(bam.second.length() > 0);
    if (bam.second.at(0) == 'n') {
      svbam_time_normal += (*bar)[bam.first]->read_time;
      num_n_reads += (*bar)[bam.first]->unique_reads;
    } else if (bam.second.at(0) == 't') {
      svbam_time_tumor += (*bar)[bam.first]->read_time;      
      num_t_reads += (*bar)[bam.first]->unique_reads;
    }
  }

#ifdef MATES
  // merge the mate regions
  GenomicRegionVector mate_normal, mate_tumor;
  for (auto& bam : opt::bam) {
    if (bam.second.at(0) == 't')
      mate_tumor.insert(mate_tumor.begin(), (*bar)[bam.first]->mate_regions.begin(), (*bar)[bam.first]->mate_regions.end());
    else if (bam.second.at(0) == 'n')
      mate_normal.insert(mate_normal.begin(), (*bar)[bam.first]->mate_regions.begin(), (*bar)[bam.first]->mate_regions.end());      
  }
  mate_normal = GenomicRegion::mergeOverlappingIntervals(mate_normal);
  mate_tumor = GenomicRegion::mergeOverlappingIntervals(mate_tumor);
  size_t mate_region_size = 0;
  // check if the tumor overlaps 
  GenomicRegionVector mate_final;
  for (auto& t : mate_tumor) {
    bool good = true;
    for (auto& n : mate_normal) {
      if (t.getOverlap(n) > 1) {
	good = false;
	break;
      }
    }
    if (good) {
      mate_final.push_back(t);
      mate_region_size += t.width();
    }
  }
  // print it out
  if (opt::verbose > 1) {
    cout << "Mate regions(" << mate_final.size() << ") covering width of " << SnowUtils::AddCommas<size_t>(mate_region_size) << endl;
    if (opt::verbose > 2)
      for (auto &i : mate_final)
	cout << "    " << i << " width " << i.width() << endl;
  }

  // grab the mate reads
  for (auto& bam : opt::bam) {
    (*bar)[bam.first]->mate_regions = mate_final;
    (*bar)[bam.first]->readMateBam();

    // add the times
    assert(bam.second.length() > 0);
    if (bam.second.at(0) == 'n') {
      svbam_time_normal += (*bar)[bam.first]->mate_read_time;
      num_n_reads += (*bar)[bam.first]->mate_unique_reads;
    } else if (bam.second.at(0) == 't') {
      svbam_time_tumor += (*bar)[bam.first]->mate_read_time;      
      num_t_reads += (*bar)[bam.first]->mate_unique_reads;
    }

    if (opt::verbose > 1) 
      cout << bam.first << " " << *(*bar)[bam.first] << endl;
  }
#endif

  // start the assembly clock
  clock_t startr = std::clock();

  // get the regions to run
  GenomicRegionVector grv_small = gr.divideWithOverlaps(5000, 500);

#ifdef MATES
  // do the discordant read clustering. Put all reads in one vector,
  // and dedupe doubles
  unordered_map<string,bool> dd;
  BamAlignmentUPVector reads_fdisc;
  for (auto& bam : opt::bam) {
    for (auto& v : (*bar)[bam.first]->arvec)
      for (auto& r : v->reads) {
	string tmp;
	r->GetTag("SR",tmp);
	if (dd.count(tmp) == 0) 
	  reads_fdisc.push_back(r);
	dd[tmp] = true;
      }
  }

  // cluster the discordant reads
  DMap this_disc = clusterDiscordantReads(reads_fdisc);
#endif

  // do the assembly
  // loop through each assembly region. There should be 1 to 1 
  // concordance between grv_small and arvec in each BamAndReads
  for (size_t i = 0; i < grv_small.size(); i++) {
    
    // be extra careful, check for dupes
    // slow-ish, but just in case for weird mate overlaps
    // If not, program crashes in SGA
    unordered_map<string, bool> dupe_check; 
    
    BamAlignmentUPVector bav_join;
    for (auto& bam : *bar) {
      if (bam.second->arvec[i]->reads.size() > 0)
	for (auto& read : bam.second->arvec[i]->reads) {
	  string tmp;
	  read->GetTag("SR", tmp);
	  if (dupe_check.count(tmp) == 1)
	    cerr << "Unexpected duplicate: " << read->RefID << ":" << read->Position << " mate: " << read->MateRefID << ":" << read->MatePosition << " sr: " << tmp << endl;
	  else {
	    dupe_check[tmp] = true;
	    bav_join.push_back(read);
	  }
	}
    }

    // make the reads tables
    ReadTable pRT;
    pRT = ReadTable(bav_join);
    
    ContigVector contigs;
    ContigVector contigs1;

    // do the first round of assembly
    string name = "c_" + to_string(grv_small[i].chr+1) + "_" + to_string(grv_small[i].pos1) + "_" + to_string(grv_small[i].pos2);
    if (opt::verbose > 2)
      cout << "Doing assembly on: " << name << " with " << bav_join.size() << " reads" << endl;

    if (bav_join.size() > 1 && bav_join.size() < 10000) {
      clock_t assembleclock = clock();
      doAssembly(&pRT, name, contigs1, 0);

      // do the second round of assembly
      ReadTable pRTc;
      pRTc = ReadTable(contigs1);
      doAssembly(&pRTc, name, contigs, 1);
      assemble_time += (double)(clock()-assembleclock);
      
      // peform alignment of contigs to reference with BWA-MEM
      clock_t bwaclock = clock();
      BWAWrapper wrap;
      BWAReadVec bwarv;
      SamRecordVec samv;
      for (auto it = contigs.begin(); it != contigs.end(); it++) 
	bwarv.push_back(BWARead(it->getID(), it->getSeq()));
      wrap.addSequences(bwarv, idx, samv);
      bwa_time += (double)(clock() - bwaclock);

      
      // make aligned contigs from SAM records
      for (auto it = samv.begin(); it != samv.end(); it++) {
	AlignedContig ac(it->record, treader_for_convert);
	if (ac.m_align.size() > 1 && ac.getMinMapq() >= 10 && ac.getMaxMapq() >= 40) { 
	  ac.alignReadsToContigs(bav_join);
	  ac.sortReads();
	  ac.splitCoverage();
	  ac.getBreakPairs();
	  if (ac.maxSplit() > 0 || true) // must have > 0 split to continue
	    cont_out->push_back(ac);
	}
      }
    } // end bav_join.size() > 1
  } // end the assembly regions loop

  // end the assembly clock
  double assembly_time = clock() - startr;
  startr = clock();

  // combine discordant reads with breaks
  combineContigsWithDiscordantClusters(this_disc, cont_out);

  // get the breakpoints
  BPVec bp_glob;
  for (auto it = cont_out->begin(); it != cont_out->end(); it++) 
    bp_glob.push_back(it->getGlobalBreak());

  // add discordant reads
  //addDiscordantPairsBreakpoints(bp_glob, dmap);
    
  ////////////////////////////////////
  // MUTEX LOCKED
  ////////////////////////////////////
  // write to the global contig out
  pthread_mutex_lock(&snow_lock);  

#ifdef MATES
  // add the clusters to the map
  for (auto& i : this_disc)
    (*dmap)[i.first] = i.second;
#endif

  // store the reads in a map of read2contig alignment
  // TODO fancy combine
  for (auto& it : (*cont_out)) 
    for (auto& r : it.m_bamreads) 
      (*r2c_reads)[to_string(r->AlignmentFlag) + r->Name] = r;
#ifdef MATES
  // add in the reads from the discordant
  for (auto& i : this_disc) {
    for (auto& r : i.second.reads)
      (*r2c_reads)[to_string(r.second->AlignmentFlag) + r.second->Name] = r.second;
    for (auto& r : i.second.mates)
      (*r2c_reads)[to_string(r.second->AlignmentFlag) + r.second->Name] = r.second;
  }
#endif
  read_counter += r2c_reads->size();

  // write the r2c
  writeR2C(false);

  // send all bps to file
  stringstream out_allbps; out_allbps << "allbps.txt";
  for (BPVec::const_iterator it = bp_glob.begin(); it != bp_glob.end(); it++) 
    (*os_allbps) << it->toString() << endl;

  // send all alignments to ASCII plot
  for (auto it = cont_out->begin(); it != cont_out->end(); it++) 
    (*all_align_stream) << it->printAlignments() << endl;

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

  pthread_mutex_unlock(&snow_lock);
  /////////////////////////////////////
  // MUTEX UNLOCKED
  /////////////////////////////////////
  
  // delete the contig 
  delete cont_out;

  // display the run time
  if (opt::verbose > 0) {
    struct timespec finish;
    clock_gettime(CLOCK_MONOTONIC, &finish);
    double elapsed = (finish.tv_sec - start.tv_sec);
    int t = clock()/CLOCKS_PER_SEC;
    int min = (int)floor(elapsed / 60.0);
    int sec = (int)(elapsed-min*60);
    char buffer[140];
    double total_time = svbam_time_normal + svbam_time_tumor + assembly_time;
    //cout << svbam_time_normal <<  " " << svbam_time_tumor << " " << assembly_time << " " << write_time << endl;
    if (total_time == 0)
      total_time = 1;
    //assembly_time = max(1, (double)assembly_time);
    if (assembly_time == 0)
      assembly_time = 1;
    int tumor_time  = static_cast<int>(floor(svbam_time_tumor /total_time * 100.0));
    int normal_time = static_cast<int>(floor(svbam_time_normal/total_time * 100.0));
    int assem_time  = static_cast<int>(floor(assembly_time    /total_time * 100.0));
    int bwa_perc_time  = static_cast<int>(floor(bwa_time    / assembly_time * 100.0));
    int assemble_perc_time  = static_cast<int>(floor(assemble_time    / assembly_time * 100.0));
    string print1 = SnowUtils::AddCommas<int>(pos1);
    string print2 = SnowUtils::AddCommas<int>(pos2);

    sprintf (buffer, "Ran chr%2d:%11s-%11s T: %5d(%2d%%) N: %5d(%2d%%) A: %2d%%[B:%2d%%,A:%2d%%] | C: %6d #R: %6d #DC: %6d-- CPU: %dm%ds Wall: %dm%ds", 
	     refID+1,print1.c_str(),print2.c_str(),num_t_reads,tumor_time,
	     num_n_reads, normal_time, 
	     assem_time,bwa_perc_time, assemble_perc_time,
	     contig_counter, read_counter, static_cast<int>(dmap->size()),
	     //static_cast<int>(reads_all_encoded->size()),
	     static_cast<int>(floor( static_cast<double>(t) /60.0)), t % 60,
	     min, sec);
    printf ("%s\n",buffer);
  }
 
 return true;
}


/*
// take a collection of reads and assemble them all together and do discordant clustering
int runAll(vector<IterPair>& ivec, string name, AlignedContigVec * cont_out, DMap * dmap, BamAlignmentVectorAncParPairVector &vec)
{

  // parse name to get genomic region
  GenomicRegion gr(name);

  // total num reads
  size_t = num_reads;
  for (auto i : ivec)
    num_reads += i.size();

  // check that this doesn't overlap a blacklist regions
  if (gr.blacklistOverlap() != 0 && num_reads > 2000) {
    if (opt::verbose > 1) 
      cout << "Blacklist region hit at: " << gr << " with num reads: " << num_reads << endl;
    return 0;
  }

  // grab discordant read sequences, figure out how to deal with it
  // grab the pairmate sequences from the BAM
  BamAlignmentVector bav_disc; // vector of heaped BamAlignment pointers
  grabPairmateReads(ivec, gr, dmap, vec);
  unsigned post_reads = bavd.size();
  if (post_reads != orig_reads && opt::verbose > 1)
    cout << " -- Disc add " << orig_reads << " to " << post_reads << " at " << name << endl;

  // de-duplicate the reads by readname and position
  BamAlignmentVector bav, bav_name, bavd_tum;
  VariantBamReader::deduplicateReads(bavd, bav);

  // make the reads tables
  ReadTable pRT;
  pRT = ReadTable(bav);

  ContigVector contigs;
  ContigVector contigs1;
  
  // do the first round of assembly
  doAssembly(&pRT, name, contigs1, 0);

  // do the second round of assembly
  ReadTable pRTc;
  pRTc = ReadTable(contigs1);
  doAssembly(&pRTc, name, contigs, 1);

  // peform alignment of contigs to reference with BWA-MEM
  BWAWrapper wrap;
  BWAReadVec bwarv;
  SamRecordVec samv;
  for (auto it = contigs.begin(); it != contigs.end(); it++) 
    bwarv.push_back(BWARead(it->getID(), it->getSeq()));
  wrap.addSequences(bwarv, idx, samv);

  // make aligned contigs from SAM records
  for (auto it = samv.begin(); it != samv.end(); it++) {
    AlignedContig ac(it->record, treader_for_convert);
    if (ac.m_align.size() > 1 && ac.getMinMapq() >= 10 && ac.getMaxMapq() >= 40) { 
      ac.alignReadsToContigs(bav);
      ac.sortReads();
      ac.splitCoverage();
      ac.getBreakPairs();
      if (ac.maxSplit() > 0) // must have > 0 split to continue
	cont_out->push_back(ac);
    }
  }

  return 0;

}
*/

void parseRunOptions(int argc, char** argv) {
  bool die = false;

  if (argc <= 2) 
    die = true;

  string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
      case 'p': arg >> opt::numThreads; break;
      case 'h': die = true; break;
    case OPT_ASQG: opt::assemb::writeASQG = true; break;
	//case 'm': arg >> opt::minOverlap; break;
	//case 'l': arg >> opt::cutoff; break;
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
      case 'i': arg >> opt::isize; break;
	//case 's': arg >> opt::sleepDelay; break;
      case 'k': arg >> opt::regionFile; break;
	//case 'w': arg >> opt::mapq; break;
	//      case 'c': arg >> opt::chunk; break;
	//case 'q': arg >> opt::qualthresh; break;
	//case 'y': arg >> opt::skip; break;
	//case 'b': arg >> opt::numBubbleRounds; break;
	//case 'x': arg >> opt::divergence; break;
	//case 'e': arg >> opt::error_rate; break;
      case 'g': arg >> opt::refgenome; break;
	//case OPT_MIN_READ_COV: arg >> opt::min_read_cov; break;
	//case OPT_MIN_CLIP: arg >> opt::min_clip; break;
	//case OPT_NODISC: opt::skip_disc = true; break;
	//case OPT_CONTIG_BAM: opt::contig_bam = true; break;
	//case OPT_ASSEMBLE_ONLY_DISC: opt::assemble_only_disc = true; break;
	//case OPT_NOSUPP: opt::skip_supp = true; break;
      case OPT_DISC_CLUSTER_ONLY: opt::disc_cluster_only = true; break;
	//case OPT_NOR2: opt::skip_r2 = true; break;
	//case OPT_NO_R2CBAM: opt::no_r2c = true; break;
	//case OPT_IGNORE_SKIP_CENT: opt::ignore_skip_cent = true; break;
	//case OPT_THROW_DISC: opt::throw_disc = true; break;
	//case OPT_MEMORY: arg >> opt::memory; break;
	//case OPT_LEARN: opt::learn_params = true; break;
	//case OPT_ASSEMBLE_NORMAL: opt::normal_assemble = true; break;
	//case OPT_MEMORY_GOAL: arg >> opt::memory_goal; break;
	//case OPT_NM: arg >> opt::nmlim; break;
	//case OPT_DEBUG: opt::debug = true; break;
	//case OPT_CONTIG_LIMIT: arg >> opt::contig_write_limit; break;
    }
  }

  // clean the outdir
  //opt::outdir = SnowUtils::getDirPath(opt::outdir);

  // check that we input something
  if (opt::bam.size() == 0) {
    cerr << "Must add a bam file " << endl;
    exit(EXIT_FAILURE);
  }

  // check file validity
  if (opt::regionFile != "x") {
    if (false) {
      cerr << "Region file does not exist: " << opt::regionFile << endl;
      exit(EXIT_FAILURE);
    }
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
      exit(1);
    }
}

// Divide up the reads and send them to the assembler
/*
void chunkReadsForAssembly(const int refID, const int pos1,
			   AlignedContigVec * cont_out, BamAlignmentVectorAncParPairVector &vec, 
			   DMap * dmap) {

  //itpairvec
  ////IterPairMap for BAM1
  ////////<c_1:333-444, Iter(start, end)>
  ////////<c_1:400-500, Iter(start, end)>
  ////IterPairMap for BAM2
  ////////<c_1:333-444, Iter(start, end)>
  ////////<c_1:400-500, Iter(start, end)>

  //keys = c(c_1:333-444,c_1:400-500, ...)

  // loop through BAMs and get iterators of bav 
  // that lie in the window
  vector<IterPairMap> itpairvec;
  for (auto vec : it) {
    if (it->size() > 0) {
      IterPairMap itmap;
      // return vectors of iterators that points to positions for chunks
      getChunkReads(it.first, refID, pos1, LITTLECHUNK, WINDOWPAD, itmap);
      itpairvec.push_back(itmap);
    }
  }

  // grab all the normal and tumor keys. A key is like c_1:12323-231232
  StringVec keys;
  for (auto& it : itpairvec) { // loop the bams
    for (auto& jt : it) { // loop the keys
      keys.push_back(jt.first);
    }
  }

  // find the unique chunks
  sort( keys.begin(), keys.end() );
  keys.erase( unique( keys.begin(), keys.end() ), keys.end());

  // send a vector of IterPairMaps to runAll for assembly
  // loop through the keys (c_1...) and then the BAMs 
  for (StringVec::const_iterator it = keys.begin(); it < keys.end(); it++) {
    vector<IterPair> this_assembly;
    for (auto v : iterpairvec) { // v is an IterPairMap
      if (v.count(*it) == 1)
	this_assembly.push_back(v[*it]); // push an IterPair on
      runAll(this_assembly, *it, cont_out, dmap, vec)
    }
  }

       BamAlignmentVector fvec;
     if (tmap.count( (*it) ) > 0)
       fvec.insert(fvec.end(), tmap[(*it)].start, tmap[(*it)].end); 
     if (nmap.count( (*it) ) > 0)
       fvec.insert(fvec.end(), nmap[(*it)].start, nmap[(*it)].end);

     bool fvec_pass = fvec.size() < opt::skip && fvec.size() >= 3;
     if ( fvec_pass || (opt::contig_bam && fvec.size() > 0)) 
       runAll(fvec, *it, cont_out, dmap);
     else if (opt::verbose > 1)
       cout << "Too many/few reads at: " << *it << " with " << fvec.size() << " reads\n";
  
  }
  
}

*/

/*
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
*/

// just get a count of how many jobs to run. Useful for limiting threads
// also set the regions
int countJobs(GenomicRegionVector &file_regions, GenomicRegionVector &run_regions) {

  // open the region file if it exists
  bool rgfile = opt::regionFile.compare("x") != 0;
  if (rgfile) 
    file_regions = GenomicRegion::regionFileToGRV(opt::regionFile, 0);
  //else if (!opt::ignore_skip_cent)
  //  file_regions = GenomicRegion::non_centromeres; // from GenomicRegion.cpp
  else {
    RefVector ref = treader_for_convert->GetReferenceData();
    for (auto& i : ref) {
      file_regions.push_back(GenomicRegion(treader_for_convert->GetReferenceID(i.RefName), 1, i.RefLength));
    }
  //  file_regions = GenomicRegion::getWholeGenome();
  }

  if (file_regions.size() == 0) {
    cerr << "ERROR: Cannot read region file: " << opt::regionFile << " or something wrong with tumor bam header" << endl;
    exit(EXIT_FAILURE);
  }

  //int threadchunk = opt::chunk;
  int threadchunk = BIGCHUNK;
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
      startr = max(1,file_regions[jj].pos1-thispad);
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
/*
void grabPairmateReads(vector<IterPair>& ivec, const GenomicRegion window, DMap * dmap, BamAlignmentVectorAncParPairVector) {

  // get count of discordant reads
  size_t disc_countr = 0;
  
  // TODO SKIP THIS
  for (auto i : ivec) {
    for (BamAlignmentVector::const_iterator it = i.start; it != i.end; it++) { // i.start etc are pointers to BamAlignments on heap
      bool good_chr = (*it)->MateRefID < 25 && (*it)->RefID < 25 && (*it)->IsMapped() && (*it)->IsMateMapped();
      if ( (abs((*it)->InsertSize) > opt::isize || ((*it)->MateRefID != (*it)->RefID)) && good_chr) 
	disc_countr++;
    }
  }

  // need at least three discordant reads to proceed
  if (disc_countr < 2) 
    return;

  // find the clusters
  GenomicRegionVector partner_regions;
  GenomicRegionVector anchor_regions;

  // loop through the reads and cluster
  for (auto i : vec) {
    for (BamAlignmentVector::const_iterator it = i.start; it != i.end; it++) {

      // add intervals to the partner
      int sign = (*it)->IsMateReverseStrand() ? -1 : 1;
      int pos = (*it)->MatePosition * sign;
      if (pos != 0)
	partner_regions.push_back(GenomicRanges((*it)->MateRefID, pos - 200, pos + 200));
      
      // add intervals to the anchor
      sign = (*it)->IsReverseStrand() ? -1 : 1;
      int pos = (*it)->Position * sign;
      if (pos != 0)
	partner_regions.push_back(GenomicRanges((*it)->RefID, pos - 200, pos + 200));
    }
  }
  
  // reduce it down to minimum set
  GenomicRegionVector partner_min = partner_regions.mergeOverlappingIntervals(partner_regions);
  sort(partner_min.begin(), partner_min.end());
  GenomicRegionVector anchor_min = partner_regions.mergeOverlappingIntervals(partner_regions);
  sort(anchor_min.begin(), anchor_min.end());

  // create an interval tree map for both anchor and partner. These are the cluster boundaries
  GenomicIntervalTreeMap treemap_p = GenomicRegion::createTreeMap(partner_min);
  GenomicIntervalTreeMap treemap_a = GenomicRegion::createTreeMap(anchor_min);

  // create a map to store clusters and associated counts
  ClusterMap cluster_map;

  // loop through and assign reads to clusters
  for (auto i : ivec) {
    for (BamAlignmentVector::const_iterator it = i.start; it != i.end; it++) { // i.start etc are pointers to BamAlignments on heap

      // only cluster reads that are in mapped pairs
      if ((*it)->IsMapped() && (*it)->IsMateMapped()) {

	// find which cluster the anchor overlaps with
	int sign = (*it)->IsReverseStrand() ? -1 : 1;
	GenomicIntervalVector anchor_overlap;
	treemap_p[(*it)->RefID].findOverlapping(sign * (*it)->Position, sign * (*it)->Position, grv_ovlp);
	assert(anchor_overlap.size() == 1);
	
	// find which cluster the partner overlaps with
	sign = (*it)->IsMateReverseStrand() ? -1 : 1;
	GenomicIntervalVector partner_overlap;
	treemap_p[(*it)->MateRefID].findOverlapping((*it)->MatePosition, (*it)->MatePosition, grv_ovlp);
	assert(partner_overlap.size() == 1);

	// add to the running count
	string cluster = to_string((*it)->RefID) + to_string(anchor_overlap.start) + "-" + to_string(anchor_overlap.stop) + 
	  "_" + to_string((*it)->RefID) + to_string(partner_overlap.start) + "-" + to_string(partner_overlap.stop);
	ClusterMap::iterator ff = cluster_map.find(cluster);
	if (ff == cluster_map.end()) { // add a new cluster
	  GenomicRegion anc((*it)->RefID, anchor_overlap.start, anchor_overlap.end);
	  GenomicRegion par((*it)->MateRefID, partner_overlap.start, partner_overlap.end);
	  cluster_map[cluster] = DiscordantCluster(cluster, anc, par);
	} else { // TODO check if tumor / normal. If tumor, add first, otherwise add second
	  ff->second.tcount++;
	}
      }
      
    }
  }
  
  // loop through the cluster map and grab the reads if the clusters are big enough
  for (auto& cluster : cluster_map) {
    if (cluster.second.tcount > 2) {

      // Open the BAMs and get the reads
      for (auto bam : opt::bam) {
	
	// grab the reads according to the mini rules
	VariantBamReader reader(bam.first, NULL, mr, opt::verbose);
	reader.setPrefix(bam.second + "d");
	reader.writeVariantBam(qc,);
	bav_vector.push_back(bav);
	
	// store the runtime data
	num_reads = bav->size();
	svbam_time = (std::clock() - startr) / (double)(CLOCKS_PER_SEC / 1000);
	varbam_data.push_back(pair<num_reads, svbam_time>);    
      }
      
    }
  }
  

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

	DMap::iterator gg = dmap->find(ff->second);
	if (gg == dmap->end()) { // make a new one
	  DiscordantCluster dcc(ff->second);
	  dcc.reads.push_back(ba);
	  (*dmap)[ff->second] = DiscordantCluster(ff->second);
	} else {
	  gg->second.reads.push_back(ba);
	}
      }
  }

  return;

}

*/

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
  int min_overlap = opt::assemb::minOverlap;
  if (pass > 0) 
    min_overlap = 60;

  //errorRate = 0;
  //if (pass > 0)
  errorRate = (double)opt::assemb::error_rate;
  //if (pass > 0)
  //  errorRate = 0.05;
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
  if (opt::assemb::writeASQG) {

    // write ASQG to file for visualization
    stringstream asqgfile;
    asqgfile << name << "pass_" << pass << ".asqg";
    ofstream ofile(asqgfile.str(), ios::out);
    ofile << asqg_stream.str();
    ofile.close();

    // write the hits stream file
    stringstream hitsfile;
    hitsfile << name << "pass_" << pass << ".hits";
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

    // PERFORM THE ASSMEBLY
    assemble(asqg_stream, min_overlap, opt::assemb::maxEdges, opt::assemb::bExact, 
	     opt::assemb::trimLengthThreshold, opt::assemb::bPerformTR, opt::assemb::bValidate, opt::assemb::numTrimRounds, 
	     opt::assemb::resolveSmallRepeatLen, opt::assemb::numBubbleRounds, opt::assemb::gap_divergence, 
	     opt::assemb::divergence, opt::assemb::maxIndelLength, 102 /*cutoff*/, name, contigs);

    delete pQueryRIT;
    asqg_stream.str("");
    hits_stream.str("");

    // print out some results
    if (opt::verbose > 2) {
      if (contigs.size() >= 1) {
	cout << "Contig Count: " << contigs.size() << " at " << name << endl;
	if (opt::verbose > 3)
	  for (ContigVector::iterator i = contigs.begin(); i != contigs.end(); i++) 
	    cout << "   " << i->getID() << " " << i->getSeq().length() << " " << i->getSeq() << std::endl;
      }
    }


}

// deduplicate stored reads combine reads across different 
/*void combineR2C(ReadMap &read_in, ReadMap &read_out) {

  // try with a hash table
  for (auto& it : read_in) { 
    string jw = it.a.Name + to_string(it.a.AlignmentFlag);
    if (read_out.count(jw) == 0) {
      read_out.insert(pair<string, BamAlignment>(jw, it));
    } 
  }
    
  return;
  
  }*/

// given a set of reads, output structure that has either all the mate sequence,
// none of the discordant reads, or empty, depending on strictness of discordant 
// assembly options
/*void handleDiscordant(BamAlignmentVector &bavd, string name, GenomicRegion gr, DMap * dmap) {

  // grab the pairmate sequences from the BAM
  unsigned orig_reads = bavd.size();
  grabPairmateReads(bavd, gr, dmap);
  unsigned post_reads = bavd.size();
  if (post_reads != orig_reads && opt::verbose > 1)
    cout << " -- Disc add " << orig_reads << " to " << post_reads << " at " << name << endl;
  
    }*/

// perform clustering of discordant reads, with supplied orientations
/*
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
*/

//
void combineContigsWithDiscordantClusters(DMap this_dmap, AlignedContigVec * cont_out) {

  int padr = 400; 
  size_t count = 0;
  
  for (auto it = cont_out->begin(); it != cont_out->end(); it++) {
    count++;

    // check the global break
    GenomicRegion bp1(it->m_farbreak.refID1, it->m_farbreak.pos1, it->m_farbreak.pos1);
    bp1.pad(padr);
    GenomicRegion bp2(it->m_farbreak.refID2, it->m_farbreak.pos2, it->m_farbreak.pos2);
    bp2.pad(padr);
 
    for (auto &kt : this_dmap) { 
    
      bool bp1reg1 = bp1.getOverlap(kt.second.reg1) != 0;
      bool bp2reg2 = bp2.getOverlap(kt.second.reg2) != 0;
      
      //debug
      bool bp1reg2 = bp1.getOverlap(kt.second.reg2) != 0;
      bool bp2reg1 = bp2.getOverlap(kt.second.reg1) != 0;

      bool pass = bp1reg1 && bp2reg2;

      //debug
      pass = pass || (bp2reg1 && bp1reg2);

      if (pass) {
	it->addDiscordantCluster(kt.second); // add disc cluster to contig
	
	cout << "FOUND AN OVERLAP ON CONTIG " << it->getContigName() << " WITH DISC " << kt.second << endl;
	// check that we haven't already added a cluster
	kt.second.contig = it->getContigName(); // add contig to disc cluster
	if (it->m_farbreak.dc.reg1.pos1 == 0) {
	  it->m_farbreak.dc = kt.second; // add cluster to global breakpoints
	} else if (it->m_farbreak.dc.ncount < kt.second.ncount) { // choose one with normal support
	  it->m_farbreak.dc = kt.second;
	} else if (it->m_farbreak.dc.tcount < kt.second.tcount) { // choose one with more tumor support
	  it->m_farbreak.dc = kt.second;
	}
      }
      
    }
  }
}


// transfer discordant pairs not mapped to a split break 
void addDiscordantPairsBreakpoints(BPVec &bp, DMap * dmap) {

  if (opt::verbose > 1)
    cout << "...transfering discordant clusters to breakpoints structure" << endl;

  for (DMap::const_iterator it = dmap->begin(); it != dmap->end(); it++) {
    if (it->second.contig == "") { // this discordant cluster isn't already associated with a break
      cout << "ADDING DISC"  << endl;
      BreakPoint tmpbp(it->second);
      bp.push_back(tmpbp);
    }
  }
}

void writeR2C(bool makeIndex /* false */) {
  
  // sort by position
  BamAlignmentUPVector r2c_vec;
  for (auto& r : (*r2c_reads))
    r2c_vec.push_back(r.second);
  sort(r2c_vec.begin(), r2c_vec.end(), ByReadPosition());

  // add the reads
  for (auto& r : r2c_vec) 
    r2c_writer->SaveAlignment(*r);
  r2c_reads->clear();
  
  //index
  if (makeIndex) {
    BamReader r;
    if (!r.Open("r2c.bam")) {
      cerr << "Cannot open r2c.bam" << endl;
    } else {
      r.CreateIndex();
    }
  }

}

// cluster the discordant reads
DMap clusterDiscordantReads(BamAlignmentUPVector &bav) {

  // remove any reads that are not present twice
  unordered_map<string, int> tmp_map;
  for (auto& i : bav)
    if (tmp_map.count(i->Name) ==0)
      tmp_map[i->Name] = 1;
    else
      tmp_map[i->Name]++;
  BamAlignmentUPVector bav_dd;
  for (auto&i : bav)
    if (tmp_map[i->Name] == 2)
      bav_dd.push_back(i);

  // sort by position
  sort(bav_dd.begin(), bav_dd.end(), ByReadPosition());

  // clear the tmp map. Now we want to use it to store if we already clustered read
  tmp_map.clear();

  vector<BamAlignmentUPVector> fwd, rev, fwdfwd, revrev, fwdrev, revfwd;
  BamAlignmentUPVector this_fwd, this_rev;

  pair<int, int> fwd_info, rev_info; // refid, pos
  fwd_info = {-1,-1};
  rev_info = {-1,-1};

  // cluster in the READ direction
  for (auto& i : bav_dd) {
    if (i->IsMapped() && i->IsMateMapped() && abs(i->InsertSize) >= opt::isize && tmp_map.count(i->Name) == 0) {
      tmp_map[i->Name] = 0;
      // forward clustering
      if (!i->IsReverseStrand()) 
	_cluster(fwd, this_fwd, fwd_info, pair<int,int>(i->RefID, i->Position), i);
      // reverse clustering 
      else 
      	_cluster(rev, this_rev, rev_info, pair<int,int>(i->RefID, i->Position), i);
    }
  }
  // finish the last clusters
  fwd.push_back(this_fwd);
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
      if (!i->IsMateReverseStrand()) 
	_cluster(fwdfwd, this_fwd, fwd_info, pair<int,int>(i->MateRefID, i->MatePosition), i);
      // reverse clustering 
      else 
	_cluster(fwdrev, this_rev, rev_info, pair<int,int>(i->MateRefID, i->MatePosition), i);
    }
    // finish the last clusters
    fwdfwd.push_back(this_fwd);
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
      if (!i->IsMateReverseStrand()) 
	_cluster(revfwd, this_fwd, fwd_info, pair<int,int>(i->MateRefID, i->MatePosition), i);
      // reverse clustering 
      else 
	_cluster(revrev, this_rev, rev_info, pair<int,int>(i->MateRefID, i->MatePosition), i);
    }
    // finish the last clusters
    revfwd.push_back(this_fwd);
    revrev.push_back(this_rev);
  }    

  DMap dd;
  for (auto& v : fwdfwd) {
    DiscordantCluster d(v, bav_dd); /// slow but works
    dd[d.id] = d;
  }
  for (auto& v : revfwd) {
    DiscordantCluster d(v, bav_dd); /// slow but works
    dd[d.id] = d;
  }
  for (auto& v : fwdrev) {
    DiscordantCluster d(v, bav_dd); /// slow but works
    dd[d.id] = d;
  }
  for (auto& v : revrev) {
    DiscordantCluster d(v, bav_dd); /// slow but works
    dd[d.id] = d;
  }


  // print the clusters
  if (opt::verbose > 1) {
    for (auto& i : dd)
      cout << i.second << endl;
  }
  
  return dd;

}

// is this a read from a tumor
bool isTumorRead(const BamAlignmentUP &a) {

  string tmp;
  if (!a->GetTag("SR", tmp)) {
    cerr << "Failed to read SR tag in isTumorRead" << endl;
    exit(EXIT_FAILURE);
  }
  
  assert(tmp.length() > 0);
  return tmp.at(0) == 't';


}


bool _cluster(vector<BamAlignmentUPVector> &main, BamAlignmentUPVector &buff, pair<int,int> &last_info, pair<int, int> this_info, BamAlignmentUP &a) {

  if ( (this_info.first == last_info.first) && (this_info.second - last_info.second) <= DISC_PAD ) {
    buff.push_back(a);
    last_info = this_info;
    return true;
  } else {
    
    if (buff.size() >= MIN_PER_CLUST) {
      main.push_back(buff);
            
      //debug
      //cout << "----------------------" << endl;
      //for (auto& i : buff)
      //	cout << "   " << i->RefID << ":" << i->Position << (i->IsReverseStrand() ? "-" : "+") << " - " << i->MateRefID << ":" << i->MatePosition <<  (i->IsMateReverseStrand() ? "-" : "+") << " " << i->Name << endl;
    }

    buff.clear();
    buff.push_back(a);
    last_info = this_info;

    return false;
  }

}


// takes an unsorted r2c bam and sorts it and combines duplicates
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

  unordered_map<string, BamAlignmentUP> map;
  BamAlignmentUPVector vec;

  BamAlignment a;
  while (br.GetNextAlignment(a)) {
    string tmp;
    a.GetTag("SR",tmp);
    if (map.count(tmp) == 0) {
      BamAlignmentUP b(new BamAlignment(a));
      map[tmp] = b;
    } else {
      string cn1, cn2, dc1, dc2, sw1, sw2, al1, al2;
      map[tmp]->GetTag("CN",cn1);
      map[tmp]->GetTag("DC",dc1);
      map[tmp]->GetTag("SW",sw1);
      map[tmp]->GetTag("AL",al1);
      a.GetTag("CN",cn2);
      a.GetTag("DC",dc2);
      a.GetTag("SW",sw2);
      a.GetTag("AL",al2);

      //debug
      cout << "duplicate found for read " << a.Name << endl;

      if (!a.HasTag("CN"))
	a.AddTag("CN", "Z", cn1 + "x" + cn2);
      else
	a.EditTag("CN", "Z", cn1 + "x" + cn2);	

      if (!a.HasTag("AL"))
	a.AddTag("AL", "Z", al1 + "x" + al2);
      else
	a.EditTag("AL", "Z", al1 + "x" + al2);	

      if (!a.HasTag("SW"))
	a.AddTag("SW", "Z", sw1 + "x" + sw2);
      else
	a.EditTag("SW", "Z", sw1 + "x" + sw2);	

      if (!a.HasTag("DC"))
	a.AddTag("DC", "Z", dc1 + "x" + dc2);
      else
	a.EditTag("DC", "Z", dc1 + "x" + dc2);	

    }
  }

  // transfer to vector
  for (auto& r : map)
    vec.push_back(r.second);

  // sort it
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

  string cmd = "rm r2c.bam";
  system(cmd.c_str());

}
