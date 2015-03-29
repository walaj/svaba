#include "bamprocess.h"
#include <string>
#include <vector>
#include <getopt.h>
#include <iostream>
#include "GenomicRegion.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "SVBamReader.h"
#include "BamQC.h"
#include "SnowUtils.h"
#include "MiniRules.h"

using namespace std;
using namespace BamTools;

static const char *PREP_USAGE_MESSAGE =
"Usage: snowman preprocess -i <input.bam> -o <output.bam> [OPTIONS] \n\n"
"  Description: Process a BAM file for use with rearrangement variant callers by removing proper pairs and bad regions\n"
"\n"
" General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h, --help                           Display this help and exit\n"
"  -i, --input-bam                      BAM file to preprocess\n"
"  -o, --output-bam                     File to write final BAM to\n"
"  -q, --qc-only                        Loop through the BAM, but only to make the QC file\n"
"  -r, --rules-file                     Rules file\n"
"  -k, --proc-regions-file              csv file of regions to proess reads from, of format chr,pos1,pos2 eg 11,2232423,2235000. Default: Whole genome\n"
  //" Read Filters\n"
  //"  -w, --min-mapq                       Minimum mapping quality a read must have to be included. Default 0\n"
  //"  -c, --min-clip                       Minimum number of bases clipped bases to considered clipped. Default: 5\n"
  //"  -n  --nm-limit                       Skip reads that have more than nm mismatches (NM tag). Default: 15\n"
  //"  -p  --perc-limit                     Exit with failure if more than this percent of reads are weird (Must see over 50m). Default: 30\n"
  //"  -s, --isize                          Insert size threshold to consider discordant. Default: 800\n"
  //"  -m, --min-readlength                 Removes reads < this cutoff. If quality score trimming is implemented, evaluates after trimming. Default: 50\n"
  //"  -z, --min-phred                      When measuring clipping and read length, dont count bases with less than this quality. Set to 0 to turn off. Default: 4\n"
  //"  -e, --exclude-with-n                 Remove reads that have N in their sequence. Default: OFF\n"
  //"      --no-pileup-check                Don't perform the low mapq pileup check. Default is to exclude reads with a pileup of > 100 mapq0 reads within 200 bp.\n"
  //"      --skip-centromeres               Don't processes centromers at all. Default OFF\n"
" Optional Input\n"
"\n";

namespace opt {

  static string bam;
  static string out;
  static int isize = 800;
  static size_t verbose = 1;
  //static int pad = 500;
  static string mutect_callstats = "";
  static string mutect2_regions = "";
  

  static int mapq = 0;
  static int minOverlap = 30;
  static int min_clip = 5;
  static bool skip_supp = true;
  static bool skip_cent = false;
  static string full_regions = "";
  static string proc_regions = "";
  static int min_length = 50;
  static int min_phred = 4;

  static int nmlim = 15;
  static int perc_limit = 30;
  static bool exclude_n = false;

  static bool no_pileup_check = false;
  static bool qc_only = false;
}

static const char* shortopts = "hv:qejb:i:o:w:n:p:s:r:m:z:c:k:u:x:";
static const struct option longopts[] = {
  { "help",                       no_argument, NULL, 'h' },
  { "verbose",                    required_argument, NULL, 'v' },
  { "input-bam",                  required_argument, NULL, 'i' },
  { "output-bam",                 required_argument, NULL, 'o' },
  { "min-mapq",                 required_argument, NULL, 'w' },
  { "min-clip",                 required_argument, NULL, 'c' },
  { "nm-limit",                 required_argument, NULL, 'n' },
  { "perc-limit",                 required_argument, NULL, 'p' },
  { "isize",                 required_argument, NULL, 's' },
  { "qc-only",                 no_argument, NULL, 'q' },
  { "full-regions-file",                 required_argument, NULL, 'r' },
  { "proc-regions-file",                 required_argument, NULL, 'k' },
  { "min-length",                 required_argument, NULL, 'm' },
  { "min-phred",                 required_argument, NULL, 'z' },
  { "exclude-with-n",                 no_argument, NULL, 'e' },
  { "no-pileup-check",                 no_argument, NULL, 'j' },
  { "skip-centromeres",                 no_argument, NULL, 'b' },
  { "mutect-callstats",               required_argument, NULL, 'u' },
  { "mutect2-regions",               required_argument, NULL, 'x' },
  { NULL, 0, NULL, 0 }
};

static struct timespec start;

void runPrep(int argc, char** argv) {

  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);

  parsePrepOptions(argc, argv);

  if (opt::verbose > 0) {
    cout << "Input BAM:  " << opt::bam << endl;
    cout << "Output BAM: " << opt::out << endl;
    cout << "Input rules file: " << opt::full_regions << endl;
    cout << "Input proc regions file: " << opt::proc_regions << endl;
    /*    cout << "MuTect callstats file:   " << opt::mutect_callstats << endl;
    cout << "MuTect2 region file:     " << opt::mutect2_regions << endl;
    cout << "Read Filters:            " << endl;
    cout << "   Min insert size (s):           " << opt::isize << endl;
    cout << "   Min clip length (c):           " << opt::min_clip << endl;
    cout << "   Min read length (m):           " << opt::min_length << endl;
    cout << "   Min phred quality (z):         " << opt::min_phred << endl;
    cout << "   Min MAPQ (w):                  " << opt::mapq  << endl;
    cout << "   Max NM (n):                    " << opt::nmlim  << endl;
    cout << "   Skip supplementary reads:      " << (opt::skip_supp ? "ON" : "OFF") << endl;
    cout << "   Skip centromere:               " << (opt::skip_cent ? "ON" : "OFF") << endl; 
    cout << "   Skip reads w/N base (e):       " << (opt::exclude_n ? "ON" : "OFF") << endl; 
    cout << "   Percent weird limit:           " << opt::perc_limit << endl;
    cout << "   Pileup-check:                  " << (!opt::no_pileup_check ? "ON" : "OFF") << endl;
    */
  }

  // open the BAM Reader
  BamReader reader;
  if (!reader.Open(opt::bam) && !opt::qc_only) {
    cerr << "Could not open BAM: " << opt::bam << endl;
    exit(EXIT_FAILURE);
  }

  // get the BAM header
  RefVector ref;
  SamHeader sam;
  if (!opt::qc_only) {
    SVBamReader::getRefVector(opt::bam, ref);
    SVBamReader::getSamHeader(opt::bam, sam);
  }

  // open the BAM Writer
  BamWriter writer;
  if (!writer.Open(opt::out, sam, ref) && !opt::qc_only) {
    cerr << "Could not open output BAM: " << opt::out << endl;
    exit(EXIT_FAILURE);
  }

  // try the MiniRules
  MiniRulesCollection * mr = new MiniRulesCollection(opt::full_regions);
  if (opt::verbose > 0)
    cout << (*mr);
  
  // make a bed file
  if (opt::verbose > 0)
    cout << "...sending merged regions to BED file" << endl;
  mr->sendToBed("merged_rules.bed");

  // parse the proc region file
  GenomicRegionVector grv_proc_regions = GenomicRegion::regionFileToGRV(opt::proc_regions, 0);
  if (grv_proc_regions.size() > 0) {
    sort(grv_proc_regions.begin(), grv_proc_regions.end());      
    if (opt::verbose > 0) {
      cout << "Number of regions to process: " << grv_proc_regions.size() << endl;
      if (opt::verbose > 1)
	for (auto it = grv_proc_regions.begin(); it != grv_proc_regions.end(); it++)
	  cout << "   Process-region: " << *it << endl;
    }
  } else if (opt::skip_cent) {
    grv_proc_regions = GenomicRegion::non_centromeres;
    if (opt::verbose > 0)
      cout << "Processing whole genome (minus centromeres)" << endl;
  } else {
    grv_proc_regions = GenomicRegion::getWholeGenome();
    if (opt::verbose > 0)
      cout << "Processing whole genome (including centromeres)" << endl;
  }

  /*
  // parse the full region file
  GenomicRegionVector mutect1 = GenomicRegion::regionFileToGRV(opt::mutect_callstats, 0);
  GenomicRegionVector mutect2 = GenomicRegion::regionFileToGRV(opt::mutect2_regions, 0);
  GenomicRegionVector regions = GenomicRegion::regionFileToGRV(opt::full_regions, 0);
  // cat them
  GenomicRegionVector grv_full_regions = regions;
  grv_full_regions.insert(grv_full_regions.begin(), mutect1.begin(), mutect1.end());
  grv_full_regions.insert(grv_full_regions.begin(), mutect2.begin(), mutect2.end());

  // merge them
  if (opt::verbose > 0)
    cout << "...merging overlapping intervals" << endl;
  GenomicRegionVector grv_m = GenomicRegion::mergeOverlappingIntervals(grv_full_regions);
  
  // send to a BED file
  if (opt::verbose > 0)
    cout << "...sending merged regions to a BED file" << endl;
  GenomicRegion::sendToBed(grv_m, "merged_regions.bed");
  */

  //GenomicRegionVector grv_full_regions = GenomicRegion::regionFileToGRV(opt::full_regions, opt::isize);
  //if (grv_m.size() > 0)
  //  sort(grv_m.begin(), grv_m.end());      

  if (grv_proc_regions.size() > 0)
    sort(grv_proc_regions.begin(), grv_proc_regions.end());      

  /*  if (opt::verbose > 0) {
    cout << "Number of full regions to include: " << grv_m.size() << endl;
    cout << "   -- Region file regions: " << regions.size() << endl;
    cout << "   -- Mutect call stats:   " << mutect1.size() << endl;
    cout << "   -- Mutect2 regions:     " << mutect2.size() << endl;
    cout << "      ** merged down from " << grv_full_regions.size() << " to " << grv_m.size() << endl;
    cout << "      ** Total span (in bp) is: " << SnowUtils::AddCommas<size_t>(span) << endl;

    int fcount = 0;
    if (opt::verbose > 1)
      for (auto it = grv_m.begin(); it != grv_m.end(); it++) {
	cout << "   Full-region: " << *it << " count " << fcount << endl;
	fcount++;
      }
  }    
  */

  BamAlignment a;
  BamQC qc; 

  //SVBamReader sv_reader(opt::bam, "", opt::isize, opt::mapq, opt::min_phred, 
  //			opt::minOverlap, opt::skip_supp, false /* skip_r2 */, opt::min_clip, opt::verbose);
  SVBamReader sv_reader(opt::bam, "", mr, opt::verbose);
  sv_reader.perclimit = opt::perc_limit;
  sv_reader.min_length = opt::min_length;
  sv_reader.exclude_n = opt::exclude_n;
  sv_reader.no_pileup_check = opt::no_pileup_check;
  if (!sv_reader.findBamIndex()) {
    cerr  << "Failed to open BAM index" << endl;
    exit(EXIT_FAILURE);
  }

  // set the iterator to loop through regions for creating trees
  //GenomicRegionVector::const_iterator grv_m_it = grv_m.begin();

  // loop through the process regions and extract files
  for (auto it = grv_proc_regions.begin(); it != grv_proc_regions.end(); it++) {
    if (!sv_reader.setBamRegion(*it)) {
      cerr << "Failed to set BAM region " << (*it) << endl;
      exit(EXIT_FAILURE);
    }

    // generate interval tree for this region only
    // assumes that everything is sorted
    // TODO make move not copy
    //GenomicRegionVector this_grv;
    //if (opt::verbose > 0)
    // cout << "...building interval tree" << endl;
    //while (grv_m_it != grv_m.end() && grv_m_it->chr == it->chr && grv_m_it->pos1 <= it->pos2) {
    //  this_grv.push_back(*grv_m_it);
    //  grv_m_it++;
    //}
    //GenomicIntervalTreeMap this_tree_map = GenomicRegion::createTreeMap(this_grv);
      
    if (opt::verbose > 0)
      cout << "Running region: "  << (*it) << endl;
    sv_reader.preprocessBam(writer, qc, opt::qc_only);
    
    if (opt::verbose > 0) {
      SnowUtils::displayRuntime(start);
      cout << endl;
    }

  }

  if (!opt::qc_only)
    writer.Close();

  // index it
  if (!reader.Open(opt::out) && !opt::qc_only) {
    cerr << "Cant read output BAM" << endl;
    exit(EXIT_FAILURE);
  }

  // create the index file. The input bam should be sorted.                                                                                                                                                                                   
  if (!reader.CreateIndex() && !opt::qc_only) {                                                                                                                                                                                                          
    cerr << "Failed to created index for "  << opt::out << endl;
    cerr << " This can happen if the input bam is not sorted. Use samtools sort inbam.bam inbam.sort; samtools index inbam.sort.bam" << endl;
    exit(EXIT_FAILURE);
  }

  // write out the stats
  ofstream stats_out;
  stats_out.open("./qcreport.txt");
  stats_out << qc;
  
}

void parsePrepOptions(int argc, char** argv) {

  bool die = false;

  if (argc < 2) 
    die = true;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'h': die = true; break;
    case 'v': arg >> opt::verbose; break;
    case 'i': arg >> opt::bam; break;
    case 'o': arg >> opt::out; break;
    case 'w': arg >> opt::mapq; break;
    case 'n': arg >> opt::nmlim; break;
    case 'u': arg >> opt::mutect_callstats; break;
    case 'x': arg >> opt::mutect2_regions; break;
    case 'p': arg >> opt::perc_limit; break;
    case 's': arg >> opt::isize; break;
    case 'q': opt::qc_only = true; break;
    case 'r': arg >> opt::full_regions; break;
    case 'k': arg >> opt::proc_regions; break;
    case 'm': arg >> opt::min_length; break;
    case 'z': arg >> opt::min_phred; break;
    case 'c': arg >> opt::min_clip; break;
    case 'e': opt::exclude_n = true; break;
    case 'b': opt::skip_cent = true; break;
    }
  }

  // dont stop the run for bad bams for quality checking only
  opt::perc_limit = opt::qc_only ? 101 : opt::perc_limit;

  // something went wrong, kill
  if (die) {
    cout << "\n" << PREP_USAGE_MESSAGE;
    exit(1);
  }

}


