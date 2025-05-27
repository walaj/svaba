//run_svaba.cpp

#include <thread>
#include <mutex>
#include "threadpool.h"

#include <getopt.h>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <map>
#include <vector>
#include <cassert>

#include "SeqLib/GenomicRegion.h"
#include "SeqLib/UnalignedSequence.h"
#include "SeqLib/ReadFilter.h"

#include "svabaLogger.h"
#include "KmerFilter.h"
#include "vcf.h"
#include "DBSnpFilter.h"
#include "svabaUtils.h"
#include "LearnBamParams.h"
#include "SeqLib/BFC.h"
#include "svaba_params.h"

#include "svabaOutputWriter.h"
#include "svabaAssemblerEngine.h"

using namespace SeqLib;
using namespace std;

// --- handy aliases ---
using BamMap              = std::map<string, string>;
using CountPair           = std::pair<int,int>;
using WalkerMap           = std::map<string, svabaBamWalker>;

//#define DEBUG 1

// 1 is normal logging
// 2 is heavy
constexpr inline int log_level = 1;

static RefGenome * ref_genome;
static unordered_map<string, BamParamsMap> params_map; // key is bam id (t000), value is map with read group as key
static BamHeader bwa_header;
static set<string> prefixes;

static int min_dscrd_size_for_variant = 0; // set a min size for what we can call with discordant reads only. 
// something like max(mean + 3*sd) for all read groups

static unordered_map<string, int> min_isize_for_disc;

//static BWAWrapper * main_bwa = nullptr;
static Filter::ReadFilterCollection * mr;
static GRC germline_svs;
static DBSnpFilter * dbsnp_filter;
static GRC file_regions;

// time
static struct timespec start;

// learned value 
static int max_mapq_possible;
static string args = "svaba"; // hold string of what the input args were
static int32_t readlen = 0;

namespace opt {

  // SGA options
  namespace sga {
    static int minOverlap = 0;
    static float error_rate = 0; 
    static int num_assembly_rounds = 3;
  }

  // error correction options
  static string ec_correct_type = "f";
  static double ec_subsample = 0.50;

  // run in single end mode?
  bool single_end = false;

  // output options
  static bool all_contigs = false;   // output all contigs

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
  static string rules = "{\"global\" : {\"duplicate\" : false, \"qcfail\" : false}, \"\" : { \"rules\" : [FRRULES,{\"rr\" : true},{\"ff\" : true}, {\"rf\" : true}, {\"ic\" : true}, {\"clip\" : 5, \"length\" : READLENLIM}, {\"ins\" : true}, {\"del\" : true}, {\"mapped\": true , \"mate_mapped\" : false}, {\"mate_mapped\" : true, \"mapped\" : false}]}}";  
  static int max_cov = 100;
  static size_t mate_lookup_min = 3;
  static size_t mate_region_lookup_limit = 400;
  static int32_t max_reads_per_assembly = -1; // set default of 50000 in parseRunOptions
  static bool no_bad_avoid = true; // if true, don't avoid previously bad regions

  // additional optional params
  static int chunk = 25000;
  static string regionFile;  // region to run on
  static string analysis_id = "no_id";

  // runtime parameters
  static int verbose = 0;
  static int numThreads = 1;
  static bool hp = false; // should run in highly-parallel mode? (no file dump til end)

  // data
  static BamMap bam;
  static string refgenome; // = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
  static string blacklist; 
  static string germline_sv_file;
  static string dbsnp; 
  static string main_bam = "-"; // the main bam

  // optimize defaults for single sample mode
  static bool germline = false; 
  
  // indel probability cutoffs
  static double lod = 8; // LOD that variant is not ref
  static double lod_db = 6; // same, but at DB snp site (want lower bc we have prior)
  static double lod_somatic = 6; // LOD that normal is REF
  static double lod_somatic_db = 10; // same, but at DBSNP (want higher bc we have prior that its germline)
  static double scale_error = 1; // how much to emphasize erorrs. 1 is standard. 0 is assume no errors

  // input options
  static bool override_reference_check = false; // allow the user, with caution, to use two different reference genomes 
                                                // (one for BAM alignment, one as ref for svaba)
}

enum { 
  OPT_HP,
  OPT_LOD,
  OPT_LOD_DB,
  OPT_LOD_SOMATIC,
  OPT_LOD_SOMATIC_DB,
  OPT_DISC_CLUSTER_ONLY,
  OPT_EC_SUBSAMPLE,
  OPT_DISCORDANT_ONLY,
  OPT_NUM_ASSEMBLY_ROUNDS,
  OPT_SCALE_ERRORS,
  OPT_OVERRIDE_REFERENCE_CHECK
};

static const char* shortopts = "hzIAt:n:p:v:r:G:e:k:c:a:m:B:D:Y:S:L:s:V:R:K:E:C:x:M:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "case-bam",               required_argument, NULL, 't' },
  { "all-contigs",             no_argument, NULL, 'A' },  
  { "id-string",               required_argument, NULL, 'a' },
  { "hp",                      no_argument, NULL, OPT_HP },
  { "control-bam",              required_argument, NULL, 'n' },
  { "threads",                 required_argument, NULL, 'p' },
  { "chunk-size",              required_argument, NULL, 'c' },
  { "region-file",             required_argument, NULL, 'k' },
  { "rules",                   required_argument, NULL, 'r' },
  { "reference-genome",        required_argument, NULL, 'G' },
  { "min-overlap",             required_argument, NULL, 'm' },
  { "dbsnp-vcf",               required_argument, NULL, 'D' },
  { "disc-sd-cutoff",          required_argument, NULL, 's' },
  { "mate-lookup-min",         required_argument, NULL, 'L' },
  { "germline-sv-database",    required_argument, NULL, 'V' },
  { "ec-subsample",            required_argument, NULL, 'E' },
  { "lod",                     required_argument, NULL, OPT_LOD },
  { "lod-dbsnp",               required_argument, NULL, OPT_LOD_DB },
  { "lod-somatic",             required_argument, NULL, OPT_LOD_SOMATIC },
  { "lod-somatic-dbsnp",       required_argument, NULL, OPT_LOD_SOMATIC_DB },
  { "scale-errors",            required_argument, NULL, OPT_SCALE_ERRORS },
  { "discordant-only",         no_argument, NULL, OPT_DISCORDANT_ONLY },
  { "ec-correct-type",         required_argument, NULL, 'K'},
  { "error-rate",              required_argument, NULL, 'e'},
  { "verbose",                 required_argument, NULL, 'v' },
  { "blacklist",               required_argument, NULL, 'B' },
  { "max-coverage",            required_argument, NULL, 'C' },
  { "max-reads",               required_argument, NULL, 'x' },
  { "max-reads-mate-region",   required_argument, NULL, 'M' },
  { "num-assembly-rounds",     required_argument, NULL, OPT_NUM_ASSEMBLY_ROUNDS },
  { "override-reference-check",no_argument, NULL, OPT_OVERRIDE_REFERENCE_CHECK},
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
"  -x, --max-reads                      Max total read count to read in from assembly region. Set 0 to turn off. [50000]\n"
"  -M, --max-reads-mate-region          Max weird reads to include from a mate lookup region. [400]\n"
"  -C, --max-coverage                   Max read coverage to send to assembler (per BAM). Subsample reads if exceeded. [500]\n"
"      --discordant-only                Only run the discordant read clustering module, skip assembly. \n"
"      --num-assembly-rounds            Run assembler multiple times. > 1 will bootstrap the assembly. [2]\n"
"      --num-to-sample                  When learning about inputs, number of reads to sample. [2,000,000]\n"
"      --hp                             Highly parallel. Don't write output until completely done. More memory, but avoids all thread-locks.\n"
"      --override-reference-check       With much caution, allows user to run svaba with different reference genomes for BAMs and -G\n"
"  Output options\n"
"  -A, --all-contigs                    Output all contigs that were assembled, regardless of mapping or length. [off]\n"
"      --read-tracking                  Track supporting reads by qname. Increases file sizes. [off]\n"
"      --write-extracted-reads          For the case BAM, write reads sent to assembly to a BAM file. [off]\n"
"  Optional external database\n"
"  -D, --dbsnp-vcf                      DBsnp database (VCF) to compare indels against\n"
"  -B, --blacklist                      BED-file with blacklisted regions to not extract any reads from.\n"
"  -V, --germline-sv-database           BED file containing sites of known germline SVs. Used as additional filter for somatic SV detection\n"
"  Assembly and EC params\n"
"  -m, --min-overlap                    Minimum read overlap, an SGA parameter. Default: 0.4* readlength\n"
"  -e, --error-rate                     Fractional difference two reads can have to overlap. See SGA. 0 is fast, but requires error correcting. [0]\n"
"  -K, --ec-correct-type                (f) Fermi-kit BFC correction, (s) Kmer-correction from SGA, (0) no correction (then suggest non-zero -e) [f]\n"
"  -E, --ec-subsample                   Learn from fraction of non-weird reads during error-correction. Lower number = faster compute [0.5]\n"
"\n";

// forward declaration, need non-static and not in anonymous namespace to keep external linkage
void  runsvaba       (int argc, char** argv);

// foward declarations
namespace {
     
  void setWalkerParams(svabaBamWalker& walk) {
    
    walk.main_bwa = main_bwa; // set the pointer
    walk.blacklist = blacklist;
    walk.do_kmer_filtering = (opt::ec_correct_type == "s" || opt::ec_correct_type == "f");
    walk.kmer_subsample = opt::ec_subsample;
    walk.max_cov = opt::max_cov;
    walk.m_mr = mr;  // set the read filter pointer
    walk.m_limit = opt::max_reads_per_assembly;
    
  }

void alignReadsToContigs(BWAWrapper& bw,
			 const UnalignedSequenceVector& usv, 
			 svabaReadVector& input_reads,
			 vector<AlignedContig>& this_alc,
			 const RefGenome*  rg) {
  
  if (!usv.size())
    return;

  // get the reference info
  GRC g;
  for (auto& a : this_alc)
    for (auto& i : a.getAsGenomicRegionVector()) {
      i.Pad(100);
      g.add(i);
    }
  g.MergeOverlappingIntervals();

  // get the reference sequence
  vector<string> ref_alleles;
  for (auto& i : g)
    //if (i.chr < 24) //1-Y
      try {
	string tmpref = rg->QueryRegion(i.ChrName(bwa_header), i.pos1, i.pos2);
	ref_alleles.push_back(tmpref); 
      } catch (...) {
	//cerr << "Caught exception for ref_allele on " << i << endl;
      }
  // make the reference allele BWAWrapper
  BWAWrapper bw_ref;
  UnalignedSequenceVector usv_ref;
  int aa = 0;
  for (auto& i : ref_alleles) {
    if (!i.empty())
      usv_ref.push_back({to_string(aa++), i, string()}); // name, seq, qual
  }
  if (!usv_ref.size())
    bw_ref.ConstructIndex(usv_ref);
  
  // set up custom alignment parameters, mean
  bw_ref.SetGapOpen(16); // default 6
  bw.SetGapOpen(16); // default 6
  bw_ref.SetMismatchPenalty(9); // default 2
  bw.SetMismatchPenalty(9); // default 4

  for (auto i : input_reads) {
    
    BamRecordVector brv, brv_ref;

    // try the corrected seq first
    //string seqr = i.GetZTag("KC");
    //  if (seqr.empty())
    //	seqr = i.QualitySequence();
    string seqr = i.Seq();
    
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
      max_as = max(max_as, thisas);
    }

    // align to the reference alleles
    if (!bw_ref.IsEmpty())
      bw_ref.AlignSequence(seqr, i.Qname(), brv_ref, hardclip, 0.60, 10);

    // get the maximum reference alignment score
    int max_as_r = 0;
    for (auto& r : brv_ref) {
      int thisas= 0;
      r.GetIntTag("AS", thisas);
      max_as_r = max(max_as_r, thisas);
    }
    
    // reject if better alignment to reference
    if (max_as_r > max_as) {
      //cerr << " Alignment Rejected for " << max_as_r << ">" << max_as << "  " << i << endl;
      //cerr << "                        " << max_as_r << ">" << max_as << "  " << brv_ref[0] << endl;
      continue;
    }

    // convert to a svabaReadVector
    svabaReadVector brv_svaba;
    for (auto& r : brv)
      brv_svaba.push_back(svabaRead(r, i.Prefix()));
    brv.clear();

    // make sure we have only one alignment per contig
    set<string> cc;

    // check which ones pass
    BamRecordVector bpass;
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

      this_r2c.AddAlignment(r);
      i.AddR2C(usv[r.ChrID()].Name, this_r2c);
      
      // add the read to the right contig (loop to check for right contig)
      for (auto& a : this_alc) {
	if (a.getContigName() != usv[r.ChrID()].Name)
	  continue;
	a.AddAlignedRead(i);
      }
      
    } // end passing bwa-aligned read loop 
  } // end main read loop
}

  
void collectAndClearReads(WalkerMap& walkers,
			     svabaReadVector& brv,
			     vector<char*>& learn_seqs,
			     unordered_set<string>& dedupe) {

  // concatenate together all the reads from the different walkers
  for (auto& w : walkers) {
    for (auto& r : w.second.reads) {
      string sr = r.SR();
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

MateRegionVector collectNormalMateRegions(WalkerMap& walkers) {
  
  MateRegionVector normal_mate_regions;
  for (auto& b : opt::bam)
    if (b.first.at(0) == 'n')
      normal_mate_regions.Concat(walkers[b.first].mate_regions);

  normal_mate_regions.MergeOverlappingIntervals();
  normal_mate_regions.CreateTreeMap();
  
  return normal_mate_regions;

}

MateRegionVector collectSomaticMateRegions(WalkerMap& walkers, MateRegionVector& bl) {


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

CountPair collectMateReads(WalkerMap& walkers,
			   const MateRegionVector& mrv,
			   int round,
			   GRC& this_bad_mate_regions) {

  CountPair counts = {0,0};  

  if (!mrv.size())
    return counts;
  
  for (auto& w : walkers) {

    int oreads = w.second.reads.size();
    w.second.m_limit = opt::mate_region_lookup_limit;

    // convert MateRegionVector to GRC
    GRC gg;
    for (auto& s : mrv) 
      gg.add(GenomicRegion(s.chr, s.pos1, s.pos2, s.strand));

    assert(w.second.SetMultipleRegions(gg));
    w.second.get_coverage = false;
    w.second.get_mate_regions = (round != MAX_MATE_ROUNDS);

    // clear out the current store of mate regions, since we 
    // already added these to the to-do pile
    w.second.mate_regions.clear();

    this_bad_mate_regions.Concat(w.second.readBam()); 
    
    // update the counts
    if (w.first.at(0) == 't') 
      counts.first += (w.second.reads.size() - oreads);
    else
      counts.second += (w.second.reads.size() - oreads);

  }
  
  return counts;
}
  
  CountPair runMateCollectionLoop(const GenomicRegion& region,
				  WalkerMap& wmap,
				  GRC& badd,
				  SvabaLogger& logger,
				  const SvabaOptions& opts) {

    GRC this_bad_mate_regions; // store the newly found bad mate regions
    
    CountPair counts = {0,0};
    
    MateRegionVector all_somatic_mate_regions;
    all_somatic_mate_regions.add(MateRegion(region.chr, region.pos1, region.pos2)); // add the origional, don't want to double back
    
    for (int jjj = 0; jjj <  MAX_MATE_ROUNDS; ++jjj) {
      
      MateRegionVector normal_mate_regions;
      if (wmap.size() > 1) // have at least one control bam
	normal_mate_regions = collectNormalMateRegions(wmap);
      
      // get the mates from somatic 3+ mate regions
      // that don't overlap with normal mate region
      MateRegionVector tmp_somatic_mate_regions = collectSomaticMateRegions(wmap, normal_mate_regions);
      
      // no more regions to check
      if (!tmp_somatic_mate_regions.size())
	break;
      
      // keep regions that haven't been visited before
      MateRegionVector somatic_mate_regions;
      for (auto& s : tmp_somatic_mate_regions) {
	
	// check if its not bad from mate region
	if (badd.size() && !opts.noBadAvoid)
	  if (badd.CountOverlaps(s))
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
      logger.log(opt::verbose > 1, log_level > 1, "...mate region ", i.ToString(bwa_header) ,
		 " case read count that triggered lookup: ", i.count, " on mate-lookup round ",
		 jjj+1);
    
    // collect the reads for this round
    pair<int,int> mate_read_counts = collectMateReads(wmap, somatic_mate_regions, jjj, this_bad_mate_regions);

    // update the counts
    counts.first += mate_read_counts.first;
    counts.second += mate_read_counts.second;

    logger.log(opts.verbose > 1, log_level > 1, "\t<case found, control found>: <",
	       AddCommas(counts.first), ",", AddCommas(counts.second), "> on round ",
	       jjj + 1);
    
    if (counts.first + counts.second == 0) // didn't get anything on first round, so quit
      break; 

    } // mate collection round loop
    
    // update this threads tally of bad mate regions
    badd.Concat(this_bad_mate_regions);
  badd.MergeOverlappingIntervals();
  badd.CreateTreeMap();
  logger.log(opts.verbose > 1, log_level > 1, "\tTotal of ", AddCommas(badd.size()), " bad mate regions for this thread");
  
  return counts;
  }

  void correctReads(vector<char*>& learn_seqs,
		    svabaReadVector& brv,
		    SvabaLogger& logger,
		    const SvabaOptions& opts) {

    if (!learn_seqs.size())
      return;
    
    if (opt::ec_correct_type == "s") {
      
      KmerFilter kmer;
      int kmer_corrected = 0;
      
      // make the index for learning correction
      kmer.makeIndex(learn_seqs);
      
      // free the training sequences
      // this memory was alloced w/strdup in collectAndClearReads
      for (size_t i = 0; i < learn_seqs.size(); ++i)
	if (learn_seqs[i]) 
	  free(learn_seqs[i]);
      
      // do the correction
      kmer_corrected = kmer.correctReads(brv); 
      //kmer_corrected = kmer.correctReads(brv); 
      
      logger.log(true/*TODO*/, false /*TODO*/, "...SGA kmer corrected ",
		 kmer_corrected, " reads of ",brv.size());
    } 
  }


void remove_hardclips(svabaReadVector& brv) {
  svabaReadVector bav_tmp;
  for (auto& i : brv)
    if (i.NumHardClip() == 0) 
      bav_tmp.push_back(i);
  brv = bav_tmp;
}

void runAssembly(const GenomicRegion& region,
		 svabaReadVector& input_reads,
		 vector<AlignedContig>& master_alc,
		 BamRecordVector& masterContigs,
		 DiscordantClusterMap& dmap,
		 unordered_map<string, CigarMap>& cigmap,
		 RefGenome* refg) {
  
  // get the local region
  string lregion;
  if (!region.IsEmpty()) {
    try {
      lregion = refg->QueryRegion(bwa_header.IDtoName(region.chr), region.pos1, region.pos2);
    } catch (...) {
      WRITELOG(" Caught exception for lregion with reg " + region.ToString(bwa_header), true, true);
      lregion = "";
    }
  }

  // make a BWA wrapper from the locally retrieved sequence
  UnalignedSequenceVector local_usv = {{"local", lregion, string()}};
  BWAWrapper local_bwa;
  if (local_usv[0].Seq.length() > 200) // have to have pulled some ref sequence
    local_bwa.ConstructIndex(local_usv);

  stringstream region_string;
  region_string << region;
  WRITELOG("...running assemblies for region " + region_string.str(), opt::verbose > 1, false);
  
  // set the contig uid
  string name = "c_" + to_string(region.chr+1) + "_" + to_string(region.pos1) + "_" + to_string(region.pos2);
  
  // where to store contigs
  UnalignedSequenceVector all_contigs_this;
  
  // setup the engine
  svabaAssemblerEngine engine(name, opt::sga::error_rate, opt::sga::minOverlap, readlen);
  engine.fillReadTable(input_reads);
  
  // do the actual assembly
  engine.performAssembly(opt::sga::num_assembly_rounds);
  
  // retrieve contigs
  all_contigs_this = engine.getContigs();
  WRITELOG("...assembled " + to_string(all_contigs_this.size()) + " contigs for " + name, opt::verbose > 1, true);

  // store the aligned contig struct
  vector<AlignedContig> this_alc;
      
  // align the contigs to the genome
  WRITELOG("...aligning " +to_string(all_contigs_this.size()) + " contigs to genome", opt::verbose > 1, false);

  UnalignedSequenceVector usv;

  // loop the unaligned contigs
  for (auto& i : all_contigs_this) {
    
    // if too short, skip
    if ((int)i.Seq.length() < (readlen * 1.2) && !opt::all_contigs)
      continue;
    
    bool hardclip = false;

    //// LOCAL REALIGNMENT
    // align to the local region
    BamRecordVector local_ct_alignments;
    if (!local_bwa.IsEmpty())
      local_bwa.AlignSequence(i.Seq, i.Name, local_ct_alignments, hardclip, SECONDARY_FRAC, SECONDARY_CAP);
    
    // check if it has a non-local alignment
    bool valid_sv = true;
    for (auto& aa : local_ct_alignments) {
      if (aa.NumClip() < MIN_CLIP_FOR_LOCAL) // || aa.GetIntTag("NM") < MAX_NM_FOR_LOCAL)
	valid_sv = false; // has a non-clipped local alignment. can't be SV. Indel only
    }
    ////////////
    
    // do the main realignment of unaligned contigs to the reference genome
    BamRecordVector human_alignments;
    main_bwa->AlignSequence(i.Seq, i.Name, human_alignments, hardclip, SECONDARY_FRAC, SECONDARY_CAP);	

    // store all contig alignment
    masterContigs.insert(masterContigs.end(),human_alignments.begin(), human_alignments.end());
    
    // if very verbose, print
    if (opt::verbose > 3)
      for (auto& ha : human_alignments)
	cerr << " aligned contig: " << ha << endl;
    
    // add in the chrosome name tag for human alignments
    if (main_bwa)
      for (auto& r : human_alignments) {
	assert(main_bwa->ChrIDToName(r.ChrID()).length());
	r.AddZTag("MC", main_bwa->ChrIDToName(r.ChrID()));
	if (!valid_sv)
	  r.AddIntTag("LA", 1); // flag as having a valid local alignment. Can't be SV
      }
    
    // 2023 addition -- sort the contigs by position
    sort(human_alignments.begin(), human_alignments.end());

    // make the aligned contigs
    AlignedContig ac(human_alignments, prefixes);
    
    // assign the local variable to each
    ac.checkLocal(region);
    
    this_alc.push_back(ac);
    usv.push_back({i.Name, i.Seq, string()});	  
  } // end loop through contigs

  assert(this_alc.size() == usv.size());

  // didnt get any contigs that made it all the way through
  if (!this_alc.size())
    return;
  
  // Align the reads to the contigs with BWA-MEM
  BWAWrapper bw;
  bw.ConstructIndex(usv);

  // align the reads to the contigs
  alignReadsToContigs(bw, usv, input_reads, this_alc, refg);

  // align the contigs to the genome
  WRITELOG("...cleaning up and callings bps", opt::verbose > 1, false);
  
  
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

  // parse the command line options
void parseRunOptions(int argc, char** argv) {
  bool die = false;

  if (argc <= 2) 
    die = true;

  bool help = false;
  stringstream ss;

  int sample_number = 0;

  string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'E' : arg >> opt::ec_subsample; break;
    case 'p': arg >> opt::numThreads; break;
    case 'm': arg >> opt::sga::minOverlap; break;
    case 'a': arg >> opt::analysis_id; break;
    case 'B': arg >> opt::blacklist; break;
    case 'V': arg >> opt::germline_sv_file; break;
    case 'L': arg >> opt::mate_lookup_min; break;
    case 'h': help = true; break;
    case OPT_OVERRIDE_REFERENCE_CHECK : opt::override_reference_check = true; break;
    case 'x' : arg >> opt::max_reads_per_assembly; break;
    case 'M' : arg >> opt::mate_region_lookup_limit; break;
    case 'A' : opt::all_contigs = true; break;
      case 'c': 
	tmp = "";
	arg >> tmp;
	if (tmp.find("chr") != string::npos) {
	  opt::chunk = 250000000; break;
	} else {
	  opt::chunk = stoi(tmp); break;
	}
    case OPT_LOD: arg >> opt::lod; break;
    case OPT_LOD_DB: arg >> opt::lod_db; break;
    case OPT_LOD_SOMATIC: arg >> opt::lod_somatic; break;
    case OPT_LOD_SOMATIC_DB: arg >> opt::lod_somatic_db; break;
    case OPT_HP: opt::hp = true; break;
    case OPT_SCALE_ERRORS: arg >> opt::scale_error; break;
    case 'C': arg >> opt::max_cov;  break;
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
    case OPT_NUM_ASSEMBLY_ROUNDS: arg >> opt::sga::num_assembly_rounds; break;
    case 'K': arg >> opt::ec_correct_type; break;
      default: die= true; 
    }
  }

  if (opt::chunk <= 0 || opt::main_bam == "-")
    opt::max_reads_per_assembly = INT_MAX;
  else if (opt::max_reads_per_assembly < 0) 
    opt::max_reads_per_assembly = 50000; //set a default

      

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
    WRITELOG("Invalid number of threads from -p flag: " + AddCommas(opt::numThreads), true, true);
    die = true;
  }

  if (die || help) 
    {
      cerr << "\n" << RUN_USAGE_MESSAGE;
      if (die)
	exit(EXIT_FAILURE);
      else 
	exit(EXIT_SUCCESS);	
    }
}

bool runWorkItem(const GenomicRegion& region,
		 svabaThreadUnit& stu,
		 size_t thread_id) {
  

#ifdef DEBUG
  WRITELOG("Creating a new BFC error corrector on region " + region.ToString(bwa_header) , false, true);
#endif  
  
  // create a new BFC read error corrector for this
  SeqPointer<BFC> bfc;
  if (opt::ec_correct_type == "f") {
    bfc = SeqPointer<BFC>(new BFC());
    for (auto& w : stu.walkers)
      w.second.bfc = bfc;
  }
  
  // setup structures to store the final data for this region
  vector<AlignedContig> alc;
  BamRecordVector all_contigs;

#ifdef DEBUG
  WRITELOG("Starting timer " + region.ToString(bwa_header) , false, true);
#endif  

  // start a timer
  svabaUtils::svabaTimer st;
  st.start();

  // setup for the BAM walkers
  CountPair read_counts = {0,0};

  // loop all of the BAMs (walkers)
  for (auto& w : stu.walkers) {

#ifdef DEBUG
    WRITELOG("Setting region -- Getting reads from BAM with prefix " + w.first + " on region " + region.ToString(bwa_header) , false, true);
#endif  
    
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

#ifdef DEBUG
    WRITELOG("Concatenating reads from BAM with prefix " + w.first + " on region " + region.ToString(bwa_header) , false, true);
#endif  
    
    // do the reading, and store the bad mate regions
    stu.badd.Concat(w.second.readBam(&log_file));

#ifdef DEBUG
    WRITELOG("...merging Overlapping Intervals from BAM with prefix " + w.first + " on region " + region.ToString(bwa_header) , false, true);
#endif  
    
    stu.badd.MergeOverlappingIntervals();
    
#ifdef DEBUG
    WRITELOG("...creating tree map from BAM with prefix " + w.first + " on region " + region.ToString(bwa_header) , false, true);
#endif  
    
    stu.badd.CreateTreeMap();
    
    // adjust the counts
    if (w.first.at(0) == 't') {
      read_counts.first += w.second.reads.size();
    } else {
      read_counts.second += w.second.reads.size();
    }

  }

#ifdef DEBUG
    WRITELOG("...collecting cigar strings on region " + region.ToString(bwa_header) , false, true);
#endif  
  
  // collect all of the cigar strings in a hash
  unordered_map<string, CigarMap> cigmap;
  for (const auto& w : stu.walkers) 
    cigmap[w.first] = w.second.cigmap;

  // setup read collectors
  vector<char*> all_seqs;
  //BamRecordVector input_reads;
  svabaReadVector input_reads;

#ifdef DEBUG
    WRITELOG("...collecting and clearing reads " + region.ToString(bwa_header) , false, true);
#endif  
  
  // collect and clear reads from main round
  unordered_set<string> dedupe;
  collectAndClearReads(stu.walkers, input_reads, all_seqs, dedupe);

  // adjust counts and timer
  st.stop("r");

  // get the mate reads, if this is local assembly and has insert-size distro
  if (!region.IsEmpty() && !opt::single_end && min_dscrd_size_for_variant) {
#ifdef DEBUG
    WRITELOG("...running mate collection loops " + region.ToString(bwa_header) , false, true);
#endif  
    runMateCollectionLoop(region, stu.walkers, stu.badd, logger);
    // collect the reads together from the mate walkers
    collectAndClearReads(stu.walkers, input_reads, all_seqs, dedupe);
    st.stop("m");
  }
  
  // do the discordant read clustering
  DiscordantClusterMap dmap, dmap_tmp;
  
  // if couldn't get insert size stats, skip discordant clustering
  if (!min_dscrd_size_for_variant || opt::single_end)
    goto afterdiscclustering;

  WRITELOG("...discordant read clustering", opt::verbose > 1, false);
  dmap = DiscordantCluster::clusterReads(input_reads, region, max_mapq_possible, &min_isize_for_disc);

  // tag FR clusters that are below min_dscrd_size_for_variant AND low support
  for (auto& d : dmap) {
    bool below_size = 	d.second.m_reg1.strand == '+' && d.second.m_reg2.strand == '-' && 
      (d.second.m_reg2.pos1 - d.second.m_reg1.pos2) < min_dscrd_size_for_variant && 
      d.second.m_reg1.chr == d.second.m_reg2.chr;

    // low support and low size, completely ditch it
    if (below_size && (d.second.tcount + d.second.ncount) < 4)
      continue;
    else
      dmap_tmp.insert(pair<string, DiscordantCluster>(d.first, d.second));
  }
  dmap = dmap_tmp;

  // print out results
  if (opt::verbose > 3)
    for (auto& i : dmap) 
      WRITELOG(i.first + " " + i.second.toFileString(bwa_header, false), true, false);

 afterdiscclustering:

  // skip all the assembly stuff?
  if (opt::disc_cluster_only)
    goto afterassembly;
  
  /////////////////////
  //// ASSEMBLY 
  /////////////////////
#ifdef DEBUG
    WRITELOG("Removing hardclips " + region.ToString(bwa_header) , false, true);
#endif    

  // remove the hardclips, don't assemble them
  remove_hardclips(input_reads);

  if (opt::disc_cluster_only)
    goto afterassembly;

  // check that we don't have too many reads
  if (input_reads.size() > (size_t)(region.Width() * 20) && region.Width() > 20000) {
    stringstream ssss;
    WRITELOG("TOO MANY READS IN REGION " + AddCommas(input_reads.size()) + "\t" + region.ToString(bwa_header), opt::verbose, false);
    goto afterassembly;
  }

  // print message about assemblies
  if (input_reads.size() > 1) {
    WRITELOG("...assembling on " + region.ToString(bwa_header), opt::verbose > 1, false);
  } else if (input_reads.size() < 3) { 
    WRITELOG("Skipping assembly (<= 2 reads) on " + region.ToString(bwa_header), opt::verbose > 1, false);
    goto afterassembly;
  }

#ifdef DEBUG
    WRITELOG("Doing kmer correction " + region.ToString(bwa_header) , false, true);
#endif    

  // do the kmer correction, in place
  if (opt::ec_correct_type == "s") {
    correctReads(all_seqs, input_reads);
  } else if (opt::ec_correct_type == "f" && input_reads.size() >= 8) {
    
    assert(bfc);
    int learn_reads_count = bfc->NumSequences();
    
#ifdef DEBUG
    WRITELOG("BFC Train " + region.ToString(bwa_header) , false, true);
#endif    
    
    bfc->Train();    
    bfc->clear();  // clear memory and reads. Keeps training data
    
#ifdef DEBUG
    WRITELOG("BFC training " + region.ToString(bwa_header) , false, true);
#endif    
    
    st.stop("t");

    // reload with the reads to be corrected
    for (auto& s : input_reads)
      bfc->AddSequence(s.Seq().c_str(), "", "");

#ifdef DEBUG
    WRITELOG("BFC Error correcting " + region.ToString(bwa_header) , false, true);
#endif    
    
    // error correct 
    bfc->ErrorCorrect();

    // retrieve the sequences
    string s, name_dum;
    for (auto& r : input_reads) {
      assert(bfc->GetSequence(s, name_dum));
      r.SetSeq(s);
    }
    
    double kcov = bfc->GetKCov();
    int kmer    = bfc->GetKMer();

    WRITELOG("...BFC attempted correct " + to_string(input_reads.size()) + " train: " + 
	     to_string(learn_reads_count) + " kcov: " + to_string(kcov) + 
	     " kmer: " + to_string(kmer), opt::verbose > 1, true);      

  }

  st.stop("k");
  
  // do the assembly, contig realignment, contig local realignment, and read realignment
  // modifes input_reads, alc, all_contigs
#ifdef DEBUG
    WRITELOG("Running assembly " + region.ToString(bwa_header) , false, true);
#endif    
  
  runAssembly(region, input_reads, alc,
	      all_contigs,
	      dmap, cigmap, stu.ref_genome.get());

afterassembly:

  // clear it out, not needed anymore
  if (bfc)
    bfc->clear();
  
  st.stop("as");
  WRITELOG("...done assembling, post processing", opt::verbose > 1, false);

  // get the breakpoints
  vector<BreakPoint> bp_glob;
  
  for (auto& i : alc) {
    vector<BreakPoint> allbreaks = i.getAllBreakPoints();
    bp_glob.insert(bp_glob.end(), allbreaks.begin(), allbreaks.end());
  }

  if (dbsnp_filter && opt::dbsnp.length()) {
    WRITELOG("...DBSNP filtering", opt::verbose > 1, false);
    for (auto & i : bp_glob) 
      dbsnp_filter->queryBreakpoint(i);
  }

#ifdef DEBUG
    WRITELOG("filtering against blacklist " + region.ToString(bwa_header) , false, true);
#endif    

  
  // filter against blacklist
  for (auto& i : bp_glob) 
    i.checkBlacklist(blacklist);

#ifdef DEBUG
    WRITELOG("Add in discordant clusters " + region.ToString(bwa_header) , false, true);
#endif    
  
  // add in the discordant clusters as breakpoints
  for (auto& i : dmap) {
    // dont send DSCRD if FR and below size
    bool below_size = 	i.second.m_reg1.strand == '+' && i.second.m_reg2.strand == '-' && 
      (i.second.m_reg2.pos1 - i.second.m_reg1.pos2) < min_dscrd_size_for_variant && 
      i.second.m_reg1.chr == i.second.m_reg2.chr;
    // DiscordantCluster not associated with assembly BP and has 2+ read support
    if (!i.second.hasAssociatedAssemblyContig() && 
	(i.second.tcount + i.second.ncount) >= MIN_DSCRD_READS_DSCRD_ONLY && i.second.valid() && !below_size) {
      BreakPoint tmpbp(i.second, main_bwa, dmap, region, bwa_header);
      bp_glob.push_back(tmpbp);
    }
  }

#ifdef DEBUG
    WRITELOG("Deduplicate breakpoints " + region.ToString(bwa_header) , false, true);
#endif    
  
    // de duplicate the breakpoints
  sort(bp_glob.begin(), bp_glob.end());
  bp_glob.erase( unique( bp_glob.begin(), bp_glob.end() ), bp_glob.end() );

#ifdef DEBUG
    WRITELOG("Adding coverage data " + region.ToString(bwa_header) , false, true);
#endif    
  
  // add the coverage data to breaks for allelic fraction computation
  unordered_map<string, STCoverage*> covs;
  for (auto& i : opt::bam) 
    covs[i.first] = &stu.walkers[i.first].cov;

  for (auto& i : bp_glob)
    i.addCovs(covs);
  
#ifdef DEBUG
    WRITELOG("Scoring breakpoints " + region.ToString(bwa_header) , false, true);
#endif    
  
  for (auto& i : bp_glob) {
    i.readlen = readlen; // set the readlength
    i.scoreBreakpoint(opt::lod, opt::lod_db, opt::lod_somatic, opt::lod_somatic_db, opt::scale_error, min_dscrd_size_for_variant);
  }

#ifdef DEBUG
    WRITELOG("Labeling somatic breakpoints " + region.ToString(bwa_header) , false, true);
#endif    

  
  // label somatic breakpoints that intersect directly with normal as NOT somatic
  unordered_set<string> norm_hash;
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
  unordered_map<string, size_t> ccc;
  for (auto& i : bp_glob) {
    if (i.evidence == "INDEL" && i.repeat_seq.length() > 6) {
      ++ccc[i.b1.hash()];
    }
  }
  for (auto& i : bp_glob) {
    if (i.evidence == "INDEL" && ccc[i.b1.hash()] > 1)
      i.confidence = "REPVAR";
  }

  // remove somatic calls if they have a germline normal SV in them or indels with 
  // 2+germline normal in same contig
  unordered_set<string> bp_hash;
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
	GenomicRegion gr1 = i.b1.gr;
	GenomicRegion gr2 = i.b2.gr;
	gr1.Pad(GERMLINE_CNV_PAD);
	gr2.Pad(GERMLINE_CNV_PAD);
	if (germline_svs.OverlapSameInterval(gr1, gr2)) {
	  i.somatic_score = -1;
	}
      }
    }
      
  }

#ifdef DEBUG
    WRITELOG("Set ref and alt tags  " + region.ToString(bwa_header) , false, true);
#endif    
  
  // add the ref and alt tags
  // WHY IS THIS NOT THREAD SAFE?
  for (auto& i : bp_glob)
    i.setRefAlt(stu.ref_genome.get()); 

  // transfer local versions to thread store
  for (const auto& a : alc)
    if (a.hasVariant()) {
      stu.m_alc.push_back(a);
      stu.m_bamreads_count += a.NumBamReads();
    }
  for (const auto& d : dmap)
    stu.m_disc_reads += d.second.reads.size();
  
  stu.m_contigs.insert(stu.m_contigs.end(), all_contigs.begin(), all_contigs.end());
  stu.m_disc.insert(dmap.begin(), dmap.end());
  for (const auto& a : alc)
    stu.m_bamreads_count += a.NumBamReads();
  for (auto& i : bp_glob) 
    if ( i.hasMinimal() && (i.confidence != "NOLOCAL" || i.complex_local ) ) 
      stu.m_bps.push_back(i);
  
  // dump if getting to much memory
  if (stu.MemoryLimit(THREAD_READ_LIMIT, THREAD_CONTIG_LIMIT) && !opt::hp) {
    WRITELOG("...writing output files on thread " + to_string(thread_id) +
	     " with limit hit of " + to_string(stu.m_bamreads_count), opt::verbose > 1, true);
    output_writer.writeUnit(stu); // mutex and flush are inside this call
  }
   
  st.stop("pp");
  
  // display the run time
  WRITELOG(svabaUtils::runTimeString(read_counts.first, read_counts.second, alc.size(), region, bwa_header, st, start), opt::verbose > 1, true);

  // clear out the reads and reset the walkers
  for (auto& w : stu.walkers) {
    w.second.clear(); 
    w.second.m_limit = opt::max_reads_per_assembly;
  }

  return true;
}

  // --- one small helper functor ---
  struct SvabaWorkItem {
    
    SeqLib::GenomicRegion region;
    int                   id;
    
    SvabaWorkItem(const SeqLib::GenomicRegion& r, int n)
      : region(r), id(n)
    {}

    /// invoked by a worker thread
    bool operator()(svabaThreadUnit& unit, size_t threadId) const {
      return runWorkItem(region, unit, threadId);
    }

  };
  

  void sendThreads(const SeqLib::GRC& regionsToRun,
		   SvabaLogger& logger,
		   const SvabaOpts& opts,
		   const BamHeader& header,
		   SvabaOutputWriter& writer,
		   BWAIndexPtr idx) {
    
    // single shared processor:
    SvabaRegionProcessor proc(logger, opts, header, writer);
    
    // our threadpool of work items:
    ThreadPool<SvabaWorkItem> pool(
				   opts.numThreads,
				   opts.refGenome,
				   opts.bamFiles,
				   idx,
				   writer
				   );
    
  int count=0;
  // submit one job per region  
  for (auto& r : regionsToRun)
    pool.submit(std::make_unique<SvabaWorkItem>(r, ++count, proc));
  // if no intervals, submit a single wholegenome job  
  if (regionsToRun.IsEmpty())
    pool.submit(std::make_unique<SvabaWorkItem>(SeqLib::GenomicRegion(), ++count, proc));

  pool.shutdown();
}
  
void makeVCFs() {

  if (opt::bam.size() == 0) {
    cerr << "makeVCFs error: must supply a BAM via -t to get header from" << endl;
    exit(EXIT_FAILURE);
  }

  // make the VCF file
  WRITELOG("...loading the bps files for conversion to VCF", opt::verbose, true);

  string file = opt::analysis_id + ".bps.txt.gz";
  //if !read_access_test(file))
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
    string fname = b.second; //bpf.filename();
    header.addSampleField(fname);
    header.colnames += "\t" + fname; 
  }


  // check if it has a matched control. If so, output "somatic / germline" vcfs
  bool case_control_run = false;
  for (auto& b : opt::bam)
    if (b.first.at(0) == 'n')
      case_control_run = true;

  // primary VCFs
  if (read_access_test(file)) {
    WRITELOG("...making the primary VCFs (unfiltered and filtered) from file " + file, opt::verbose, true);
    VCFFile snowvcf(file, opt::analysis_id, bwa_header, header, !false,
		    opt::verbose > 0);

    string basename = opt::analysis_id + ".svaba.unfiltered.";
    snowvcf.include_nonpass = true;
    WRITELOG("...writing unfiltered VCFs", opt::verbose, true);
    snowvcf.writeIndels(basename, false, !case_control_run);
    snowvcf.writeSVs(basename, false,    !case_control_run);

    WRITELOG("...writing filtered VCFs", opt::verbose, true);
    basename = opt::analysis_id + ".svaba.";
    snowvcf.include_nonpass = false;
    snowvcf.writeIndels(basename, false, !case_control_run);
    snowvcf.writeSVs(basename, false,    !case_control_run);

  } else {
    WRITELOG("ERROR: Failed to make VCF. Could not file bps file " + file, true, true);
  }

}  
  
} //anonymous namespace

void runsvaba(int argc, char** argv) {

  // shared information to be passed around thread units
  SvabaShardConfig sc;

  // parse command line
  try {
    sc.opts = SvabaOptions::parse(argc, argv);
    if (opts.help) {
      SvabaOptions::printUsage();
      return;
    }
  }
  catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return;
  }

  // open message logger
  sc.logger.init(opts.analysis_id + ".log");
  
  // write to log
  sc.logger.log(true, true, 
	     "-----------------------------------------------------------\n",
	     "---  Running svaba SV and indel detection on ",
	     " threads ---", (sc.opts.numThreads >= 10 ? "\n" : "-\n"),
	     "---  Version: ", SVABA_VERSION, " - ", SVABA_DATE, "                           ---\n",
	     "---    (inspect *.log for real-time progress updates)   ---\n",
	     "-----------------------------------------------------------");
  
  // print the options
  sc.opts.printToLogger(sc.logger);

  // open the human reference
  sc.logger.log(true, true, "...loading the human reference sequence for BWA");
  sc.bwa_idx = std::make_shared<SeqLib::BWAIndex>();
  sc.bwa_idx->loadIndex(sc.opt.refGenome);
  
  /*  main_bwa = new BWAWrapper();
  main_bwa->SetAScore(opts.bwa.sequence_match_score);
  main_bwa->SetGapOpen(opts.bwa.gap_open_penalty);
  main_bwa->SetGapExtension(opts.bwa.gap_extension_penalty);
  main_bwa->SetMismatchPenalty(opts.bwa.mismatch_penalty);
  main_bwa->SetZDropoff(opts.bwa.zdrop);
  main_bwa->SetBandwidth(opts.bwa.bandwidth);
  main_bwa->SetReseedTrigger(opts.bwa.reseed_trigger);
  main_bwa->Set3primeClippingPenalty(opts.bwa.clip3_pen);
  main_bwa->Set5primeClippingPenalty(opts.bwa.clip5_pen);
  */
  
  // load main BWA aligner reference genome
  if (!sc.bwa_index->LoadIndex(opts.refgenome)) {
    sc.logger.log(true, true, "Failed to load bwa index ",sc.opts.refgenome,"May need to run bwa index on reference fasta");
    return;
  }
  
  // get the dictionary from reference
  sc.bwa_header = sc.bwa_idx->HeaderFromIndex();
  
  // load the same index, but for querying seq from ref
  ref_genome = new RefGenome;  
  if (!ref_genome->LoadIndex(opts.refgenome) || ref_genome->IsEmpty()) {
    logger.log("Failed to load samtools index ",opts.refgenome, "May need to run samtools faidx on reference fasta");
    return;
  }

  // open the main bam to get header info
  SeqLib::BamReader first_tumor_bam_reader;  
  if (!first_tumor_bam_reader.Open(sc.opts.main_bam)) {
    if (sc.opts.main_bam == "-")
      cerr << "ERROR: Cannot read from stdin" << endl;
    else
      cerr << "ERROR: Cannot open main bam file: " << opt::main_bam << endl;
    exit(EXIT_FAILURE);
  }
  SeqLib::BamHeader first_header = first_tumor_bam_reader.Header();
  
  // then open the main header
  if (sc.header.isEmpty()) {
    cerr << "ERROR: empty header in main bam file" << endl;
    exit(EXIT_FAILURE);
  }
  
  // check that the two headers are equivalant
  if (sc.opts.override_reference_check) {
    sc.logger.log(true, true, "!!! Will NOT perform check of reference compatability with BAM");
    sc.logger.log(true, true, "!!! Only if really *sure* that reference and BAM have same chr in same order");
  } else {
    svabaUtils::checkHeaderCompatibility(first_header, sc.header, sc.logger);
  }

  // open file loader
  SvabaFileLoader loader(sc.logger, sc.opts, sc.header);
  
  // open the blacklist
  loader.loadBedRegions(sc.opts.blacklist, sc.blacklist);
  
  // open the germline sv database
  loader.loadBedRegions(opts.germline_sv_file, germline_svs);
  
  // open the DBSnpFilter
  if (opts.dbsnp.length()) {
    dbsnp_filter = new DBSnpFilter(sc.opts.dbsnp, sc.header);
  }

  // needed for aligned contig
  for (auto& b : opts.bam)
    prefixes.insert(b.first);

  // parse the region file, count number of jobs
  SeqLib::GRC regions_torun;
  loader.countJobs(file_regions, regions_torun);  

  // --- learn the insert-sizes ---
  sc.logger.log(true, true,"... learning insert size distribution across all BAMs; this may take a while");
  
  // learn from the BAM files
  BamLearningResult bamStats = LearnBamParams::learnAll(
    opts.bamFiles,         // map<sampleName,bamPath>
    PER_RG_LEARN_LIMIT,
    opts.sdDiscCutoff,
    logger
  );

  logger.log(true, true,"... done learning insert size distribution");

  // --- set the SGA min overlap if user didn't ---
  if (opts.sga.minOverlap == 0) {
    opts.sga.minOverlap = std::max(30, int(0.6 * bamStats.globalReadLen));
  }

  // --- report what we learned ---
  logger.log(true, true,
    "... found max read length = ", bamStats.globalReadLen,
    "; SGA minOverlap = ", opts.sga.minOverlap,
    "; max MAPQ = ", bamStats.globalMaxMapQ
  );

  // --- compute seed parameters ---
  int seedLength, seedStride;
  {
    svabaAssemblerEngine tester(
      "test", opts.sga.error_rate, opts.sga.minOverlap, readlen
    );
    tester.calculateSeedParameters(
      readlen, opts.sga.minOverlap, seedLength, seedStride
    );
  }
  logger.log(false, true,
    "... calculated seedLength = ", seedLength,
    " (error_rate=", opts.sga.error_rate,
    ", readlen=", bamStats.globalReadLen, ")"
  );

  // --- build per-RG discordant size rules ---
  std::stringstream ssRules;
  std::unordered_set<std::string> seenRG;
  std::map<std::string,int> minIsizeForDisc;
  for (auto const& [sample, pm] : bamStats.perFile) {
    for (auto const& [rg, bp] : pm) {
      if (!seenRG.count(rg)) {
        int cutoff = int(
          std::floor(bp.mean_isize + bp.sd_isize * opts.sd_disc_cutoff)
        );
        ssRules << R"({"isize":[)" << cutoff << R"(,0],"rg":")"
                << rg << R"("})" << ",";
        seenRG.insert(rg);
        minIsizeForDisc[rg] = cutoff;
      }
    }
  }
  if (!ssRules.str().empty())
    ssRules.seekp(-1, ssRules.cur);  // drop trailing comma

  logger.log(false, true, "[INFO] Learned discordant size cutoffs by RG:");
  for (auto const& kv : minIsizeForDisc) {
    logger.log(false, true,
      "   RG=", kv.first, " -> cutoff=", kv.second
    );
  }

  // plug into your JSON rules template
  opts.rules = svabaUtils::myreplace(
    opts.rules, "FRRULES", ssRules.str().empty()
      ? "{}"
      : ssRules.str()
  );

  // similarly replace any READLENLIM token
  if (opts.rules.find("READLENLIM") != std::string::npos) {
    int clipLen = std::max(30, int(0.3 * readlen));
    opts.rules = svabaUtils::myreplace(
      opts.rules, "READLENLIM", std::to_string(clipLen)
    );
  }
  
  /*  
  // learn the insert sizes
#ifndef TCGA_ASSUME
  logger.log(true,true,
	     "...learning insert size distribution across all BAMs; this may take a while");
  std::map<std::string,BamParamsMap> paramsMap;
  int maxReadLen, maxMapQ, minDiscordantSize;
  LearnBamParams::learnAll(
    opts.bamFiles,
    opts.sdDiscCutoff,
    logger,
    paramsMap,
    maxReadLen,
    maxMapQ,
    minDiscordantSize
  );
  logger.log(false,true,
    "...done learning insert size distribution");
#else
  logger.log(true,true,
    "...skipping insert size learning for speed, assuming cutoff=2000");
  int minDiscordantSize = 2000;
#endif
  
  // set the min overlap
  if (!opts.sga.minOverlap) 
    opts.sga.minOverlap = (0.6 * readlen) < 30 ? 30 : 0.6 * readlen;

  ss << "...found read length of " << readlen << ". Min Overlap is " << opt::sga::minOverlap << endl;
  ss << "...max read MAPQ detected: " << max_mapq_possible << endl;
  WRITELOG(ss.str(), opt::verbose > 1, true);
  ss.str(string());

  // if didn't learn anything, then make sure we set an overlap
  if (!readlen && !opt::sga::minOverlap) {
    cerr << "ERROR: Didn't learn from reads (stdin assembly?). Need to explicitly set readlen with --readlen" << endl;
    exit(EXIT_FAILURE);
  }
  
  // get the seed length for printing
  int seedLength, seedStride;
  svabaAssemblerEngine enginetest("test", opt::sga::error_rate, opt::sga::minOverlap, readlen);
  enginetest.calculateSeedParameters(readlen, opt::sga::minOverlap, seedLength, seedStride);
  WRITELOG("...calculated seed size for error rate of " + to_string(opt::sga::error_rate) + " and read length " +
	   to_string(readlen) + " is " + to_string(seedLength), opt::verbose, true);

  
  // loop through and construct the readgroup rules
  stringstream ss_rules;
  unordered_set<string> rg_seen;

  // set the insert size bounds for each read group
  for (auto& a : params_map) {
    for (auto& i : a.second) {
      if (!rg_seen.count(i.second.read_group)) {
	int mi = floor(i.second.mean_isize + i.second.sd_isize * opt::sd_disc_cutoff);
	ss_rules << "{\"isize\" : [ " << mi << ",0], \"rg\" : \"" << i.second.read_group << "\"},";
	rg_seen.insert(i.second.read_group);
	min_isize_for_disc.insert(pair<string, int>(i.second.read_group, mi));
      }
    } 
  }

  // print learned isize cutoffs
  ss << "[INFO] Learned isize cutoffs by read group: " << endl;
  for (const auto& kv : min_isize_for_disc) {
    ss << "  RG: \"" << kv.first << "\" -> cutoff: " << kv.second << endl;
  }
  WRITELOG(ss.str(), opt::verbose, true);

  // format the rules JSON from the above string
  if (opt::rules.find("FRRULES") != string::npos) {
    string string_rules = ss_rules.str();
    if (!string_rules.empty()) // cut last comma
      string_rules = string_rules.substr(0, string_rules.length() - 1);
    else
      string_rules = "{}";
    opt::rules = svabaUtils::myreplace(opt::rules, "FRRULES", string_rules);
  }

  // set min length for clips
  if (opt::rules.find("READLENLIM") != string::npos) {
    if (readlen == 0)
      readlen = 30; // set some small default, in case we didn't learn the bam
    opt::rules = svabaUtils::myreplace(opt::rules, "READLENLIM", to_string((int) (readlen * 0.3)));
  }
  */

  // set the ReadFilterCollection to be applied to each region
  logger.log(opts.verbose > 1, true, opts.rules);
  mr = new Filter::ReadFilterCollection(opts.rules, bwa_header);
  logger.log(opts.verbose > 1, true, *mr);

  // put args into string for VCF later
  args += "(v" + string(SVABA_VERSION) + ") ";
  for (int i = 0; i < argc; ++i)
    args += string(argv[i]) + " ";

  // start the timer
#ifndef __APPLE__
  clock_gettime(CLOCK_MONOTONIC, &start);
#endif

  // send the jobs to the queue
  logger.log(true, true, "Starting detection pipeline");
  sendThreads(regions_torun, logger, opts, bwa_header, output_writer, bwaIndex);

  // more clean up 
  if (ref_genome)
    delete ref_genome;
  
  // make the VCF file
  makeVCFs();
  
#ifndef __APPLE__
    cerr << displayRuntime(start) << endl;
#endif
}


