 #include "benchmark.h"

 #include <getopt.h>
 #include <string>
 #include <sstream>
 #include <fstream>
 #include <iostream>
 #include <cstdlib>
 #include <algorithm>

 #include "vcf.h"
 #include "Fractions.h"
 #include "SeqLib/BWAWrapper.h"
 #include "SeqLib/SeqLibCommon.h"
 #include "SeqLib/GenomicRegion.h"

 #include "ReadSim.h"
 #include "SeqFrag.h"
 #include "SimGenome.h"
 #include "KmerFilter.h"
 #include "PowerLawSim.h"
 #include "BamSplitter.h"
 #include "SnowmanUtils.h"
 #include "AlignedContig2.h"
 #include "SimTrainerWalker.h"
 #include "SnowmanAssemblerEngine.h"


 static std::vector<double> snv_error_rates;
 static std::vector<double> del_error_rates;
 static std::vector<double> ins_error_rates;
 static std::vector<double> coverages;
 static std::vector<double> fractions;

 static SeqLib::GRC regions;
 static Fractions fractions_bed;
 static SeqLib::BamReader bwalker;
 static faidx_t * findex;

 #define DEFAULT_SNV_RATE 0.01
 #define DEFAULT_DEL_RATE 0.05
 #define DEFAULT_INS_RATE 0.05
 #define DEFAULT_COV 10

 namespace opt {

   static std::string refgenome = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
   static int mode = -1;
   static size_t readlen = 101;
   static int num_runs = 100;
   static uint32_t seed = 0;
   static std::string regionFile = "";
   static std::string bam = "";
   static int isize_mean = 250;
   static int isize_sd = 50;

   static int nbreaks = 10;
   static int nindels = 10;

   static std::string string_id = "noid";

   static std::string frac_bed_file;

   static bool scramble = false;

   static int viral_count = 0;

   static std::string blacklist;

   static std::set<std::string> prefixes;
 }

 enum { 
   OPT_ASSEMBLY,
   OPT_SIMBREAKS,
   OPT_ISIZE_MEAN,
   OPT_ISIZE_SD,
   OPT_SPLITBAM,
   OPT_SCRAMBLE,
   OPT_POWERSIM,
   OPT_REALIGN,
   OPT_BLACKLIST, 
   OPT_REALIGN_SV
 };


 static const char *BENCHMARK_USAGE_MESSAGE =
 "Usage: snowman benchmark\n\n"
 "  Description: Various benchmarking tests for Snowman\n"
 "\n"
 "  General options\n"
 "  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
 "  -G, --reference-genome               Indexed ref genome for BWA-MEM. Default (Broad): /seq/reference/...)\n"
 "  -s, --seed                           Seed for the random number generator\n"
 "  -A, --string-id                      String to name output files with (e.g. <string-id>_0_01.bam\n"
 "  Choose one of the following:\n"     
 "      --test-assembly                  Generate single-end reads from small contigs to test assembly/remapping\n"
 "      --sim-breaks-power               Simulate rearrangements and indels and output paired-end reads\n"
 "      --split-bam                      Divide up a BAM file into smaller sub-sampled files, with no read overlaps between files. Preserves read-pairs\n"
 "      --realign-test                   Randomly sample the reference genome and test ability of BWA-MEM to realign to reference for different sizes / error rates\n"
 "      --realign-sv-test                Make an SV (rearrangment) and simulate contigs from it. Test size of contigs vs alignment accuracy for SVs\n"
 "  Shared Options for Test and Simulate:\n"
 "  -c, --read-covearge                  Desired coverage. Input as comma-separated to test multiple (test assembly)\n"
 "  -b, --bam                            BAM file to train the simulation with\n"
 "  -k, --regions                        Regions to simulate breaks or test assembly\n"
 "  -E, --snv-error-rate                 The random SNV error rate per base. Input as comma-separated to test multiple (test assembly)\n"
 "  -I, --ins-error-rate                 The random insertion error rate per read. Input as comma-separated to test multiple (test assembly)\n"
 "  -D, --del-error-rate                 The random deletion error rate per read. Input as comma-separated to test multiple (test assembly)\n"
 "  Test Assembly (--test-assembly) Options:\n"
 "  -n, --num-runs                       Number of random trials to run\n"
   // "  Simulate Breaks (--sim-breaks)  Options:\n"
   //"      --isize-mean                     Desired mean insert size for the simulated reads\n"
   //"      --isize-sd                       Desired std. dev. forinsert size for the simulated reads\n"
   //"  -M, --viral-integration              Number of segments to integrate viruses into. (Reduces amount of space for indels)\n"
   //"      --add-scrambled-inserts          Add scrambled inserts at junctions, randomly between 0 and 100 bp with p(0) = 50% and p(1-100) = 50%;\n"
 "  Simulate Breaks with Power Law(--sim-breaks-power)  Options:\n"
 "      --blacklist                      BED file specifying blacklist regions not to put breaks in\n"
 "  -R, --num-rearrangements             Number of rearrangements to simulate\n"
 "  -X, --num-indels                     Number of indels to simulate\n"
 "  Split Bam (--split-bam)  Options:\n"
 "  -f, --fractions                      Fractions to split the bam into\n"
 "\n";


 static const char* shortopts = "haG:c:n:s:k:b:E:I:D:R:X:A:f:M:";
 static const struct option longopts[] = {
   { "help",                 no_argument, NULL, 'h' },
   { "reference-genome",     required_argument, NULL, 'G' },
   { "string-id",            required_argument, NULL, 'A' },
   { "seed",                 required_argument, NULL, 's' },
   { "num-runs",             required_argument, NULL, 'n' },
   { "regions",              required_argument, NULL, 'k' },
   { "bam",                  required_argument, NULL, 'b' },
   { "read-coverage",        required_argument, NULL, 'c' },
   { "snv-error-rate",       required_argument, NULL, 'E' },
   { "del-error-rate",       required_argument, NULL, 'D' },
   { "ins-error-rate",       required_argument, NULL, 'I' },
   { "num-rearrangements",   required_argument, NULL, 'R' },
   { "fractions",            required_argument, NULL, 'f' },
   { "num-indels",           required_argument, NULL, 'X' },
   { "viral-integration",    required_argument, NULL, 'M' },
   { "isize-mean",           required_argument, NULL, OPT_ISIZE_MEAN},
   { "isize-sd",             required_argument, NULL, OPT_ISIZE_SD},
   { "test-assembly",        no_argument, NULL, OPT_ASSEMBLY},
   { "blacklist",            required_argument, NULL, OPT_BLACKLIST},
   { "sim-breaks-power",     no_argument, NULL, OPT_POWERSIM},
   { "sim-breaks",           no_argument, NULL, OPT_SIMBREAKS},
   { "split-bam",            no_argument, NULL, OPT_SPLITBAM},
   { "realign-sv-test",      no_argument, NULL, OPT_REALIGN_SV},
   { "add-scrambled-inserts",no_argument, NULL, OPT_SCRAMBLE},
   { "realign-test",no_argument, NULL, OPT_REALIGN},
   { NULL, 0, NULL, 0 }
 };

 void runBenchmark(int argc, char** argv) {

   parseBenchmarkOptions(argc, argv);

   opt::prefixes.insert(opt::bam);

   std::cerr << 
     "-----------------------------------------" << std::endl << 
     "--- Running Snowman Benchmarking Test ---" << std::endl <<
     "-----------------------------------------" << std::endl;
   if (opt::mode == OPT_ASSEMBLY)
     std::cerr << "********* RUNNING ASSEMBLY TEST ***********" << std::endl;
   else if (opt::mode == OPT_SIMBREAKS)
     std::cerr << "********* RUNNING SIMULATE BREAKS ***********" << std::endl;
   else if (opt::mode == OPT_SPLITBAM)
     std::cerr << "********* RUNNING SPLIT BAM ***********" << std::endl;
   else if (opt::mode == OPT_REALIGN)
     std::cerr << "********* RUNNING REALIGN TEST ***********" << std::endl;    
   else if (opt::mode == OPT_REALIGN_SV)
     std::cerr << "********* RUNNING REALIGN SV TEST ***********" << std::endl;    

   if (opt::mode == OPT_ASSEMBLY || opt::mode == OPT_SIMBREAKS) {
     std::cerr << "    Error rates:" << std::endl;
     std::cerr << errorRateString(snv_error_rates, "SNV") << std::endl;
     std::cerr << errorRateString(ins_error_rates, "Del") << std::endl;
     std::cerr << errorRateString(del_error_rates, "Ins") << std::endl;
     std::cerr << errorRateString(coverages, "Coverages") << std::endl;
     std::cerr << "    Insert size: " << opt::isize_mean << "(" << opt::isize_sd << ")" << std::endl;
   } else if (opt::mode == OPT_SPLITBAM) {
     std::cerr << errorRateString(fractions, "Fractions") << std::endl;
   }

   // open the BAM
   if (opt::bam.length() && SeqLib::read_access_test(opt::bam)) {
     assert(bwalker.Open(opt::bam));
   }
   else if (opt::mode == OPT_REALIGN || opt::mode == OPT_SPLITBAM || opt::mode == OPT_REALIGN_SV) {
     std::cerr << "NEED TO INPUT VALID BAM (perhaps just to get header info)" << std::endl;
     exit(EXIT_FAILURE);
  }
    // parse the region file
  if (opt::regionFile.length()) {
    if (SeqLib::read_access_test(opt::regionFile)) {
      regions = SeqLib::GRC(opt::regionFile, bwalker.Header());
      regions.MergeOverlappingIntervals();
    }
    // samtools format
    else if (opt::regionFile.find(":") != std::string::npos && opt::regionFile.find("-") != std::string::npos) {
      if (bwalker.Header().isEmpty()) {
	std::cerr << "Error: To parse a samtools style string, need a BAM header. Input bam with -b" << std::endl;
	exit(EXIT_FAILURE);
      }
      regions.add(SeqLib::GenomicRegion(opt::regionFile, bwalker.Header()));    
    } else {
      std::cerr << "Can't parse the regions. Input as BED file or Samtools style string (requires BAM with -b to for header info)" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (!regions.size()) {
      std::cerr << "ERROR: Must input a region to run on " << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // seed the RNG
  if (opt::mode != OPT_SPLITBAM) {
    if (opt::seed == 0)
      opt::seed = (unsigned)time(NULL);
    srand(opt::seed);
    std::cerr << "   Seed: " << opt::seed << std::endl;
  }

  // read the fractions file
  if (opt::frac_bed_file.length() && opt::mode == OPT_SPLITBAM) {
    fractions_bed.readFromBed(opt::frac_bed_file, bwalker.Header());
    //  fractions_bed.readBEDfile(opt::frac_bed_file.length(), 0, bwalker.header());
  }

  findex = fai_load(opt::refgenome.c_str());  // load the reference

  SeqLib::GRC blacklist;
  if (opt::mode == OPT_POWERSIM && !opt::blacklist.empty()) {
    std::cerr << "...reading blacklist file " << opt::blacklist << std::endl;
    blacklist = SeqLib::GRC(opt::blacklist, bwalker.Header());
    blacklist.CreateTreeMap();
    std::cerr << "...read in " << blacklist.size() << " blacklist regions " << std::endl;
  }

  // 
  if (opt::mode == OPT_ASSEMBLY)
    assemblyTest();
  else if (opt::mode == OPT_SIMBREAKS)
    genBreaks();
  else if (opt::mode == OPT_SPLITBAM)
    splitBam();
  else if (opt::mode == OPT_REALIGN)
    realignRandomSegments();
  else if (opt::mode == OPT_REALIGN_SV)
    realignBreaks();
  else if (opt::mode == OPT_POWERSIM)  {

    std::cerr << "...opening output" << std::endl;
    std::ofstream outfasta;
    SnowmanUtils::fopen("tumor_seq.fa", outfasta);
    
    std::ofstream events;
    SnowmanUtils::fopen("events.txt", events);
    
    PowerLawSim(findex, opt::nbreaks, -1.0001, blacklist, outfasta, events);
    outfasta.close();
    events.close();
  }
  else 
    std::cerr << "Mode not recognized. Chose from: --test-assembly, --sim-breaks, --split-bam" << std::endl;

}

std::string genBreaks() {

  // train on the input BAM
  std::vector<std::string> quality_scores;
  SeqLib::GenomicRegion v(0, 1000000,2000000);
    /*  std::vector<SeqLib::GenomicRegion> v = {
    SeqLib::GenomicRegion(0, 1000000, 2000000)
    SeqLib::GenomicRegion(0, 60000000,70000000),
    SeqLib::GenomicRegion(1, 1000000, 10000000), 
    SeqLib::GenomicRegion(1, 60000000,70000000),
    SeqLib::GenomicRegion(2, 1000000, 10000000), 
    SeqLib::GenomicRegion(3, 60000000,70000000),
    SeqLib::GenomicRegion(16, 1000000,1110000),
    SeqLib::GenomicRegion(17, 1000000,1110000),
    SeqLib::GenomicRegion(21, 1000000,1110000) 
    };*/
  
  /*
  SimTrainerWalker stw(opt::bam);
  stw.setBamWalkerRegions(v);

  stw.train();
  std::ofstream bamstats("bam_stats.txt");
  bamstats << stw.printBamStats() << std::endl;
  bamstats.close();
  */

  bwalker.SetRegion(v);
  SeqLib::BamRecord r; 
  std::cerr << "...sampling reads to learn quality scores" << std::endl;
  while (bwalker.GetNextRecord(r)) {
    std::string ss = r.Sequence();
    if (ss.find("AAAAAAAA") == std::string::npos && ss.find("TTTTTTTT") == std::string::npos) // already handle homopolymers
      quality_scores.push_back(r.Qualities());
  }
  
  std::cerr << "...loading the reference genome" << std::endl;


  SeqLib::GenomicRegion gg = regions[0];
  std::cerr << "--Generating breaks on: " << gg << std::endl; 
  std::cerr << "--Total number of rearrangement breaks: " << opt::nbreaks << std::endl; 
  std::cerr << "--Total (approx) number of indels: " << opt::nindels << std::endl; 

  SimGenome sg(gg, opt::nbreaks, opt::nindels, findex, opt::scramble, opt::viral_count);
  
  std::string final_seq = sg.getSequence();
  
  // write the final seq to a file
  std::ofstream fseq;
  fseq.open("tumor_seq.fa", std::ios::out);
  fseq << ">tumor\n" << final_seq << std::endl;
  fseq.close();
  //exit(1);

  ReadSim rs;

  std::ofstream ind;
  ind.open("indels.tsv", std::ios::out);
  for (auto& i : sg.m_indels)
    ind << i << std::endl;
  ind.close();


  std::ofstream con;
  con.open("connections.tsv", std::ios::out);
  con << sg.printBreaks();
  con.close();

  std::ofstream mic;
  mic.open("microbe_spikes.tsv", std::ios::out);
  mic << sg.printMicrobeSpikes();
  mic.close();
  
  exit(0);

  rs.addAllele(final_seq, 1);

  // sample paired reads
  std::vector<std::string> reads1;
  std::vector<std::string> reads2;
  std::vector<std::string> qual1;
  std::vector<std::string> qual2;
  std::cerr << "Simulating reads at coverage of " << coverages[0] << " del rate " << del_error_rates[0] << 
    " ins rate " << ins_error_rates[0] << " snv-rate " << snv_error_rates[0] << 
    " isize " << opt::isize_mean << "(" << opt::isize_sd << ")" << std::endl;
  rs.samplePairedEndReadsToCoverage(reads1, reads2, qual1, qual2, coverages[0], snv_error_rates[0], ins_error_rates[0], del_error_rates[0], 
				    opt::readlen, opt::isize_mean, opt::isize_sd, quality_scores);
  assert(reads1.size() == reads2.size());
  
  // write the paired end fastq. Give random errors
  std::cerr << "...writing/errorring reads 1" << std::endl;
  std::ofstream pe1;
  pe1.open("paired_end1.fastq", std::ios::out);
  size_t ccc= 0;
  for (size_t i = 0; i < reads1.size(); ++i) {
    rs.baseQualityRelevantErrors(reads1[i], qual1[i]);
    pe1 << "@r" << ccc++ << std::endl << reads1[i] << std::endl << "+\n" << qual1[i] << std::endl;
  }
  pe1.close();

  std::cerr << "...writing/erroring reads 2" << std::endl;  
  std::ofstream pe2;
  pe2.open("paired_end2.fastq", std::ios::out);
  ccc= 0;
  for (size_t i = 0; i < reads2.size(); ++i) { 
    //std::string qs = quality_scores[rand() % quality_scores.size()];
    rs.baseQualityRelevantErrors(reads2[i], qual2[i]);
    pe2 << "@r" << ccc++ << std::endl << reads2[i] << std::endl << "+\n" << qual2[i] << std::endl;
  }
  pe2.close();

  std::cerr << "********************************" << std::endl;
  std::cerr << "Suggest running: " << std::endl;
  std::cerr << "bwa mem -t 10 $REFHG19 paired_end1.fastq paired_end2.fastq | samtools sort -O bam -T /tmp -l 9 -m 16G > sim.bam && samtools index sim.bam" << std::endl;
  std::cerr << "snowman benchmark --split-bam -f 0.3,0.6 -b $n11 -k <region> -A norm" << std::endl;
  std::cerr << "samtools merge sim_wnormal.bam sim.bam norm_0.30_subsampled.bam" << std::endl;
  std::cerr << "samtools index norm_0.60_subsampled.bam && samtools index sim_wnormal.bam" << std::endl;
  //std::cerr << "bwa mem $REFHG19 paired_end1.fastq paired_end2.fastq | samtools sort -O bam -T /tmp -l 9 -m 16G > sim.bam && samtools sort sim.bam" << std::endl;
  //std::cerr << "bwa mem $REFHG19 paired_end1.fastq paired_end2.fastq > sim.sam && samtools view sim.sam -Sb > tmp.bam && " << 
  //  "samtools sort -m 4G tmp.bam sim && rm sim.sam tmp.bam && samtools index sim.bam" << std::endl;
  std::cerr << "********************************" << std::endl;

  return "";
}

void splitBam() {

  BamSplitter bs(opt::seed);
  assert(bs.Open(opt::bam));

  // set the regions to split on
  if (regions.size()) 
    bs.SetMultipleRegions(regions);

  std::cerr << "...set " << regions.size() << " walker regions " << std::endl;

  if (fractions_bed.size()) {
    bs.fractionateBam(opt::string_id + ".fractioned.bam", fractions_bed);
  } else {

    // set the output bams
    std::vector<std::string> fnames;
    for (auto& i : fractions) {
      std::stringstream ss;
      ss.precision(2);
      ss << opt::string_id << "_" << i << "_subsampled.bam";
      fnames.push_back(ss.str());
    }
    
    bs.setWriters(fnames, fractions);
    bs.splitBam();
  }

}

void assemblyTest() {

  SeqLib::GenomicRegion gr(16, 7565720, 7575000); //"chr17:7,569,720-7,592,868");
 
  std::cerr << "...loading the reference genome" << std::endl;
  findex = fai_load(opt::refgenome.c_str());  // load the reference
  std::string local_ref = ""; // debug // getRefSequence("16", gr, findex);

  size_t seqlen = local_ref.length();
  if (seqlen * 2 <= opt::readlen) {
    std::cerr << "**** Read length must be > 2 * sequence length" << std::endl;
    exit(EXIT_FAILURE);
  }

  // make the BWA Wrapper
  std::cerr << "...constructing local_seq index" << std::endl;
  SeqLib::BWAWrapper local_bwa; 
  local_bwa.ConstructIndex({{"local_ref", local_ref, std::string()}});

  // align local_seq to itself
  SeqLib::BamRecordVector self_align;
  local_bwa.AlignSequence(local_ref, "local_ref", self_align, false, false, 0);

  // write out the index
  local_bwa.WriteIndex("local_ref.fa");
  std::ofstream fa;
  fa.open("local_ref.fa");
  fa << ">local_ref" << std::endl << local_ref << std::endl;
  fa.close();

  std::cout << "coverage\tnumreads\tnumcontigs\tnumfinal\tcontig_coverage\tkmer_corr\terror_rate" << std::endl;
  for (int rep = 0; rep < opt::num_runs; ++rep) {
    std::cerr << "...assembly test. Working on iteration " << rep << " of " << opt::num_runs << std::endl;
    for (int k = 1; k <= 1; ++k) {
      for (auto& c : coverages) {
	for (auto& E : snv_error_rates) {
	  for (auto& D : del_error_rates) {
	    for (auto& I : ins_error_rates) {
	  	  
	      // make the read vector
	      ReadSim rs;
	      rs.addAllele(local_ref, 1);
	      
	      // sample reads randomly
	      //std::cerr << "...making read vector" << std::endl;
	      std::vector<std::string> reads;  
	      rs.sampleReadsToCoverage(reads, c, E, I, D, opt::readlen);
	      
	      // sample paired reads
	      //std::cerr << "...making paired-end read vector" << std::endl;
	      std::vector<std::string> reads1;
	      std::vector<std::string> reads2;
	      std::vector<std::string> qual1;
	      std::vector<std::string> qual2;
	      std::vector<std::string> qual = {std::string('I',opt::readlen)};
	      rs.samplePairedEndReadsToCoverage(qual1, qual2, reads1, reads2, c, E, I, D, opt::readlen, 350, 50, qual);
	      assert(reads1.size() == reads2.size());
	      
	      // align these reads to the local_seq
	      //std::cerr << "...realigned to local_seq" << std::endl;
	      SeqLib::BamRecordVector reads_to_local;
	      int count = 0;
	      for (auto& i : reads) {
		if (i.find("N") == std::string::npos) {
		  SeqLib::BamRecordVector read_hits;
		  local_bwa.AlignSequence(i, "read_" + std::to_string(++count), read_hits, false,false, 0);
		  if (read_hits.size())
		    reads_to_local.push_back(read_hits[0]);
		}
	      }
	      
	      // kmer filter the reads
	      KmerFilter kmer;
	      if (k == 1)
		kmer.correctReads(reads_to_local, reads_to_local);
	      
	      //std::cerr << " Attempted align of " << reads.size() << " to local_seq. Got hits on " << reads_to_local.size() << std::endl;
	  
	  // make plot of reads to contig
	  //AlignedContig sa(self_align);  
	  //sa.alignReads(reads_to_local);
	  
	  // assemble them
	  //std::cerr << "...assembling" << std::endl;
	  double error_rate = 0;
	  if (k == 0)
	    error_rate = 0.05;
	  int min_overlap = 35;
	  SnowmanAssemblerEngine engine("test", error_rate, min_overlap, opt::readlen);
	  engine.fillReadTable(reads_to_local);
	  engine.performAssembly(2);
	  
	  // align them back
	  SeqLib::BamRecordVector contigs_to_local;
	  for (auto& i : engine.getContigs()) {
	    SeqLib::BamRecordVector ct_alignments;
	    local_bwa.AlignSequence(i.Seq, i.Name, ct_alignments, false, false, 0);
	    AlignedContig ac(ct_alignments, opt::prefixes);
	    //ac.alignReads(reads_to_local);
	    //std::cout << ac;
	    contigs_to_local.insert(contigs_to_local.begin(), ct_alignments.begin(), ct_alignments.end());
	  }
	  
	  // write the results
	  SeqLib::GRC grc(contigs_to_local);
	  grc.MergeOverlappingIntervals();
	  double width = 0;
	  for (auto& i : grc)
	    width = std::max(width, (double)i.Width());
	  width = width / local_ref.length();
	  std::cout << c << "\t" << reads_to_local.size() << "\t" << engine.getContigs().size() << 
	    "\t" << grc.size() << "\t" << width << "\t" << k << "\t" << E << std::endl;
	  
	  if (k == 1 && c == 20 && E == 0.01) {
	    // write out the contig to local ref bam
	    SeqLib::BamWriter bw2;
	    bw2.Open("contigs_to_ref.bam");
	    bw2.SetHeader(local_bwa.HeaderFromIndex());
	    bw2.WriteHeader();
	    for (auto& i : contigs_to_local)
	      bw2.WriteRecord(i);
	    
	    // write the paired end fasta
	    std::ofstream pe1;
	    pe1.open("paired_end1.fa", std::ios::out);
	    size_t ccc= 0;
	    for (auto& i : reads1)
	      pe1 << ">r" << ccc++ << std::endl << i << std::endl;
	    pe1.close();

	    std::ofstream pe2;
	    pe2.open("paired_end2.fa", std::ios::out);
	    ccc= 0;
	    for (auto& i : reads2)
	      pe2 << ">r" << ccc++ << std::endl << i << std::endl;
	    pe2.close();
	    
	    // write out the read to local ref aligned bam
	    SeqLib::BamWriter bw;
	    bw.Open("reads_to_ref_" + std::to_string(c) + ".bam");
	    bw.SetHeader(local_bwa.HeaderFromIndex());
	    bw.WriteHeader();
	    for (auto& i : reads_to_local)
	      bw.WriteRecord(i);
	    
	    SeqLib::BamWriter bwk;
	    bwk.Open("k.bam");
	    bwk.SetHeader(local_bwa.HeaderFromIndex());
	    bwk.WriteHeader();
	    for (auto& i : reads_to_local) {
	      std::string kc = i.GetZTag("KC");
	      if (kc.length())
		i.SetSequence(kc);
	      bwk.WriteRecord(i);
	    }
	    
	  }
	    } // end kmer loop
	  } // end error loop
	} // end coverage
      }
    }
  }
}

  
    void parseBenchmarkOptions(int argc, char** argv) {
    
    bool die = false;

  if (argc < 2) 
    die = true;

  std::string del_er;
  std::string snv_er;
  std::string ins_er;
  std::string covs;
  std::string frac;

  std::string t;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'h': die = true; break;
    case OPT_ASSEMBLY: opt::mode = OPT_ASSEMBLY; break;
    case OPT_REALIGN: opt::mode = OPT_REALIGN; break;
    case OPT_REALIGN_SV: opt::mode = OPT_REALIGN_SV; break;
    case OPT_POWERSIM: opt::mode = OPT_POWERSIM; break;
    case OPT_BLACKLIST: arg >> opt::blacklist; break;
    case 'G': arg >> opt::refgenome; break;
    case 'n': arg >> opt::num_runs; break;
    case 's': arg >> opt::seed; break;
    case 'c': arg >> covs; break;
    case 'k': arg >> opt::regionFile; break;
    case 'b': arg >> opt::bam; break;
    case 'E': arg >> snv_er; break;
    case 'R': arg >> opt::nbreaks; break;
    case 'X': arg >> opt::nindels; break;
    case 'D': arg >> del_er; break;
    case 'I': arg >> ins_er; break;
    case 'A': arg >> opt::string_id; break;
    case 'f': arg >> frac; break;
    case 'M': arg >> opt::viral_count; break;
    case OPT_SIMBREAKS: opt::mode = OPT_SIMBREAKS; break;
    case OPT_SPLITBAM: opt::mode = OPT_SPLITBAM; break;
    case OPT_ISIZE_MEAN: arg >> opt::isize_mean; break;
    case OPT_ISIZE_SD: arg >> opt::isize_sd; break;
    case OPT_SCRAMBLE: opt::scramble = true; break;
    default: die= true; 
    }
  }

  if (die) {
    std::cerr << "\n" << BENCHMARK_USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }

  // parse the error rates
  snv_error_rates = parseErrorRates(snv_er);
  del_error_rates = parseErrorRates(del_er);
  ins_error_rates = parseErrorRates(ins_er);
  coverages = parseErrorRates(covs);

  // parse the fractions string or read file
  if (!SeqLib::read_access_test(frac))
    fractions = parseErrorRates(frac);
  else
    opt::frac_bed_file = frac;

  // check that the bam is valid
  if (opt::mode == OPT_SIMBREAKS && !SeqLib::read_access_test(opt::bam)) {
    std::cerr << "ERROR: Input BAM required for --sim-breaks" << std::endl;
    exit(EXIT_FAILURE);
  }

  // set the default error rates
  if (!snv_error_rates.size())
    snv_error_rates.push_back(DEFAULT_SNV_RATE);
  if (!del_error_rates.size())
    del_error_rates.push_back(DEFAULT_DEL_RATE);
  if (!ins_error_rates.size())
    ins_error_rates.push_back(DEFAULT_INS_RATE);
  if (!coverages.size())
    coverages.push_back(DEFAULT_COV);
  if ((!fractions.size() && !opt::frac_bed_file.length()) && opt::mode == OPT_SPLITBAM) {
    std::cerr << "Error: Must specify fractions to split into with -f (e.g. -f 0.1,0.8), or as BED file" << std::endl;
    exit(EXIT_FAILURE);
  }
  
    }

// parse the error rates string from the options
  std::vector<double> parseErrorRates(const std::string& s) {
    
    std::vector<double> out;
  
  std::istringstream is(s);
  std::string val;
  while(std::getline(is, val, ',')) {
    try {
      out.push_back(std::stod(val));
    } catch (...) {
      std::cerr << "Could not convert " << val << " to number. If you're inputting a file not CSV, then check file exists" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
    
  return out;
}
 
 std::string errorRateString(const std::vector<double>& v, const std::string& name) {
  
  std::stringstream ss;
  ss << "        " << name << ":";
  for (auto& i : v)
    ss << " " << i << ",";
  
  std::string out = ss.str();
  out.pop_back();
  return out;
  
}

void realignBreaks() {
  
  ReadSim rs;
  
  std::cerr << "...loading reference genome" << std::endl;
  SeqLib::BWAWrapper bwa; 
  bwa.LoadIndex(opt::refgenome);

  std::ofstream results;
  SnowmanUtils::fopen("realign.sv.results.csv", results);
  
  results << "id\tsize\tAlign1\tAlign2\tAlignBoth" << std::endl;
  const int RUN_NUM = 1000;
  //size_t count = 0;

  std::string chrstring1, chrstring2;
  for (int k = 0; k < RUN_NUM; ++k) {
  for (int i = 30; i < 500; i += 10) {

    SeqLib::GenomicRegion gr1, gr2;
    gr1.Random();
    gr2.Random();
    gr1.pos2 = gr1.pos1 + i/2; ; //+ k - 1 + (ins == 0 ? iii : 0); // add sequence to deletion ones, because it gets removed later
    gr2.pos2 = gr2.pos2 + i/2;
    chrstring1 = bwalker.Header().IDtoName(gr1.chr); //std::to_string(gr.chr+1);
    chrstring2 = bwalker.Header().IDtoName(gr2.chr); //std::to_string(gr.chr+1);
	
    int len;
    char * seq1 = faidx_fetch_seq(findex, const_cast<char*>(chrstring1.c_str()), gr1.pos1, gr1.pos2 - 1, &len);
    char * seq2 = faidx_fetch_seq(findex, const_cast<char*>(chrstring2.c_str()), gr2.pos1, gr2.pos2 - 1, &len);
    
    std::string s1, s2;
    if (seq1)
      s1 = std::string(seq1);
    else 
      continue;
    if (s1.find("N") != std::string::npos)
      continue;

    if (seq2)
      s2 = std::string(seq2);
    else 
      continue;
    if (s2.find("N") != std::string::npos)
      continue;

    
    if (rand() % 2)
      SeqLib::rcomplement(s2);
    std::string ss = s1 + s2;

    SeqLib::BamRecordVector aligns;
    bwa.AlignSequence(ss, std::to_string(i), aligns, false, 0.90, 2);
    
    bool a1 = false; bool a2 = false;
    for (auto& jj : aligns) {
      if (gr1.GetOverlap(jj.asGenomicRegion())) 
	a1 = true;
      if (gr2.GetOverlap(jj.asGenomicRegion()))       
	a2 = true;
    }

    results << k << "\t" << i << "\t" << a1 << "\t" << a2 << "\t" << (a1 && a2) << std::endl;
	    
  }
  }
    
}

void realignRandomSegments() {
  
  ReadSim rs;

  std::cerr << "...loading reference genome" << std::endl;
  SeqLib::BWAWrapper bwa; 
  bwa.LoadIndex(opt::refgenome);

  std::ofstream results;
  SnowmanUtils::fopen("realign.results.csv", results);

  results << "ID" << "\t" << "chr" << "\t" << "pos1" << "\t" << "pos2" << "\t" << "width" << "\t" << "num_aligns" << "\t" << "correct_hit_num" << "\t" << "snv_rate" << "\t" << "del_size" << "\t" << "ins_size" << std::endl;
  
  const int RUN_NUM = 1000;

  size_t tcount = 0;
    

  std::string chrstring;

  std::vector<int> widths = {70, 101, 250};
  std::vector<double> snv_rate = {0, 0.01, 0.05};
  std::vector<int> indel_size = {0, 1, 5, 20, 50};

  for (size_t ins = 0; ins < 2; ++ins) { //insertion or del
    for (auto& iii : indel_size) {
      for (auto& snv : snv_rate) {
	for (auto& k : widths) {
	  ++tcount;
	  for (size_t i = 0; i < RUN_NUM; ++i) {
	    if (i == 0)
	      std::cerr << "...working on width " << k << " and SNV rate " << snv << "\t" << (ins == 1 ? "INS" : "DEL") << "\t" << iii << ". " << tcount << " of " << (2 * widths.size() * snv_rate.size() * indel_size.size()) << std::endl;
	    
	    SeqLib::GenomicRegion gr;
	    gr.Random();
	    gr.pos2 = gr.pos1 + k - 1 + (ins == 0 ? iii : 0); // add sequence to deletion ones, because it gets removed later
	    chrstring = bwalker.Header().IDtoName(gr.chr); //bwalker.header()->target_name[gr.chr]; //std::to_string(gr.chr+1);
	    
	    int len;
	    char * seq = faidx_fetch_seq(findex, const_cast<char*>(chrstring.c_str()), gr.pos1, gr.pos2 - 1, &len);
	    
	    std::string s;
	    if (seq)
	      s = std::string(seq);
	    else 
	      continue;
	    
	    if (s.find("N") != std::string::npos)
	      continue;
	    
	    rs.makeSNVErrors(s, snv);
	    
	    if (gr.Width() < iii + 10)
	      continue;

	    if (ins == 1 && iii > 0)
	      rs.makeInsErrors(s, true, iii);
	    else if (iii > 0)
	      rs.makeDelErrors(s, iii);      
	    
	    SeqLib::BamRecordVector aligns;
	    bwa.AlignSequence(s, std::to_string(i), aligns, false, 0.90, 50);

	    int align_num = -1;
	    for (size_t j = 0; j < aligns.size(); ++j) {
	      
	      SeqLib::GenomicRegion gr_this = aligns[j].asGenomicRegion();
	      gr_this.Pad(100);
	      if (gr.GetOverlap(gr_this)) {
		
		if (iii == 0) { 
		  align_num = j;
		  break;
		} else {
		  if (ins == 1) {
		    int g = aligns[j].MaxInsertionBases();
		    if (g && g - 2 < iii && g + 2 > iii) {
		      align_num = j;
		      break;
		    }
		  } else {
		    int d = aligns[j].MaxDeletionBases();
		    if (d && d - 2 < iii && d + 2 > iii) {
		      align_num = j;
		      break;
		    }
		  }
		}
	      }
	      
	    }
	    
	    results << i << "\t" << chrstring << "\t" << gr.pos1 << "\t" << gr.pos2 << "\t" << k << "\t" << aligns.size() << "\t" << align_num << "\t" << snv << "\t" << (ins == 0 ? iii : 0) << "\t" << (ins == 1 ? iii : 0) << std::endl;
	    
	  }
	}
      }
    }
  }
  results.close();
  
}
