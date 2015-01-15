#include "SVBamReader.h"
#include "grabReads.h"
#include <regex>
#include <seqan/align.h>
#include <seqan/graph_msa.h>
#include "SVBamReader.h"
#include "unistd.h"
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <time.h>
#include <stdlib.h>

using namespace std;
using namespace seqan;
using std::ifstream;

typedef String<Dna5> TString;
typedef StringSet<TString> TStringSet;
typedef Index<StringSet<TString>, FMIndex<> > TIndex;
typedef Iterator<TIndex, TopDown<ParentLinks<> > >::Type TIter;
typedef String<char> TSequence;
typedef Align<TSequence,ArrayGaps> TAlign; 

typedef std::unordered_map<std::string, BamTools::RefData> RefMap;
typedef std::vector<BamTools::BamAlignment> BAVec;
typedef std::unordered_map<std::string, BamTools::BamAlignment> BAMap;

static struct timespec start;

struct IterPair {
  IterPair(BamAlignmentVector::const_iterator s, BamAlignmentVector::const_iterator e) : start(s), end(e) {}
  IterPair() {}
  ~IterPair() {}
  BamAlignmentVector::const_iterator start;
  BamAlignmentVector::const_iterator end;  
};

inline bool existTest (const std::string& name) {
  ifstream f(name.c_str());
  if (f.good()) {
    f.close();
    return true;
  } else {
    f.close();
    return false;
  }   
}

typedef std::map<std::string, IterPair> IterPairMap;
typedef std::vector<std::string> StringVec;

static const int CHR_LEN [25] = {249250621, 243199373, 198022430,191154276,180915260,171115067,
    159138663, 146364022, 141213431, 135534747, 135006516, 133851895,                                                                                                                                        
115169878, 107349540, 102531392, 90354753,  81195210,  78077248,                                                                                                                                         
59128983,  63025520,  48129895,  51304566,  155270560, 59373566,                                                                                                                                         
16571}; 

static const uint32_t CHR_CLEN [25] = {0, 249250621, 492449994,  690472424, 881626700, 1062541960, 1233657027,
				       1392795690,1539159712,1680373143,1815907890,1950914406,2084766301,
				       2199936179, 2307285719, 2409817111, 2500171864, 2581367074, 2659444322,                                                                                                                   
				       2718573305, 2781598825, 2829728720, 2881033286, 3036303846, 3095677412};

static const StringVec CHR_NAME {"1", "2", "3", "4", "5", "6", "7", "8", "9",
				  "10", "11", "12", "13", "14", "15", "16", "17", 
                                  "18", "19", "20", "21", "22", "X", "Y", "M"};

const int MAX_CHARS_PER_LINE = 512;
const int MAX_TOKENS_PER_LINE = 20;
const char* const DELIMITER = ",";

namespace opt
{
  static unsigned verbose = 1;
  static unsigned numThreads = 1;
  static unsigned minOverlap = 40;
  static int cutoff = 100;
  static string tbam;
  static string nbam;
  static unsigned isize = 2000;
  static unsigned sleepDelay = 0;
  static std::string regionFile = "x";
  static std::string outdir = "./";
  static unsigned mapq = 10;
  static unsigned chunk = 50e6;
  static int qualthresh = 4;
  static unsigned skip = 2000;
  static unsigned numbubble = 3;
  static float divergence = 0.05;
  static float gap_divergence = 0.05;
  static float error_rate = 0.05;
  static bool writeASQG = false;
}

static const char* shortopts = "t:n:p:m:l:v:i:s:r:ho:w:c:q:y:g:b:r:ae:";
static const struct option longopts[] = {
  { "tumor-bam",     no_argument, NULL, 't' },
  { "help",     no_argument, NULL, 'h' },
  { "normal-bam",     no_argument, NULL, 'n' },
  { "output-directory",     no_argument, NULL, 'o' },
  { "threads (\"processes\")",     no_argument, NULL, 'p' },
  { "min-overlap",     no_argument, NULL, 'm' },
  { "length-cutoff",     no_argument, NULL, 'l' },
  { "insert-size",     no_argument, NULL, 'i' },
  { "sleep-delay",     no_argument, NULL, 'i' },
  { "mapq",     no_argument, NULL, 'w'},
  { "skip-num",     no_argument, NULL, 'y'},
  { "qual-thresh",     no_argument, NULL, 'q'},
  { "region-file",     no_argument, NULL, 'r' },
  { "chunk-size",     no_argument, NULL, 'c' },
  { "error-overlap",     no_argument, NULL, 'e' },
  { "verbose",     no_argument, NULL, 'v' },
  { NULL, 0, NULL, 0 }
};

static const char *RUN_USAGE_MESSAGE =
"Usage: taiga run [OPTION] -b Tumor BAM\n\n"
"  Description: Grab SV reads from the SV bam and perform assembly with SGA\n"
"\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h  --help                           Display this help and exit\n"
"  -t, --tumor-bam                      Tumor SV bam\n"
"  -n, --normal-bam                     Normal SV bam\n"
"  -p, --threads                        use NUM threads to run Taiga DDefault: 1\n"
"  -m, --min-overlap                    Minimum read overlap, an SGA parameter. Default: 0.4* readlength\n"
"  -s, --sleep-delay                    Delay (in seconds) between starting processes to avoid BAM IO stampede. Default: 0\n"
"  -c, --chunk-size                     Size of region to read-in at once. Lower numbers require less memory, but more I/O. Default: 50e6\n"
"  -l  --length-cutoff                  Minimum length of the contig (discared otherwwise). Default: 1.5*readlength\n"
"  -r  --region-file                    Set a region txt file. Format: one region per line, Ex: 1,10000000,11000000\n"
"  -o  --output-directory               Where to place all of the intermediate and final files. Defaul: ./\n"
"  Read filters\n"
"  -i, --insert-size                    Minimum insert size to consider reads discordant. Default: 1000\n"
"  -w  --mapq                           Minimum MAPQ a read must have to move on to assembly. Default: 30\n"
"  -q  --qual-thresh                    Clip all bases with optical quality less than qual-thresh. Default: 4\n"
"  -y  --skip-num                       Do not assemble regions with > skip-num reads within 5000 bp window. Default: 2000\n"
"  Assembly params\n"
"  -g                                   Gap divergence. Default: 0.05\n"
"  -x                                   Divergence. Default: 0.05\n"
"  -b                                   Num rounds bubble. Default: 3\n"
"  -a                                   Output an ASQG graph file for each 5000bp window. Default: false\n"
"  -e                                   What is the error tolerance on read overlaps. Default: 0.05\n"
"\n";

// declare methods
void chunkReadsForAssembly(const int refID, const int pos1, const int chunk, const int pad, 
			   ContigVector * cont_out, const BamAlignmentVector * tbav, const BamAlignmentVector * nbav);
bool parseRegionFile(std::vector<int> &chr, std::vector<int> &pos1, std::vector<int> &pos2);
int sendWholeGenomeJobs(wqueue<TaigaWorkItem*> &mqueue);
int sendRegionJobs(wqueue<TaigaWorkItem*> &mqueue);
void getChunkReads(const BamAlignmentVector * srv, const unsigned refID, const unsigned pos1, const unsigned chunk, const unsigned pad, IterPairMap &mmap);
int alignReadsToContigs(TSequence &cseq, TSequence &qseq, TAlign &align);
void rcomplement(std::string &a);
uint32_t convertPos(unsigned refid, unsigned pos);

bool runTaiga(int argc, char** argv) {
  
  parseTaigaOptions(argc, argv);

  // check if Snowman is already started

  //create marker file to let know started, terminate if already made
  std::stringstream odir;
  odir << opt::outdir << "/started.txt";
  if (existTest(odir.str())) {
    std::cerr << "Snowman already started. Terminating...\n";
    return true;
  }
  std::ofstream marker(odir.str(), ios::out);
  marker << "Snowman started";
  marker.close();

  if (opt::verbose > 0) {
    std::cout << "Num threads: " << opt::numThreads << std::endl;
    std::cout << "Min overlap: " << opt::minOverlap << std::endl;
    std::cout << "Length cutoff: " << opt::cutoff << std::endl;
    std::cout << "Output directory: " << opt::outdir << std::endl;
    std::cout << "Chunk size: " << opt::chunk << std::endl;
    std::cout << "Read filters:\n"; 
    std::cout << "   Min MAPQ: " << opt::mapq << std::endl;
    std::cout << "   Insert size: " << opt::isize << std::endl;
  }

  // TODO: Adjust number of threads if region file is small

  // Create the queue and consumer (worker) threads
  wqueue<TaigaWorkItem*>  queue;
  std::vector<ConsumerThread<TaigaWorkItem>*> threadqueue;
  //std::vector<ConsumerThread<TaigaWorkItem, Contig>> threadqueue;
  for (unsigned i = 0; i < opt::numThreads; i++) {
    //ContigVector * tmp = new ContigVector();
    //ConsumerThread<TaigaWorkItem, Contig>* threadr = new ConsumerThread<TaigaWorkItem, Contig>(queue, opt::verbose > 0, tmp);
    ConsumerThread<TaigaWorkItem>* threadr = new ConsumerThread<TaigaWorkItem>(queue, opt::verbose > 0);
    threadr->start();
    threadqueue.push_back(threadr);
  }

  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);

  int num_jobs; 
  if (opt::regionFile.size() == 1) 
    num_jobs = sendWholeGenomeJobs(queue);
  else 
    num_jobs = sendRegionJobs(queue);
  
  // wait for the threads to finish
  for (unsigned i = 0; i < opt::numThreads; i++) 
    threadqueue[i]->join();

  if (opt::verbose > 0)  
    std::cout << "All threads done" << std::endl;

  // concatenate the files
  /*std::stringstream allread;
  allread << opt::outdir << "/all_reads.tab";
  if (opt::verbose > 0)
    std::cout << "...concatenating .fa and .tab files\n";
  allread.str(std::string());
  allread << "cat " << opt::outdir << "/*tmp.tab > " << opt::outdir << "/all_reads.tab";
  if (opt::verbose > 1)
    std::cerr << allread.str() << std::endl;;
  std::system(allread.str().c_str());
  */

  std::stringstream ss;
  ss << opt::outdir << "/all_contigs.fa";
  ss.str(std::string());
  ss << "cat " << opt::outdir << "/*tmp.fa > " << opt::outdir << "/all_contigs.fa";
  if (opt::verbose > 1)
    std::cout << ss.str() << std::endl;
  std::system(ss.str().c_str());

  // merge the BAM files
  sleep(2);
  std::stringstream bam_merge;
  if (num_jobs == 1) 
    bam_merge << "samtools sort " << opt::outdir << "/*reads2contig.bam r2c; samtools index " << opt::outdir << "/r2c.bam; rm " << opt::outdir << "/*reads2contig.bam";
  else
    bam_merge << "samtools merge " << opt::outdir << "/r2c_tmp.bam " << opt::outdir << "/*reads2contig.bam;" << 
      " samtools sort " << opt::outdir << "r2c_tmp.bam r2c; rm " << opt::outdir << "/r2c_tmp.bam; samtools index " << opt::outdir << "/r2c.bam; " <<
      "rm " << opt::outdir << "/*reads2contig.bam;";

  std::system(bam_merge.str().c_str());
  if (opt::verbose > 0)
    std::cerr << bam_merge.str() << std::endl;

  // remove tmp files
  std::stringstream rms;
  rms << " rm -f " << opt::outdir << "/*tmp.fa";
  std::system(rms.str().c_str());

  // set the ending marker
  std::stringstream odir2;
  odir2 << opt::outdir << "/ended.txt";
  std::ofstream marker2(odir2.str(), ios::out);
  marker2 << "Snowman ended";
  marker2.close();

  return 0;

}

//TODO turn to const
int runAll(BamAlignmentVector &bavd, std::string name, ContigVector * cont_out)
{

  // de-duplicate the reads 
  BamAlignmentVector bav;
  std::unordered_map<std::string, int> name_map;
  for (BamAlignmentVector::const_iterator it = bavd.begin(); it != bavd.end(); it++) {
    std::string uname = it->Name + "_" + std::to_string(it->AlignmentFlag); //.substr(it->Name.length() - 20, 20); // just get the last 20 chars
    //uname = uname + "FLAG_" + std::to_string(it->AlignmentFlag); // 10 for base 10
    std::string tmp; 
    it->GetTag("JW", tmp);
    unsigned tmp2;
    it->GetTag("HP", tmp2);
    // debug
    //if (name.compare("contig_1:10007000-10012000_") == 0)
    //  std::cerr << "UNAME: " << uname << " Pos: " << it->Position << " Tag " << tmp << " HP: " << tmp2 << std::endl;

    if (name_map.find(uname) == name_map.end()) {  // its not already added
      name_map.insert(std::pair<std::string, int>(uname, 0));
      bav.push_back(*it);
    } else { // it's already here
      //std::cerr << "De-duping " << uname << std::endl;
    }
  }

  // grab the pairmates
  BamAlignmentVector bav_disc, bav_disc_keep;
  for (BamAlignmentVector::iterator it = bav.begin(); it != bav.end(); it++) {
    if (std::abs(it->InsertSize) > opt::isize || (it->MateRefID != it->RefID)) {
      uint32_t val = convertPos(it->MateRefID+1, it->MatePosition);
      //uint32_t cscale = 1e8;
      //uint32_t lin_pos = ((uint32_t)it->MateRefID + 1)  * cscale + (uint32_t)(it->MatePosition/100);
      it->AddTag("RP", "I", val);
      bav_disc.push_back(*it);
    }
    //std::cout << it->Name << " " << it->MateRefID+1 << ":" << it->MatePosition << std::endl;
  }

  // create the keep vector
  if (bav_disc.size() > 1) { // need at least two disc reads to proceed
    std::sort( bav_disc.begin(), bav_disc.end(), BamTools::Algorithms::Sort::ByTag<uint32_t>("RP", BamTools::Algorithms::Sort::AscendingOrder));

    // get the differences in mate-pair position
    // c(1,2,3,10,21,22) -> 1,1,1,7,11,1 
    std::vector<uint32_t> diff_vec;
    uint32_t last_pos = 0;  
    for (BamAlignmentVector::const_iterator it = bav_disc.begin(); it != bav_disc.end(); it++) {

      uint32_t tmpr;
      it->GetTag("RP", tmpr);

      uint32_t diff_back = std::abs((int)(tmpr - (int)last_pos));
      diff_vec.push_back(diff_back);
      last_pos = tmpr;
    }

    // keep only some of the alignments
    for (int i = 0; i < bav_disc.size() - 1; i++) {
      if (diff_vec[i] <= 10 || diff_vec[i+1] <= 1000)
	bav_disc_keep.push_back(bav_disc[i]);
    }

    // sort by mate position
    if (bav_disc_keep.size() > 0) {
      std::sort( bav_disc_keep.begin(), bav_disc_keep.end(), BamTools::Algorithms::Sort::ByTag<uint32_t>("RP", BamTools::Algorithms::Sort::AscendingOrder));

      std::cerr << "Original number discordant: " << bav_disc.size() << " Final number: " << bav_disc_keep.size() << std::endl;
      for (BamAlignmentVector::const_iterator it = bav_disc_keep.begin(); it != bav_disc_keep.end(); it++)
        std::cerr << it->Name << " " << it->RefID+1 << ":" << it->Position << " Mate: " << it->MateRefID+1 << ":" << it->MatePosition << std::endl;

      // get the regions
      std::vector<int> disc_refID, disc_pos1, disc_pos2, disc_count;
      disc_refID.push_back(bav_disc_keep[0].MateRefID);
      disc_pos1.push_back(bav_disc_keep[0].MatePosition);
      int dref  = disc_refID.back();
      int dpos1 = disc_pos1.back();
      int dpos2 = dpos1;
      int dcount = 1;

      for (BamAlignmentVector::const_iterator it = bav_disc_keep.begin() + 1; it != bav_disc_keep.end(); it++) {
        std::cerr << "MateRefID: " << it->MateRefID << " dref: " << dref << " MatePos: " << (int)it->MatePosition << " dpso2: " << dpos2 << " DIST: " << std::abs((int)it->MatePosition - dpos2) << std::endl;
	if (it->MateRefID != dref || std::abs((int)it->MatePosition - dpos2) > 1000) { // set new region'
	  dref = it->MateRefID;
	  disc_pos2.push_back(dpos2); // finish the old one
	  disc_refID.push_back(it->MateRefID); // start a new one
	  disc_pos1.push_back(it->MatePosition);
	  disc_count.push_back(dcount);
	  dcount = 0;
	}
	dcount++;
        dpos2 = it->MatePosition;
      }
      disc_pos2.push_back(bav_disc_keep.back().MatePosition); // finish the last one
      disc_count.push_back(dcount);

      // print the regions for debugging
      std::cerr << "SIZES: " << disc_pos2.size() << " " << disc_pos1.size() << " " << disc_refID.size() << " " << disc_count.size() << std::endl;
      for (int i = 0; i < disc_pos2.size(); i++)
	std::cerr << "Region: " << (disc_refID[i]+1) << ":" << disc_pos1[i] << "-" << disc_pos2[i] << " Count: " << disc_count[i] << std::endl;

      // open the BAM and get the reads
      SVBamReader t_disc_reader(opt::tbam, "tumor", opt::isize, opt::mapq, opt::qualthresh, opt::minOverlap, opt::verbose);    
      SVBamReader n_disc_reader(opt::nbam, "normal", opt::isize, opt::mapq, opt::qualthresh, opt::minOverlap, opt::verbose);    
      for (int i = 0; i < disc_pos2.size(); i++) {
	if (disc_count[i] > 1) {
	  BamAlignmentVector tbav_disc, nbav_disc;

	  int buff = 100;
	  int p1 = disc_pos1[i] - buff;
	  int p2 = disc_pos2[i] + buff;
	  if (!t_disc_reader.setBamRegion(disc_refID[i], disc_refID[i], p1, p2))
	    std::cerr << "Failed to set BAM position in Tumor for Disc Read on position: " << disc_refID[i] << ":" << p1 << "-" << p2 << std::endl;
	  if (!n_disc_reader.setBamRegion(disc_refID[i], disc_refID[i], disc_pos1[i], disc_pos2[i]))
	    std::cerr << "Failed to set BAM position in Normal for Disc Read" << disc_refID[i] << ":" << p1 << "-" << p2 << std::endl;
	  
	  // read the tumor
          std::cerr << "SUCCESSFULLY OPENED INDICIES FOR " << disc_refID[i] << ":" << p1 << "-" << p2 << std::endl;
	  if (!t_disc_reader.bamToBAVec(tbav_disc)) {
	    std::cerr << "Failed to get BAM DISC reads for Tumor" << std::endl;
	    std::exit(EXIT_FAILURE);
	    return false;
	  }
	  std::cerr << "got " << tbav_disc.size() << " tum disc-region reads" << std::endl;
	  
	  // read the normal
	  if (!n_disc_reader.bamToBAVec(nbav_disc)) {
	    std::cerr << "Failed to get BAM DISC reads for Normal" << std::endl;
	    std::exit(EXIT_FAILURE);
	    return false;
	  }
	  std::cerr << "got " << nbav_disc.size() << " norm disc-region reads" << std::endl;
	}
      }

    }
  }


  //    if (diffr < 3) {// 3 is bp diff of 300
  //  bav_disc_keep.push_back(*it);
  //}
  //std::cerr << it->Name << " " << it->RefID+1 << ":" << it->Position << " Mate: " << it->MateRefID+1 << ":" << it->MatePosition << " RP: " << tmpr << " DIFF: " << diffr << std::endl;
    
  //DEBUG
  //if (name.compare("contig_1:10007000-10012000_") == 0)
  //  for (BAVec::const_iterator jj = bav.begin(); jj != bav.end(); jj++)
  //    std::cerr << "BAV: " << jj->Name << " " << jj->QueryBases << std::endl;

  // DEBUG
  //std::cerr << name << std::endl;
  //if (name.compare("contig_1:10011000-10016000_") == 0)
//   for (BamAlignmentVector::const_iterator it = bavd.begin(); it != bavd.end(); it++) {
//       std::string tmpr;
//       it->GetTag("JW", tmpr);
//       std::string seqr;
//       it->GetTag("TS", seqr);
//        std::cerr << it->AlignmentFlag << " " << it->Name << " JW: " << tmpr << " TS " << seqr <<std::endl;
//   }
//   std::cerr <<"BAVSIZE: " << bav.size() << std::endl;

  ReadTable* pRT = new ReadTable(bav);

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

  double errorRate = (double)opt::error_rate;
  int seedLength = opt::minOverlap;
  int seedStride = seedLength;
  bool bIrreducibleOnly = true; // default

  OverlapAlgorithm* pOverlapper = new OverlapAlgorithm(pBWT, pRBWT, 
                                                       errorRate, seedLength, 
                                                       seedStride, bIrreducibleOnly);
  pOverlapper->setExactModeOverlap(false);
  pOverlapper->setExactModeIrreducible(false);

  std::stringstream hits_stream;
  std::stringstream asqg_stream;

  ASQG::HeaderRecord headerRecord;
  headerRecord.setOverlapTag(opt::minOverlap);
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
    OverlapResult rr = pOverlapper->overlapReadInexact(read, opt::minOverlap, &obl);

    pOverlapper->writeOverlapBlocks(hits_stream, workid, rr.isSubstring, &obl);

    ASQG::VertexRecord record(read.id, read.seq.toString());
    record.setSubstringTag(rr.isSubstring);
    record.write(asqg_stream);

    workid++;
  }
  std::string line;
  bool bIsSelfCompare = true;
  ReadInfoTable* pQueryRIT = new ReadInfoTable(pRT);
  delete pRT;

  while(getline(hits_stream, line)) {
    size_t readIdx;
    size_t totalEntries;
    bool isSubstring; 
    OverlapVector ov;
    OverlapCommon::parseHitsString(line, pQueryRIT, pQueryRIT, pSAf, pSAr, bIsSelfCompare, readIdx, totalEntries, ov, isSubstring);
    for(OverlapVector::iterator iter = ov.begin(); iter != ov.end(); ++iter)
    {
       ASQG::EdgeRecord edgeRecord(*iter);
       edgeRecord.write(asqg_stream);
    }


  }
  //std::cout << "Hits Stream:\n" << hits_stream.str() << std::endl;
  //std::cout << "ASQG Stream:\n" << asqg_stream.str() << std::endl;

  // optionally output the graph structure
  if (opt::writeASQG) {
    std::stringstream asqgfile;
    asqgfile << opt::outdir << "/" << name << ".asqg";
    // write ASQG to file for visualization
    std::ofstream ofile(asqgfile.str(), ios::out);
    ofile << asqg_stream.str();
    ofile.close();
    }

    // Get the number of strings in the BWT, this is used to pre-allocated the read table
    delete pOverlapper;
    delete pBWT; 
    delete pRBWT;
    delete pSAf;
    delete pSAr;

    int maxEdges = 128;
    int trimLengthThreshold = opt::cutoff; //useless for numTrimRounds = 0
    int numTrimRounds = 0; //
    bool bPerformTR = false; // transitivie edge reducetion
    bool bValidate = false;
    int resolveSmallRepeatLen = -1; 
    int numBubbleRounds = opt::numbubble;
    double maxBubbleGapDivergence = opt::gap_divergence;
    double maxBubbleDivergence = opt::divergence;
    int maxIndelLength = 20;
    bool bExact = true;
    ContigVector contigs;

    assemble(asqg_stream, opt::minOverlap, maxEdges, bExact, 
	      trimLengthThreshold, bPerformTR, bValidate, numTrimRounds, 
              resolveSmallRepeatLen, numBubbleRounds, maxBubbleGapDivergence, 
	     maxBubbleDivergence, maxIndelLength, opt::cutoff, name, contigs);
    if (contigs.size() >= 1 && opt::verbose > 1) 
      std::cout << "Contig Count: " << contigs.size() << " at " << name << std::endl;

    // local alignment parameters
    //unsigned counter = 0;

    // keep only the top 5 contigs
    //std::sort(contigs.rbegin(), contigs.rend());
    //contigs.erase(std::min(contigs.begin() + 5, contigs.end()), contigs.end());    
    //if (contigs.size() > 1)
    //  if ( (contigs[0].getLength() - contigs[1].getLength()) > 300)
	//contigs.erase(contigs.begin() + 1, contigs.end());

    //std::cerr << "Size: " << contigs.size() << std::endl;
    //for (ContigVector::const_iterator it = contigs.begin(); it != contigs.end(); it++)
    //  std::cerr << it->getLength() << " ";
    //std::cerr << std::endl;

    // MATCHING BY FIND
    int buff = 10;
    for (ContigVector::iterator i = contigs.begin(); i != contigs.end(); i++) {
      for (BamAlignmentVector::iterator j = bav.begin(); j != bav.end(); j++) {
	std::string QB;
	j->GetTag("TS", QB);
        //int seqlen = j->QueryBases.length();
	int seqlen = QB.length();
        size_t posa = i->getSeq().find(QB.substr(0,buff)); 
	
        // PROCEED IF ALIGNS TO FORWARD
        if (posa != std::string::npos) {
	  size_t posb = i->getSeq().find(QB.substr(std::max(seqlen-buff,0),buff)); 
          if (posb != std::string::npos) {
            //i->addRead(j->Name, posa, QB);
	    i->addRead(*j, posa);
	  }
	} 	 
 
        //OTHERWISE TRY REVERSE
        else {
	  std::string rstring = QB;
	  rcomplement(rstring); 
          posa = i->getSeq().find(rstring.substr(0,buff)); 
          if (posa != std::string::npos) {
	    size_t posb = i->getSeq().find(rstring.substr(std::max(seqlen-buff,0),buff)); 
            if (posb != std::string::npos) {
	      //i->addRead(j->Name, posa, rstring);
	      j->EditTag("TS", "Z", rstring); // edit the tag to be reverseComplemented
	      i->addRead(*j, posa);
	    }
	  }
	}

	//debug
	//	if (j->Name.compare("H7AJUADXX130905:2:2111:6629:64887")==0) {
	//  std::cerr << "POS: " << posa << " Rseq: " << QB << " CSeq: " << i->getSeq() << std::endl;
	//	}

      } // end read loop

      // add the alignment if there are 3+ reads
      if (i->getReadCount() > 2)
	cont_out->push_back(*i); 

      //if (i->getID().compare("contig_1:10915000-10920000_4")==0)
      //  i->printBamAlignments();

      // print verbose message
      if (opt::verbose > 1) 
	std::cerr << "Contig Size: " << i->getSeq().length() << " Read Support: " << i->getReadCount() << std::endl;

    }
    
    delete pQueryRIT;
    asqg_stream.str("");
    hits_stream.str("");

    if (opt::verbose > 1) {
      std::cerr << "Num contigs: " << contigs.size() << std::endl;
    }
    return 0;
}

bool grabReads(int refID, int pos1, int pos2, ContigVector * cont_out) {
 
  //  SeqRecordVector tsrv, nsrv;
  BamAlignmentVector tbav, nbav;

  ////// TUMOR
  // open the reads
  SVBamReader t_reader(opt::tbam, "tumor", opt::isize, opt::mapq, opt::qualthresh, opt::minOverlap, opt::verbose);
  // find the bai
  if (!t_reader.findBamIndex())
    std::cerr << "Failed to open BAM index in Tumor" << std::endl;
  // set the BAM region
  if (!t_reader.setBamRegion(refID, refID, pos1, pos2))
    std::cerr << "Failed to set BAM position in Tumor" << std::endl;
  // add the reads to the SeqRecordVector
  if (!t_reader.bamToBAVec(tbav)) {
    std::cerr << "Failed to get DFDF BAM reads for Tumor" << std::endl;
    std::exit(EXIT_FAILURE);
    return false;
  }
  unsigned num_t_reads = tbav.size();

  ////// NORMAL
  // open the reads
  SVBamReader n_reader(opt::nbam, "normal", opt::isize, opt::mapq, opt::qualthresh, opt::minOverlap, opt::verbose);
  // find the bai
  n_reader.findBamIndex();
  // set the BAM region
  if (!n_reader.setBamRegion(refID, refID, pos1, pos2))
    std::cerr << "Failed to set BAM position in Normal" << std::endl;
  // add the reads to the SeqRecordVector
  n_reader.bamToBAVec(nbav);
  unsigned num_n_reads = nbav.size(); 

  // verbose
  if (opt::verbose > 0) {
    std::string print1 = AddCommas<int>(pos1);
    std::string print2 = AddCommas<int>(pos2);
    //char buffer[200];
    
    //sprintf(buffer, "Running region- %2d:%11s", )
    std::cout << "Running region- " << refID + 1 << ":" << print1 << "-" << print2 << 
      "\tTumor: " << num_t_reads << "\tNormal: " << num_n_reads << std::endl;
  }

  // Divide up the reads and send them to the assembler
  const int chunk = 4000;
  const int pad = 1000;
  chunkReadsForAssembly(refID, pos1, chunk, pad, cont_out, &tbav, &nbav);

  // write the cont out to a file
  std::ofstream ostream;
  std::stringstream ofile;
  ofile << opt::outdir << "/" << refID + 1 << "_" << pos1 << ".tmp.fa";
  ostream.open(ofile.str());

  // write the reads out
  //std::ofstream rstream; 
  //std::stringstream rfile;
  //rfile << opt::outdir << "/" << refID + 1 << "_" << pos1 << ".tmp.tab";
  //rstream.open(rfile.str());

  std::vector<Contig>::iterator it = cont_out->begin();
  for(std::vector<Contig>::iterator it = cont_out->begin(); it != cont_out->end(); ++it) {
    ostream << ">" << it->getID() << "\n";
    ostream << it->getSeq() << "\n";
    //rstream << it->printAlignments();
  }
  ostream.close();
  //rstream.close();

  // set the reference for the BAM
  BamTools::RefVector ref;  
  for (int i = 0; i < 25; i++) {
    BamTools::RefData rf(CHR_NAME[i], CHR_LEN[i]);
    ref.push_back(rf);      
  }

  // open the BAM for writing
  BamTools::SamHeader sam("none");
  std::stringstream r2cfile;
  r2cfile << opt::outdir << "/" << refID+1 << "_" << pos1 << "_reads2contig.bam";
  std::string outbam = r2cfile.str();
  //std::string rmr = "rm " + outbam;
  //std::system(rmr.c_str());
  if (opt::verbose > 0)
    std::cerr << "Writing the read2contig BAM: " << outbam << std::endl;
  BamTools::BamWriter writer;
  if (!writer.Open(outbam, sam, ref))  
    std::cerr << "Error initializing the BAM for: " << outbam << std::endl;

  // transfer the alignments to a vector
  BAVec read_vec;
  for(std::vector<Contig>::const_iterator it = cont_out->begin(); it != cont_out->end(); ++it) {
    BAVec tmpvec = it->getBamAlignments();
    BAVec::const_iterator jt = tmpvec.begin();
    for (; jt != tmpvec.end(); jt++)
      read_vec.push_back(*jt);  
  }

  BAVec read_vec_reduced;
  // sort by read-name to add tags
  std::sort(read_vec.begin(), read_vec.end(), BamTools::Algorithms::Sort::ByTag<std::string>("JW", BamTools::Algorithms::Sort::AscendingOrder) );  
  unsigned aa = 0;
  unsigned bb = 1;
  while (aa < read_vec.size()) {
    std::string curr_name; // get the curent read name
    read_vec[aa].GetTag("JW", curr_name);
    std::string newtag; // current contig tag
    read_vec[aa].GetTag("CN", newtag);
    while (aa+bb < read_vec.size()) { // keep looping until dont hit repeat name
      std::string tmpname; // get the current read name
      read_vec[aa+bb].GetTag("JW", tmpname);
      //std::cerr << tmpname << " " << curr_name << std::endl;
      if (tmpname.compare(curr_name)==0) { // have the same read name, add CN tags
	std::string this_cn; 
	read_vec[aa+bb].GetTag("CN", this_cn);
	newtag = newtag + "x" + this_cn;
	bb++;
	//std::cerr << tmpname << " " << newtag << std::endl;
      } else { // diff name. Set the tags {
	// add all the tags in
	bb--; // if no hits above, move bb back to 0.
	if (bb > 0) // only edit if it needs to be changed
	  read_vec[aa].EditTag("CN", "Z", newtag);

	// add to the final output
	read_vec_reduced.push_back(read_vec[aa]);
	//std::string ccc;
	//read_vec_reduced[read_vec_reduced.size() - 1].GetTag("CN", ccc);
	//std::cerr << "CCC: " << ccc << std::endl;
	// reset everything
	aa = aa + bb; // skip over all the stuff we already did
	break; // break the inner while
      }
    } // end inner while

    bb = 1;
    aa++;
  } 

  // sort it
  std::sort(read_vec_reduced.begin(), read_vec_reduced.end(), BamTools::Algorithms::Sort::ByPosition(BamTools::Algorithms::Sort::AscendingOrder) );

  // write the alignments to the BAM
  for (BAVec::const_iterator it = read_vec_reduced.begin(); it != read_vec_reduced.end(); it++) 
    writer.SaveAlignment(*it);

  writer.Close();
  
  // read the bam to index it
  //BamTools::BamReader reader_ind;
  //if(!reader_ind.Open(outbam))
  //  std::cerr << "Failed to open BAM which was just made: " << outbam << std::endl;
  //if(!reader_ind.CreateIndex())
  //  std::cerr << "Failed to create BAM index for BAM: " << outbam << std::endl;    

  // sort it and index it
  //std::string sort_call = "samtools sort " + outbam + " " + 
  // opt::outdir + "/r2c";
  //std::string rm_call = "rm " + outbam;
  //std::string index_call = "samtools index " + opt::outdir + "r2c.bam";

  //std::cerr << sort_call << std::endl;
  //std::system(sort_call.c_str());
  //std::cerr << rm_call << std::endl;
  //std::system(rm_call.c_str());
  //std::cerr << index_call << std::endl;
  //std::system(index_call.c_str());

  // clear the read structure
  tbav.clear();
  nbav.clear();

  // delete the contig 
  delete cont_out;

  // display the run time
  struct timespec finish;
  double elapsed;
  clock_gettime(CLOCK_MONOTONIC, &finish);
  elapsed = (finish.tv_sec - start.tv_sec);
  int t = std::clock()/CLOCKS_PER_SEC;
  int min = (int)std::floor(elapsed / 60.0);
  int sec = (int)(elapsed-min*60);
  char buffer[100];
  sprintf (buffer, "CPU time: %dm%ds    Wall time: %dm%ds", 
            (int)std::floor( ((double)t) /60.0), t % 60, min, sec);
  printf ("%s\n",buffer);

  return true;
}

void parseTaigaOptions(int argc, char** argv) {
  bool die = false;

  if (argc < 2) 
    die = true;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
      case 'p': arg >> opt::numThreads; break;
      case 'h': die = true; break;
      case 'a': opt::writeASQG = true; break;
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
      case 'c': arg >> opt::chunk; break;
      case 'q': arg >> opt::qualthresh; break;
      case 'y': arg >> opt::skip; break;
      case 'g': arg >> opt::gap_divergence; break;
      case 'b': arg >> opt::numbubble; break;
      case 'x': arg >> opt::divergence; break;
      case 'e': arg >> opt::error_rate; break;
    }
  }

  // clean the outdir
  

  if(opt::numThreads <= 0)
    {
      std::cout << "run: invalid number of threads: " << opt::numThreads << "\n";
      die = true;
    }

  if (opt::cutoff < 0)
      opt::cutoff = opt::minOverlap * 1.5;

  if (die) 
    {
      std::cout << "\n" << RUN_USAGE_MESSAGE;
      exit(1);
    }

}

// Divide up the reads and send them to the assembler
void chunkReadsForAssembly(const int refID, const int pos1, const int chunk, const int pad, 
			   ContigVector * cont_out, const BamAlignmentVector * tbav, const BamAlignmentVector *nbav) {

  //if (opt::verbose > 0)
  //  std::cerr << "...preparing read chunks\n";
  IterPairMap tmap, nmap;

  // return vectors of iterators that points to positions for chunks
  getChunkReads(tbav, refID, pos1, chunk, pad, tmap);
  getChunkReads(nbav, refID, pos1, chunk, pad, nmap);

  /*std::cerr << "TMAP SIZE: " << tmap.size() << std::endl;
  std::cerr << "NMAP SIZE: " << nmap.size() << std::endl;
  for(IterPairMap::const_iterator it = nmap.begin(); it != nmap.end(); it++) {
    std::cerr << it->first << "size: " << it->second.end - it->second.start << std::endl;
    for (BamAlignmentVector::const_iterator jt = it->second.start; jt != it->second.end; jt++)
      std::cerr << jt->Name << " " << jt->Position << std::endl; 
    std::cerr << "\n\n";
    }*/

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

     /*     std::cerr << *it << std::endl;
     std::cerr << tmap[*it].end - tmap[*it].start << std::endl;
     std::cerr << nmap[*it].end - nmap[*it].start << std::endl;
     std::cerr << fvec.size() << std::endl;*/

     if (fvec.size() < opt::skip && fvec.size() >= 5) 
       runAll(fvec, *it, cont_out);
     else if (fvec.size() >= 500 && opt::verbose > 1)
       std::cerr << "Too many reads at: " << *it << " with " << fvec.size() << " reads\n";
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

    unsigned myk_pos;
    myk->GetTag("HP", myk_pos);
    //std::cerr << myk_pos << " MYK " << bav->end() - myk << std::endl;

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
	std::stringstream sstm;
	sstm << "contig_" << refID + 1 << ":" << max(nextstart - chunk, pos1) << "-" << currend << "_";
	//std::cerr << sstm.str() << std::endl;
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
    std::stringstream sstm;
    sstm << "contig_" << refID + 1 << ":" << max(nextstart - chunk, pos1) << "-" << currend << "_";
    IterPair p(myk_curr_start, myk-1); //myk-1 is the same as bav.end()
    mmap.insert(pair<string, IterPair>(sstm.str(), p));
  } 

}


// Parse the region file specified by the -r flag. Must be comma delimited
bool parseRegionFile(std::vector<int> &chr, std::vector<int> &pos1, std::vector<int> &pos2) {

  ifstream fin;
  fin.open(opt::regionFile); // open a file
  if (!fin.good()) 
    return false; // exit if file not found

  // process (print) the tokens
  if (opt::verbose > 0)
     std::cout << "Regions:\n";

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
	{
	  for (n = 1; n < MAX_TOKENS_PER_LINE; n++)
	    {
	      token[n] = strtok(0, DELIMITER); // subsequent tokens
	      if (!token[n]) break; // no more tokens
	    }
	}

      for (int i = 0; i < n; i += 3) {
	std::string mchr  = token[i];
	std::string mpos1 = token[i+1];
	std::string mpos2 = token[i+2];
        chr.push_back( atoi(mchr.c_str()));
        pos1.push_back(atoi(mpos1.c_str()));
        pos2.push_back(atoi(mpos2.c_str()));
        if (opt::verbose > 0) 
          std::cout << "   " << token[i] << ":" << token[i+1] << "-" << token[i+2] << std::endl;

      }
  } // end while

  return true;
  
}

// fill the worker queue with Taiga jobs tiling the whole genome
int sendWholeGenomeJobs(wqueue<TaigaWorkItem*>  &mqueue) {

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
}


// fill the worker queue with jobs from the region file
int sendRegionJobs(wqueue<TaigaWorkItem*>  &mqueue) {

  // set the regions, if they exist
  std::vector<int> chr;
  std::vector<int> pos1;
  std::vector<int> pos2;

  // open the region file if it exists
  parseRegionFile(chr, pos1, pos2);

  for (int it = 0; it != pos1.size(); it++) {
    if (pos2[it] > CHR_LEN[chr[it]-1])
      pos2[it] = CHR_LEN[chr[it]-1];
  }

  // reduce the number of threads if not needed
  opt::numThreads = min((int)opt::numThreads, (int)chr.size());

  int threadchunk = opt::chunk;
  unsigned jj = 0; 
  int startr, endr;;
  int kk = 0;

  //for (int it= 0; it != chr.size(); it++)
  //  std::cerr << chr[it] << " " << pos1[it] << " " << pos2[it] << std::endl;

  int num_jobs = 0;

  // loop through each region
  bool stoploop = false;
  while (jj < chr.size()) {
    // if regions are greater than chunk, breakup
    if ( (pos2[jj] - pos1[jj]) > threadchunk) {
      startr = pos1[jj]; 
      endr = startr + threadchunk;

      do {
	kk++;
	num_jobs++;
        ContigVector * cont_out = new ContigVector();
        TaigaWorkItem * item = new TaigaWorkItem(opt::tbam, opt::nbam, chr[jj]-1, startr, endr, kk*1e6 + jj, cont_out);
        mqueue.add(item);

	if (endr == pos2[jj])
	  stoploop = true;

        endr   = min(pos2[jj], (kk+1)*threadchunk + pos1[jj]);
        startr = min(pos2[jj],  kk*threadchunk + pos1[jj]);
        sleep(opt::sleepDelay);

      } while (!stoploop);

    } else { // region size is below chunk
        num_jobs++;
        ContigVector * cont_out = new ContigVector();
        TaigaWorkItem * item = new TaigaWorkItem(opt::tbam, opt::nbam, chr[jj]-1, pos1[jj], pos2[jj], jj, cont_out);
        mqueue.add(item);
        sleep(opt::sleepDelay);
    }
    jj++;
    kk = 0;
    stoploop = false;
  } // end big while

  return num_jobs;
}

// 
int alignReadsToContigs(TSequence &cseq, TSequence &qseq, TAlign &align) {

    int match = 4;
    int mismatch = -2;
    int gapopen = -4;
    int gapextend = -2;
    int score;
 
    resize(rows(align), 2);
    assignSource(row(align,0),cseq);
    assignSource(row(align,1),qseq);
    score = localAlignment(align, Score<int,Simple>(match, mismatch, gapopen, gapextend));
    return score;

}

void rcomplement(std::string &a) {

  std::reverse(&a[0], &a[a.size()]);
  std::string::iterator it = a.begin();
  for (; it != a.end(); it++)
    if (*it == 'A')
      *it = 'T';
    else if (*it == 'T')
      *it = 'A';
    else if (*it == 'C')
      *it = 'G';
    else
      *it = 'C';

}

uint32_t convertPos(unsigned refid, unsigned pos) {

  if (refid > 24)
    refid = 24;
  uint32_t out = CHR_CLEN[refid] + pos;
  return out;
}
