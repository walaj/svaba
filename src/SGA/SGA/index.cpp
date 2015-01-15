//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// index - Build a BWT/FM-index for a set of reads
//
#include <iostream>
#include <fstream>
#include <algorithm>
#include "SGACommon.h"
#include "Util.h"
#include "index.h"
#include "SuffixArray.h"
#include "SeqReader.h"
#include "SACAInducedCopying.h"
#include "BWTDiskConstruction.h"
#include "BWT.h"
#include "Timer.h"
#include "BWTCABauerCoxRosone.h"
#include "BWTCARopebwt.h"
#include "SampledSuffixArray.h"

//
// Getopt
//
#define SUBPROGRAM "index"

static const char *INDEX_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
  "Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *INDEX_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Index the reads in READSFILE using a suffixarray/bwt\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"  -a, --algorithm=STR                  BWT construction algorithm. STR can be:\n"
"                                       sais - induced sort algorithm, slower but works for very long sequences (default)\n"
"                                       ropebwt - very fast and memory efficient. use this for short (<200bp) reads\n"
"  -d, --disk=NUM                       use disk-based BWT construction algorithm. The suffix array/BWT will be constructed\n"
"                                       for batchs of NUM reads at a time. To construct the suffix array of 200 megabases of sequence\n"
"                                       requires ~2GB of memory, set this parameter accordingly.\n"
"  -t, --threads=NUM                    use NUM threads to construct the index (default: 1)\n"
"  -c, --check                          validate that the suffix array/bwt is correct\n"
"  -p, --prefix=PREFIX                  write index to file using PREFIX instead of prefix of READSFILE\n"
"      --no-reverse                     suppress construction of the reverse BWT. Use this option when building the index\n"
"                                       for reads that will be error corrected using the k-mer corrector, which only needs the forward index\n"
"      --no-forward                     suppress construction of the forward BWT. Use this option when building the forward and reverse index separately\n"
"      --no-sai                         suppress construction of the SAI file. This option only applies to -a ropebwt\n"
"  -g, --gap-array=N                    use N bits of storage for each element of the gap array. Acceptable values are 4,8,16 or 32. Lower\n"
"                                       values can substantially reduce the amount of memory required at the cost of less predictable memory usage.\n"
"                                       When this value is set to 32, the memory requirement is essentially deterministic and requires ~5N bytes where\n"
"                                       N is the size of the FM-index of READS2.\n"
"                                       The default value is 8.\n"
  "\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
  static unsigned int verbose;
  static std::string readsFile;
  static std::string prefix;
  static std::string algorithm = "sais";
  static int numReadsPerBatch = 2000000;
  static int numThreads = 1;
  static bool bDiskAlgo = false;
  static bool bBuildReverse = true;
  static bool bBuildForward = true;
  static bool bBuildSAI = true;
  static bool validate;
  static int gapArrayStorage = 4;
}

static const char* shortopts = "p:a:m:t:d:g:cv";

enum { OPT_HELP = 1, OPT_VERSION, OPT_NO_REVERSE, OPT_NO_FWD, OPT_NO_SAI };

static const struct option longopts[] = {
  { "verbose",     no_argument,       NULL, 'v' },
  { "check",       no_argument,       NULL, 'c' },
  { "prefix",      required_argument, NULL, 'p' },
  { "threads",     required_argument, NULL, 't' },
  { "disk",        required_argument, NULL, 'd' },
  { "gap-array",   required_argument, NULL, 'g' },
  { "algorithm",   required_argument, NULL, 'a' },
  { "no-reverse",  no_argument,       NULL, OPT_NO_REVERSE },
  { "no-forward",  no_argument,       NULL, OPT_NO_FWD },
  { "no-sai",      no_argument,       NULL, OPT_NO_SAI },
  { "help",        no_argument,       NULL, OPT_HELP },
  { "version",     no_argument,       NULL, OPT_VERSION },
  { NULL, 0, NULL, 0 }
};

int indexMain(int argc, char** argv)
{
  Timer t("sga index");
  parseIndexOptions(argc, argv);
  if(!opt::bDiskAlgo)
    {
      if(opt::algorithm == "sais")
	indexInMemorySAIS();
      else if(opt::algorithm == "bcr")
	indexInMemoryBCR();
      else if(opt::algorithm == "ropebwt")
	indexInMemoryRopebwt();
    }
  else
    {
      indexOnDisk();
    }
  return 0;
}

//
void indexInMemoryBCR()
{
  std::cout << "Building index for " << opt::readsFile << " in memory using BCR\n";

  if(opt::bBuildForward || opt::bBuildReverse)
    {
      // Parse the initial read table
      std::vector<DNAEncodedString> readSequences;
      SeqReader reader(opt::readsFile);
      SeqRecord sr;
      while(reader.get(sr))
	readSequences.push_back(sr.seq.toString());

      if(opt::bBuildForward)
	{
	  BWTCA::runBauerCoxRosone(&readSequences, opt::prefix + BWT_EXT, opt::prefix + SAI_EXT);
	}


      if(opt::bBuildReverse)
	{
	  // Reverse all the reads
	  for(size_t i = 0; i < readSequences.size(); ++i)
	    readSequences[i] = reverse(readSequences[i].toString());
	  BWTCA::runBauerCoxRosone(&readSequences, opt::prefix + RBWT_EXT, opt::prefix + RSAI_EXT);
	}
    }
}

//
void indexInMemoryRopebwt()
{
  std::cout << "Building index for " << opt::readsFile << " in memory using ropebwt\n";

  bool use_threads = opt::numThreads >= 4;

  if(opt::bBuildForward)
    {
      std::string bwt_filename = opt::prefix + BWT_EXT;
      std::string sai_filename = opt::prefix + SAI_EXT;
      BWTCA::runRopebwt(opt::readsFile, bwt_filename, use_threads, false);

      if(opt::bBuildSAI)
        {
	  std::cout << "\t done bwt construction, generating .sai file\n";
	  BWT* pBWT = new BWT(bwt_filename);
	  SampledSuffixArray ssa;
	  ssa.buildLexicoIndex(pBWT, opt::numThreads);
	  ssa.writeLexicoIndex(sai_filename);
	  delete pBWT;
        }
    }

  if(opt::bBuildReverse)
    {
      std::string rbwt_filename = opt::prefix + RBWT_EXT;
      std::string rsai_filename = opt::prefix + RSAI_EXT;
      BWTCA::runRopebwt(opt::readsFile, rbwt_filename, use_threads, true);

      if(opt::bBuildSAI)
        {
	  std::cout << "\t done rbwt construction, generating .rsai file\n";
	  BWT* pRBWT = new BWT(rbwt_filename);
	  SampledSuffixArray ssa;
	  ssa.buildLexicoIndex(pRBWT, opt::numThreads);
	  ssa.writeLexicoIndex(rsai_filename);
	  delete pRBWT;
        }
    }
}
//
void indexInMemorySAIS()
{
  std::cout << "Building index for " << opt::readsFile << " in memory using SAIS\n";

  if(opt::bBuildForward || opt::bBuildReverse)
    {
      // Parse the initial read table
      ReadTable* pRT = new ReadTable(opt::readsFile);

      // Create and write the suffix array for the forward reads
      if(opt::bBuildForward)
	{
	  buildIndexForTable(opt::prefix, pRT, false);
	}

      if(opt::bBuildReverse)
	{
	  // Reverse all the reads
	  pRT->reverseAll();

	  // Build the reverse suffix array
	  buildIndexForTable(opt::prefix, pRT, true);
	}

      delete pRT;
    }
}

//
void indexOnDisk()
{
  std::cout << "Building index for " << opt::readsFile << " on disk\n";
  BWTDiskParameters parameters;
  parameters.inFile = opt::readsFile;
  parameters.outPrefix = opt::prefix;
  parameters.bwtExtension = BWT_EXT;
  parameters.saiExtension = SAI_EXT;
  parameters.numReadsPerBatch = opt::numReadsPerBatch;
  parameters.numThreads = opt::numThreads;
  parameters.storageLevel = opt::gapArrayStorage;
  parameters.bBuildReverse = false;
  parameters.bUseBCR = (opt::algorithm == "bcr");
  
  if(opt::bBuildForward)
    {
      buildBWTDisk(parameters);
    }
  
  if(opt::bBuildReverse)
    {
      parameters.bwtExtension = RBWT_EXT;
      parameters.saiExtension = RSAI_EXT;
      parameters.bBuildReverse = true;
      buildBWTDisk(parameters);
    }
}

//
void buildIndexForTable(std::string prefix, const ReadTable* pRT, bool isReverse)
{
  // Create suffix array from read table
  SuffixArray* pSA = new SuffixArray(pRT, opt::numThreads);

  if(opt::validate)
    {
      std::cout << "Validating suffix array\n";
      pSA->validate(pRT);
    }

  std::string bwt_filename = prefix + (!isReverse ? BWT_EXT : RBWT_EXT);
  pSA->writeBWT(bwt_filename, pRT);

  std::string sufidx_filename = prefix + (!isReverse ? SAI_EXT : RSAI_EXT);
  pSA->writeIndex(sufidx_filename);

  delete pSA;
  pSA = NULL;
}

// 
// Handle command line arguments
//
void parseIndexOptions(int argc, char** argv)
{
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
      std::istringstream arg(optarg != NULL ? optarg : "");
      switch (c) 
        {
	case 'p': arg >> opt::prefix; break;
	case '?': die = true; break;
	case 'c': opt::validate = true; break;
	case 'd': opt::bDiskAlgo = true; arg >> opt::numReadsPerBatch; break;
	case 't': arg >> opt::numThreads; break;
	case 'g': arg >> opt::gapArrayStorage; break;
	case 'a': arg >> opt::algorithm; break;
	case 'v': opt::verbose++; break;
	case OPT_NO_REVERSE: opt::bBuildReverse = false; break;
	case OPT_NO_FWD: opt::bBuildForward = false; break;
	case OPT_NO_SAI: opt::bBuildSAI = false; break;
	case OPT_HELP:
	  std::cout << INDEX_USAGE_MESSAGE;
	  exit(EXIT_SUCCESS);
	case OPT_VERSION:
	  std::cout << INDEX_VERSION_MESSAGE;
	  exit(EXIT_SUCCESS);
        }
    }

  // Transform algorithm parameter to lower case
  std::transform(opt::algorithm.begin(), opt::algorithm.end(), opt::algorithm.begin(), ::tolower);
    
  if (argc - optind < 1) 
    {
      std::cerr << SUBPROGRAM ": missing arguments\n";
      die = true;
    } 
  else if (argc - optind > 1) 
    {
      std::cerr << SUBPROGRAM ": too many arguments\n";
      die = true;
    }

  if(opt::gapArrayStorage != 4 && opt::gapArrayStorage != 8 &&
     opt::gapArrayStorage != 16 && opt::gapArrayStorage != 32)
    {
      std::cerr << SUBPROGRAM ": invalid argument, --gap-array,-g must be one of 4,8,16,32 (found: " << opt::gapArrayStorage << ")\n";
      die = true;
    }

  if(opt::numThreads <= 0)
    {
      std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
      die = true;
    }

  if(opt::algorithm != "sais" && opt::algorithm != "bcr" && opt::algorithm != "ropebwt")
    {
      std::cerr << SUBPROGRAM ": unrecognized algorithm string " << opt::algorithm << ". --algorithm must be sais, bcr or ropebwt\n";
      die = true;
    }

  if(opt::algorithm == "ropebwt" && opt::bDiskAlgo)
    {
      std::cerr << SUBPROGRAM ": the options -a ropebwt and -d are not compatible, please only use one.\n";
      die = true;
    }

  if (die) 
    {
      std::cout << "\n" << INDEX_USAGE_MESSAGE;
      exit(EXIT_FAILURE);
    }

  // Parse the input filenames
  opt::readsFile = argv[optind++];
  if(opt::prefix.empty())
    opt::prefix = stripFilename(opt::readsFile);

  // Check if input file is empty
  size_t filesize = getFilesize(opt::readsFile);
  if(filesize == 0)
    {
      std::cerr << SUBPROGRAM ": input file is empty\n";
      exit(EXIT_FAILURE);
    }
}
