#include "benchmark.h"
#include "assembly2vcf.h"
#include "splitcounter.h"

#define AUTHOR "Jeremiah Wala <jwala@broadinstitute.org>"

static const char *SUTILS_USAGE_MESSAGE =
"-----------------------------------------------------------------------\n"
"--- Snowman Utils - Utililty methods for analyzing/benchmarking SVs ---\n"
"-----------------------------------------------------------------------\n"
"Program: SnowmanUtils \n"
"Contact: Jeremiah Wala [ jwala@broadinstitute.org ]\n"
"Usage: snowmanutils <command> [options]\n\n"
"Commands:\n"
"           benchmark      Run benchmarking tests for Snowman\n"
"           assembly2vcf   Run Snowman filtering on BWA-MEM aligned assembly (eg Discovar, SGA, etc) and aligned reads\n"
"           splitcounter   Simple utility to count locations of split-alignment breakpoints (eg for PacBio)\n"
"           splitfasta     Simple utility to split a fasta into smaller subsequences, splitting seq in the middle\n"
"\nReport bugs to jwala@broadinstitute.org \n\n";

int main(int argc, char** argv) {
  
  if (argc <= 1) {
    std::cerr << SUTILS_USAGE_MESSAGE;
    return 0;
  } else {
    std::string command(argv[1]);
    if (command == "help" || command == "--help") {
      std::cerr << SUTILS_USAGE_MESSAGE;
      return 0;
    } else if (command == "benchmark") {
      runBenchmark(argc-1, argv+1);
    } else if (command == "splitfasta") {
      runSplitFasta(argc-1, argv+1);
    }
    else if (command == "assembly2vcf") {
      runAssembly2VCF(argc-1, argv+1);
    } else if (command == "splitcounter") {
      runSplitCounter(argc-1, argv+1);
    }
    else {
      std::cerr << SUTILS_USAGE_MESSAGE;
      return 0;
    }
  } 
  
  std::cerr << "Done with snowman utils" << std::endl;
  return 0;

}
