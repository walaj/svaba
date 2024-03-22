/* SvABA - Somatic Structural Variation Dectection
 * Copyright 2014 Broad Institute 
 * Written by Jeremiah Wala (jwala@broadinstitute.org)
 * Released under the included license detailed below.
 *
 * SvABA incorportes the core of String Graph Assembler, 
 * -- String Graph Assembler -- Copyright 2009 Wellcome Trust Sanger Institute
 * -- Written by Jared Simpson
 * -- Released under the GPL
 */

#include "refilter.h"
#include "run_svaba.h"

void runToVCF(int argc, char** argv);

#define AUTHOR "Jeremiah Wala <jeremiah.wala@gmail.com"

static const char *SVABA_USAGE_MESSAGE =
"------------------------------------------------------------\n"
"-------- SvABA - SV and indel detection by assembly --------\n"
"------------------------------------------------------------\n"
"Program: SvABA \n"
"Version: " SVABA_VERSION " - March 2024\n"
"Contact: Jeremiah Wala [ jeremiah.wala@gmail.org ]\n"
"Usage: svaba <command> [options]\n\n"
"Commands:\n"
"           run            Run SvABA SV and Indel detection on BAM(s)\n"
"           refilter       Refilter the SvABA breakpoints with additional/different criteria to created filtered VCF and breakpoints file.\n"
"           tovcf          Convert the bps.txt.gz file to a VCF\n"
"\nReport bugs to jwala@broadinstitute.org \n\n";

int main(int argc, char** argv) {

  if (argc <= 1) {
    std::cerr << SVABA_USAGE_MESSAGE;
    return 0;
  } else {
    std::string command(argv[1]);
    if (command == "help" || command == "--help") {
      std::cerr << SVABA_USAGE_MESSAGE;
      return 0;
    } else if (command == "run") {
      runsvaba(argc -1, argv + 1);
    } else if (command == "refilter") {
      runRefilterBreakpoints(argc-1, argv+1);
    } else if (command == "tovcf") {
      runToVCF(argc-1, argv+1);
    }
    else {
      std::cerr << SVABA_USAGE_MESSAGE;
      return 0;
    }
  } 
  
  std::cerr << "Done with SvABA" << std::endl;
  return 0;

}
