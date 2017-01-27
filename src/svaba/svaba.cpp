/* SvAbA - Somatic Structural Variation Dectection
 * Copyright 2014 Broad Institute 
 * Written by Jeremiah Wala (jwala@broadinstitute.org)
 * Released under the included license detailed below.
 *
 * SvAbA incorportes the core of String Graph Assembler, 
 * -- String Graph Assembler -- Copyright 2009 Wellcome Trust Sanger Institute
 * -- Written by Jared Simpson
 * -- Released under the GPL
 */

#include "refilter.h"
#include "run_svaba.h"

#define AUTHOR "Jeremiah Wala <jwala@broadinstitute.org>"

static const char *SVABA_USAGE_MESSAGE =
"------------------------------------------------------------\n"
"--- SvAbA (sah-bah) - SV and indel detection by assembly ---\n"
"------------------------------------------------------------\n"
"Program: SvAbA \n"
"FH Version: 134 \n"
"Contact: Jeremiah Wala [ jwala@broadinstitute.org ]\n"
"Usage: svaba <command> [options]\n\n"
"Commands:\n"
"           run            Run SvAbA SV and Indel detection on BAM(s)\n"
"           refilter       Refilter the SvAbA breakpoints with additional/different criteria to created filtered VCF and breakpoints file.\n"
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
    }
    else {
      std::cerr << SVABA_USAGE_MESSAGE;
      return 0;
    }
  } 
  
  std::cerr << "Done with SvAbA" << std::endl;
  return 0;

}
