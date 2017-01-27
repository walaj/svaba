/* SVaBA - Somatic Structural Variation Dectection
 * Copyright 2014 Broad Institute 
 * Written by Jeremiah Wala (jwala@broadinstitute.org)
 * Released under the included license detailed below.
 *
 * SVaBA incorportes the core of String Graph Assembler, 
 * -- String Graph Assembler -- Copyright 2009 Wellcome Trust Sanger Institute
 * -- Written by Jared Simpson
 * -- Released under the GPL
 */

#include "refilter.h"
#include "run_svaba.h"

#define AUTHOR "Jeremiah Wala <jwala@broadinstitute.org>"

static const char *SVABA_USAGE_MESSAGE =
"--------------------------------------------------------\n"
"--- SVaBA - Structural Variant and Indel Detection ---\n"
"--------------------------------------------------------\n"
"Program: SVaBA \n"
"FH Version: 134 \n"
"Contact: Jeremiah Wala [ jwala@broadinstitute.org ]\n"
"Usage: svaba <command> [options]\n\n"
"Commands:\n"
"           run            Run SVaBA SV and Indel detection on BAM(s)\n"
"           refilter       Refilter the SVaBA breakpoints with additional/different criteria to created filtered VCF and breakpoints file.\n"
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
      runSVaBA(argc -1, argv + 1);
    } else if (command == "refilter") {
      runRefilterBreakpoints(argc-1, argv+1);
    }
    else {
      std::cerr << SVABA_USAGE_MESSAGE;
      return 0;
    }
  } 
  
  std::cerr << "Done with SVaBA" << std::endl;
  return 0;

}
