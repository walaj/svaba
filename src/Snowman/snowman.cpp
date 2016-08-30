/* Snowman - Somatic Structural Variation Dectection
 * Copyright 2014 Broad Institute 
 * Written by Jeremiah Wala (jwala@broadinstitute.org)
 * Released under the included license detailed below.
 *
 * Snowman incorportes the core of String Graph Assembler, 
 * -- String Graph Assembler -- Copyright 2009 Wellcome Trust Sanger Institute
 * -- Written by Jared Simpson
 * -- Released under the GPL
 */

#include "refilter.h"
#include "run_snowman.h"

#define AUTHOR "Jeremiah Wala <jwala@broadinstitute.org>"

static const char *SNOWMAN_USAGE_MESSAGE =
"--------------------------------------------------------\n"
"--- Snowman - Structural Variant and Indel Detection ---\n"
"--------------------------------------------------------\n"
"Program: SnowmanSV \n"
"FH Version: 110 \n"
"Contact: Jeremiah Wala [ jwala@broadinstitute.org ]\n"
"Usage: snowman <command> [options]\n\n"
"Commands:\n"
"           run            Run Snowman SV and Indel detection on BAM(s)\n"
"           refilter       Refilter the Snowman breakpoints with additional/different criteria to created filtered VCF and breakpoints file.\n"
"\nReport bugs to jwala@broadinstitute.org \n\n";

int main(int argc, char** argv) {

  if (argc <= 1) {
    std::cerr << SNOWMAN_USAGE_MESSAGE;
    return 0;
  } else {
    std::string command(argv[1]);
    if (command == "help" || command == "--help") {
      std::cerr << SNOWMAN_USAGE_MESSAGE;
      return 0;
    } else if (command == "run") {
      runSnowman(argc -1, argv + 1);
    } else if (command == "refilter") {
      runRefilterBreakpoints(argc-1, argv+1);
    }
    else {
      std::cerr << SNOWMAN_USAGE_MESSAGE;
      return 0;
    }
  } 
  
  std::cerr << "Done with snowman" << std::endl;
  return 0;

}
