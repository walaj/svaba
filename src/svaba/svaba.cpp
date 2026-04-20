/* SvABA - Somatic Structural Variation Dectection
 * Copyright 2014 Broad Institute 
 * Written by Jeremiah Wala (jeremiah.wala@gmail.com)
 * Released under the included license detailed below.
 *
 * SvABA incorportes the core of String Graph Assembler, 
 * -- String Graph Assembler -- Copyright 2009 Wellcome Trust Sanger Institute
 * -- Written by Jared Simpson
 * -- Released under the GPL
 */

// svaba.cpp
#include <iostream>
#include <string_view>
#include "SvabaOptions.h"

void runToVCF(int argc, char** argv);
void runsvaba(int argc, char** argv);
void runRefilterBreakpoints(int argc, char** argv);
void runPostprocess(int argc, char** argv);

static void printUsage() {
    constexpr std::string_view header =
        "------------------------------------------------------------\n"
        "-------- SvABA - SV and indel detection by assembly --------\n"
        "------------------------------------------------------------\n";
    std::cout << header
              << "Program: SvABA\n"
              << "Version: " << SVABA_VERSION << " - " << SVABA_DATE << "\n\n"
              << "Usage: svaba <command> [options]\n\n"
              << "Commands:\n"
              << "  run          Run SV and indel detection on BAM(s)\n"
              << "  refilter     Re-run LOD/PASS filtering on an existing bps.txt.gz\n"
              << "  postprocess  Sort + streaming-dedup per-suffix output BAMs\n"
              << "  tovcf        Convert a deduped bps.txt.gz into VCFv4.5 output\n\n"
              << "Report issues at https://github.com/walaj/svaba/issues\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printUsage();
        return EXIT_SUCCESS;
    }

    std::string_view cmd = argv[1];

    if (cmd == "help" || cmd == "--help") {
        printUsage();
        return EXIT_SUCCESS;
    }
    else if (cmd == "run") {
        // strip off the run before passing to parser
        return (runsvaba(argc - 1, argv + 1), EXIT_SUCCESS);
    }
    else if (cmd == "refilter") {
        return (runRefilterBreakpoints(argc - 1, argv + 1), EXIT_SUCCESS);
    }
    else if (cmd == "postprocess") {
        return (runPostprocess(argc - 1, argv + 1), EXIT_SUCCESS);
    }
    else if (cmd == "tovcf") {
        return (runToVCF(argc - 1, argv + 1), EXIT_SUCCESS);
    }
    else {
        std::cerr << "Unknown command: " << cmd << "\n\n";
        printUsage();
        return EXIT_FAILURE;
    }
    
    std::cerr << "Done with SvABA" << std::endl;
    return 0;

}
