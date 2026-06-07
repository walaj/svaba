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
#include <string>
#include <string_view>
#include "SvabaOptions.h"
#include "SvabaGitVersion.h"   // generated: SVABA_GIT_HASH / DESCRIBE / DIRTY

void runToVCF(int argc, char** argv);
void runsvaba(int argc, char** argv);
void runRefilterBreakpoints(int argc, char** argv);
void runPostprocess(int argc, char** argv);
void runExtractPairs(int argc, char** argv);
// "Secret" benchmark subcommand: deliberately NOT listed in printUsage()
// below so it doesn't show up in the user-facing --help. Invoke with
// `svaba test -h` to see its own usage. Used to isolate the per-phase
// CPU cost of the main hot paths (walk / correct / assemble / align)
// independently of the real pipeline's bookkeeping and scoring.
void runTest(int argc, char** argv);

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
              << "  run            Run SV and indel detection on BAM(s)\n"
              << "  refilter       Re-run LOD/PASS filtering on an existing bps.txt.gz\n"
              << "  postprocess    Sort + streaming-dedup per-suffix output BAMs\n"
              << "  tovcf          Convert a deduped bps.txt.gz into VCFv4.5 output\n"
              << "  extract-pairs  Extract read pairs from a BAM by SEQ match (+ rev-comp)\n\n"
              << "Run 'svaba --version' to print the version and build commit.\n"
              << "Report issues at https://github.com/walaj/svaba/issues\n";
}

// Print the semantic version (from SvabaOptions.h) plus the git provenance
// of this build (from the CMake-generated SvabaGitVersion.h). The commit
// lines are omitted when svaba wasn't built from a git checkout.
static void printVersion() {
    std::cout << "svaba " << SVABA_VERSION << " (" << SVABA_DATE << ")\n";
    if (std::string_view(SVABA_GIT_HASH) != "unknown") {
        std::cout << "  commit:   " << SVABA_GIT_HASH
                  << (SVABA_GIT_DIRTY ? " (dirty)" : "") << "\n"
                  << "  describe: " << SVABA_GIT_DESCRIBE << "\n";
    }
    std::cout << "  built:    " << __DATE__ << " " << __TIME__ << "\n"
              << "  source:   https://github.com/walaj/svaba\n";
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
    else if (cmd == "--version" || cmd == "-v" || cmd == "version") {
        printVersion();
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
    else if (cmd == "extract-pairs") {
        return (runExtractPairs(argc - 1, argv + 1), EXIT_SUCCESS);
    }
    // Secret benchmark subcommand — intentionally not in the user-facing
    // command list. Entry point: `svaba test -h` for its own usage.
    else if (cmd == "test") {
        return (runTest(argc - 1, argv + 1), EXIT_SUCCESS);
    }
    else {
        std::cerr << "Unknown command: " << cmd << "\n\n";
        printUsage();
        return EXIT_FAILURE;
    }
    
    std::cerr << "Done with SvABA" << std::endl;
    return 0;

}
