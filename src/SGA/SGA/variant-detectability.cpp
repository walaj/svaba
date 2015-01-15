//-----------------------------------------------
//-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Detect the power to detect variants at specific
// locations of a reference genome using k-mers
//
#include <iostream>
#include <fstream>
#include <algorithm>
#include "SGACommon.h"
#include "Util.h"
#include "variant-detectability.h"
#include "BWTAlgorithms.h"
#include "BWTIndexSet.h"
#include "Timer.h"

// Structs
struct KmerCounts
{
    size_t total;
    size_t zero;
};

// Local functions
void computeDetectableAll(StringVector& ref_sequences, const BWTIndexSet& ref_index);
void computeDetectableSampling(StringVector& ref_sequences, const BWTIndexSet& ref_index);
KmerCounts computeChangeCounts(const BWTIndexSet& ref_index, std::string& sequence, size_t base_idx, char new_base);

//
// Getopt
//
#define SUBPROGRAM "index"

static const char *VARIANT_DETECTABILITY_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *VARIANT_DETECTABILITY_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Index the reads in READSFILE using a suffixarray/bwt\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"  -k, --kmer=K                         set the k-mer length\n"
"  -n, --num-samples=N                  perform the calculation by randomly sampling N mutations\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string referenceFile;
    static size_t kmer = 31;
    static size_t num_samples = 10000;
}

static const char* shortopts = "k:n:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_NO_REVERSE };

static const struct option longopts[] = {
    { "kmer",        required_argument, NULL, 'k' },
    { "num-samples", required_argument, NULL, 'n' },
    { "verbose",     no_argument,       NULL, 'v' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

int variantDetectabilityMain(int argc, char** argv)
{
    parseVarDetectOptions(argc, argv);

    // Load the reference BWT
    std::string bwt_name = stripExtension(opt::referenceFile) + BWT_EXT;
    BWTIndexSet ref_index;
    ref_index.pBWT = new BWT(bwt_name);
    ref_index.pCache = new BWTIntervalCache(11, ref_index.pBWT);

    // Read reference
    ReadTable ref_table(opt::referenceFile);
        
    // Convert to string vector
    StringVector ref_sequences;
    for(size_t i = 0; i < ref_table.getCount(); ++i) {
        ref_sequences.push_back(ref_table.getRead(i).seq.toString());   
    }
    computeDetectableSampling(ref_sequences, ref_index);

    delete ref_index.pBWT;
    delete ref_index.pCache;
    return 0;
}

//
void computeDetectableAll(StringVector& ref_sequences, const BWTIndexSet& ref_index)
{
    size_t total_tested = 0;
    size_t total_detected = 0;

    for(size_t ri = 0; ri < ref_sequences.size(); ++ri) {
        std::string& sequence = ref_sequences[ri];
        size_t l = sequence.length();
        for(size_t i = 0; i < l; ++i) {

            char b = sequence[i];
            for(size_t j = 0; j < 4; ++j) {
                char m = "ACGT"[j];
                if(b == m)
                    continue;
                KmerCounts counts = computeChangeCounts(ref_index, sequence, i, m);
                total_tested += 1;
                total_detected += (counts.zero > 0) ? 1 : 0;
            }
        }

        if(total_tested % 1000 == 0)
            printf("Tested %zu\n", total_tested);
    }
    printf("Tested: %zu\n", total_tested);
    printf("Detectable: %zu\n", total_detected);
}

void computeDetectableSampling(StringVector& ref_sequences, const BWTIndexSet& ref_index)
{
    size_t total_tested = 0;
    size_t total_detected = 0;
    size_t total_complete = 0;
    size_t num_ref = ref_sequences.size();

    // seed rng
    srandom( time(NULL) );

    while(total_tested < opt::num_samples) {

        // Choose a random chromosome
        size_t chr_idx = random() % num_ref;
        std::string& sequence = ref_sequences[chr_idx];
        size_t l = sequence.length();

        // Choose a random position
        size_t base_idx = random() % l;
        char b = sequence[base_idx];

        // Choose a random mutation
        char m = b;
        do {
            int j = random() % 4;
            m = "ACGT"[j];
        } while(m == b);

        KmerCounts counts = computeChangeCounts(ref_index, sequence, base_idx, m);
        total_tested += 1;
        total_detected += (counts.zero > 0) ? 1 : 0;
        total_complete += (counts.zero == counts.total) ? 1 : 0;
        if(opt::verbose > 0)
            printf("tt: %zu td: %zu tc: %zu\n", total_tested, total_detected, total_complete);
    }

    printf("k=%zu tested=%zu detected=%zu complete=%zu\n", opt::kmer, total_tested, total_detected, total_complete);
}


// Returns true if changing the given reference base is detectable with kmers
KmerCounts computeChangeCounts(const BWTIndexSet& ref_index, std::string& sequence, size_t base_idx, char new_base)
{
    // Introduce the change
    char old_base = sequence[base_idx];
    sequence[base_idx] = new_base;
    size_t l = sequence.length();

    // Iterate over kmers covering this position
    size_t start_k_idx = (base_idx + 1) > opt::kmer ? base_idx + 1 - opt::kmer : 0;
    size_t end_k_idx = (base_idx + opt::kmer) < l ? base_idx : l - opt::kmer;
    assert(end_k_idx - start_k_idx <= opt::kmer);

    KmerCounts counts;
    counts.total = 0;
    counts.zero = 0;
    for(size_t ki = start_k_idx; ki <= end_k_idx; ++ki) {
        std::string ks = sequence.substr(ki, opt::kmer);
        size_t occ = BWTAlgorithms::countSequenceOccurrences(ks, ref_index) + 
                     BWTAlgorithms::countSequenceOccurrences(reverseComplement(ks), ref_index);

        counts.total += 1;
        counts.zero += (occ == 0) ? 1 : 0;
    }

    // Reset the base
    sequence[base_idx] = old_base;
    return counts;
}

// 
// Handle command line arguments
//
void parseVarDetectOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case '?': die = true; break;
            case 'k': arg >> opt::kmer; break;
            case 'n': arg >> opt::num_samples; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << VARIANT_DETECTABILITY_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << VARIANT_DETECTABILITY_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
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

    if(die) 
    {
        std::cerr << "Try `" << SUBPROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::referenceFile = argv[optind++];
}
