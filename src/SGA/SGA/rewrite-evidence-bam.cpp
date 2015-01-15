//-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// rewrite-evidence-bam - fill in read name and quality
// information in an SGA graph-diff evidence bam file
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include "Util.h"
#include "SeqReader.h"
#include "rewrite-evidence-bam.h"
#include "api/algorithms/Sort.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"

//
// Getopt
//
#define SUBPROGRAM "rewrite-evidence-bam"
static const char *REWRITE_EVIDENCE_BAM_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2012 Wellcome Trust Sanger Institute\n";

static const char *REWRITE_EVIDENCE_BAM_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... EVIDENCE_BAM_FILE\n"
"Discard mate-pair alignments from a BAM file that are potentially erroneous\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -f, --fastq=FILE                 parse the read names and sequences from the fastq file in F (required)\n"
"      -m, --merge-bam=FILE             merge the evidence BAM alignments with the alignments in F\n"
"      -o, --outfile=FILE               write the new BAM file to F\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string outFile;
    static std::string inputFastqFile;
    static std::string inputEvidenceBAMFile;
    static std::string mergeWithBAM;
}

static const char* shortopts = "o:m:f:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",            no_argument,         NULL, 'v' },
    { "outfile",            required_argument,   NULL, 'o' },
    { "merge-bam",          required_argument,   NULL, 'm' },
    { "fastq",              required_argument,   NULL, 'f' },
    { "help",               no_argument,         NULL, OPT_HELP },
    { "version",            no_argument,         NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int rewriteEvidenceBAMMain(int argc, char** argv)
{
    parseRewriteEvidenceBAMOptions(argc, argv);

    // Open the bam files for reading/writing
    BamTools::BamReader* pBamReader = new BamTools::BamReader;
    pBamReader->Open(opt::inputEvidenceBAMFile);


    if(!opt::mergeWithBAM.empty())
        std::cout << "Should merge alignments into: " << opt::mergeWithBAM << "\n";

    // Load alignments into a vector
    std::vector<BamTools::BamAlignment> evidence_alignments;
    BamTools::BamAlignment record;
    while(pBamReader->GetNextAlignment(record))
        evidence_alignments.push_back(record);
    printf("Read %zu evidence alignments\n", evidence_alignments.size());

    // Build a read index->alignment index map
    std::map<size_t, size_t> read_index_to_alignment_map;
    for(size_t i = 0; i < evidence_alignments.size(); ++i)
    {
        // We assume the name of the read that SGA writes is prefixed with "idx-"
        assert(evidence_alignments[i].Name.find("idx-") == 0);
        std::stringstream parser(evidence_alignments[i].Name.substr(4));
        size_t read_index;
        parser >> read_index;
        read_index_to_alignment_map.insert(std::make_pair(read_index, i));
    }

    // Read the original FASTQ file and file in read name/quality for the alignments
    SeqReader seqReader(opt::inputFastqFile);
    SeqRecord seqRecord;
    size_t read_index = 0;
    while(seqReader.get(seqRecord))
    {
        // Check if this read is in the evidence alignment set
        std::map<size_t, size_t>::iterator iter = read_index_to_alignment_map.find(read_index);
        if(iter != read_index_to_alignment_map.end())
        {
            BamTools::BamAlignment& alignment = evidence_alignments[iter->second];

            // Flip the quality string if the read sequence is flipped
            bool isRC = seqRecord.seq.toString() == reverseComplement(alignment.QueryBases);
            if(isRC)
                std::reverse(seqRecord.qual.begin(), seqRecord.qual.end());

            // Sanity check that the sequences match
            if(!isRC && seqRecord.seq.toString() != alignment.QueryBases)
            {
                std::cerr << "Read at index " << read_index << " in the FASTQ file does not match alignment with idx-" << read_index << "\n";
                exit(EXIT_FAILURE);
            }
            
            // Rewrite name
            alignment.Name = seqRecord.id;
            if(seqRecord.hasQuality())
                alignment.Qualities = seqRecord.qual;
        }
        read_index += 1;
    }

    // Resort the alignment vector by reference position
    std::sort(evidence_alignments.begin(), evidence_alignments.end(), BamTools::Algorithms::Sort::ByPosition());

    // Write the result
    BamTools::BamWriter* pBamWriter = new BamTools::BamWriter;
    pBamWriter->Open(opt::outFile, pBamReader->GetHeaderText(), pBamReader->GetReferenceData());
    for(size_t i = 0; i < evidence_alignments.size(); ++i)
    {
        pBamWriter->SaveAlignment(evidence_alignments[i]);
    }
   

    pBamWriter->Close();
    pBamReader->Close();

    delete pBamReader;
    delete pBamWriter;
    return 0;
}

// 
// Handle command line arguments
//
void parseRewriteEvidenceBAMOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'o': arg >> opt::outFile; break;
            case 'f': arg >> opt::inputFastqFile; break;
            case 'm': arg >> opt::mergeWithBAM; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << REWRITE_EVIDENCE_BAM_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << REWRITE_EVIDENCE_BAM_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 1) 
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } 
    else if(argc - optind > 1)
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if(opt::outFile.empty())
    {
        std::cerr << SUBPROGRAM ": --outfile is required\n";
        die = true;
    }

    if(opt::inputFastqFile.empty())
    {
        std::cerr << SUBPROGRAM ": --fastq is required\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << REWRITE_EVIDENCE_BAM_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input evidence bam filename
    opt::inputEvidenceBAMFile = argv[optind++];
}
