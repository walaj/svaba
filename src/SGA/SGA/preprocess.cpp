//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// preprocess - prepare data files for assembly
//
#include <iostream>
#include <fstream>
#include "Util.h"
#include "preprocess.h"
#include "Timer.h"
#include "SeqReader.h"
#include "PrimerScreen.h"
#include "Alphabet.h"
#include "Quality.h"

static unsigned int DEFAULT_MIN_LENGTH = 40;
static int LOW_QUALITY_PHRED_SCORE = 3;

//
// Getopt
//
#define SUBPROGRAM "preprocess"
static const char *PREPROCESS_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *PREPROCESS_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] READS1 READS2 ...\n"
"Prepare READS1, READS2, ... data files for assembly\n"
"If pe-mode is turned on (pe-mode=1) then if a read is discarded its pair will be discarded as well.\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"\nInput/Output options:\n"
"      -o, --out=FILE                   write the reads to FILE (default: stdout)\n"
"      -p, --pe-mode=INT                0 - do not treat reads as paired (default)\n"
"                                       1 - reads are paired with the first read in READS1 and the second\n"
"                                       read in READS2. The paired reads will be interleaved in the output file\n"
"                                       2 - reads are paired and the records are interleaved within a single file.\n"
"          --pe-orphans=FILE            if one half of a read pair fails filtering, write the passed half to FILE\n"
"\nConversions/Filtering:\n"
"          --phred64                    convert quality values from phred-64 to phred-33.\n"
"          --discard-quality            do not output quality scores\n"
"      -q, --quality-trim=INT           perform Heng Li's BWA quality trim algorithm. \n"
"                                       Reads are trimmed according to the formula:\n"
"                                       argmax_x{\\sum_{i=x+1}^l(INT-q_i)} if q_l<INT\n"
"                                       where l is the original read length.\n"
"      -f, --quality-filter=INT         discard the read if it contains more than INT low-quality bases.\n"
"                                       Bases with phred score <= 3 are considered low quality. Default: no filtering.\n"
"                                       The filtering is applied after trimming so bases removed are not counted.\n"
"                                       Do not use this option if you are planning to use the BCR algorithm for indexing.\n"
"      -m, --min-length=INT             discard sequences that are shorter than INT\n"
"                                       this is most useful when used in conjunction with --quality-trim. Default: 40\n"
"      -h, --hard-clip=INT              clip all reads to be length INT. In most cases it is better to use\n"
"                                       the soft clip (quality-trim) option.\n"
"      --permute-ambiguous              Randomly change ambiguous base calls to one of possible bases.\n"
"                                       If this option is not specified, the entire read will be discarded.\n"
"      -s, --sample=FLOAT               Randomly sample reads or pairs with acceptance probability FLOAT.\n"
"      --dust                           Perform dust-style filtering of low complexity reads.\n"
"      --dust-threshold=FLOAT           filter out reads that have a dust score higher than FLOAT (default: 4.0).\n"
"      --suffix=SUFFIX                  append SUFFIX to each read ID\n"
"\nAdapter/Primer checks:\n"
"          --no-primer-check            disable the default check for primer sequences\n"
"      -r, --remove-adapter-fwd=STRING\n"
"      -c, --remove-adapter-rev=STRING  Remove the adapter STRING from input reads.\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum QualityScaling
{
    QS_UNDEFINED,
    QS_NONE,
    QS_SANGER,
    QS_PHRED64
};

namespace opt
{
    static unsigned int verbose;
    static std::string outFile;
    static unsigned int qualityTrim = 0;
    static unsigned int hardClip = 0;
    static unsigned int minLength = DEFAULT_MIN_LENGTH;
    static int qualityFilter = -1;
    static unsigned int peMode = 0;
    static double sampleFreq = 1.0f;

    static bool bDiscardAmbiguous = true;
    static bool bDiscardQuality = false;
    static QualityScaling qualityScale = QS_SANGER;

    static bool bFilterGC = false;
    static bool bDustFilter = false;
    static double dustThreshold = 4.0f;
    static std::string suffix;
    static double minGC = 0.0f;
    static double maxGC = 1.0;
    static bool bIlluminaScaling = false;
    static bool bDisablePrimerCheck = false;
    static std::string orphanFile;
    static std::string adapterF;  // adapter sequence forward
    static std::string adapterR; // adapter sequence reverse
}

static const char* shortopts = "o:q:m:h:p:r:c:s:f:vi";

enum { OPT_HELP = 1, OPT_VERSION, OPT_PERMUTE, 
       OPT_QSCALE, OPT_MINGC, OPT_MAXGC, 
       OPT_DUST, OPT_DUST_THRESHOLD, OPT_SUFFIX, 
       OPT_PHRED64, OPT_OUTPUTORPHANS, OPT_DISABLE_PRIMER,
       OPT_DISCARD_QUALITY };

static const struct option longopts[] = {
    { "verbose",                no_argument,       NULL, 'v' },
    { "out",                    required_argument, NULL, 'o' },
    { "quality-trim",           required_argument, NULL, 'q' },
    { "quality-filter",         required_argument, NULL, 'f' },
    { "pe-mode",                required_argument, NULL, 'p' },
    { "hard-clip",              required_argument, NULL, 'h' },
    { "min-length",             required_argument, NULL, 'm' },
    { "sample",                 required_argument, NULL, 's' },
    { "remove-adapter-fwd",     required_argument, NULL, 'r' },
    { "remove-adapter-rev",     required_argument, NULL, 'c' },
    { "dust",                   no_argument,       NULL, OPT_DUST},
    { "dust-threshold",         required_argument, NULL, OPT_DUST_THRESHOLD },
    { "suffix",                 required_argument, NULL, OPT_SUFFIX },
    { "phred64",                no_argument,       NULL, OPT_PHRED64 },
    { "pe-orphans",             required_argument, NULL, OPT_OUTPUTORPHANS },
    { "min-gc",                 required_argument, NULL, OPT_MINGC},
    { "max-gc",                 required_argument, NULL, OPT_MAXGC},
    { "help",                   no_argument,       NULL, OPT_HELP },
    { "version",                no_argument,       NULL, OPT_VERSION },
    { "permute-ambiguous",      no_argument,       NULL, OPT_PERMUTE },
    { "discard-quality",        no_argument,       NULL, OPT_DISCARD_QUALITY },
    { "no-primer-check",        no_argument,       NULL, OPT_DISABLE_PRIMER },
    { NULL, 0, NULL, 0 }
};

static int64_t s_numReadsRead = 0;
static int64_t s_numReadsKept = 0;
static int64_t s_numBasesRead = 0;
static int64_t s_numBasesKept = 0;
static int64_t s_numReadsPrimer = 0;
static int64_t s_numInvalidPE = 0;
static int64_t s_numFailedDust = 0;

//
// Main
//
int preprocessMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("sga preprocess");
    parsePreprocessOptions(argc, argv);

    std::cerr << "Parameters:\n";
    std::cerr << "QualTrim: " << opt::qualityTrim << "\n";

    if(opt::qualityFilter >= 0)
        std::cerr << "QualFilter: at most " << opt::qualityFilter << " low quality bases\n";
    else
        std::cerr << "QualFilter: no filtering\n";

    std::cerr << "HardClip: " << opt::hardClip << "\n";
    std::cerr << "Min length: " << opt::minLength << "\n";
    std::cerr << "Sample freq: " << opt::sampleFreq << "\n";
    std::cerr << "PE Mode: " << opt::peMode << "\n";
    std::cerr << "Quality scaling: " << opt::qualityScale << "\n";
    std::cerr << "MinGC: " << opt::minGC << "\n";
    std::cerr << "MaxGC: " << opt::maxGC << "\n";
    std::cerr << "Outfile: " << (opt::outFile.empty() ? "stdout" : opt::outFile) << "\n";
    std::cerr << "Orphan file: " << (opt::orphanFile.empty() ? "none" : opt::orphanFile) << "\n";
    if(opt::bDiscardAmbiguous)
        std::cerr << "Discarding sequences with ambiguous bases\n";
    if(opt::bDiscardQuality)
        std::cerr << "Discarding quality scores\n";
    if(opt::bDustFilter)
        std::cerr << "Dust threshold: " << opt::dustThreshold << "\n";
    if(!opt::suffix.empty())
        std::cerr << "Suffix: " << opt::suffix << "\n";

    if(opt::adapterF.length() && opt::adapterR.length())
    {
        std::cerr << "Adapter sequence fwd: " << opt::adapterF << "\n";
        std::cerr << "Adapter sequence rev: " << opt::adapterR << "\n";
    }

    // Seed the RNG
    srand(time(NULL));

    std::ostream* pWriter;
    if(opt::outFile.empty())
    {
        pWriter = &std::cout;
    }
    else
    {
        std::ostream* pFile = createWriter(opt::outFile);
        pWriter = pFile;
    }

    // Create a filehandle to write orphaned reads to, if necessary
    std::ostream* pOrphanWriter = NULL;
    if(!opt::orphanFile.empty())
        pOrphanWriter = createWriter(opt::orphanFile);

    if(opt::peMode == 0)
    {
        // Treat files as SE data
        while(optind < argc)
        {
            std::string filename = argv[optind++];
            std::cerr << "Processing " << filename << "\n\n";
            SeqReader reader(filename, SRF_NO_VALIDATION);
            SeqRecord record;

            while(reader.get(record))
            {
                bool passed = processRead(record);
                if(passed && samplePass())
                {
                    if(!opt::suffix.empty())
                        record.id.append(opt::suffix);

                    record.write(*pWriter);
                    ++s_numReadsKept;
                    s_numBasesKept += record.seq.length();
                }
            }
        }
    }
    else
    {
        assert(opt::peMode == 1 || opt::peMode == 2);
        int numFiles = argc - optind;
        if(opt::peMode == 1 && numFiles % 2 == 1)
        {
            std::cerr << "Error: An even number of files must be given for pe-mode 1\n";
            exit(EXIT_FAILURE);
        }

        while(optind < argc)
        {
            SeqReader* pReader1;
            SeqReader* pReader2;

            if(opt::peMode == 1)
            {
                // Read from separate files
                std::string filename1 = argv[optind++];
                std::string filename2 = argv[optind++];
                
                if(filename1 == "-" || filename2 == "-")
                {
                    std::cerr << "Reading from stdin is not supported in --pe-mode 1\n";
                    std::cerr << "Maybe you meant --pe-mode 2 (interleaved pairs?)\n";
                    exit(EXIT_FAILURE);
                }

                pReader1 = new SeqReader(filename1, SRF_NO_VALIDATION);
                pReader2 = new SeqReader(filename2, SRF_NO_VALIDATION);

                std::cerr << "Processing pe files " << filename1 << ", " << filename2 << "\n";

            }
            else
            {
                // Read from a single file
                std::string filename = argv[optind++];
                pReader1 = new SeqReader(filename, SRF_NO_VALIDATION);
                pReader2 = pReader1;
                std::cerr << "Processing interleaved pe file " << filename << "\n";
            }

            SeqRecord record1;
            SeqRecord record2;
            while(pReader1->get(record1) && pReader2->get(record2))
            {
                // If the names of the records are the same, append a /1 and /2 to them
                if(record1.id == record2.id)
                {
                    if(!opt::suffix.empty())
                    {
                        record1.id.append(opt::suffix);
                        record2.id.append(opt::suffix);
                    }

                    record1.id.append("/1");
                    record2.id.append("/2");
                }

                // Ensure the read names are sensible
                std::string expectedID2 = getPairID(record1.id);
                std::string expectedID1 = getPairID(record2.id);

                if(expectedID1 != record1.id || expectedID2 != record2.id)
                {
                    std::cerr << "Warning: Pair IDs do not match (expected format /1,/2 or /A,/B)\n";
                    std::cerr << "Read1 ID: " << record1.id << "\n";
                    std::cerr << "Read2 ID: " << record2.id << "\n";
                    s_numInvalidPE += 2;
                }

                bool passed1 = processRead(record1);
                bool passed2 = processRead(record2);

                if(!samplePass())
                    continue;

                if(passed1 && passed2)
                {
                    record1.write(*pWriter);
                    record2.write(*pWriter);
                    s_numReadsKept += 2;
                    s_numBasesKept += record1.seq.length();
                    s_numBasesKept += record2.seq.length();
                }
                else if(passed1 && pOrphanWriter != NULL)
                {
                    record1.write(*pOrphanWriter);
                }
                else if(passed2 && pOrphanWriter != NULL)
                {
                    record2.write(*pOrphanWriter);
                }
            }

            if(pReader2 != pReader1)
            {
                // only delete reader2 if it is a distinct pointer
                delete pReader2;
                pReader2 = NULL;
            }
            delete pReader1;
            pReader1 = NULL;
        }

    }

    if(pWriter != &std::cout)
        delete pWriter;
    if(pOrphanWriter != NULL)
        delete pOrphanWriter;

    std::cerr << "\nPreprocess stats:\n";
    std::cerr << "Reads parsed:\t" << s_numReadsRead << "\n";
    std::cerr << "Reads kept:\t" << s_numReadsKept << " (" << (double)s_numReadsKept / (double)s_numReadsRead << ")\n";
    std::cerr << "Reads failed primer screen:\t" << s_numReadsPrimer << " (" << (double)s_numReadsPrimer / (double)s_numReadsRead << ")\n";
    std::cerr << "Bases parsed:\t" << s_numBasesRead << "\n";
    std::cerr << "Bases kept:\t" << s_numBasesKept << " (" << (double)s_numBasesKept / (double)s_numBasesRead << ")\n";
    std::cerr << "Number of incorrectly paired reads that were discarded: " << s_numInvalidPE << "\n";
    if(opt::bDustFilter)
        std::cerr << "Number of reads failed dust filter: " << s_numFailedDust << "\n";
    delete pTimer;
    return 0;
}

// Process a single read by quality trimming, filtering
// returns true if the read should be kept
bool processRead(SeqRecord& record)
{

    // let's remove the adapter if the user has requested so
    // before doing any filtering
    if(!opt::adapterF.empty())
    {
        std::string _tmp(record.seq.toString());
        size_t found = _tmp.find(opt::adapterF);
        int _length;

        if(found != std::string::npos)
        {
            _length = opt::adapterF.length();
        }
        else
        { 
            // Couldn't find the fwd adapter; Try the reverse version
            found = _tmp.find(opt::adapterR);
           _length = opt::adapterR.length();
        }

        if(found != std::string::npos) // found the adapter
        {
            _tmp.erase(found, _length);
            record.seq = _tmp;

            // We have to remove the qualities of the adapter
            if(!record.qual.empty())
            {
                _tmp = record.qual;
                _tmp.erase(found, _length);
                record.qual = _tmp;
            }
        }
    }

    // Check if the sequence has uncalled bases
    std::string seqStr = record.seq.toString();
    std::string qualStr = record.qual;

    ++s_numReadsRead;
    s_numBasesRead += seqStr.size();

    // If ambiguity codes are present in the sequence
    // and the user wants to keep them, we randomly
    // select one of the DNA symbols from the set of
    // possible bases
    if(!opt::bDiscardAmbiguous) 
    {
        for(size_t i = 0; i < seqStr.size(); ++i)
        {
            // Convert '.' to 'N'
            if(seqStr[i] == '.')
                seqStr[i] = 'N';

            if(!IUPAC::isAmbiguous(seqStr[i]))
                continue;

            // Get the string of possible bases for this ambiguity code
            std::string possibles = IUPAC::getPossibleSymbols(seqStr[i]);

            // select one of the bases at random
            int j = rand() % possibles.size();
            seqStr[i] = possibles[j];
        }
    }

    // Ensure sequence is entirely ACGT
    size_t pos = seqStr.find_first_not_of("ACGT");
    if(pos != std::string::npos)
        return false;

    // Validate the quality string (if present) and
    // perform any necessary transformations
    if(!qualStr.empty())
    {
        // Calculate the range of phred scores for validation
        bool allValid = true;
        for(size_t i = 0; i < qualStr.size(); ++i)
        {
	  if(opt::qualityScale == QS_PHRED64)
                qualStr[i] = Quality::phred64toPhred33(qualStr[i]);
            allValid = Quality::isValidPhred33(qualStr[i]) && allValid;
        }

        if(!allValid)
        {
            std::cerr << "Error: read " << record.id << " has out of range quality values.\n";
            std::cerr << "Expected phred" << (opt::qualityScale == QS_SANGER ? "33" : "64") << ".\n";
            std::cerr << "Quality string: "  << qualStr << "\n";
            std::cerr << "Check your data and re-run preprocess with the correct quality scaling flag.\n";
            exit(EXIT_FAILURE);
        }
    }

    // Hard clip
    if(opt::hardClip > 0)
    {
        seqStr = seqStr.substr(0, opt::hardClip);
        if(!qualStr.empty())
            qualStr = qualStr.substr(0, opt::hardClip);
    }

    // Quality trim
    if(opt::qualityTrim > 0 && !qualStr.empty())
        softClip(opt::qualityTrim, seqStr, qualStr);

    // Quality filter
    if(opt::qualityFilter >= 0 && !qualStr.empty())
    {
        int numLowQuality = countLowQuality(seqStr, qualStr);
        if(numLowQuality > opt::qualityFilter)
            return false;
    }

    // Dust filter
    if(opt::bDustFilter)
    {
        double dustScore = calculateDustScore(seqStr);
        bool bAcceptDust = dustScore < opt::dustThreshold;

        if(!bAcceptDust)
        {
            s_numFailedDust += 1;
            if(opt::verbose >= 1)
            {
                printf("Failed dust: %s %s %lf\n", record.id.c_str(),
                                                   seqStr.c_str(),
                                                   dustScore);
            }
            return false;
        }
    }

    // Filter by GC content
    if(opt::bFilterGC)
    {
        double gc = calcGC(seqStr);
        if(gc < opt::minGC || gc > opt::maxGC)
            return false;
    }

    // Primer screen
    if(!opt::bDisablePrimerCheck)
    {
        bool containsPrimer = PrimerScreen::containsPrimer(seqStr);
        if(containsPrimer)
        {
            ++s_numReadsPrimer;
            return false;
        }
    }

    record.seq = seqStr;

    if(opt::bDiscardQuality)
        record.qual.clear();
    else
        record.qual = qualStr;

    if(record.seq.length() == 0 || record.seq.length() < opt::minLength)
        return false;

    return true;
}

// return true if the random value is lower than the acceptance value
bool samplePass()
{
    if(opt::sampleFreq >= 1.0f)
        return true; // no sampling

    double r = rand() / (RAND_MAX + 1.0f);
    return r < opt::sampleFreq;
}

// Perform a soft-clipping of the sequence by removing low quality bases from the
// 3' end using Heng Li's algorithm from bwa
void softClip(int qualTrim, std::string& seq, std::string& qual)
{
    assert(seq.size() == qual.size());

    int endpoint = 0; // not inclusive
    int max = 0;
    int i = seq.length() - 1;
    int terminalScore = Quality::char2phred(qual[i]);
    // Only perform soft-clipping if the last base has qual less than qualTrim
    if(terminalScore >= qualTrim)
        return;

    int subSum = 0;
    while(i >= 0)
    {
        int ps = Quality::char2phred(qual[i]);
        int score = qualTrim - ps;
        subSum += score;
        if(subSum > max)
        {
            max = subSum;
            endpoint = i;
        }
        --i;
    }

    // Clip the read
    seq = seq.substr(0, endpoint);
    qual = qual.substr(0, endpoint);
}

// Count the number of low quality bases in the read
int countLowQuality(const std::string& seq, const std::string& qual)
{
    assert(seq.size() == qual.size());

    int sum = 0;
    for(size_t i = 0; i < seq.length(); ++i)
    {
        int ps = Quality::char2phred(qual[i]);
        if(ps <= LOW_QUALITY_PHRED_SCORE)
            ++sum;
    }
    return sum;
}

double calcGC(const std::string& seq)
{
    double num_gc = 0.0f;
    double num_total = 0.0f;
    for(size_t i = 0; i < seq.size(); ++i)
    {
        if(seq[i] == 'C' || seq[i] == 'G')
            ++num_gc;
        ++num_total;
    }
    return num_gc / num_total;
}

//
// Handle command line arguments
//
void parsePreprocessOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'o': arg >> opt::outFile; break;
            case 'q': arg >> opt::qualityTrim; break;
            case 'f': arg >> opt::qualityFilter; break;
            case 'i': arg >> opt::bIlluminaScaling; break;
            case 'm': arg >> opt::minLength; break;
            case 'h': arg >> opt::hardClip; break;
            case 'p': arg >> opt::peMode; break;
            case 'r': arg >> opt::adapterF; break;
            case 'c': arg >> opt::adapterR; break;
            case 's': arg >> opt::sampleFreq; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_DUST_THRESHOLD: arg >> opt::dustThreshold; opt::bDustFilter = true; break;
            case OPT_SUFFIX: arg >> opt::suffix; break;
            case OPT_MINGC: arg >> opt::minGC; opt::bFilterGC = true; break;
            case OPT_OUTPUTORPHANS: arg >> opt::orphanFile; break;
            case OPT_MAXGC: arg >> opt::maxGC; opt::bFilterGC = true; break;
            case OPT_PHRED64: opt::qualityScale = QS_PHRED64; break;
            case OPT_PERMUTE: opt::bDiscardAmbiguous = false; break;
            case OPT_DUST: opt::bDustFilter = true; break;
            case OPT_DISABLE_PRIMER: opt::bDisablePrimerCheck = true; break;
            case OPT_DISCARD_QUALITY: opt::bDiscardQuality = true; break;
            case OPT_HELP:
                std::cout << PREPROCESS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << PREPROCESS_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 1)
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    }

    if (die)
    {
        std::cout << "\n" << PREPROCESS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    if(opt::peMode > 2)
    {
        std::cerr << SUBPROGRAM ": error pe-mode must be 0,1 or 2 (found: " << opt::peMode << ")\n";
        exit(EXIT_FAILURE);
    }

    if(opt::minLength < DEFAULT_MIN_LENGTH)
    {
        std::cerr << SUBPROGRAM ": WARNING - it is suggested that the min read length is " << DEFAULT_MIN_LENGTH << "\n";
        std::cerr << SUBPROGRAM ": Using very short reads may considerably impact the performance\n";
    }

    if(opt::adapterF.empty() != opt::adapterR.empty())
    {
        std::cerr << SUBPROGRAM ": Forward and Reverse sequence is necessary to perform adapter removal.\n";
        exit(EXIT_FAILURE);
    }

    if(!opt::outFile.empty() && opt::outFile == opt::orphanFile)
    {
        std::cerr << SUBPROGRAM ": Output file and orphan file must be different\n";
        exit(EXIT_FAILURE);
    }
}
