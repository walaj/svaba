//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// graph-concordance - check variants in a VCF file
// against an assembly graph
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "graph-concordance.h"
#include "BWTAlgorithms.h"
#include "BWTIndexSet.h"
#include "SGACommon.h"
#include "Timer.h"
#include "SequenceProcessFramework.h"
#include "HapgenUtil.h"
#include "SGAStats.h"
#include "HashMap.h"
#include "DindelRealignWindow.h"
#include "VariantIndex.h"
#include "HapgenUtil.h"

//
// Getopt
//
#define SUBPROGRAM "graph-concordance"
static const char *GRAPH_CONCORDANCE_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2013 Wellcome Trust Sanger Institute\n";

static const char *GRAPH_CONCORDANCE_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... -r variant.fastq -b base.fastq VCF_FILE\n"
"Count read support for variants in a vcf file\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"          --reference=STR              load the reference genome from FILE\n"
"      -t, --threads=NUM                use NUM threads to compute the overlaps (default: 1)\n"
"      -g, --germline=FILE              load germline variants from FILE\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int k = 51;
    static std::string variantFile;
    static std::string baseFile;
    static std::string haplotypeFile;
    static std::string vcfFile;
    static std::string referenceFile;
    static std::string germlineFile;
}

static const char* shortopts = "b:r:o:k:t:g:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_REFERENCE, OPT_HAPLOID };

static const struct option longopts[] = {
    { "verbose",               no_argument,       NULL, 'v' },
    { "threads",               required_argument, NULL, 't' },
    { "outfile",               required_argument, NULL, 'o' },
    { "reads",                 required_argument, NULL, 'r' },
    { "base",                  required_argument, NULL, 'b' },
    { "variant",               required_argument, NULL, 'r' },
    { "germline",              required_argument, NULL, 'g' },
    { "reference",             required_argument, NULL, OPT_REFERENCE },
    { "haploid",               no_argument,       NULL, OPT_HAPLOID },
    { "help",                  no_argument,       NULL, OPT_HELP },
    { "version",               no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

std::string applyVariant(const std::string& in, int pos,
                         const std::string& ref, const std::string& alt,
                         std::vector<int>& coordinate_map, 
                         bool update_map, int variant_tag = -1)
{
    assert(variant_tag < 0);
    std::string out = in;

    // Search the coordinate map for the index of the ref string to be substituted
    size_t idx = 0;
    while(idx < coordinate_map.size())
    {
        if(coordinate_map[idx] == pos)
            break;
        idx += 1;
    }

    // If we could not find the position in the coordinate map
    // then the input variant is for a reference base that has been removed
    // already. Skip it.
    if(idx == coordinate_map.size())
        return "";

    /*
    std::cout << "In: " << in << "-\n";
    std::cout << "Current map: \n";
    for(size_t i = 0; i < coordinate_map.size(); ++i)
        std::cout << i << " " << coordinate_map[i] << "\n";
    std::cout << "\n";
    std::cout << "Replace idx: " << idx << "\n";
    std::cout << "Replace coordinate: " << coordinate_map[idx] << "\n";
    std::cout << "ref: " << ref << "\n";
    std::cout << "Alt: " << alt << "\n";
    */

    // Ensure that the reference string at the variant matches the expected
    if(out.substr(idx, ref.length()) != ref)
    {
        std::cerr << "Warning: REF does not match reference\n";
        return "";
    }

    // Replace the reference sequence with the alternative
    out.replace(idx, ref.length(), alt);

    if(update_map)
    {
        // shift coordinate map to account for indels

        // Make an iterator to the first position that is replaced
        std::vector<int>::iterator fi = coordinate_map.begin() + idx;

        // Iterator which is one past the last position to remove
        std::vector<int>::iterator li = fi + ref.size();
        
        // erase returns an iterator pointing to the position
        // following the last position removed, this is where
        // we want to insert 
        std::vector<int>::iterator ii = coordinate_map.erase(fi, li);

        // Insert new elements for the incoming bases
        // These elements can be tagged with an integer
        // by the calling function to indicate the sample
        // the new bases came from
        coordinate_map.insert(ii, alt.size(), variant_tag);
        
        /*
        std::cout << "Updated map: \n";
        for(size_t i = 0; i < coordinate_map.size(); ++i)
            std::cout << coordinate_map[i] << " ";
        std::cout << "\n";
        */
        assert(coordinate_map.size() == out.size());
    }
    return out;
}

// Extract somatic k-mers from the haplotype accordining to the tagged
// entries in the coordinate map
StringVector extractSomaticKmers(const std::string& haplotype,
                                 const std::vector<int> coordinate_map,
                                 int somatic_tag)
{
    StringVector kmers;

    assert(haplotype.size() == coordinate_map.size());

    // Get the index in the vector of the first somatic variant
    int first_idx = 0;
    while(coordinate_map[first_idx] != somatic_tag && 
          first_idx != (int)coordinate_map.size()) {
        first_idx++;
    }
    
    // no somatic variants?
    if(first_idx == (int)coordinate_map.size())
        return kmers;
    
    // Get the index of the last somatic variant
    int last_idx = first_idx;
    while(coordinate_map[last_idx] == somatic_tag)
        last_idx++;
    last_idx -= 1;

    // Start k - 1 bases before the first index
    first_idx -= (opt::k - 1);
    if(first_idx < 0)
        first_idx = 0;

    for(size_t i = first_idx; i <= (size_t)last_idx; ++i)
    {
        std::string t = haplotype.substr(i, opt::k);
        if(t.size() == (size_t)opt::k)
            kmers.push_back(t);
    }

    return kmers;
}

std::string reconstructHaplotype(const VCFRecord& somatic_record, 
                                 const VariantRecordVector nearby_germline, 
                                 const ReadTable& refTable,
                                 const BWTIndexSet& variantIndex,
                                 size_t flanking_size)
{
    // Extract the reference haplotype
    int eventLength = somatic_record.varStr.length();
    int zeroBasedPos = somatic_record.refPosition - 1;
    int start = zeroBasedPos - flanking_size - 1;
    int end = zeroBasedPos + eventLength + 2 * flanking_size;
    const SeqItem& chr = refTable.getRead(somatic_record.refName);
    if(end > chr.seq.length())
        return "";

    std::string reference_haplotype = chr.seq.substr(start, end - start);
    
    // Set up a vector mapping positions on the current haplotype
    // to original reference coordinates. This allows use to apply
    // the germline variants to the haplotype when the coordinates
    // shift due to indels
    std::vector<int> coordinate_map(reference_haplotype.size());
    for(int i = 0; i < (int)reference_haplotype.size(); ++i)
        coordinate_map[i] = i + start;

    // Do not attempt haplotype reconstruction if there are too many variants
    if(nearby_germline.size() > 4)
        return "";
    
    int SOMATIC_TAG = -2;
    int GERMLINE_TAG = -3;

    // Iterate over all possible combinations of germline variants
    size_t total_haplotypes = 1 << nearby_germline.size(); // 2 ^ n

    for(size_t haplotype_id = 0; haplotype_id < total_haplotypes; haplotype_id++)
    {
        std::string haplotype = reference_haplotype;
        std::vector<int> haplotype_coordinates = coordinate_map;

        // Iterate over the germline variants and check if the i-th bit
        // is set in the haplotype id. If so, apply the variant to the haplotype
        for(size_t var_id = 0; var_id < nearby_germline.size(); ++var_id)
        {
            if( (haplotype_id >> var_id) & 1)
            {
                std::string new_haplotype = 
                    applyVariant(haplotype, 
                                 nearby_germline[var_id].position - 1,
                                 nearby_germline[var_id].ref_sequence, 
                                 nearby_germline[var_id].alt_sequence,
                                 haplotype_coordinates, 
                                 true, 
                                 GERMLINE_TAG);
                haplotype = new_haplotype;

                if(haplotype.empty())
                    break;
            }
        }

        if(haplotype.empty())
            continue;

        // Apply the somatic variant to this haplotype
        std::string somatic_haplotype = 
            applyVariant(haplotype, 
                         somatic_record.refPosition - 1,
                         somatic_record.refStr, 
                         somatic_record.varStr,
                         haplotype_coordinates, 
                         true, 
                         SOMATIC_TAG);
    
        if(somatic_haplotype.empty())
            continue;

        StringVector kmers = extractSomaticKmers(somatic_haplotype, haplotype_coordinates, SOMATIC_TAG);
        
        // Count the number of somatic kmers that are present in the variant reads
        size_t read_kmers = 0;
        for(size_t ki = 0; ki < kmers.size(); ++ki)
        {
            if(BWTAlgorithms::countSequenceOccurrences(kmers[ki], variantIndex) > 0)
                read_kmers += 1;
        }

        // Check if all somatic kmers are present in the reads
        if(read_kmers == kmers.size())
            return somatic_haplotype;
    }

    return "";
}

//
// Main
//
int graphConcordanceMain(int argc, char** argv)
{
    parseGraphConcordanceOptions(argc, argv);

    Timer* pTimer = new Timer(PROGRAM_IDENT);

    // Load the reference
    ReadTable refTable(opt::referenceFile, SRF_NO_VALIDATION);
    refTable.indexReadsByID();

    // Index the germline variants so we can make germline haplotypes
    std::cerr << "Loading germline variant index..." << std::flush;
    VariantIndex germlineIndex(opt::germlineFile, refTable);
    std::cerr << "done\n";

    // Load FM-index of the reads
    BWTIndexSet variantIndex;
    std::cerr << "Loading variant read index... " << std::flush;
    std::string variantPrefix = stripGzippedExtension(opt::variantFile);
    variantIndex.pBWT = new BWT(variantPrefix + BWT_EXT, 256);
    variantIndex.pCache = new BWTIntervalCache(11, variantIndex.pBWT);
    std::cerr << "done\n";
    
    //
    BWTIndexSet baseIndex;
    std::cerr << "Loading base read index index... " << std::flush;
    std::string basePrefix = stripGzippedExtension(opt::baseFile);
    baseIndex.pBWT = new BWT(basePrefix + BWT_EXT, 256);
    baseIndex.pCache = new BWTIntervalCache(11, baseIndex.pBWT);
    std::cerr << "done\n";

    std::ifstream input(opt::vcfFile.c_str());
    std::string line;
    while(getline(input, line))
    {
        if(line.empty())
            continue;

        if(line[0] == '#')
        {
            std::cout << line << "\n";
            continue;
        }
        
        // parse record
        VCFRecord record(line);
        if(record.isMultiAllelic())
        {
            std::cerr << "Error: multi-allelic VCF found, please run vcfbreakmulti\n";
            exit(EXIT_FAILURE);
        }

        // Grab nearby variants
        size_t flanking_size = opt::k;
        VariantRecordVector nearby_vector = germlineIndex.getNearVariants(record.refName, 
                                                                          record.refPosition, 
                                                                          flanking_size);
        
        if(opt::verbose > 0)
        {
            std::cerr << "\n\n Variant: " << record << "\n";
            std::cerr << "Nearby: " << nearby_vector.size() << "\n";
        }

        std::string classification = "UNKNOWN";

        std::string haplotype;
        if(nearby_vector.size() <= 4)
            haplotype = reconstructHaplotype(record, nearby_vector, refTable, variantIndex, flanking_size);
        else
            classification = "TOO_COMPLEX";

        size_t max_unique_variant_kmers = 0;
        if(!haplotype.empty())
        {    
            IntVector v_cp = HapgenUtil::makeCountProfile(haplotype, opt::k, 9, variantIndex);
            IntVector b_cp = HapgenUtil::makeCountProfile(haplotype, opt::k, 9, baseIndex);
            
        
            size_t unique_variant_kmers = 0;
            for(size_t i = 0; i < v_cp.size(); ++i)
            {
                if(v_cp[i] > 0 && b_cp[i] == 0)
                    unique_variant_kmers++;
            }

            if(unique_variant_kmers > max_unique_variant_kmers)
                max_unique_variant_kmers = unique_variant_kmers;
            
            if(opt::verbose > 0)
            {
                std::cerr << "k-mer profile somatic (variant):  ";
                for(size_t j = 0; j < v_cp.size(); ++j)
                    std::cerr << v_cp[j];
                std::cerr << "\n";

                std::cerr << "k-mer profile somatic (base):     ";
                for(size_t j = 0; j < b_cp.size(); ++j)
                    std::cerr << b_cp[j];
                std::cerr << "\n";
                fprintf(stderr, "unique variant kmers: %zu\n\n", unique_variant_kmers);
                fprintf(stderr, "max unique variant kmers: %zu\n\n", max_unique_variant_kmers);
            }
        }
        else
        {
            classification = "CANNOT_RECONSTRUCT";
        }

        if(max_unique_variant_kmers > 5)
            classification = "SOMATIC";
        else
            classification = "GERMLINE";
        
        // write out the record
        record.addComment("MaxUniqueVariantKmers", (int)max_unique_variant_kmers);
        record.addComment("KmerClassification", classification);
        std::cout << record << "\n";
    }
    
    // Cleanup
    delete variantIndex.pBWT;
    delete variantIndex.pCache;
    delete baseIndex.pBWT;
    delete baseIndex.pCache;
    delete pTimer;

    return 0;
}

// 
// Handle command line arguments
//
void parseGraphConcordanceOptions(int argc, char** argv)
{
    std::string algo_str;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'b': arg >> opt::baseFile; break;
            case 'r': arg >> opt::variantFile; break;
            case 'g': arg >> opt::germlineFile; break;
            case OPT_REFERENCE: arg >> opt::referenceFile; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << GRAPH_CONCORDANCE_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << GRAPH_CONCORDANCE_VERSION_MESSAGE;
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

    if(opt::referenceFile.empty())
    {
        std::cerr << SUBPROGRAM ": a reference file must be given.\n";
        die = true;
    }

    if(opt::germlineFile.empty())
    {
        std::cerr << SUBPROGRAM ": a germline variants file must be given.\n";
        die = true;
    }

    opt::vcfFile = argv[optind++];

    if (die) 
    {
        std::cout << "\n" << GRAPH_CONCORDANCE_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}
