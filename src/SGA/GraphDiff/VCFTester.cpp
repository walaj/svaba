///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// VCFTester - Run dindel on each variant in a vcf file
//
#include "VCFTester.h"
#include "StdAlnTools.h"
#include "DindelUtil.h"
#include "BWTAlgorithms.h"
#include "HaplotypeBuilder.h"
#include <sstream>

//
//
//
VCFTester::VCFTester(const GraphCompareParameters& params) : m_parameters(params),
                                                             m_graphComparer(params),
                                                             m_baseVCFFile("vtest.base.vcf","w"),
                                                             m_variantVCFFile("vtest.variant.vcf","w"),
                                                             m_callsVCFFile("vtest.calls.vcf","w")
{
    m_baseVCFFile.outputHeader("stub", "stub");
    m_variantVCFFile.outputHeader("stub", "stub");
    m_callsVCFFile.outputHeader("stub", "stub");
    DindelUtil::initializeCodeCounts(m_returnCodes);
}

//
VCFTester::~VCFTester()
{
    std::cout << "Done testing variants\n";
    DindelUtil::printReturnReport(m_returnCodes);
}

void VCFTester::process(const VCFFile::VCFEntry& record)
{
    // Extract the reference string for this record
    // We want a full kmer to the left and right of the variant
    int eventLength = record.alt.length();

    int zeroBasedPos = record.pos - 1;
    int start = zeroBasedPos - m_parameters.kmer - 1;
    int end = zeroBasedPos + eventLength + m_parameters.kmer;

    const SeqItem& chr = m_parameters.pRefTable->getRead(record.chrom);

    std::string refStr = chr.seq.substr(start, end - start);

    // Apply the variant to the reference string
    // Update the coordinate to be wrt the start of the reference substring
    int relative_pos = zeroBasedPos - start;
    std::string varStr = applyVariant(refStr, relative_pos, record.ref, record.alt);

    std::cout << "\n**************************************\n";
    std::cout << "Debugging variant (" << record.chrom << "-" << record.pos << " " << record.ref << "/" << record.alt << ")\n";

    std::cout << "Running dindel on variant";
    record.write(std::cout);

    StdAlnTools::globalAlignment(refStr, varStr, true);

    // Run dindel
    std::stringstream baseSS;
    std::stringstream variantSS;
    std::stringstream callsSS;
    std::string id = ".";

    StringVector ref_haplotypes;
    ref_haplotypes.push_back(refStr);

    StringVector var_haplotypes;
    var_haplotypes.push_back(varStr);
    DindelReturnCode code = DindelUtil::runDindelPairMatePair(id,
                                                              ref_haplotypes,
                                                              var_haplotypes,
                                                              m_parameters,
                                                              baseSS,
                                                              variantSS,
                                                              callsSS,
                                                              NULL);

    m_returnCodes[code] += 1;

    // Try to find out why we didn't find it in the graph
    std::vector<int> bv;
    std::vector<int> vv;

    size_t k = m_parameters.kmer;
    std::string variantKmer;
    for(size_t i = 0; i < varStr.size() - k + 1; ++i)
    {
        std::string ks = varStr.substr(i, k);

        size_t baseCount = BWTAlgorithms::countSequenceOccurrences(ks, m_parameters.baseIndex);
        size_t varCount = BWTAlgorithms::countSequenceOccurrences(ks, m_parameters.variantIndex);

        bv.push_back(baseCount > 9 ? 9 : baseCount);
        vv.push_back(varCount > 9 ? 9 : varCount);

        if(baseCount == 0 && varCount > 0)
            variantKmer = ks;
    }



    /*
    std::cout << "base: " << baseSS.str() << "\n";
    std::cout << "variant: " << variantSS.str() << "\n";
    */
    std::cout<< " B:\t";
    
    for(size_t i = 0; i < bv.size(); ++i)
        std::cout << bv[i];
    std::cout << "\n";
    
    std::cout<< " V:\t";
    for(size_t i = 0; i < vv.size(); ++i)
        std::cout << vv[i];
    std::cout << "\n";

    bool canAssemble = !variantKmer.empty();
    std::cout << "Can assemble? " << canAssemble << "\n";
    
    //
    // Attempt HaplotypeBuilder assembly
    //
    AnchorSequence startAnchor;
    startAnchor.sequence = varStr.substr(0,k);
    startAnchor.count = 0;
    startAnchor.position = 0;

    AnchorSequence endAnchor;
    endAnchor.sequence = varStr.substr(varStr.size() - k, k);
    endAnchor.count = 0;
    endAnchor.position = 0;

    // Attempt the actual GraphCompare assembly process
    if(canAssemble)
        m_graphComparer.testKmer(variantKmer);
}

std::string VCFTester::applyVariant(const std::string& in, int pos,
                                    const std::string& ref, const std::string& alt)
{
    std::string out = in;
    assert(pos > 0 && pos < (int)out.length());

    // Ensure that the reference string at the variant matches the expected
    assert(out.substr(pos, ref.length()) == ref);
    out.replace(pos, ref.length(), alt);
    return out;
}
