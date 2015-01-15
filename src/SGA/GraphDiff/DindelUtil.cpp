///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
//
// DindelUtil - Wrappers and utility functions
// for the dindel haplotype scoring functions
//
#include "DindelUtil.h"
#include "HapgenUtil.h"
#include "VCFUtil.h"
#include "Profiler.h"
#include "Verbosity.h"
#include "overlapper.h"
#include "multiple_alignment.h"

// Run dindel on a pair of samples
DindelReturnCode DindelUtil::runDindelPairMatePair(const std::string& id,
                                                   const StringVector& base_haplotypes,
                                                   const StringVector& variant_haplotypes,
                                                   const GraphCompareParameters& parameters,
                                                   std::ostream& baseOut,
                                                   std::ostream& variantOut,
                                                   std::ostream& callsOut,
                                                   DindelReadReferenceAlignmentVector* pReadAlignments)
{
    PROFILE_FUNC("runDindelPairMatePair")

    StringVector inHaplotypes;
    inHaplotypes.insert(inHaplotypes.end(), base_haplotypes.begin(), base_haplotypes.end());
    inHaplotypes.insert(inHaplotypes.end(), variant_haplotypes.begin(), variant_haplotypes.end());

    //
    // First, extract the reads from the normal and variant data sets that match each haplotype
    //
    assert(inHaplotypes.size() > 0);

    // Get canidate alignments for the input haplotypes
    HapgenAlignmentVector candidateAlignments;

    // Choose the kmer size for alignment
    size_t align_kmer = 31;
    for(size_t i = 0; i < inHaplotypes.size(); ++i)
    {
        HapgenAlignmentVector thisCandidateAlignments;
        HapgenUtil::alignHaplotypeToReferenceKmer(align_kmer,
                                                  inHaplotypes[i],
                                                  parameters.referenceIndex,
                                                  parameters.pRefTable,
                                                  thisCandidateAlignments);

        candidateAlignments.insert(candidateAlignments.end(), thisCandidateAlignments.begin(), thisCandidateAlignments.end());
    }
   
    // Remove duplicate or bad alignment pairs
    HapgenUtil::coalesceAlignments(candidateAlignments);

    if(Verbosity::Instance().getPrintLevel() > 3)
        printf("runDindel -- %zu candidate alignments found\n", candidateAlignments.size());
    
    // Check if this group of haplotypes is a SNP-only event
    bool has_complex_variant = false;
    for(size_t i = 0; i < candidateAlignments.size(); ++i)
    {
        if(candidateAlignments[i].numEvents > 0 && !candidateAlignments[i].isSNP)
            has_complex_variant = true;
    }

    size_t MAX_ALIGNMENTS = 10;
    if(candidateAlignments.size() > MAX_ALIGNMENTS)
        return DRC_AMBIGUOUS_ALIGNMENT;

    // Join each haplotype with flanking sequence from the reference genome for each alignment
    // This function also adds a haplotype (with flanking sequence) for the piece of the reference
    int FLANKING_SIZE = 0;
    if (parameters.dindelRealignParameters.realignMatePairs)
        FLANKING_SIZE = 1000;
    StringVector flankingHaplotypes;

    // This vector contains the internal portion of the haplotypes, without the flanking sequence
    // It is used to extract reads
    StringVector candidateHaplotypes;
    for(size_t i = 0; i < candidateAlignments.size(); ++i)
    {
        HapgenUtil::makeFlankingHaplotypes(candidateAlignments[i],
                                           parameters.pRefTable,
                                           FLANKING_SIZE,
                                           inHaplotypes,
                                           flankingHaplotypes,
                                           candidateHaplotypes);
    
    }

    if(Verbosity::Instance().getPrintLevel() > 3)
        printf("runDindel -- made %zu flanking haplotypes\n", candidateHaplotypes.size());

    // Normal reads
    SeqRecordVector normalReads;
    SeqRecordVector normalRCReads;

    // Remove non-unique candidate haplotypes
    std::sort(candidateHaplotypes.begin(), candidateHaplotypes.end());
    StringVector::iterator haplotype_iterator = std::unique(candidateHaplotypes.begin(), candidateHaplotypes.end());
    candidateHaplotypes.resize(haplotype_iterator - candidateHaplotypes.begin());

    // Set the value to use for extracting reads that potentially match the haplotype
    // Do not use a kmer for extraction greater than this value
    size_t KMER_CEILING = 31;
    size_t extractionKmer = parameters.kmer < KMER_CEILING ? parameters.kmer : KMER_CEILING;
    
    bool extractOK = true;
    if(!parameters.bReferenceMode)
    {
        // Reads on the same strand as the haplotype
        extractOK = HapgenUtil::extractHaplotypeReads(candidateHaplotypes, parameters.baseIndex, extractionKmer, 
                                                      false, parameters.maxReads, parameters.maxExtractionIntervalSize, &normalReads, NULL);

        if(!extractOK)
            return DRC_OVER_DEPTH;

        // Reads on the reverse strand
        extractOK = HapgenUtil::extractHaplotypeReads(candidateHaplotypes, parameters.baseIndex, extractionKmer, 
                                                      true, parameters.maxReads, parameters.maxExtractionIntervalSize, &normalRCReads, NULL);

        if(!extractOK)
            return DRC_OVER_DEPTH;
    }

    // Variant reads
    SeqRecordVector variantReads;
    SeqRecordVector variantRCReads;

    extractOK = HapgenUtil::extractHaplotypeReads(candidateHaplotypes, parameters.variantIndex, extractionKmer, 
                                                  false, parameters.maxReads, parameters.maxExtractionIntervalSize, &variantReads, NULL);

    if(!extractOK)
        return DRC_OVER_DEPTH;

    extractOK = HapgenUtil::extractHaplotypeReads(candidateHaplotypes, parameters.variantIndex, extractionKmer, 
                                                  true, parameters.maxReads, parameters.maxExtractionIntervalSize, &variantRCReads, NULL);

    if(!extractOK)
        return DRC_OVER_DEPTH;

    size_t normal_reads = normalReads.size() + normalRCReads.size();
    size_t variant_reads = variantReads.size() + variantRCReads.size();
    size_t total_reads = normal_reads + variant_reads;
    
    if(Verbosity::Instance().getPrintLevel() > 3)
        printf("Extracted %zu normal reads, %zu variant reads\n", normal_reads, variant_reads);

    if(total_reads > parameters.maxReads)
        return DRC_OVER_DEPTH;

    if (total_reads == 0)
        return DRC_UNDER_DEPTH;

    // Generate the input haplotypes for dindel
    // We need at least 2 haplotypes (one is the reference)
    size_t totFlankingHaplotypes = flankingHaplotypes.size();

    if(totFlankingHaplotypes < 2)
        return DRC_NO_ALIGNMENT;

    // Ensure the reference haplotype is a non-empty string
    if(flankingHaplotypes[0].size() == 0)
        return DRC_NO_ALIGNMENT;

    // Make Dindel referenceMappings
    StringVector dindelHaplotypes;
    std::set<DindelReferenceMapping> refMappings;

    //
    for(size_t i = 0; i < candidateAlignments.size(); ++i)
    {
        std::string upstream, defined, downstream;
        std::string refName = parameters.pRefTable->getRead(candidateAlignments[i].referenceID).id;

        HapgenUtil::extractReferenceSubstrings(candidateAlignments[i],parameters.pRefTable, 
                                               FLANKING_SIZE, upstream, defined, downstream);

        std::string refSeq = upstream + defined + downstream;
     
        int refStart = candidateAlignments[i].position - int(upstream.size()) + 1;

        // Here the score is used as an estimate of how unique "defined" is in the reference sequence.
        // "defined" is not the reference sequence but a candidate haplotype.
        // It is conservative because the flanking sequence is not used in this estimation.
    	DindelReferenceMapping rm(refName, 
                                  refSeq, 
                                  refStart, 
                                  candidateAlignments[i].score+2*FLANKING_SIZE, 
                                  candidateAlignments[i].isRC);

        std::set<DindelReferenceMapping>::iterator rmit = refMappings.find(rm);
        if(rmit == refMappings.end())
        {
            refMappings.insert(rm);
        }
	    else
        {
            if(rm.referenceAlignmentScore > rmit->referenceAlignmentScore) 
                rmit->referenceAlignmentScore = rm.referenceAlignmentScore;
        }
    }
    
    // RESET MAPPING SCORES
    for(std::set<DindelReferenceMapping>::iterator it = refMappings.begin(); it != refMappings.end(); it++) 
        it->referenceAlignmentScore = 1000;
    
    // make flankingHaplotypes unique
    std::set< std::string > setFlanking(flankingHaplotypes.begin(), flankingHaplotypes.end());

    for(std::set< std::string >::const_iterator it = setFlanking.begin(); it != setFlanking.end(); it++)
    {
        dindelHaplotypes.push_back(*it);
        //dindelRefMappings[i] = std::vector<DindelReferenceMapping>(refMappings.begin(),refMappings.end());
    }

    std::vector<DindelReferenceMapping> dRefMappings(refMappings.begin(),refMappings.end());
    DindelWindow dWindow(dindelHaplotypes, dRefMappings);

    //
    // Run Dindel
    //
    
    // Initialize VCF collections
    VCFCollection vcfCollections[2];

    // If in multisample mode, load the sample names into the VCFCollection
    if(parameters.variantIndex.pPopIdx != NULL)
    {
        for(size_t i = 0; i <= 1; ++i)
            vcfCollections[i].samples = parameters.variantIndex.pPopIdx->getSamples();
    }

    size_t start_i = parameters.bReferenceMode ? 1 : 0;

    DindelRealignWindowResult *pThisResult = NULL;
    DindelRealignWindowResult *pPreviousResult = NULL;

    for(size_t i = start_i; i <= 1; ++i)
    {
        SeqRecordVector& fwdReads = (i == 0) ? normalReads : variantReads;
        SeqRecordVector& rcReads = (i == 0) ? normalRCReads : variantRCReads;
        const BWTIndexSet* indices = &parameters.variantIndex;

        // Create dindel reads
        // Mates must be at the end of the array.
        std::vector<DindelRead> dReads;
        for(size_t j = 0; j < fwdReads.size(); ++j)
            dReads.push_back(convertToDindelRead(indices, fwdReads[j], true));

        for(size_t j = 0; j < rcReads.size(); ++j)
        {
            rcReads[j].seq.reverseComplement();
            std::reverse(rcReads[j].qual.begin(), rcReads[j].qual.end());
            dReads.push_back(convertToDindelRead(indices, rcReads[j], false));
        }

        pThisResult = new DindelRealignWindowResult();

        std::stringstream out_ss;
        try
        {
            DindelRealignWindow dRealignWindow(&dWindow, dReads, parameters.dindelRealignParameters);
            dRealignWindow.run("hmm", vcfCollections[i], pReadAlignments, id, pThisResult, pPreviousResult, parameters.pRefTable);
        }
        catch(std::string e)
        {
            std::cerr << "Dindel Exception: " << e << "\n";
            exit(DRC_EXCEPTION);
        }


        if(i == 0)
            pPreviousResult = pThisResult;
    }

    // Copy raw VCFRecords to output
    for(size_t i = 0; i <= 1; ++i)
    {
        std::ostream& curr_out = i == 0 ? baseOut : variantOut;
        for(size_t j = 0; j < vcfCollections[i].records.size(); ++j)
            curr_out << vcfCollections[i].records[j] << "\n";
    }

    // Make comparative calls
    size_t VARIANT_IDX = 1;
    size_t BASE_IDX = 0;
    bool has_base_calls = !vcfCollections[BASE_IDX].records.empty();
    for(size_t i = 0; i < vcfCollections[1].records.size(); ++i)
    {
        bool not_called_in_base = true;
        if(has_base_calls)
            not_called_in_base = vcfCollections[BASE_IDX].records[i].passStr == "NoCall" ||
                                 vcfCollections[BASE_IDX].records[i].passStr == "NoSupp";

        bool called_in_variant = vcfCollections[VARIANT_IDX].records[i].passStr == "PASS";
        if(called_in_variant && not_called_in_base)
            callsOut << vcfCollections[VARIANT_IDX].records[i] << "\n";
    }

    baseOut.flush();
    variantOut.flush();

    delete pThisResult;
    delete pPreviousResult;

    return DRC_OK;
}

DindelRead DindelUtil::convertToDindelRead(const BWTIndexSet* indices, const SeqRecord& record, bool is_forward)
{
    double MAP_QUAL = 40.0;

    std::string sample = "SAMPLE";
    if(indices->pPopIdx != NULL)
    {
        // Parse the read index from the name of the SeqItem
        if(record.id.find("idx-") == std::string::npos)
        {
            std::cerr << "Unexpected sequence id: " << record.id << "\n";
            assert(false);
        }
        std::stringstream parser(record.id.substr(4));
        size_t index;
        parser >> index;

        sample = indices->pPopIdx->getName(index);
    }

    return DindelRead(record, sample, MAP_QUAL, is_forward);
}

//
void DindelUtil::doMultipleReadHaplotypeAlignment(const std::vector<DindelRead> & dReads,
                                                  const StringVector & haplotypes)
{

    // globally align haplotypes to the first haplotype (arbitrary)
    assert(haplotypes.size()>0);
    for (size_t h = 0; h < haplotypes.size(); ++h)
    {
        std::cout << "ALIGNING EVERYTHING AGAINST HAPLOTYPE " << h << "\n";
        MultipleAlignment ma;
        const std::string rootSequence = haplotypes[h];
        ma.addBaseSequence("root", rootSequence, "");

        std::string hid;
        for(size_t j = 0; j < haplotypes.size(); j++)
        {
            std::stringstream ss;
            if (j!=h)
                ss << "haplotype-" << j;
            else
                ss << "HAPLOTYPE-" << j;
            SequenceOverlap o = Overlapper::computeOverlap(rootSequence, haplotypes[j]);
            ma.addOverlap(ss.str(), haplotypes[j], "", o);
        }
    

        for(size_t r = 0; r < dReads.size(); ++r)
        {
            std::stringstream ss;
            if (r<dReads.size()/2) 
                ss << "read-" << r << "("  << dReads[r].getID() << ")"; 
            else 
                ss << "MATE read-" << r;
            
            SequenceOverlap o = Overlapper::computeOverlap(rootSequence, dReads[r].getSequence());
            ma.addOverlap(ss.str(), dReads[r].getSequence(), "", o);
        }
        ma.print(100000);
    }
}

// Initialize an array of code counts
void DindelUtil::initializeCodeCounts(int counts[DRC_NUM_CODES])
{
    for(int i = 0; i < DRC_NUM_CODES; ++i)
        counts[i] = 0;
}

// Print a report of the dindel return codes
void DindelUtil::printReturnReport(const int counts[DRC_NUM_CODES])
{
    // Sum the total number of haplotype sets processed
    int sum = 0;
    for(int i = 0; i < DRC_NUM_CODES; ++i)
        sum += counts[i];

    printf("Total variants processed: %d\n", sum);
    printf("    number failed due to depth check: %d\n", counts[DRC_OVER_DEPTH]);
    printf("    number failed due to no alignment: %d\n", counts[DRC_NO_ALIGNMENT]);
    printf("    number failed due to poor alignment: %d\n", counts[DRC_POOR_ALIGNMENT]);
    printf("    number failed due to ambiguous alignment: %d\n", counts[DRC_AMBIGUOUS_ALIGNMENT]);
    printf("    number failed due to dindel exception: %d\n", counts[DRC_EXCEPTION]);
    printf("    number passed to dindel: %d\n", counts[DRC_OK]);
}

