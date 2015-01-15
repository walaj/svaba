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
#ifndef DINDEL_UTIL_H
#define DINDEL_UTIL_H

#include "GraphCompare.h"
#include "DindelRealignWindow.h"
#include "HapgenUtil.h"


enum DindelReturnCode
{
    DRC_OK,
    DRC_NO_ALIGNMENT,
    DRC_POOR_ALIGNMENT,
    DRC_AMBIGUOUS_ALIGNMENT,
    DRC_OVER_DEPTH,
    DRC_UNDER_DEPTH,
    DRC_EXCEPTION,
    DRC_NUM_CODES
};

namespace DindelUtil
{
    // Run dindel on a pair of samples
    DindelReturnCode runDindelPairMatePair(const std::string& id,
                                           const StringVector& base_haplotypes,
                                           const StringVector& variant_haplotypes,
                                           const GraphCompareParameters& parameters,
                                           std::ostream& baseOut,
                                           std::ostream& variantOut,
                                           std::ostream& callsOut,
                                           DindelReadReferenceAlignmentVector* pReadAlignments);

    // Run a naive caller
    DindelReturnCode runNaiveCaller(const std::string& normalString,
                                    const std::string& variantString,
                                    const GraphCompareParameters& parameters,
                                    std::ostream& baseOut,
                                    std::ostream& variantOut);

    // show multiple alignment of all reads and all haplotypes
    void doMultipleReadHaplotypeAlignment(const std::vector<DindelRead> & dReads,
                                          const StringVector & haplotypes);

    // Convert a read-haplotype alignment to a reference alignment in BAM format
    void readHaplotypeAlignmentToBAM(const std::string& readID,
                                     const std::string& readSequence,
                                     const std::string& readHaplotypeCigar,
                                     const int haplotype_reference_position,
                                     const std::string& haplotypeReferenceCigar);


    // Convert a SeqItem to a dindel read
    DindelRead convertToDindelRead(const BWTIndexSet* indices, const SeqRecord& record, bool is_forward);

    void initializeCodeCounts(int counts[DRC_NUM_CODES]);
    void printReturnReport(const int counts[DRC_NUM_CODES]);
};

#endif
