//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// VCFUtil - collection of utility functions
// for producing VCF records.
//
#ifndef VCFUTIL_H
#define VCFUTIL_H

#include <string>
#include <vector>

// Enums
enum VCFClassification
{
    VCF_SUB,
    VCF_DEL,
    VCF_INS,
    VCF_COMPLEX,
    VCF_NUM_CLASSIFICATIONS
};

// 
enum VCFReturnCode
{
    VCF_EXACT_MATCH_FAILED, // failed the flanking-base match test
    VCF_MAP_QUALITY_FAILED, // mapping quality check failed
    VCF_BASE_MULTIMAP_FAILED, // mapping quality check failed
    VCF_BASE_PARTIALMAP_FAILED, // the base sequence did not fully align to the reference
    VCF_INVALID_MULTIALIGNMENT, // failed because the multiple alignment was not valid
    VCF_OK, // VCF generation succeeded
    VCF_NUM_RETURN_CODES
};

// Structs 
struct VCFRecord
{
    // functions

    // initialize an empty record
    VCFRecord();

    // parse a record from a tab-delimited string
    VCFRecord(const std::string& line);

    void addComment(const std::string& key, const std::string& value);
    void addComment(const std::string& key, const int& value);
    void addComment(const std::string& key, const double& value);

    void setPassStr(const std::string str);
    void setQuality(const double q);

    void printVerbose() const;
    VCFClassification classify() const;
    friend std::ostream& operator<<(std::ostream& o, const VCFRecord& record);

    static bool sort(const VCFRecord& a, const VCFRecord& b) 
    { 
        if(a.refName != b.refName)
            return a.refName < b.refName;
        else
            return a.refPosition < b.refPosition;
    }

    bool isMultiAllelic() const { return varStr.find(",") != std::string::npos; }

    // data
    std::string refName;
    size_t refPosition;
    std::string id;
    std::string refStr;
    std::string varStr;
    double quality;
    std::string passStr;
    std::string formatStr;
    std::vector<std::string> comments;
    std::vector<std::string> sampleStr;
};

// Typedefs
typedef std::vector<VCFRecord> VCFVector;
typedef std::vector<std::string> StringVector;

// A set of VCF records, with sample names in a predefined order
struct VCFCollection
{
    StringVector samples;
    VCFVector records;
};

namespace VCFUtil
{
    // Generate a VCF record from a cancer mutation
    // Base is the sequence in the individuals normal genome
    // Variant is the mutated copy of base
    // Ref is the sequence of the reference that base represents.
    // Ref will typically be very close to base, but may have some
    // differences due to SNPs or indels, which we do not want
    // to call as variants.
    // The minExactMatch parameter constrains base and variant to have
    // at least a perfect alignment of minExactMatch at the start 
    // and end of the sequences. This is used as a sanity check for the alignments
    // If the variants were produced by graph-diff, base and variant
    // should share exactly k bases at the ends.
    VCFReturnCode generateVCFFromCancerVariant(const std::string& ref,
                                               const std::string& base,
                                               const std::string& variant,
                                               const std::string& refName,
                                               size_t refPosition,
                                               const std::string& varName,
                                               int minExactMatch,
                                               double dustThreshold,
                                               int verbose,
                                               VCFVector& outRecords);

    // Write a VCF header to the given file handle
    void writeHeader(std::ostream* pWriter, const std::string& program, 
                     const std::string& bamFile, const std::string& reference);

    // constants
    static const int HOMOPOLYMER_ANNOTATE_THRESHOLD = 5;

};

#endif
