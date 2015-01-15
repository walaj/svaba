//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// VCFUtil - collection of utility functions
// for producing VCF records.
//
#include <iostream>
#include <stdlib.h>
#include <iterator>
#include <sstream>
#include "VCFUtil.h"
#include "StdAlnTools.h"
#include "MultiAlignment.h"

// 
// VCFRecord
//
VCFRecord::VCFRecord()
{
    passStr = "PASS";
    quality = 20.0f;
}

// parse from a line
VCFRecord::VCFRecord(const std::string& line)
{
    std::stringstream parser(line);

    parser >> refName;
    parser >> refPosition;
    parser >> id;
    parser >> refStr;
    parser >> varStr;

    // Non-conformant VCF might represent unknown qualities with '.'
    std::string quality_str;
    parser >> quality_str;
    if(quality_str == ".")
    {
        quality = 255;
    }
    else
    {
        std::stringstream qparser(quality_str);
        qparser >> quality;
    }

    parser >> passStr;
    
    std::string tmp;
    parser >> tmp;

    // Split the comments on the semicolon delimiter
    comments = split(tmp, ';');
    
    parser >> formatStr;
    if(!formatStr.empty())
    {
        // parse genotypes
        while(parser >> tmp)
            sampleStr.push_back(tmp);
    }
}

void VCFRecord::addComment(const std::string& key, const std::string& value)
{
    std::string out = key;
    out.append(1, '=');
    out.append(value);
    comments.push_back(out);
}

void VCFRecord::addComment(const std::string& key, const int& value)
{
    std::stringstream v_str;
    v_str << value;
    addComment(key, v_str.str());
}

void VCFRecord::addComment(const std::string& key, const double& value)
{
    std::stringstream v_str;
    v_str.precision(5);
    v_str.setf(std::ios::fixed,std::ios::floatfield);
    v_str << value;
    addComment(key, v_str.str());
}


void VCFRecord::setPassStr(const std::string str)
{
    passStr = str;
}

void VCFRecord::setQuality(const double q)
{
    quality = q;
}

void VCFRecord::printVerbose() const
{
    std::cout << "RefName: " << refName << "\n";
    std::cout << "Position: " << refPosition << "\n";
    std::cout << "RefStr: " << refStr << "\n";
    std::cout << "VarStr: " << varStr << "\n";
}

// Write the VCF record to the stream
std::ostream& operator<<(std::ostream& o, const VCFRecord& record)
{
    o.precision(5);
    o.setf(std::ios::fixed,std::ios::floatfield);

    o << record.refName << "\t";
    o << record.refPosition << "\t",
    o << (record.id.empty() ? "." : record.id) << "\t";
    o << record.refStr << "\t";
    o << record.varStr << "\t";
    o << (int)record.quality << "\t";
    o << record.passStr << "\t";

    for(size_t i = 0; i < record.comments.size(); ++i)
    {
        if(i > 0)
            o << ";";
        o << record.comments[i];
    }

    if(record.sampleStr.size())
    {
        o << "\t" << record.formatStr;
        for(size_t i = 0; i < record.sampleStr.size(); ++i)
            o << "\t" << record.sampleStr[i];
    }

    return o;
}


//
VCFClassification VCFRecord::classify() const
{
    if(refStr.size() == varStr.size())
        return VCF_SUB;
    
    if(refStr.size() > varStr.size())
    {
        // If a deletion, varStr should be a prefix of refStr
        if(refStr.find(varStr) == 0)
            return VCF_DEL;
        else
            return VCF_COMPLEX;
    }

    if(refStr.size() < varStr.size())
    {
        // If an insertion, refStr should be a prefix of refStr
        if(varStr.find(refStr) == 0)
            return VCF_INS;
        else
            return VCF_COMPLEX;
    }
    return VCF_COMPLEX;
}

//
// VCFUtil
//

// See header for description
VCFReturnCode VCFUtil::generateVCFFromCancerVariant(const std::string& ref,
                                                    const std::string& base,
                                                    const std::string& variant,
                                                    const std::string& refName,
                                                    size_t refPosition,
                                                    const std::string& varName,
                                                    int minExactMatch,
                                                    double dustThreshold,
                                                    int verbose,
                                                    VCFVector& outRecords)
{
    VCFVector tempRecords;

    // Construct a multiple alignment from pairwise alignments
    // of the base and variant strings to the reference

    // We start by constructing a global alignment between the sequences
    // and representing the alignment as a Cigar string
    std::string cigarRefBase = StdAlnTools::globalAlignmentCigar(base, ref);
    std::string cigarRefVariant = StdAlnTools::globalAlignmentCigar(variant, ref);
    
    // Initialize the multiple alignment data
    
    MAlignData baseMA;
    baseMA.position = 0;
    baseMA.str = base;
    baseMA.expandedCigar = StdAlnTools::expandCigar(cigarRefBase);
    baseMA.name = "base";

    MAlignData varMA;
    varMA.position = 0;
    varMA.str = variant;
    varMA.expandedCigar = StdAlnTools::expandCigar(cigarRefVariant);
    varMA.name = "variant";

    MAlignDataVector maVector;
    maVector.push_back(baseMA);
    maVector.push_back(varMA);
    
    MultiAlignment ma(ref, maVector);

    /*
    std::cout << "\nMultipleAlignment:\n";
    ma.print();
    */

    // Get the row index of the variant and base strings
    size_t baseRow = ma.getIdxByName("base");
    size_t varRow = ma.getIdxByName("variant");
    size_t refRow = ma.getRootIdx();

    size_t numCols = ma.getNumColumns();

    // Iterate over every column of the multiple alignment and detect variants
    int leftExactMatch = 0;
    int rightExactMatch = 0;

    int eventStart = -1;
    int eventEnd = -1;
    bool isIndel = false;

    bool inLeftExact = true;

    for(size_t i = 0; i < numCols; ++i)
    {
        char baseSymbol = ma.getSymbol(baseRow, i);
        char varSymbol = ma.getSymbol(varRow, i);
        bool isVariant = (varSymbol != baseSymbol);

        // Update the counter of the number of exact matches for the leftmost and rightmost bases
        if(isVariant)
        {
            inLeftExact = false; // stop the count of leftmost exact matches
            rightExactMatch  = 0; // reset the rightmost count
        }
        else // this is a good match
        {
            rightExactMatch += 1;
            if(inLeftExact)
                leftExactMatch += 1;
        }
        
        // Process variants indicated by this portion of the alignment
        if(isVariant)
        {
            if(eventStart == -1)
            {
                // Not currently in a variant event, start a new one
                eventStart = i;
                eventEnd = i;
            }
            else
            {
                // Extend the event
                assert(eventEnd != -1);
                eventEnd += 1;
            }

            isIndel = isIndel || varSymbol == '-' || baseSymbol == '-';
        }
        else
        {
            // Check if this is the end of a variant
            if(eventStart != -1)
            {
                if(verbose > 0)
                {
                    std::cout << "\n\nProcessing variant\n";
                    ma.print();
                }

                // Extract the substrings of the three sequences describing
                // this event.
                // If the event is an indel, we move the start point
                // of the event backwards until all 3 columns are not padding characters.
                // This is to avoid the following situation:
                // AATATTT--CGT ref
                // AATATTTT-CGT base
                // AATATTTTTCGT var
                // or 
                // TTAGGTTTTTTT ref 
                // TTAGG-TTTTTT base 
                // TTAGG--TTTTT var
                //
                // In the first case, it is a 1 base insertion wrt base but a 2 base insertion wrt to the ref.
                // So we need to move the event to the first column with all Ts to properly capture what is happening.
                // In the second case it is a 1 base deletion wrt base and a 2 base deletion wrt to the ref. 
                // Again, we need to move back to the last full column.
                if(isIndel)
                {
                    // It is impossible to extract a reference string if the indel is at the beginning
                    // or end of th multiple alignment so we just fail out here
                    if(eventStart == 0 || eventEnd == (int)numCols - 1)
                        return VCF_INVALID_MULTIALIGNMENT;
                    
                    while(eventStart >= 0)
                    {
                        char rc = ma.getSymbol(refRow, eventStart);
                        char bc = ma.getSymbol(baseRow, eventStart);
                        char vc = ma.getSymbol(varRow, eventStart);
                        if(rc == '-' || bc == '-' || vc == '-')
                            eventStart -= 1;
                        else
                            break;
                    }
                    assert(eventStart >= 0);
                }

                size_t eventLength = eventEnd - eventStart + 1;
                std::string refString = StdAlnTools::unpad(ma.getPaddedSubstr(refRow, eventStart, eventLength));
                std::string baseString = StdAlnTools::unpad(ma.getPaddedSubstr(baseRow, eventStart, eventLength));
                std::string varString = StdAlnTools::unpad(ma.getPaddedSubstr(varRow, eventStart, eventLength));

                if(verbose > 0)
                {
                    printf("Ref start: %d col: %d eventStart: %d eventEnd: %d\n", (int)refPosition, (int)i, eventStart, eventEnd);
                    std::cout << "RefString " << refString << "\n";
                    std::cout << "baseString " << baseString << "\n";
                    std::cout << "varString " << varString << "\n";
                }
                assert(!refString.empty());
                //assert(!baseString.empty());
                assert(!varString.empty());
                if(baseString.empty())
                    std::cerr << "Warning: Empty base string\n";

                // Get the base position in the reference string. This is not necessarily the 
                // same as the eventStart column as the reference may be padded
                int refBaseOffset = ma.getBaseIdx(refRow, eventStart);
                
                // Generate VCF Record
                VCFRecord record;
                record.refName = refName;
                record.refPosition = refPosition + refBaseOffset + 1; // VCF uses 1-based coordinates

                record.refStr = refString;
                record.varStr = varString;

                // Add a comment with the variant name
                record.addComment("id", varName);

                // If the base string is not the same as the reference string, not it in the comment
                // This can happen when a het SNP is mutated or we have some other weird site
                if(refString != baseString)
                    record.addComment("NS", baseString);

                // Check if the reference context is a homopolymer
                int preHP = ma.countHomopolymer(refRow, eventStart, 0);
                int sufHP = ma.countHomopolymer(refRow, eventEnd);
                int maxHP = std::max(preHP, sufHP);
                if(maxHP >= HOMOPOLYMER_ANNOTATE_THRESHOLD)
                {
                    std::stringstream hps;
                    hps << maxHP;
                    record.addComment("HP", hps.str());
                }

                // Calculate the highest dust score on the reference sequence
                // using 64-base windows. 
                double globalDustScore = maxDustWindow(ref);
                if(globalDustScore > dustThreshold)
                {
                    std::stringstream dts;
                    dts.precision(2);
                    dts << globalDustScore;
                    record.addComment("DT", dts.str());
                }

                if(verbose > 0)
                    std::cout << record << "\n";

                tempRecords.push_back(record);

                // Reset state
                eventStart = -1;
                eventEnd = -1;
                isIndel = false;
            }
        }
    }

    // Perform the start/end exact match sanity check
    if(leftExactMatch >= minExactMatch && rightExactMatch >= minExactMatch)
    {
        outRecords.insert(outRecords.end(), tempRecords.begin(), tempRecords.end());
        return VCF_OK;
    }
    else
    {
        return VCF_EXACT_MATCH_FAILED;
    }
}

// Write a VCF header to the given file handle
void VCFUtil::writeHeader(std::ostream* pWriter, const std::string& program, const std::string& inputBAM, const std::string& reference)
{
    // Get and format the current date
    time_t rawtime;
    struct tm * timeinfo;
    static int MAX_DATE_CHARS = 80;
    char dateBuffer [MAX_DATE_CHARS];
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    int copied = strftime(dateBuffer,MAX_DATE_CHARS,"%Y%m%d",timeinfo);
    assert(copied < MAX_DATE_CHARS);

    *pWriter << "##fileformat=VCFv4.0\n";
    *pWriter << "##fileDate=" << dateBuffer << "\n";
    *pWriter << "##command=" << program << "\n";
    *pWriter << "##inputBAM=" <<  inputBAM << "\n";
    *pWriter << "##reference=" << reference << "\n";
    *pWriter << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
}
