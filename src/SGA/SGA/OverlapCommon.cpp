//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapCommon - Common wrapper used for finding overlaps
// for a set of reads
//
#include "OverlapCommon.h"

// Convert a line from a hits file into a vector of overlaps and sets the flag
// indicating whether the read was found to be a substring of other reads
// Only the forward read table is used since we only care about the IDs and length
// of the read, not the sequence, so that we don't need an explicit reverse read table
void OverlapCommon::parseHitsString(const std::string& hitString, 
                                    const ReadInfoTable* pQueryRIT, 
                                    const ReadInfoTable* pTargetRIT, 
                                    const SuffixArray* pFwdSAI, 
                                    const SuffixArray* pRevSAI, 
                                    bool bCheckIDs,
                                    size_t& readIdx,
                                    size_t& sumBlockSize,
                                    OverlapVector& outVector, 
                                    bool& isSubstring)
{
    OverlapVector outvec;
    std::istringstream convertor(hitString);

    sumBlockSize = 0;
    // Read the overlap blocks for a read
    size_t numBlocks;
    convertor >> readIdx >> isSubstring >> numBlocks;

    for(size_t i = 0; i < numBlocks; ++i)
    {
        // Read the block
        OverlapBlock record;
        convertor >> record;
        //std::cout << "\t" << record << "\n";

        // Iterate through the range and write the overlaps
        for(int64_t j = record.ranges.interval[0].lower; j <= record.ranges.interval[0].upper; ++j)
        {
            sumBlockSize += 1;
            const SuffixArray* pCurrSAI = (record.flags.isTargetRev()) ? pRevSAI : pFwdSAI;
            const ReadInfo& queryInfo = pQueryRIT->getReadInfo(readIdx);

            int64_t saIdx = j;

            // The index of the second read is given as the position in the SuffixArray index
            const ReadInfo& targetInfo = pTargetRIT->getReadInfo(pCurrSAI->get(saIdx).getID());

	    //std::cout << " HERE RER " << queryInfo.id << " " << targetInfo.id << " " << queryInfo.length << " " << targetInfo.length << "\n";
            // Skip self alignments and non-canonical (where the query read has a lexo. higher name)
            if(queryInfo.id != targetInfo.id)
            {    
                Overlap o = record.toOverlap(queryInfo.id, targetInfo.id, queryInfo.length, targetInfo.length);
		//std::cout << queryInfo.id << " " << targetInfo.id << " " << queryInfo.length << " " << targetInfo.length << "\n";

                // The alignment logic above has the potential to produce duplicate alignments
                // To avoid this, we skip overlaps where the id of the first coord is lexo. lower than 
                // the second or the match is a containment and the query is reversed (containments can be 
                // output up to 4 times total).
                if(bCheckIDs && (o.id[0] < o.id[1] || (o.match.isContainment() && record.flags.isQueryRev())))
                    continue;

                outVector.push_back(o);
            }
        }
    }
}
