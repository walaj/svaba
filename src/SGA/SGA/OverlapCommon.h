//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapCommon - Common functions used in overlap methods
//
#ifndef OVERLAPCOMMON_H
#define OVERLAPCOMMON_H

#include "Util.h"
#include "overlap.h"
#include "SuffixArray.h"
#include "SGACommon.h"
#include "Timer.h"
#include "ReadInfoTable.h"

namespace OverlapCommon
{

// Parse a line from a .hits file into a variety of overlap-related data
void parseHitsString(const std::string& hitString, 
                     const ReadInfoTable* pQueryRIT, 
                     const ReadInfoTable* pTargetRIT, 
                     const SuffixArray* pFwdSAI, 
                     const SuffixArray* pRevSAI,
                     bool bCheckIDs,
                     size_t& readIdx, 
                     size_t& sumBlockSize,
                     OverlapVector& outVector, 
                     bool& isSubstring);
};

#endif
