//-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// BWTIndexSet.h - Struct holding the different
// types of index for a set of sequences. Used
// to simplify the interface to some functions.
//
#ifndef BWTINDEXSET_H
#define BWTINDEXSET_H

#include "BWT.h"
#include "BWTIntervalCache.h"
#include "SampledSuffixArray.h"
#include "PopulationIndex.h"
#include "QualityTable.h"

// A collection of indices. For some algorithms
// all indices are not necessary so some of these
// can be NULL. The algorithms will check as preconditions
// which indices they need.
struct BWTIndexSet
{
    // Constructor
    BWTIndexSet() : pBWT(NULL), pRBWT(NULL), pCache(NULL), pSSA(NULL), pPopIdx(NULL), pQualityTable(NULL) {}

    // Data
    const BWT* pBWT;
    const BWT* pRBWT;
    const BWTIntervalCache* pCache;
    const SampledSuffixArray* pSSA;
    const PopulationIndex* pPopIdx;
    const QualityTable* pQualityTable;
};

#endif
