///----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapExtractorWithCorrection - Extract reads from an 
// FM-index that share a kmer with a query sequence. Performs
// error correction on the extracted reads, such that exact
// post-correction matches can be found.
//
#ifndef OVERLAP_EXTRACTOR_WITH_CORRECTION_H
#define OVERLAP_EXTRACTOR_WITH_CORRECTION_H

#include <string>
#include <map>
#include "overlapper.h"
#include "BWTIndexSet.h"
#include "KmerOverlaps.h"
#include "IOverlapExtractor.h"
#include "ErrorCorrectProcess.h"
#include "HashMap.h"

class OverlapExtractorWithCorrection : public IOverlapExtractor
{
    public:
        OverlapExtractorWithCorrection(size_t k, const BWTIndexSet& index_set, 
                                       ErrorCorrectProcess* corrector, int minOverlap, double minIdentity);
        
        // Extract reads from the FM-index that share a k-mer with the query sequence,
        // correct them, then return the subset that matches the query exactly
        SequenceOverlapPairVector getExactOverlaps(const std::string& query);

    private:

        // Get inexact overlaps between the query and reads by direct lookup in the FM-index
        void getRawOverlapsDirect(const std::string& query, SequenceOverlapPairVector* out_vector);

        // Get inexact overlaps between the query and reads using the kmer cache
        void getRawOverlapsCached(const std::string& query, bool is_reverse, SequenceOverlapPairVector* out_vector);
    
        // Convert raw overlaps to corrected, exact overlaps
        void getCorrectedExactOverlaps(const std::string& query, 
                                       const SequenceOverlapPairVector* raw_vector, 
                                       SequenceOverlapPairVector* out_vector);

        //
        void populateOutput(const std::string& query, bool is_reverse, SequenceOverlapPairVector* out_vector);

        // Extract reads with k-mers matching the query from the FM-index
        // and update the cache
        void extractAndUpdate(const std::string& query);
        
        // Update the collection of cached sequences
        size_t updateStrings(const std::string& incoming);

        // Update the cache map
        void updateCache(const std::string& str, size_t index);
        
        // Data
        size_t m_k;
        BWTIndexSet m_index_set;
        ErrorCorrectProcess* m_corrector;
        int m_minOverlap;
        double m_minIdentity;

        // The raw reads extracted from the FM-index
        std::vector<std::string> m_strings;

        // A map from raw read to error corrected read
        HashMap<std::string, std::string> m_correction_cache;
        
        // A cache from k-mers that have been used to extract reads to their index in m_strings
        std::map<std::string, std::vector<size_t> > m_cache_map;
};

#endif

