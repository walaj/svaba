///----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// StringHaplotypeBuilder - Build read coherent haplotypes
// using read overlaps.
//
#ifndef STRING_HAPLOTYPE_BUILDER_H
#define STRING_HAPLOTYPE_BUILDER_H
#include "BWT.h"
#include "SampledSuffixArray.h"
#include "BWTInterval.h"
#include "BWTIntervalCache.h"
#include "SGUtil.h"
#include "SGWalk.h"
#include "VariationBuilderCommon.h"
#include "HaplotypeBuilder.h"
#include "multiple_alignment.h"
#include "GraphCompare.h"
#include "ErrorCorrectProcess.h"
#include "SGWalk.h"
#include "KmerOverlaps.h"
#include "OverlapExtractorWithCorrection.h"
#include <queue>

// Build haplotypes starting from a given sequence.
class StringHaplotypeBuilder
{
    public:

        StringHaplotypeBuilder(const GraphCompareParameters& parameters);
        ~StringHaplotypeBuilder();

        void setInitialHaplotype(const std::string& str);

        // Run the bubble construction process
        // Returns true if the graph was successfully built between the two sequences
        HaplotypeBuilderReturnCode run(StringVector& out_haplotypes);
        
    private:
    
        // Get all the reads that contain at least one of the kmers in the vector
        void getReadsForKmers(const StringVector& kmers, size_t k, StringVector* reads);

        // Correct sequencing errors in the set of reads
        void correctReads(StringVector* reads);

        //
        SequenceOverlapPairVector computeExtensions(std::string sequence, EdgeDir direction);
        
        // Get overlaps for the given sequence
        SequenceOverlapPairVector getCorrectedOverlaps(const std::string& sequence, EdgeDir direction);
        SequenceOverlapPairVector getCorrectedOverlaps2(const std::string& sequence, EdgeDir direction);
        
        // Add an edge to the graph between the described vertices
        void addEdge(Vertex* vertex1, Vertex* vertex2, SequenceOverlap overlap);

        // Build the graph in the given direction, starting from a seed vertex
        Vertex* extendSeed(Vertex* seed_vertex, EdgeDir direction);

        // Get the label of an edge to s2, represented by the overlap object
        std::string getEdgeLabel(const std::string& s2, const SequenceOverlap& overlap);

        // Remove tips from the graph
        void trimTip(Vertex* x, EdgeDir direction);

        // Returns true if the sequence represents a junction in the variation graph
        bool isJoinSequence(const std::string& sequence, EdgeDir dir);

        //
        // Data
        //
        GraphCompareParameters m_parameters;
        ErrorCorrectProcess* m_corrector;
        IOverlapExtractor* m_extractor;

        StringGraph* m_graph;
        size_t m_numReads;
        double m_overlapLengthFrac;

        std::string m_initial_kmer_string;

        // Cache the sequences of reads in the graph to avoid adding redundant reads
        HashSet<std::string> m_used_reads;

        // A map from kmer to vertex ids. We use this when inserting new reads into the graph to avoid
        // comparing incoming reads against everything
        HashMap<std::string, StringVector> m_kmer_vertex_id_cache;
        static const size_t m_vertex_map_kmer = 31;
};

#endif
