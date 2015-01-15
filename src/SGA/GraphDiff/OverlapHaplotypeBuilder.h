///----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapHaplotypeBuilder - Build read coherent haplotypes
// using read overlaps.
//
#ifndef OVERLAP_HAPLOTYPE_BUILDER_H
#define OVERLAP_HAPLOTYPE_BUILDER_H
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
#include <queue>

// A struct holding a kmer that is shared between two reads
struct SharedVertexKmer
{
    size_t kmer_index;
    Vertex* vertex;

    static bool sortByVertex(const SharedVertexKmer& a, const SharedVertexKmer& b) { return a.vertex < b.vertex; }
    static bool equalByVertex(const SharedVertexKmer& a, const SharedVertexKmer& b) { return a.vertex == b.vertex; }
};
typedef std::vector<SharedVertexKmer> SharedVertexKmerVector;

// A tip in the graph that can be extended
struct ExtendableTip
{
    std::string id;
    EdgeDir direction;
};
typedef std::vector<ExtendableTip> ExtendableTipVector;

// Build haplotypes starting from a given sequence.
class OverlapHaplotypeBuilder
{
    public:

        OverlapHaplotypeBuilder(const GraphCompareParameters& parameters);
        ~OverlapHaplotypeBuilder();

        void setInitialHaplotype(const std::string& str);

        // Run the bubble construction process
        // Returns true if the graph was successfully built between the two sequences
        HaplotypeBuilderReturnCode run(StringVector& out_haplotypes);
        
    private:
    
        // Check for walks through the graph that cover all seed vertices
        // If there is a valid covering assembly of the reads, it will be
        // written to the vector
        void checkWalks(StringVector* out_strings);
        
        // Check if the walk set is a valid variant path
        bool areWalksValid(const SGWalkVector& walks);
        
        // Get all the reads that contain at least one of the kmers in the vector
        void getReadsForKmers(const StringVector& kmers, size_t k, StringVector* reads);

        // Correct sequencing errors in the set of reads
        void correctReads(StringVector* reads);

        // Initialize the string graph using the corrected reads
        // Returns true if the graph was successfully initialized
        bool buildInitialGraph(const StringVector& reads);

        // Extend the graph by finding overlaps for the tip vertices
        void extendGraph();

        // Insert a new vertex into the graph with the specified sequence
        void insertVertexIntoGraph(const std::string& prefix, const std::string& sequence);

        // Use the kmer to vertex map to find vertices in the graph that possibly
        // overlap the incoming vertex
        SharedVertexKmerVector getCandidateOverlaps(const Vertex* incoming_vertex);

        // Update the kmer to vertex map to include kmers for the new vertex
        void updateKmerVertexMap(const Vertex* incoming_vertex);

        // Find corrected reads that share a perfect overlap to the input sequence
        StringVector getCorrectedOverlaps(const std::string& sequence, EdgeDir direction);

        // Find tip vertices in the graph
        ExtendableTipVector findTips() const;

        // Trim a tip off the graph.
        void trimTip(Vertex* x, EdgeDir direction);

        // Returns true if the sequence represents a junction in the variation graph
        bool isJoinSequence(const std::string& sequence, EdgeDir dir);

        // Returns true if the vertex x is at a bifurcation of the graph, a dead-end tip
        // or unambiguously connected to a non-join vertex. Used to avoid generating
        // redudant covering walks
        bool isUniqueJoin(Vertex* x, EdgeDir direction);

        // Returns the number of vertices in the vector that are contained in any walk
        size_t countCoveredVertices(const VertexPtrVec& vertices, const SGWalkVector& walks);

        // Coerece the set of reads into an ordered overlapping sequences
        // This function relies on the fact that all reads share the same kmer
        void orderReadsInitial(const std::string& initial_kmer, const StringVector& reads, StringVector* ordered_reads);

        // Build a multiple alignment from an ordered set of reads
        MultipleAlignment buildMultipleAlignment(const StringVector& ordered_vector);

        //
        // Data
        //
        GraphCompareParameters m_parameters;
        ErrorCorrectProcess* m_corrector;
        StringGraph* m_graph;
        size_t m_numReads;

        std::string m_initial_kmer_string;

        // Cache the sequences of reads in the graph to avoid adding redundant reads
        HashSet<std::string> m_used_reads;

        // Cache the corrected sequences for reads in the graph. Since we may visit the same
        // raw sequence multiple times, this lets us avoid redundantly correcting the reads.
        HashMap<std::string, std::string> m_correction_cache;

        // A map from kmer to vertex ids. We use this when inserting new reads into the graph to avoid
        // comparing incoming reads against everything
        HashMap<std::string, StringVector> m_kmer_vertex_id_cache;
        static const size_t m_vertex_map_kmer = 31;
};

#endif
