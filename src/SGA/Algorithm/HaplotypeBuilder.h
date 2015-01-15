///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// HaplotypeBuilder - Construct candidate
// haplotypes from a pair of k-mer seeds.
//
#ifndef HAPLOTYPE_BUILDER_H
#define HAPLOTYPE_BUILDER_H
#include "BWT.h"
#include "BWTInterval.h"
#include "SGUtil.h"
#include "SGWalk.h"
#include "VariationBuilderCommon.h"
#include <queue>

// Structs
struct AnchorSequence
{
    std::string sequence;
    int position;
    int count;

    friend std::ostream& operator<<(std::ostream& out, const AnchorSequence& a) 
    { 
        out << a.position << " " << a.sequence << " " << a.count; 
        return out; 
    }
};
typedef std::vector<AnchorSequence> AnchorVector;

// The result, a set of haplotypes
struct HaplotypeBuilderResult
{
    StringVector haplotypes;
};

// Return codes to indicate why the process may have failed
enum HaplotypeBuilderReturnCode
{
    HBRC_OK,
    HBRC_TOO_MANY_VERTICES,
    HBRC_NO_PATH,
    HBRC_WALK_FAILED,
};

//
// Class to build a variant bubble starting at a particular sequence
//
class HaplotypeBuilder
{
    public:

        HaplotypeBuilder();
        ~HaplotypeBuilder();

        void setTerminals(const AnchorSequence& leftAnchor, const AnchorSequence& rightAnchor);
        void setIndex(const BWT* pBWT, const BWT* pRBWT);

        // Set the threshold of kmer occurrences to use it as an edge
        void setKmerParameters(size_t k, size_t t);
    
        // Run the bubble construction process
        // Returns true if the graph was successfully built between the two sequences
        HaplotypeBuilderReturnCode run();
        
        // Parse walks from the constructed graph
        HaplotypeBuilderReturnCode parseWalks(HaplotypeBuilderResult& results) const;

    private:
        
        // Add a vertex to the graph and record the sequence coverage value
        void addVertex(Vertex* pVertex, int coverage);

        // Make the sequence of a new deBruijn vertex using the edge details
        std::string makeDeBruijnVertex(const std::string& v, char edgeBase, EdgeDir direction);
        void addDeBruijnEdges(const Vertex* pX, const Vertex* pY, EdgeDir direction);

        // Count the number of extensions of a de Bruijn node that are above
        // the required k-mer coverage
        size_t countValidExtensions(const AlphaCount64& ac) const;
        
        //
        // Data
        //
        const BWT* m_pBWT;
        const BWT* m_pRevBWT;

        StringGraph* m_pGraph;
        StrIntMap m_vertexCoverageMap;

        BuilderExtensionQueue m_queue;
        Vertex* m_pStartVertex;
        Vertex* m_pJoinVertex;
        
        //
        size_t m_kmerThreshold;
        size_t m_kmerSize;
};


#endif
