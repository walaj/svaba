//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldVisitors - Functors that perform
// some operation on a ScaffoldVertex/Graph
//
#ifndef SCAFFOLDVISITORS_H
#define SCAFFOLDVISITORS_H

#include "ScaffoldGraph.h"
#include "SGUtil.h"

//
namespace ScaffoldAlgorithms
{
    bool areEdgesAmbiguous(ScaffoldEdge* pXY, ScaffoldEdge* pXZ);
    void inferScaffoldEdgeYZ(ScaffoldEdge* pXY, ScaffoldEdge* pXZ,
                                             int& dist, double& sd, EdgeDir& dir_yz, 
                                             EdgeDir& dir_zy, EdgeComp& comp);

};

// Write summary statistics of the scaffold graph to stdout
class ScaffoldStatsVisitor
{
    public:
        
        void previsit(ScaffoldGraph* pGraph);
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* pGraph);

    private:
        size_t m_numVertices;
        size_t m_numEdges;
};

// Classify vertices as unique or repetitive based on their
// a-statistic values
class ScaffoldAStatisticVisitor
{
    public:
        
        ScaffoldAStatisticVisitor(double uniqueThreshold, double minCopyNumber);
        void previsit(ScaffoldGraph* /*pGraph*/);
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* /*pGraph*/);

    private:
        double m_uniqueThreshold;
        double m_minCopyNumber;

        size_t m_numUnique;
        size_t m_numRepeat;
        size_t m_numLowCN;
        size_t m_sumUnique;
        size_t m_sumRepeat;
        size_t m_sumLowCN;
};

// Walk the graph validating that the links present are valid
class ScaffoldLinkValidator
{
    public:
        
        ScaffoldLinkValidator(int maxOverlap, double threshold, int verbose);
        void previsit(ScaffoldGraph* /*pGraph*/);
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* /*pGraph*/);

    private:

        int m_maxOverlap;
        double m_threshold;
        size_t m_numUnique;
        size_t m_numRepeat;
        size_t m_numCut;
        int m_verbose;
};

// Walk through the scaffold graph making chains of vertices
class ScaffoldChainVisitor
{
    public:
        
        ScaffoldChainVisitor(int maxOverlap);
        void previsit(ScaffoldGraph* /*pGraph*/);
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* /*pGraph*/);

    private:
        int m_maxOverlap;
};

// Remove transitive edges from the graph
class ScaffoldTransitiveReductionVisitor
{
    public:
        ScaffoldTransitiveReductionVisitor();

        void previsit(ScaffoldGraph* /*pGraph*/);
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* /*pGraph*/);

};

// Detect and remove polymorphmic vertices in the scaffold
class ScaffoldPolymorphismVisitor
{
    public:
        ScaffoldPolymorphismVisitor(int maxOverlap);
        
        void previsit(ScaffoldGraph* /*pGraph*/);
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* /*pGraph*/);

    private:
        int m_maxOverlap;
        int m_numMarked;        
};

// Detect structural variation in the scaffold
class ScaffoldSVVisitor
{
    public:
        ScaffoldSVVisitor(int maxSize);
        
        void previsit(ScaffoldGraph* /*pGraph*/);
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* /*pGraph*/);

    private:
        int m_maxSVSize;
        int m_numMarked;
};

// Remove the vertices with conflicting distance estimates
class ScaffoldConflictingVisitor
{
    public:
        ScaffoldConflictingVisitor() {}
        void previsit(ScaffoldGraph*);
        bool visit(ScaffoldGraph*, ScaffoldVertex*);
        void postvisit(ScaffoldGraph*);

    private:
        int m_numMarked;
};

// Compute a layout of the contigs linked to each vertex
class ScaffoldLayoutVisitor
{
    public:
        ScaffoldLayoutVisitor();
        
        void previsit(ScaffoldGraph* /*pGraph*/);
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* /*pGraph*/);    
};

class ScaffoldDistanceRefinementVisitor
{
    public:
        ScaffoldDistanceRefinementVisitor(const StringGraph* pStringGraph);
        
        void previsit(ScaffoldGraph* /*pGraph*/);
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* /*pGraph*/);

    private:
        const StringGraph* m_pStringGraph;

};

// Walk through the graph breaking up components whenever a vertex has more
// than one edge in a particular direction
class ScaffoldMultiEdgeRemoveVisitor
{
    public:
        
        ScaffoldMultiEdgeRemoveVisitor() {}
        void previsit(ScaffoldGraph* /*pGraph*/) {}
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* /*pGraph*/) {}

};

struct ScaffoldStats
{
    int bases;
    int numContigs;
    int span; // size + gaps

    static bool sortSpanDesc(const ScaffoldStats& a, const ScaffoldStats& b)
    {
        return a.span > b.span;
    }
};

// Output scaffolds
class ScaffoldWriterVisitor
{
    public:
        ScaffoldWriterVisitor(const std::string& filename);
        ~ScaffoldWriterVisitor();

        void previsit(ScaffoldGraph* /*pGraph*/);
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* /*pGraph*/);

    private:
        std::ostream* m_pWriter;
        std::vector<ScaffoldStats> m_statsVector;
};

#endif
