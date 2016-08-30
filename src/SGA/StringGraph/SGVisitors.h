//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGVisitors - Algorithms that visit
// each vertex in the graph and perform some
// operation
//
#include "SGAlgorithms.h"
#include "SGUtil.h"
#include "Util.h"
//#include "contigs.h"
#include "SeqLib/UnalignedSequence.h"

#ifndef SGVISITORS_H
#define SGVISITORS_H

// Visit each node, writing it to a file as a fasta record
struct SGFastaVisitor
{
    // constructor
    SGFastaVisitor(std::string filename) : m_fileHandle(filename.c_str()) {}
    ~SGFastaVisitor() { m_fileHandle.close(); }

    // functions
    void previsit(StringGraph* /*pGraph*/) {}
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph* /*pGraph*/) {}

    // data
    std::ofstream m_fileHandle;
};

// Visit each node, writing it to a file as a fasta record
// jwala modified
struct SGVisitorContig
{
  // constructor
  SGVisitorContig() {}
  ~SGVisitorContig() {}
  
  // functions
  void previsit(StringGraph* /*pGraph*/) {}
  bool visit(StringGraph* pGraph, Vertex* pVertex);
  void postvisit(StringGraph* /*pGraph*/) {}
  
  // data
  SeqLib::UnalignedSequenceVector m_ct;
};

// Visit each node, writing it to a file as a fasta record
// Run the Myers transitive reduction algorithm on each node
struct SGTransitiveReductionVisitor
{
    SGTransitiveReductionVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    int marked_verts;
    int marked_edges;
};

// Remove identical vertices from the graph
struct SGIdenticalRemoveVisitor
{
    SGIdenticalRemoveVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph* pGraph);
    int count;
};

// Remove contained vertices from the graph
struct SGContainRemoveVisitor
{
    SGContainRemoveVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph* pGraph);
};

// Validate that the graph does not contain
// any extra edges or missing irreducible edges
struct SGValidateStructureVisitor
{
    SGValidateStructureVisitor() {}
    void previsit(StringGraph*) {}
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*) {}
};

// Remodel the graph to infer missing edges or remove erroneous edges
struct SGSmallRepeatResolveVisitor
{
    SGSmallRepeatResolveVisitor(int minDiff) : m_minDiff(minDiff) {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    int m_minDiff;
};

// Remove edges from the graph when the ratio between an edge's overlap length
// and the longest overlap for the vertex is less than the given parameter
struct SGOverlapRatioVisitor
{
    SGOverlapRatioVisitor(double minRatio) : m_minRatio(minRatio) {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    double m_minRatio;
};

// Detects and removes small "tip" vertices from the graph
// when they are less than minLength in size
struct SGTrimVisitor
{
    SGTrimVisitor(size_t minLength) : m_minLength(minLength) {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    size_t m_minLength;
    int num_island;
    int num_terminal;
};

// Detect and remove duplicate edges
struct SGDuplicateVisitor
{
    SGDuplicateVisitor(bool silent = false) : m_bSilent(silent) {}
    void previsit(StringGraph*);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    bool m_hasDuplicate;
    bool m_bSilent;
};

// Remove the edges of super-repetitive vertices in the graph
struct SGSuperRepeatVisitor
{
    SGSuperRepeatVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    size_t m_num_superrepeats;
};

// Smooth out variation in the graph
struct SGSmoothingVisitor
{
    SGSmoothingVisitor(std::string filename, 
                       double maxGapDiv, 
                       double maxTotalDiv, 
                       int maxIndelLength) : m_numRemovedTotal(0), 
                                             m_maxGapDivergence(maxGapDiv),
                                             m_maxTotalDivergence(maxTotalDiv),
                                             m_maxIndelLength(maxIndelLength),
                                             m_outFile(filename.c_str()) {}

    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    int m_simpleBubblesRemoved;
    int m_complexBubblesRemoved;
    int m_numRemovedTotal;

    double m_maxGapDivergence;
    double m_maxTotalDivergence;
    int m_maxIndelLength;
    std::ofstream m_outFile;
};

// Compile summary statistics for the graph
struct SGGraphStatsVisitor
{
    SGGraphStatsVisitor() {}
    void previsit(StringGraph* pGraph);
    bool visit(StringGraph* pGraph, Vertex* pVertex);
    void postvisit(StringGraph*);

    int num_terminal;
    int num_island;
    int num_monobranch;
    int num_dibranch;
    int num_simple;
    int num_edges;
    int num_vertex;
    size_t sum_edgeLen;
};

#endif
