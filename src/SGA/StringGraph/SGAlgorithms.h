//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGAlgorithms - Algorithms for manipulating
// the string graph
//
#ifndef SGALGORITHMS_H
#define SGALGORITHMS_H

#include "Bigraph.h"
#include "SGUtil.h"
#include <queue>


namespace SGAlgorithms
{

//
// Helper data structures
//

typedef std::set<VertexID> VertexIDSet;
typedef std::set<EdgeDesc> EdgeDescSet;
typedef std::map<EdgeDesc, Overlap> EdgeDescOverlapMap;

// Remodel the graph by finding new edges for the given vertex to avoid
// causing a disconnection when removing pDeleteEdge
void remodelVertexForExcision(StringGraph* pGraph, Vertex* pVertex, Edge* pDeleteEdge);

// Create the edges in pGraph described by the overlap
// A pointer to the first edge of the edge/edge twin is returned or NULL
// if the edges cannot be added
Edge* createEdgesFromOverlap(StringGraph* pGraph, const Overlap& o, bool allowContained, size_t maxEdges = -1);

// Calculate the error rate between the two vertex sequences
double calcErrorRate(const Vertex* pX, const Vertex* pY, const Overlap& ovrXY);

// Update the containment flags in the graph using the described overlaps
void updateContainFlags(StringGraph* pGraph, Vertex* pVertex, EdgeDescOverlapMap& containMap);
void updateContainFlags(StringGraph* pGraph, Vertex* pVertex, const EdgeDesc& ed, const Overlap& ovr);

//
// Overlap inference algorithms
//

// Infer an overlap from two edges
// The input overlaps are between X->Y Y->Z
// and the returned overlap is X->Z
Overlap inferTransitiveOverlap(const Overlap& ovrXY, const Overlap& ovrYZ);
EdgeDesc overlapToEdgeDesc(Vertex* pY, const Overlap& ovrXY);

// Returns true if given overlaps X->Y, X->Z, the overlap X->Z is transitive
bool isOverlapTransitive(const Vertex* pY, const Vertex* pZ, const Overlap& ovrXY, 
                         const Overlap& ovrXZ, const double maxER, const int minOverlap);

// Returns true if XZ has a non-zero length overlap
bool hasTransitiveOverlap(const Overlap& ovrXY, const Overlap& ovrYZ);

// Partition a set of overlaps into transitive and irreducible edges
void partitionTransitiveOverlaps(EdgeDescOverlapMap* pOverlapMap, 
                                 EdgeDescOverlapMap* pTransitive,
                                 double maxER, int minLength);

// Each read can have at most one edge to any other read in a given direction
// This function removes any duplicates
void removeSubmaximalOverlaps(EdgeDescOverlapMap* pOverlapMap);

// Simple getters for std::transform
EdgeDesc getEdgeDescFromEdge(Edge* pEdge);

void printOverlapMap(const EdgeDescOverlapMap& overlapMap);

};

#endif
