//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// RemovalAlgorithm - Determine the set of edges that
// must be added to a given vertex to keep the structure of the
// graph intact when a vertex that it shares an edge with
// is going to be removed from the graph.
// This is similar to the CompleteOverlapSet but does
// not perform an exhaustive exploration of the edges
//
#ifndef REMOVALALGORITHM_H
#define REMOVALALGORITHM_H

#include "Bigraph.h"
#include "SGAlgorithms.h"
#include "CompleteOverlapSet.h"

namespace RemovalAlgorithm
{
    
// Returns the set of overlaps that must be added to the graph
// if the vertex at the end of pRemovalEdge is going to be deleted
SGAlgorithms::EdgeDescOverlapMap computeRequiredOverlaps(const Vertex* pVertex, const Edge* pRemovalEdge, double maxER, int minLength);
void findPotentialOverlaps(const Vertex* pX, const Edge* pRemovalEdge, double maxER, int minLength, SGAlgorithms::EdgeDescOverlapMap& outMap);
void eliminateReachableEdges(const Vertex* pVertex, const Edge* pRemovalEdge, double maxER, int minLength, SGAlgorithms::EdgeDescOverlapMap& outMap);
void enqueueEdges(const Vertex* pY, EdgeDir dirY, const Overlap& ovrXY, const EdgeDesc& edXY, int minOverlap, ExploreQueue& outQueue, EdgeDescSet* pSeenSet);

};

#endif
