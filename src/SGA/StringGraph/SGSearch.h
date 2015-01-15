//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGSearch - Algorithms and data structures
// for searching a string graph
//
#ifndef SGSEARCH_H
#define SGSEARCH_H

#include "Bigraph.h"
#include "SGWalk.h"
#include "GraphSearchTree.h"
#include <deque>

// Returns the extension distance indicated
// by the given edge
struct SGDistanceFunction
{
    int operator()(const Edge* pEdge) const
    {
        return pEdge->getSeqLen();
    }
};

// 
typedef GraphSearchTree<Vertex, Edge, SGDistanceFunction> SGSearchTree;

//
struct SGWalkBuilder
{
    public:
        SGWalkBuilder(SGWalkVector& outWalks, bool bIndexWalk);
        ~SGWalkBuilder();

        // These three functions must be provided by the builder object
        // the generic graph code calls these to describe the walks through
        // the graph
        void startNewWalk(Vertex* pStartVertex);
        void addEdge(Edge* pEdge);
        void finishCurrentWalk();

    private:
        SGWalkVector& m_outWalks;
        SGWalk* m_pCurrWalk;
        bool m_bIndexWalk;

};

// String Graph searching algorithms
namespace SGSearch
{
    //
    bool findWalks(Vertex* pX, 
                   Vertex* pY, 
                   EdgeDir initialDir,
                   int maxDistance, 
                   size_t maxNodes, 
                   bool exhaustive,
                   SGWalkVector& outWalks);

    void findVariantWalks(Vertex* pX, 
                          EdgeDir initialDir, 
                          int maxDistance,
                          size_t maxWalks, 
                          SGWalkVector& outWalks);

    void findCollapsedWalks(Vertex* pX, EdgeDir initialDir, 
                            int maxDistance, size_t maxNodes,
                            SGWalkVector& outWalks);

    // Count the number of vertices that span the sequence junction
    // described by edge XY. Returns -1 if the search was not completed
    int countSpanningCoverage(Edge* pXY, size_t maxQueue);

    // Returns true if all the endpoints of the edges in epv are in vertexSet
    bool checkEndpointsInSet(EdgePtrVec& epv, std::set<Vertex*>& vertexSet);
};

#endif
