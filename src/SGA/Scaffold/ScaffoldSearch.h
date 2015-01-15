//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldSearch - Find walks through a scaffold graph
//
#ifndef SCAFFOLDSEARCH_H
#define SCAFFOLDSEARCH_H

#include "GraphSearchTree.h"
#include "ScaffoldWalk.h"
#include "ScaffoldGraph.h"

//
struct ScaffoldDistanceFunction
{
    int operator()(const ScaffoldEdge* pEdge) const
    {
        return pEdge->getDistance();
    }
};

//
typedef GraphSearchTree<ScaffoldVertex, ScaffoldEdge, ScaffoldDistanceFunction> ScaffoldSearchTree;

//
struct ScaffoldWalkBuilder
{
    public:
        ScaffoldWalkBuilder(ScaffoldWalkVector& outWalks);
        ~ScaffoldWalkBuilder();

        // These three functions must be provided by the builder object
        // the generic graph code calls these to describe the walks through
        // the graph
        void startNewWalk(ScaffoldVertex* pStartVertex);
        void addEdge(ScaffoldEdge* pEdge);
        void finishCurrentWalk();

    private:
        ScaffoldWalkVector& m_outWalks;
        ScaffoldWalk* m_pCurrWalk;
};

//
namespace ScaffoldSearch
{
    void findVariantWalks(ScaffoldVertex* pX, 
                          EdgeDir initialDir, 
                          int maxDistance,
                          size_t maxWalks, 
                          ScaffoldWalkVector& outWalks);

    
    // Return the index of a walk in allWalks that every vertex in coverVector
    // If no covering walk can be found, returns -1
    int findCoveringWalk(const ScaffoldWalkVector& allWalks, 
                         const ScaffoldVertexPtrVector& coverVector);

    // Returns walks from pX to pY that exclusively
    // travel along primary links
    void findPrimaryWalks(ScaffoldVertex* pX,
                          ScaffoldVertex* pY,
                          EdgeDir intialDir,
                          int maxDistance,
                          size_t maxNodes, 
                          ScaffoldWalkVector& outWalks);

    void printWalks(const ScaffoldWalkVector& walkVector);

};

#endif
