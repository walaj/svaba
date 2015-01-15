//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldWalk - a walk through a scaffold graph
//
#ifndef SCAFFOLDWALK_H
#define SCAFFOLDWALK_H

#include <vector>
#include "ScaffoldEdge.h"

typedef std::vector<ScaffoldEdge*> ScaffoldEdgePtrVector;
typedef std::vector<ScaffoldVertex*> ScaffoldVertexPtrVector;

// A walk through the scaffold graph
class ScaffoldWalk
{
    public:
        ScaffoldWalk(ScaffoldVertex* pStartVertex);
        void addEdge(ScaffoldEdge* pEdge);

        // Find the index of the first instance of pvertex in the walk
        // if pvertex is not found, returns -1
        int findVertex(ScaffoldVertex* pVertex) const;

        // Find the orientation of the vertex with respect
        // to the beginning of the walk. This function assumes
        // pVertex is in the walk.
        EdgeComp findOrientation(ScaffoldVertex* pVertex) const;

        // Returns all the ordered list of vertices in the walk
        ScaffoldVertexPtrVector getVertices() const;
        ScaffoldEdgePtrVector getEdges() const;

        // returns the sum of the gaps plus contig lengths from
        // the starting vertex to pVertex.
        // Precondition: pVertex is required to be in the walk
        int getDistanceToVertex(ScaffoldVertex* pVertex) const;

        // Returns the sum of gap distances in the walk
        int64_t getGapSum() const;

        // Returns the sum of contig lengths
        int64_t getContigLengthSum() const;

        // Return the last vertex in the walk
        ScaffoldVertex* getStartVertex() const;
        ScaffoldVertex* getLastVertex() const;

        void print() const;
        void printDot(std::ostream& out) const;

    private:
        ScaffoldVertex* m_pStartVertex;
        ScaffoldEdgePtrVector m_edges;

};
typedef std::vector<ScaffoldWalk> ScaffoldWalkVector;

#endif
