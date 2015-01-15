//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGWalk - Data structure holding a walk through
// the graph
//
#ifndef SGWALK_H
#define SGWALK_H

#include "Bigraph.h"
#include <deque>

// A walk on the string graph is given by the starting vertex
// then a vector of edges used in the walk
enum SGWalkType
{
    SGWT_START_TO_END, // include every vertex of the walk
    SGWT_EXTENSION, // do not include the start vertex
    SGWT_INTERNAL // do not include the start and end vertex
};

// Placement data for the position of a vertex sequence in a walk
struct SGWalkVertexPlacement
{
    Vertex* pVertex;
    int position; // 0-based
    bool isRC; 
};
typedef std::vector<SGWalkVertexPlacement> SGWalkVertexPlacementVector;

class SGWalk
{
    public:
        
        SGWalk(Vertex* pStartVertex, bool bIndexWalk = false);
        SGWalk(const SGWalk& other);
        SGWalk(const EdgePtrVec& edgeVec, bool bIndexWalk = false);

        ~SGWalk();

        SGWalk& operator=(const SGWalk& other);

        void addEdge(Edge* pEdge);
        void popLast();
        void setFinished(bool b);
        bool isFinished() const;

        VertexPtrVec getVertices() const;
        Vertex* getStartVertex() const;
        Vertex* getLastVertex() const;

        // Returns a vector of EdgeComps of the orientation of each
        // vertex in the path with respect to the start of the walk
        std::vector<EdgeComp> getOrientationsToStart() const;

        Edge* getFirstEdge() const;
        Edge* getLastEdge() const;
        Edge* getEdge(size_t idx) const;
        Vertex* getVertex(size_t idx) const;

        size_t getNumVertices() const;
        size_t getNumEdges() const;

        // Returns true if the walk contains the specified vertex
        // If the walk is not indexed, this will assert
        bool isIndexed() const;
        bool containsVertex(const VertexID& id) const;

        // Truncate the walk after the first instance of id
        void truncate(const VertexID& id);

        // Get the string this walk represents. This function supports multiple different walk types
        // (see the enum description). If the pointer to the VertexPlacementVector is not NULL,
        // the position of the vertices within the walk are written there.
        std::string getString(SGWalkType type, SGWalkVertexPlacementVector* pPlacementVector = NULL) const;
        
        // Get the substring of the full path string starting from position fromX
        // to position toY on the first and last vertices, respectively
        std::string getFragmentString(const Vertex* pX, const Vertex* pY,
                                      int fromX, int toY, 
                                      EdgeDir dirX, EdgeDir dirY) const;

        // Return a string identifying the path through the graph this walk represents
        std::string pathSignature() const;

        // distance calculations
        int getExtensionDistance() const;
        int getEndToEndDistance() const;
        int getStartToEndDistance() const;
        int getEndToStartDistance() const;

        // Compare a pair of walks by the total string length
        // Returns true if the walk of a is longer than the walk of b
        static bool compareByTotalLength(const SGWalk& a, const SGWalk& b);

        void print() const;
        void printSimple() const;

    private:
        
        Vertex* m_pStartVertex;
        EdgePtrVec m_edges;
        
        typedef std::set<VertexID> WalkIndex;
        WalkIndex* m_pWalkIndex;

        // The distance from the end of pStart to the last vertex in the walk
        // This equals the length of the extension string
        // x -----------
        // y     ------------
        // z              ----------
        // distance     ************
        int m_extensionDistance;
        bool m_extensionFinished;
};

typedef std::vector<SGWalk> SGWalkVector;
typedef std::deque<SGWalk> WalkQueue;

#endif
