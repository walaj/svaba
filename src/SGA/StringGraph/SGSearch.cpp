//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGSearch - Algorithms and data structures
// for searching a string graph
//
#include "SGSearch.h"
#include <queue>

//
SGWalkBuilder::SGWalkBuilder(SGWalkVector& outWalks, bool bIndexWalk) : m_outWalks(outWalks), m_pCurrWalk(NULL), m_bIndexWalk(bIndexWalk)
{
    
}

//
SGWalkBuilder::~SGWalkBuilder()
{
    // The pointer to the current walk should be NULL
    // or else finishCurrentWalk() was not called for the last
    // walk and the graph search algorithm has a bug
    assert(m_pCurrWalk == NULL);
}

//
void SGWalkBuilder::startNewWalk(Vertex* pStartVertex)
{
    m_pCurrWalk = new SGWalk(pStartVertex, m_bIndexWalk);
}

//
void SGWalkBuilder::addEdge(Edge* pEdge)
{
    m_pCurrWalk->addEdge(pEdge);
}

//
void SGWalkBuilder::finishCurrentWalk()
{
    m_outWalks.push_back(*m_pCurrWalk);
    delete m_pCurrWalk;
    m_pCurrWalk = NULL;
}

// Find all the walks between pX and pY that are within maxDistance
// If the exhaustive flag is set, only return walks if all the possible
// solutions have been found. If exhaustive is false, any walks found will be
// returned in outWalks even if the search is aborted.
// Returns true if all the possible walks were found.
bool SGSearch::findWalks(Vertex* pX, Vertex* pY, EdgeDir initialDir,
                         int maxDistance, size_t maxNodes, bool exhaustive, SGWalkVector& outWalks)
{
    SGSearchTree searchTree(pX, pY, initialDir, maxDistance, maxNodes);

    // Iteravively perform the BFS using the search tree.
    while(searchTree.stepOnce()) { }

    // If the search was aborted, do not return any walks
    // because we do not know if there are more valid paths from pX
    // to pY that we could not find because the search space was too large
    if(!searchTree.wasSearchAborted() || !exhaustive)
    {
        // Extract the walks from the graph as a vector of edges
        SGWalkBuilder builder(outWalks, false);
        searchTree.buildWalksToGoal(builder);
    }
    return !searchTree.wasSearchAborted();
}

// Search the graph for a set of walks that represent alternate
// versions of the same sequence. These walks are found by searching
// the graph for a set of walks that start/end at a common vertex and cover
// every vertex inbetween the start/end points. Additionally, the internal
// vertices cannot have an edge to any vertex that is not in the set of vertices
// indicated by the walks. Finally, all walks must have the same orientation for the
// end vertex.
// If these conditions are met, we can retain one of the walks and cleanly remove
// the others.
void SGSearch::findVariantWalks(Vertex* pX, 
                                EdgeDir initialDir, 
                                int maxDistance,
                                size_t maxWalks, 
                                SGWalkVector& outWalks)
{
    findCollapsedWalks(pX, initialDir, maxDistance, 500, outWalks);

    if(outWalks.size() <= 1 || outWalks.size() > maxWalks)
    {
        outWalks.clear();
        return;
    }

    Edge* pLastEdge = outWalks.front().getLastEdge();
    Vertex* pLastVertex = pLastEdge->getEnd();
    EdgeDir lastDir = pLastEdge->getTwinDir();

    // Validate that any of the returned walks can be removed cleanly.
    // This means that for all the internal vertices on each walk (between
    // the common endpoints) the only links are to other vertices in the set,
    // and that the last vertex has the same orientation for all walks.

    // Construct the set of vertex IDs
    std::set<Vertex*> completeVertexSet;
    bool sameOrientation = true;
    for(size_t i = 0; i < outWalks.size(); ++i)
    {
        if (outWalks[i].getLastEdge()->getTwinDir() != lastDir)
        {
            sameOrientation = false;
            break;
        }

        VertexPtrVec verts = outWalks[i].getVertices();
        for(size_t j = 0; j < verts.size(); ++j)
        {
            completeVertexSet.insert(verts[j]);
        }
    }

    // Return if the walks do not have the same orientation for the last vertex
    if(!sameOrientation)
    {
        outWalks.clear();
        return;
    }

    // Check that each vertex in the internal nodes only has
    // edges to the vertices in the set

    bool cleanlyRemovable = true;

    // Ensure that all the vertices linked to the start vertex
    // in the specified dir are present in the set.
    EdgePtrVec epv = pX->getEdges(initialDir);
    cleanlyRemovable = checkEndpointsInSet(epv, completeVertexSet);

    // Ensure that all the vertex linked to the last vertex
    // in the incoming direction are preset
    epv = pLastVertex->getEdges(lastDir);
    cleanlyRemovable = cleanlyRemovable && checkEndpointsInSet(epv, completeVertexSet);

    // Check that each vertex connected to an interval vertex is also present
    for(std::set<Vertex*>::iterator iter = completeVertexSet.begin(); iter != completeVertexSet.end(); ++iter)
    {
        Vertex* pY = *iter;
        if(pY == pX || pY == pLastVertex)
            continue;
        epv = pY->getEdges();
        cleanlyRemovable = cleanlyRemovable && checkEndpointsInSet(epv, completeVertexSet);
    }

    if(!cleanlyRemovable)
    {
        outWalks.clear();
    }
    return;
}

// Check that all the endpoints of the edges in the edge pointer vector
// are members of the set
bool SGSearch::checkEndpointsInSet(EdgePtrVec& epv, std::set<Vertex*>& vertexSet)
{
    for(size_t i = 0; i < epv.size(); ++i)
    {
        if(vertexSet.find(epv[i]->getEnd()) == vertexSet.end())
            return false;
    }
    return true;
}

// Return a set of walks that all start from pX and join together at some later vertex
// If no such walk exists, an empty set is returned
void SGSearch::findCollapsedWalks(Vertex* pX, EdgeDir initialDir, 
                                  int maxDistance, size_t maxNodes, 
                                  SGWalkVector& outWalks)
{
    SGSearchTree searchTree(pX, NULL, initialDir, maxDistance, maxNodes);

    // Iteravively perform the BFS using the search tree. After each step
    // we check if the search has collapsed to a single vertex.
    bool done = false;
    while(!done)
    {
        done = !searchTree.stepOnce();
        if(searchTree.wasSearchAborted())
            break;

        Vertex* pCollapsedVertex;
        bool isCollapsed = searchTree.hasSearchConverged(pCollapsedVertex);
        if(isCollapsed)
        {
            assert(pCollapsedVertex != NULL);
            
            // pCollapsedVertex is common between all walks.
            // Check that the extension distance for any walk is no longer than maxDistance.
            // If this is the case, we truncate all walks so they end at pCollapsedVertex and return 
            // all the non-redundant walks in outWalks
            SGWalkBuilder builder(outWalks, true);
            searchTree.buildWalksContainingVertex(pCollapsedVertex, builder);
            return;
        }
    }

    // no collapsed walk found, return empty set
    outWalks.clear();
}

// Count the number of reads that span the junction described by edge XY
// X --------------
// Y        ------------
// Z            -----------
// W                ----------
// In this case Z spans the junction but W does not. 
int SGSearch::countSpanningCoverage(Edge* pXY, size_t maxQueue)
{
    (void)pXY;
    (void)maxQueue;
    
    assert(false && "deprecated");
    return 0;
}
