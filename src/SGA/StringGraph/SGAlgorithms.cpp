//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGAlgorithms - Collection of algorithms for operating on string graphs
//
#include "SGAlgorithms.h"
#include "SGUtil.h"
#include "CompleteOverlapSet.h"
#include "RemovalAlgorithm.h"
#include <iterator>

// add edges to the graph for the given overlap
Edge* SGAlgorithms::createEdgesFromOverlap(StringGraph* pGraph, const Overlap& o, bool allowContained, size_t maxEdges)
{
    // Initialize data and perform checks
    Vertex* pVerts[2];
    //Vertex* pVerts = new Vertex[2]; //JEREMIAH PUT ON HEAP
    EdgeComp comp = (o.match.isRC()) ? EC_REVERSE : EC_SAME;

    bool isContainment = o.match.isContainment();
    assert(allowContained || !isContainment);
    (void)allowContained;
    for(size_t idx = 0; idx < 2; ++idx)
    {
        pVerts[idx] = pGraph->getVertex(o.id[idx]);

        // If one of the vertices is not in the graph, skip this edge
        // This can occur if one of the verts is a strict substring of some other vertex so it will
        // never be added to the graph
        if(pVerts[idx] == NULL)
            return NULL;
    }

    // Check if this is a substring containment, if so mark the contained read
    // but do not create edges
    for(size_t idx = 0; idx < 2; ++idx)
    {
        if(!o.match.coord[idx].isExtreme())
        {
	  std::cerr << "contained "; //JEREMIAH

            size_t containedIdx = 1 - idx;
            assert(o.match.coord[containedIdx].isExtreme());
            pVerts[containedIdx]->setColor(GC_RED);
            pGraph->setContainmentFlag(true);
            return NULL;
        }
    }

    // If either vertex has the maximum number of edges,
    // do not add any more. This is to protect against ultra-dense
    // regions of the graph inflating memory usage. The nodes that reach
    // this limit, and nodes connected to them are marked as super repeats.
    // After loading the graph, all edges to super repeats are cut to prevent
    // misassembly.
    size_t num_edges_0 = pVerts[0]->countEdges();
    size_t num_edges_1 = pVerts[1]->countEdges();
    if(num_edges_0 > maxEdges || num_edges_1 > maxEdges)
    {
        //WARN_ONCE("Edge limit reached for vertex when loading graph");
        pVerts[0]->setSuperRepeat(true);
        pVerts[1]->setSuperRepeat(true);
        return NULL;
    }

    if(!isContainment)
    {
        Edge* pEdges[2];
        for(size_t idx = 0; idx < 2; ++idx)
        {
            EdgeDir dir = o.match.coord[idx].isLeftExtreme() ? ED_ANTISENSE : ED_SENSE;
            const SeqCoord& coord = o.match.coord[idx];
            //pEdges[idx] = new(pGraph->getEdgeAllocator()) Edge(pVerts[1 - idx], dir, comp, coord); // JEREMIAH
            pEdges[idx] = new Edge(pVerts[1 - idx], dir, comp, coord); // JEREMIAH
        }

        pEdges[0]->setTwin(pEdges[1]);
        pEdges[1]->setTwin(pEdges[0]);
         
        pGraph->addEdge(pVerts[0], pEdges[0]);
        pGraph->addEdge(pVerts[1], pEdges[1]);
        return pEdges[0];
    }
    else
    {
        // Contained edges don't have a direction, they can be travelled from
        // one vertex to the other in either direction. Hence, we 
        // add two edges per vertex. Later during the contain removal
        // algorithm this is important to determine transitivity
        Edge* pEdges[4];
        for(size_t idx = 0; idx < 2; ++idx)
        {
            const SeqCoord& coord = o.match.coord[idx];
	    //            pEdges[idx] = new(pGraph->getEdgeAllocator()) Edge(pVerts[1 - idx], ED_SENSE, comp, coord);
            //            pEdges[idx + 2] = new(pGraph->getEdgeAllocator()) Edge(pVerts[1 - idx], ED_ANTISENSE, comp, coord);
            pEdges[idx] = new Edge(pVerts[1 - idx], ED_SENSE, comp, coord);
            pEdges[idx + 2] = new Edge(pVerts[1 - idx], ED_ANTISENSE, comp, coord);

        }
        
        // Twin the edges and add them to the graph
        pEdges[0]->setTwin(pEdges[1]);
        pEdges[1]->setTwin(pEdges[0]);

        pEdges[2]->setTwin(pEdges[3]);
        pEdges[3]->setTwin(pEdges[2]);
    
        pGraph->addEdge(pVerts[0], pEdges[0]);
        pGraph->addEdge(pVerts[0], pEdges[2]);

        pGraph->addEdge(pVerts[1], pEdges[1]);
        pGraph->addEdge(pVerts[1], pEdges[3]);
     
        // Set containment flags
        updateContainFlags(pGraph, pVerts[0], pEdges[0]->getDesc(), o);
        return pEdges[0];
    }
}

// Find new edges for pVertex that are required if pDeleteEdge is removed from the graph
void SGAlgorithms::remodelVertexForExcision(StringGraph* pGraph, Vertex* pVertex, Edge* pDeleteEdge)
{
    assert(pVertex == pDeleteEdge->getStart());
    // If the edge is a containment edge, nothing needs to be done. No edges can be transitive
    // through containments
    if(pDeleteEdge->getOverlap().isContainment())
        return;

    double maxER = pGraph->getErrorRate();
    int minLength = pGraph->getMinOverlap();
    
    EdgeDescOverlapMap addMap = RemovalAlgorithm::computeRequiredOverlaps(pVertex, pDeleteEdge, maxER, minLength);
    for(EdgeDescOverlapMap::iterator iter = addMap.begin();
        iter != addMap.end(); ++iter)
    {
        //std::cout << "Adding edge " << iter->second << " during removal of " << pDeleteEdge->getEndID() << "\n";
        createEdgesFromOverlap(pGraph, iter->second, false);
    }

    /*
    // Set the contain flags based on newly discovered edges
    updateContainFlags(pGraph, pVertex, containMap);
    */
}

// Set containment flags in the graph based on the overlap map
void SGAlgorithms::updateContainFlags(StringGraph* pGraph, Vertex* pVertex, EdgeDescOverlapMap& containMap)
{
    for(EdgeDescOverlapMap::iterator iter = containMap.begin(); iter != containMap.end(); ++iter)
    {
        updateContainFlags(pGraph, pVertex, iter->first, iter->second);
    }
}

//
void SGAlgorithms::updateContainFlags(StringGraph* pGraph, Vertex* pVertex, const EdgeDesc& ed, const Overlap& ovr)
{
    assert(ovr.isContainment());
    // Determine which of the two vertices is contained
    Vertex* pOther = ed.pVertex;
    if(ovr.getContainedIdx() == 0)
        pVertex->setContained(true);
    else
        pOther->setContained(true);
    pGraph->setContainmentFlag(true);
}

// Calculate the error rate between the two vertices
double SGAlgorithms::calcErrorRate(const Vertex* pX, const Vertex* pY, const Overlap& ovrXY)
{
    int num_diffs = ovrXY.match.countDifferences(pX->getSeq().toString(), pY->getSeq().toString());
    return static_cast<double>(num_diffs) / static_cast<double>(ovrXY.match.getMinOverlapLength());
}

// Infer an overlap from two edges
// The input edges are between X->Y Y->Z
// and the returned overlap is X->Z
Overlap SGAlgorithms::inferTransitiveOverlap(const Overlap& ovrXY, const Overlap& ovrYZ)
{
    // Construct the match
    Match match_yx = ovrXY.match;
    match_yx.swap(); 
    Match match_yz = ovrYZ.match;

    // Infer the match_ij based match_i and match_j
    Match match_xz = Match::infer(match_yx, match_yz);
    match_xz.expand();

    // Convert the match to an overlap
    Overlap ovr(ovrXY.id[0], ovrYZ.id[1], match_xz);
    return ovr;
}

// Returns true if, given overlaps X->Y, X->Z, the overlap X->Z is transitive
bool SGAlgorithms::isOverlapTransitive(const Vertex* pY, const Vertex* pZ, const Overlap& ovrXY, 
                                       const Overlap& ovrXZ, const double maxER, const int minOverlap)
{
    // Ensure the overlaps are in the correct order, the overlap with Y should be at least
    // as large as the overlap with Z
    assert(ovrXY.getOverlapLength(0) >= ovrXZ.getOverlapLength(0));
    assert(pY != pZ);

    // Compute the overlap YZ
    Overlap ovrYX = ovrXY;
    ovrYX.swap();
    Overlap ovrYZ = SGAlgorithms::inferTransitiveOverlap(ovrYX, ovrXZ);

    // If ovrYZ is a containment, then Z is not transitive wrt Y
    if(ovrYZ.match.isContainment())
        return false;

    // If the overlap between Y and Z is not long enough then Z is not transitive
    if(ovrYZ.getOverlapLength(0) < minOverlap)
        return false;

    // Finally, check that the error rate is below the threshold
    double error_rate = SGAlgorithms::calcErrorRate(pY, pZ, ovrYZ);
    if(isErrorRateAcceptable(error_rate, maxER))
        return true;
    else
        return false;
}

// The following algorithms use these local types
// defining an edge/overlap pair.
typedef std::pair<EdgeDesc, Overlap> EdgeDescOverlapPair;

// Compare two edges by their overlap length
struct EDOPairCompare
{
    bool operator()(const EdgeDescOverlapPair& edpXY, const EdgeDescOverlapPair& edpXZ) {
        return edpXY.second.match.coord[0].length() < edpXZ.second.match.coord[0].length();
    }
};

// typedefs
typedef std::priority_queue<EdgeDescOverlapPair, 
                            std::vector<EdgeDescOverlapPair>,
                            EDOPairCompare> EDOPairQueue;

// Move the transitive edges from pOverlapMap to pTransitive
void SGAlgorithms::partitionTransitiveOverlaps(EdgeDescOverlapMap* pOverlapMap, 
                                               EdgeDescOverlapMap* pTransitive,
                                               double maxER, int minLength)
{
    EDOPairQueue overlapQueue;
    for(SGAlgorithms::EdgeDescOverlapMap::iterator iter = pOverlapMap->begin();
        iter != pOverlapMap->end(); ++iter)
    {
        overlapQueue.push(*iter);
    }

    // Traverse the list of overlaps in order of length and move elements from
    // the irreducible map to the transitive map
    while(!overlapQueue.empty())
    {
        EdgeDescOverlapPair edoPair = overlapQueue.top();
        overlapQueue.pop();

        EdgeDesc& edXY = edoPair.first;
        Overlap& ovrXY = edoPair.second;

        assert(!ovrXY.match.isContainment());

        SGAlgorithms::EdgeDescOverlapMap::iterator iter = pOverlapMap->begin();
        while(iter != pOverlapMap->end())
        {
            bool move = false;
            const EdgeDesc& edXZ = iter->first;
            const Overlap& ovrXZ = iter->second;

            // Four conditions must be met to mark an edge X->Z transitive through X->Y
            // 1) The overlaps must be in the same direction
            // 2) The overlap X->Y must be strictly longer than X->Z
            // 3) The overlap between Y->Z must not be a containment
            // 4) The overlap between Y->Z must be within the error and length thresholds
            if(!(edXZ == edXY) && edXY.dir == edXZ.dir && ovrXY.getOverlapLength(0) > ovrXZ.getOverlapLength(0))
            {
                move = SGAlgorithms::isOverlapTransitive(edXY.pVertex, edXZ.pVertex, ovrXY, ovrXZ, maxER, minLength);
            }
            
            if(move)
            {
                //std::cout << "Marking overlap: " << iter->second << " as trans via " << ovrXY << "\n";
                if(pTransitive != NULL)
                    pTransitive->insert(*iter);
                pOverlapMap->erase(iter++);
            }
            else
            {
                ++iter;
            }
        }
    }
}

void SGAlgorithms::removeSubmaximalOverlaps(EdgeDescOverlapMap* pOverlapMap)
{
    EDOPairQueue overlapQueue;
    for(SGAlgorithms::EdgeDescOverlapMap::iterator iter = pOverlapMap->begin();
        iter != pOverlapMap->end(); ++iter)
    {
        overlapQueue.push(*iter);
    }

    // Traverse the list of overlaps in order of length
    // Only add the first seen overlap for each vertex
    pOverlapMap->clear();
    VertexIDSet seenVerts;
    // the irreducible map to the transitive map
    while(!overlapQueue.empty())
    {
        EdgeDescOverlapPair edoPair = overlapQueue.top();
        overlapQueue.pop();

        EdgeDesc& edXY = edoPair.first;
        if(seenVerts.count(edXY.pVertex->getID()) == 0)
        {
            seenVerts.insert(edXY.pVertex->getID());
            pOverlapMap->insert(edoPair);
        }
    }
}

// Return a descriptor of the edge describing ovrXY
EdgeDesc SGAlgorithms::overlapToEdgeDesc(Vertex* pY, const Overlap& ovrXY)
{
    EdgeDesc edXY;
    edXY.pVertex = pY;
    edXY.comp = (ovrXY.match.isRC()) ? EC_REVERSE : EC_SAME;
    edXY.dir = ovrXY.match.coord[0].isLeftExtreme() ? ED_ANTISENSE : ED_SENSE; // X -> Y
    return edXY;
}

// Return true if XZ has an overlap
bool SGAlgorithms::hasTransitiveOverlap(const Overlap& ovrXY, const Overlap& ovrYZ)
{
    Match match_yx = ovrXY.match;
    match_yx.swap(); 
    Match match_yz = ovrYZ.match;
    return Match::doMatchesIntersect(match_yx, match_yz);
}

//
EdgeDesc SGAlgorithms::getEdgeDescFromEdge(Edge* pEdge)
{
    return pEdge->getDesc();
}

void SGAlgorithms::printOverlapMap(const EdgeDescOverlapMap& overlapMap)
{
    for(EdgeDescOverlapMap::const_iterator iter = overlapMap.begin(); iter != overlapMap.end(); ++iter)
    {
        std::cout << "  Overlap:" << iter->second << "\n";
    }
}

