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
#include "RemovalAlgorithm.h"

//
SGAlgorithms::EdgeDescOverlapMap RemovalAlgorithm::computeRequiredOverlaps(const Vertex* pVertex, const Edge* pRemovalEdge, double maxER, int minLength)
{
    // this procedure is only valid for proper overlaps, the edge cannot be a containment
    assert(!pRemovalEdge->getOverlap().isContainment());
    SGAlgorithms::EdgeDescOverlapMap outMap;
    findPotentialOverlaps(pVertex, pRemovalEdge, maxER, minLength, outMap);
    eliminateReachableEdges(pVertex, pRemovalEdge, maxER, minLength, outMap);   
    return outMap;
}

// Find edges of pRemovalEdge->getEnd() that are potentially irreducible edges of pX once pRemovalEdge
// has been eliminated from the graph
void RemovalAlgorithm::findPotentialOverlaps(const Vertex* pX, const Edge* pRemovalEdge, 
                                             double maxER, int minLength, SGAlgorithms::EdgeDescOverlapMap& outMap)
{
    ExploreQueue queue;
    Vertex* pY = pRemovalEdge->getEnd();
    Overlap ovrXY = pRemovalEdge->getOverlap();
    EdgeDesc edXY = pRemovalEdge->getDesc();

    // Get the edges of pY that face in the opposite direction of edge X <-- Y
    // New edges of pX must be in this direction. Edges in the direction of X
    // must have a longer overlap with X than Y does and therefore cannot be transitive
    // wrt to Y.
    EdgeDescSet visitedSet;
    EdgeDir dirY = !pRemovalEdge->getTwin()->getDir();
    enqueueEdges(pY, dirY, ovrXY, edXY, minLength, queue, &visitedSet);

    int extra_extension = 0;
    while(!queue.empty())
    {
        ExploreElement ee = queue.front();
        queue.pop();
        
        // If the overlap between this element and X is valid, add it to the output map
        // and proceed no further. If the overlap is above minLength but the error rate
        // with X is higher than the threshold, do not add it to the map but add its neighbors
        // as they may be valid overlaps to X
        EdgeDesc& edXZ = ee.ed;
        Vertex* pZ = edXZ.pVertex;
        Overlap& ovrXZ = ee.ovr;
        int overlapLen = ovrXZ.getOverlapLength(0);
        
        if(pZ == pX)
            continue;

        if(overlapLen >= minLength)
        {
            double error_rate = SGAlgorithms::calcErrorRate(pX, pZ, ovrXZ);
            if(isErrorRateAcceptable(error_rate, maxER))
            {
                outMap.insert(std::make_pair(edXZ, ovrXZ));
            }
            else
            {
                // Enqueue neighbors of pZ
                EdgeDir dirZ = edXZ.getTransitiveDir();
                enqueueEdges(pZ, dirZ, ovrXZ, edXZ, minLength, queue, &visitedSet);
                ++extra_extension;
            }
        }
    }

    // Remove non-maximal overlaps from the collection
    SGAlgorithms::removeSubmaximalOverlaps(&outMap);
}

// Using the edges of pVertex, eliminate edges from outMap if they are transitive
void RemovalAlgorithm::eliminateReachableEdges(const Vertex* pVertex, const Edge* pRemovalEdge,
                                               double maxER, int minLength, SGAlgorithms::EdgeDescOverlapMap& outMap)
{
    // During the enqueue process, edges may have been added to outMap that are transitive to another
    // edge in outMap. We get rid of these first.
    SGAlgorithms::partitionTransitiveOverlaps(&outMap, NULL, maxER, minLength);

    // Now, eliminate any edges in outMap that are transitive wrt some edge reachable from pVertex
    // We control the depth of search using the shortest overlap in outMap
    // Get the length of the shortest overlap on pX in outMap
    int shortestOverlap = -1;
    for(SGAlgorithms::EdgeDescOverlapMap::iterator iter = outMap.begin();
        iter != outMap.end(); ++iter)
    {
        int currOverlap = iter->second.getOverlapLength(0);
        if(shortestOverlap == -1 || currOverlap < shortestOverlap)
        {
            shortestOverlap = currOverlap;
      //      std::cout << "SHORT: " << iter->second << "\n";
        }
    }
    
    // Avoid loops by only enqueuing vertices that have not been visited before
    EdgeDescSet visitedSet;
    ExploreQueue queue;

    // Enqueue the initial overlaps of pX to the queue if they are longer than the shortest overlap
    EdgeDir dirX = pRemovalEdge->getDir();
    EdgePtrVec edges = pVertex->getEdges(dirX);
    for(size_t i = 0; i < edges.size(); ++i)
    {
        Edge* pEdge = edges[i];

        // Skip the edge that is being removed
        if(pEdge == pRemovalEdge)
            continue;

        EdgeDesc ed = pEdge->getDesc();
        Overlap ovr = pEdge->getOverlap();
        if(ovr.getOverlapLength(0) >= shortestOverlap && !ovr.isContainment())
        {
            queue.push(ExploreElement(ed, ovr));
            visitedSet.insert(ed);
        }
    }

    // The elements in the queue represent the irreducible overlaps of pX and any overlaps
    // that are transitive wrt the irreducible edges. If the edges in outMap are transitive
    // wrt to this edge set, remove them from the map. Any remaining edges are new irreducible edges.
    int num_ext = 0;
    while(!queue.empty())
    {
        if(outMap.empty())
            return;

        ExploreElement ee = queue.front();
        queue.pop();
        
        // If the overlap between this element and X is valid, add it to the output map
        // and proceed no further. If the overlap is above minLength but the error rate
        // with X is higher than the threshold, do not add it to the map but add its neighbors
        // as they may be valid overlaps to X
        EdgeDesc& edXY = ee.ed;
        Vertex* pY = edXY.pVertex;
        Overlap& ovrXY = ee.ovr;
        //std::cout << "OVRXY: " << ovrXY << "\n";
        int overlapLen = ovrXY.getOverlapLength(0);
        assert(overlapLen >= shortestOverlap);
        (void)overlapLen;

        // Check all the overlaps in the map to see if they are transitive wrt pY
        SGAlgorithms::EdgeDescOverlapMap::iterator iter = outMap.begin();
        while(iter != outMap.end())
        {
            bool eliminate = false;
            const EdgeDesc& edXZ = iter->first;
            const Overlap& ovrXZ = iter->second;
            
            if(ovrXY.getOverlapLength(0) >= ovrXZ.getOverlapLength(0))
            {
                if(edXY.pVertex == edXZ.pVertex)
                    eliminate = true;
                else
                    eliminate = SGAlgorithms::isOverlapTransitive(edXY.pVertex, edXZ.pVertex, ovrXY, ovrXZ, maxER, minLength);
            }

            if(eliminate)
            {
                outMap.erase(iter++);
            }
            else
                ++iter;
        }

        // If any overlaps remain in the map, enqueue the neighbors of this vertex
        // that have an overlap with X that is at least as long as shortest overlap
        if(!outMap.empty())
        {
            num_ext++;
            enqueueEdges(pY, edXY.getTransitiveDir(), ovrXY, edXY, shortestOverlap, queue, &visitedSet);
        }
    }
}

// Add the edges of pY in direction dirY to the explore queue if they have a valid overlap with X
void RemovalAlgorithm::enqueueEdges(const Vertex* pY, EdgeDir dirY, const Overlap& ovrXY, const EdgeDesc& /*edXY*/, int minOverlap, 
                                    ExploreQueue& outQueue, EdgeDescSet* pSeenSet)
{
    EdgePtrVec edges = pY->getEdges(dirY);
    for(size_t i = 0; i < edges.size(); ++i)
    {
        Edge* pEdge = edges[i];
        Vertex* pZ = pEdge->getEnd();

        // Compute the edgeDesc and overlap on pX for this edge
        Overlap ovrYZ = pEdge->getOverlap();

        if(!ovrYZ.isContainment() && SGAlgorithms::hasTransitiveOverlap(ovrXY, ovrYZ))
        {
            Overlap ovrXZ = SGAlgorithms::inferTransitiveOverlap(ovrXY, ovrYZ);
            EdgeDesc edXZ = SGAlgorithms::overlapToEdgeDesc(pZ, ovrXZ);

            if((pSeenSet == NULL || pSeenSet->count(edXZ) == 0) && ovrXZ.getOverlapLength(0) >= minOverlap)
            {
                outQueue.push(ExploreElement(edXZ, ovrXZ));
                if(pSeenSet != NULL)
                    pSeenSet->insert(edXZ);
            }
        }
    }
}

