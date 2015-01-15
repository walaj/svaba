//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldAlgorithms - Collection of algorithms
// for computing linear scaffolds from a scaffold graph
//
#include "ScaffoldAlgorithms.h"

// Build the scaffolds.
// The algorithm first finds the connected components of the graph then
// all the terminal vertices of each component. It then searches the
// components starting from the terminals for a walk that maximizes
// the size of the contigs in the primary scaffold.
void ScaffoldAlgorithms::makeScaffolds(ScaffoldGraph* pGraph)
{
    // Set every edge to be colored black
    // Edges that are kept in the scaffold will be colored white
    pGraph->setEdgeColors(GC_BLACK);

    ScaffoldConnectedComponents connectedComponents;
    ScaffoldAlgorithms::connectedComponents(pGraph, connectedComponents);
 
    // Iterate over each connected component and linearize the scaffolds   
    for(size_t i = 0; i < connectedComponents.size(); ++i)
    {
        // Discover the terminal vertices in the component
        ScaffoldVertexPtrVector& component = connectedComponents[i];

        if(component.size() == 1) 
            continue; //nothing to do

        //
        ScaffoldVertexPtrVector terminalVertices;
        computeTerminalsForConnectedComponent(component, terminalVertices);

        if(terminalVertices.empty())
        {
            std::cerr << "Warning: scaffold component of size " << component.size() << 
                         " does not have a terminal vertex. Skipping\n";
            continue;
        }
        // Construct a scaffold layout for each terminal vertex of the component
        // We select the layout that contains the greatest amount of sequence as
        // the initial layout of the scaffold
        size_t bestLayoutBases = 0;
        ScaffoldWalk bestWalk(NULL);

        for(size_t j = 0; j < terminalVertices.size(); j++)
        {
            ScaffoldWalkVector layouts;
            computeLayout(terminalVertices[j], layouts);

            for(size_t k = 0; k < layouts.size(); ++k)
            {
                size_t layoutSum = layouts[k].getContigLengthSum();

                if(layoutSum > bestLayoutBases)
                {
                    bestLayoutBases = layoutSum;
                    bestWalk = layouts[k];
                }
            }
        }

        // Ensure a valid layout is chosen
        assert(bestLayoutBases > 0 && bestWalk.getStartVertex() != NULL);

        // Remove every edge that is not a part of the chosen walk
        // Since every edge in the graph was initially colored black,
        // we color the kept edges white then remove every black edge
        // in a single pass later
        ScaffoldEdgePtrVector keptEdges = bestWalk.getEdges();
        for(size_t j = 0; j < keptEdges.size(); ++j)
        {
            keptEdges[j]->setColor(GC_WHITE);
            keptEdges[j]->getTwin()->setColor(GC_WHITE);
        }
    }

    // Remove all edges that aren't on a chosen layouts
    pGraph->deleteEdgesByColor(GC_BLACK);
}

// Remove internal cycles from the connected components of the graph
void ScaffoldAlgorithms::removeInternalCycles(ScaffoldGraph* pGraph)
{
    bool done = false;
    while(!done)
    {
        bool cycleFound = false;
        pGraph->setVertexColors(GC_WHITE);

        // Compute the connected components of the graph
        ScaffoldConnectedComponents connectedComponents;
        ScaffoldAlgorithms::connectedComponents(pGraph, connectedComponents);

        // Check each CC for a cycle
        for(size_t i = 0; i < connectedComponents.size(); ++i)
        {
            // Discover the terminal vertices in the component
            ScaffoldVertexPtrVector& component = connectedComponents[i];

            if(component.size() == 1) 
                continue; //nothing to do

            // Compute the terminal vertices for the CC
            // We start the cycle search from these positions
            ScaffoldVertexPtrVector terminalVertices;
            computeTerminalsForConnectedComponent(component, terminalVertices);     

            for(size_t j = 0; j < terminalVertices.size(); j++)
            {
                ScaffoldEdge* pBackEdge = checkForInternalCycle(terminalVertices[j]);
                if(pBackEdge != NULL)
                {
                    std::cout << "Internal cycle found between: " << pBackEdge->getStartID() << " and " << pBackEdge->getEndID() << "\n";

                    // Mark the endpoints of the cycle as repeats
                    pBackEdge->getStart()->setClassification(SVC_REPEAT);
                    pBackEdge->getEnd()->setClassification(SVC_REPEAT);
                    cycleFound = true;
                    break;
                }
            }
        }

        if(cycleFound)
        {
            pGraph->deleteVertices(SVC_REPEAT);
        }
        else
        {
            done = true;
        }
    }
}

// Destroy simple cycles in the graph
void ScaffoldAlgorithms::destroyStrictCycles(ScaffoldGraph* pGraph, std::string out_filename)
{
    std::ofstream cycle_writer(out_filename.c_str());

    bool done = false;
    while(!done)
    {
        done = true;
        pGraph->setVertexColors(GC_WHITE);
        pGraph->setEdgeColors(GC_WHITE);

        // Compute the connected components of the graph
        ScaffoldConnectedComponents connectedComponents;
        ScaffoldAlgorithms::connectedComponents(pGraph, connectedComponents);

        // Check each CC for a cycle. 
        for(size_t i = 0; i < connectedComponents.size(); ++i)
        {
            ScaffoldVertexPtrVector& component = connectedComponents[i];

            // Trivially not a cycle
            if(component.size() == 1)
                continue;
            
            assert(!component.empty());

            // The strict cycle check will find the first cycle in the connected component reachable from this vertex
            ScaffoldVertex* test_vertex = component.front();
            ScaffoldVertexVector cycle_vertices = checkForStrictCycle(test_vertex);
            if(!cycle_vertices.empty())
            {
                //std::cout << "\tcycle starting at " << cycle_vertices[0]->getID() << " length " << cycle_vertices.size() << " found\n";

                // Delete the edges for this vertex
                for(size_t j = 0; j < cycle_vertices.size(); ++j)
                {
                    cycle_vertices[j]->deleteEdgesAndTwins();
                    cycle_writer << cycle_vertices[j]->getID() << "\n";
                }
                
                // iterate the cycle detection
                done = false;
            }
        }
    }
}

// Compute the connected components of the graph
void ScaffoldAlgorithms::connectedComponents(ScaffoldGraph* pGraph, ScaffoldConnectedComponents& outComponents)
{
    ScaffoldVertexPtrVector allVertices = pGraph->getAllVertices();
    ScaffoldSearchTree::connectedComponents(allVertices, outComponents);
}

// Compute the terminal vertices of a connected component
void ScaffoldAlgorithms::computeTerminalsForConnectedComponent(const ScaffoldVertexPtrVector& component, 
                                                               ScaffoldVertexPtrVector& terminals)
{
    for(size_t i = 0; i < component.size(); ++i)
    {
        ScaffoldVertex* pVertex = component[i];
        size_t asCount = pVertex->getEdges(ED_ANTISENSE).size();
        size_t sCount = pVertex->getEdges(ED_SENSE).size();

        if(asCount == 0 || sCount == 0)
            terminals.push_back(pVertex);
    }
}


// Compute a layout of a component of the scaffold graph starting
// from a particular vertex. Writes a path each reachable terminal vertex
// into outWalks.
// Precondition: The connected component containing pVertex does not contain
// a cycle.
void ScaffoldAlgorithms::computeLayout(ScaffoldVertex* pStartVertex, ScaffoldWalkVector& outWalks)
{
    // Detect the direction to begin constructing the component
    size_t asCount = pStartVertex->getEdges(ED_ANTISENSE).size();
    size_t sCount = pStartVertex->getEdges(ED_SENSE).size();
    assert(asCount == 0 || sCount == 0);
    if(asCount == 0 && sCount == 0)
        return; // nothing to do
    EdgeDir dir = (asCount > 0) ? ED_ANTISENSE : ED_SENSE;

    LayoutDistanceMap distanceMap;
    LayoutEdgeMap edgeMap;
    LayoutTerminalSet terminalSet;

    distanceMap[pStartVertex] = 0;
    edgeMap[pStartVertex] = NULL;

    // Enqueue the edges of this node to the search structure
    LayoutQueue visitQueue;
    ScaffoldEdgePtrVector startEdges = pStartVertex->getEdges(dir);
    for(size_t i = 0; i < startEdges.size(); ++i)
    {
        ScaffoldVertex* pZ = startEdges[i]->getEnd();
        LayoutNode initial;
        initial.pEdge = startEdges[i];
        initial.distance = startEdges[i]->getDistance();

        distanceMap[pZ] = initial.distance;
        edgeMap[pZ] = startEdges[i];

        visitQueue.push(initial);
    }

    while(!visitQueue.empty())
    {
        LayoutNode node = visitQueue.front();
        visitQueue.pop();
        ScaffoldEdge* pXY = node.pEdge;
        ScaffoldVertex* pY = pXY->getEnd();

        // Get the edges of the endpoint and add them to the queue
        // to explore they are closer than the current estimate
        EdgeDir yDir = !pXY->getTwin()->getDir();
        ScaffoldEdgePtrVector yEdges = pY->getEdges(yDir);
        if(yEdges.empty())
        {
            terminalSet.insert(pY);
        }

        // Add the edges of y to the queue if a new or shorter path has been found
        for(size_t i = 0; i < yEdges.size(); ++i)
        {
            ScaffoldEdge* pYZ = yEdges[i];
            ScaffoldVertex* pZ = pYZ->getEnd();

            int zDistance = node.distance + pYZ->getDistance();

            // Check if the gap distance to z is less than the current
            // If so, update the edge to z to be the current edge and re-queue
            // the edge
            LayoutDistanceMap::iterator distIter = distanceMap.find(pZ);
            if(distIter == distanceMap.end() || distIter->second > zDistance)
            {
                distanceMap[pZ] = zDistance;
                edgeMap[pZ] = pYZ;
                LayoutNode updateNode;
                updateNode.pEdge = pYZ;
                updateNode.distance = zDistance;
                visitQueue.push(updateNode);
            }
        }
    }

    // Construct paths from the starting node to the terminal vertices
    LayoutTerminalSet::iterator terminalIter = terminalSet.begin();
    for(; terminalIter != terminalSet.end(); ++terminalIter)
    {
        ScaffoldEdgePtrVector reverseEdges;
        ScaffoldVertex* pCurrent = *terminalIter;
        while(pCurrent != pStartVertex)
        {
            ScaffoldEdge* pEdge = edgeMap[pCurrent];
            reverseEdges.push_back(pEdge);
            assert(pEdge != NULL);
            pCurrent = pEdge->getStart();
        }

        // Build the real walk
        ScaffoldWalk walk(pStartVertex);
        ScaffoldEdgePtrVector::reverse_iterator rIter = reverseEdges.rbegin();
        for(; rIter != reverseEdges.rend(); ++rIter)
            walk.addEdge(*rIter);
        outWalks.push_back(walk);
    }
}

// Check for internal cycles in the connected component reachable from pVertex.
// If a cycle is found, the back edge pointer making the cycle is returned. Otherwise
// NULL is returned.
ScaffoldEdge* ScaffoldAlgorithms::checkForInternalCycle(ScaffoldVertex* pVertex)
{
    // Detect the direction to search for cycle
    size_t asCount = pVertex->getEdges(ED_ANTISENSE).size();
    size_t sCount = pVertex->getEdges(ED_SENSE).size();
    assert(asCount == 0 || sCount == 0);
    if(asCount == 0 && sCount == 0)
        return NULL; // nothing to do
    EdgeDir dir = (asCount > 0) ? ED_ANTISENSE : ED_SENSE;

    assert(pVertex->getColor() == GC_WHITE);
    LayoutEdgeMap predMap;
    predMap[pVertex] = NULL;
    ScaffoldEdge* pBackEdge = _cycleDFS(pVertex, dir, predMap);

    // reset colors
    for(LayoutEdgeMap::iterator iter = predMap.begin(); iter != predMap.end(); ++iter)
        iter->first->setColor(GC_WHITE);
    return pBackEdge;
}

// Recursive function for finding a cycle using a directional. If a cycle is found, a pointer to the back edge
// is returned. If NULL is returned, there is no cycle.
ScaffoldEdge* ScaffoldAlgorithms::_cycleDFS(ScaffoldVertex* pVertex, EdgeDir dir, LayoutEdgeMap& predMap)
{
    pVertex->setColor(GC_GRAY);
    ScaffoldEdgePtrVector edges = pVertex->getEdges(dir);
    for(size_t i = 0; i < edges.size(); ++i)
    {
        ScaffoldEdge* pEdge = edges[i];
        ScaffoldVertex* pEnd = pEdge->getEnd();
        if(pEnd->getColor() == GC_GRAY)
        {
            return pEdge;
        }

        if(pEnd->getColor() == GC_WHITE)
        {
            predMap[pEnd] = pEdge;
            EdgeDir continueDir = !pEdge->getTwin()->getDir();
            ScaffoldEdge* pBackEdge = _cycleDFS(pEnd, continueDir, predMap);
            if(pBackEdge != NULL)
                return pBackEdge;
        }
    }
    pVertex->setColor(GC_BLACK);
    return NULL;
}

// Check for strict cycles in the connected component reachable from pVertex.
// If a cycle is found, the back edge pointer making the cycle is returned. Otherwise
// NULL is returned.
ScaffoldVertexVector ScaffoldAlgorithms::checkForStrictCycle(ScaffoldVertex* pVertex)
{
    assert(pVertex->getColor() == GC_WHITE);
    
    LayoutEdgeMap predMap;
    predMap[pVertex] = NULL;
    ScaffoldEdge* pCycleEdge = _cycleDFSBothDir(pVertex, predMap);

    ScaffoldVertexVector out;
    if(pCycleEdge != NULL)
    {
        // Find all vertices contained in this cycle
        ScaffoldEdge* current = pCycleEdge;
        while(current != NULL)
        {
            out.push_back(current->getStart());
            current = predMap[current->getStart()];
        }
    }

    // reset colors
    for(LayoutEdgeMap::iterator iter = predMap.begin(); iter != predMap.end(); ++iter)
    {
        iter->first->setColor(GC_WHITE);

        if(iter->second != NULL)
        {
            iter->second->setColor(GC_WHITE);
            iter->second->getTwin()->setColor(GC_WHITE);
        }
    }

    return out;
}

// Recursive function for finding a cycle by searching both directions from a vertex. 
// If a cycle is found, a pointer to the back edge is returned. If NULL is returned, there is no cycle.
ScaffoldEdge* ScaffoldAlgorithms::_cycleDFSBothDir(ScaffoldVertex* pVertex, LayoutEdgeMap& predMap)
{
    pVertex->setColor(GC_GRAY);
    ScaffoldEdgePtrVector edges = pVertex->getEdges();
    for(size_t i = 0; i < edges.size(); ++i)
    {
        ScaffoldEdge* pEdge = edges[i];

        // If an edge's twin has already been used, it is
        // marked red. skip these edges
        if(pEdge->getColor() == GC_RED)
            continue;
        
        // Check if this is a back edge
        ScaffoldVertex* pEnd = pEdge->getEnd();
        if(pEnd->getColor() == GC_GRAY)
            return pEdge;

        if(pEnd->getColor() == GC_WHITE)
        {
            pEdge->setColor(GC_RED);
            pEdge->getTwin()->setColor(GC_RED);
            predMap[pEnd] = pEdge;
            
            // Extend the DFS to the next node
            ScaffoldEdge* pBackEdge = _cycleDFSBothDir(pEnd, predMap);
            if(pBackEdge != NULL)
                return pBackEdge;
        }
    }
    pVertex->setColor(GC_BLACK);
    return NULL;
}
