///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// StringGraphGenerator - Iteratively construct 
// a local String Graph using the FM-index
//
#include "StringGraphGenerator.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"

//#define DEBUGGENERATE 1

StringGraphGenerator::StringGraphGenerator(const OverlapAlgorithm* pOverlapper, 
                                           const SeqRecord& startRead, 
                                           const SeqRecord& endRead, 
                                           int minOverlap,
                                           EdgeDir startDir,
                                           int maxDistance) : m_pOverlapper(pOverlapper), 
                                                              m_minOverlap(minOverlap), 
                                                              m_pGraph(NULL), 
                                                              m_startDir(startDir), 
                                                              m_maxDistance(maxDistance)
{
    m_pGraph = new StringGraph;

    // Add the start and end vertices to the graph
    m_pStartVertex = addTerminalVertex(startRead);
    m_pEndVertex = addTerminalVertex(endRead);

    // Set the color of the starting node to be UNEXPLORED
    // and the color of the end node to be EXPLORED.
    // This indicates to the subsequent search which vertices
    // should be expanded
    m_pStartVertex->setColor(UNEXPLORED_COLOR);
    m_pEndVertex->setColor(EXPLORED_COLOR);

    // Set up the expansion frontier
    FrontierQueue queue;
    GraphFrontier node;
    node.pVertex = m_pStartVertex;
    node.dir = m_startDir;
    node.distance = m_pStartVertex->getSeqLen();
    queue.push(node);

    buildGraph(queue);

    //m_pGraph->writeDot("local.dot");

    SGDuplicateVisitor dupVisit(true);
    m_pGraph->visit(dupVisit);

    // If the terminal vertices are marked as contained, reset the containment flags so they will not be removed
    resetContainmentFlags(m_pStartVertex);
    resetContainmentFlags(m_pEndVertex);

    SGContainRemoveVisitor containVisit;
    m_pGraph->visit(containVisit);
    //m_pGraph->writeDot("local-final.dot");
}

//
StringGraphGenerator::~StringGraphGenerator()
{
    delete m_pGraph;
    m_pGraph = NULL;

    // m_pStartVertex and m_pEndVertex are deleted by the graph destructor
    // so they do not need to be explicitly freed here.
}

// Build the graph by expanding nodes on the frontier
void StringGraphGenerator::buildGraph(FrontierQueue& queue)
{
    while(!queue.empty())
    {
        if(queue.size() > 200)
            break;

        GraphFrontier node = queue.front();
        queue.pop();
        if(node.pVertex->getColor() == EXPLORED_COLOR)
            continue; // node has been visited already
        
        // Search the FM-index for the current vertex
        SeqRecord record;
        record.id = node.pVertex->getID();
        record.seq = node.pVertex->getSeq().toString();
        
        OverlapBlockList blockList;
        assert(blockList.empty());
        m_pOverlapper->overlapRead(record, m_minOverlap, &blockList);

        // Update the graph and the frontier queue with newly found vertices
        updateGraphAndQueue(node, queue, blockList);
        node.pVertex->setColor(EXPLORED_COLOR);
    }

    m_pGraph->setColors(GC_WHITE);
}

// Search for walks between the start and end vertex
SGWalkVector StringGraphGenerator::searchWalks()
{
    SGWalkVector walks;
    SGSearch::findWalks(m_pStartVertex, m_pEndVertex, m_startDir, m_maxDistance, true, 10000, walks);
    return walks;
}

// 
void StringGraphGenerator::updateGraphAndQueue(GraphFrontier& currNode, FrontierQueue& queue, OverlapBlockList& blockList)
{
    // Partition the block list into containment blocks and extension (valid) blocks
    // We do not add containment edges to the graph so the containments are discarded
    OverlapBlockList containList;
    OverlapBlockList overlapList;

    Vertex* pX = currNode.pVertex;

    //partitionBlockList(pX->getSeqLen(), &blockList, &overlapList, &containList);

    // Process the overlap blocks, adding new vertices and edges where necessary
    for(OverlapBlockList::iterator iter = blockList.begin(); iter != blockList.end(); ++iter)
    {
        if(iter->getEdgeDir() != currNode.dir)
            continue;

        std::string vertexID = iter->toCanonicalID();
        if(vertexID == pX->getID())
            continue; // skip self-edges


        std::string vertexSeq = iter->getFullString(pX->getSeq().toString());
        Overlap o = iter->toOverlap(pX->getID(), vertexID, pX->getSeqLen(), vertexSeq.length());

/*
#if DEBUGGENERATE
        std::cout << "has overlap to: " << vertexID << " len: " << iter->overlapLen << " flags: " << iter->flags << "\n";
        std::cout << "Overlap string: " << iter->getOverlapString(pX->getSeq().toString()) << "\n";
#endif
*/      
        // Check if a vertex with endVertexID exists in the graph
        Vertex* pVertex = m_pGraph->getVertex(vertexID);
        if(pVertex == NULL)
        {

#if DEBUGGENERATE
            std::cout << "Vertex with ID: " << vertexID << " does not exist, creating\n";
            std::cout << "Vertex sequence: " << vertexSeq << "\n";
#endif
            // Generate the new vertex
            vertexSeq = iter->getFullString(pX->getSeq().toString());
            //pVertex = new(m_pGraph->getVertexAllocator()) Vertex(vertexID, vertexSeq);
            pVertex = new Vertex(vertexID, vertexSeq);
            pVertex->setColor(UNEXPLORED_COLOR);
            m_pGraph->addVertex(pVertex);
        }

        // Construct the found edge
        Edge* pXY = SGAlgorithms::createEdgesFromOverlap(m_pGraph, o, true);

        // If the endpoint vertex is unexplored, queue it
        if(pVertex->getColor() == UNEXPLORED_COLOR)
        {
            GraphFrontier node;
            node.pVertex = pVertex;
            node.dir = !pXY->getTwin()->getDir(); // continuation direction
            node.distance = currNode.distance + pXY->getSeqLen();
            if(node.distance <= m_maxDistance)
                queue.push(node);
        }
    }
}

//
Vertex* StringGraphGenerator::addTerminalVertex(const SeqRecord& record)
{
    assert(m_pGraph != NULL);

    // Build the vertex by performing a full-length search for the
    // sequence in the FM-index. We set the ID of the vertex to be the 
    // lowest index in the returned block list.
    OverlapBlockList endBlockList;
    m_pOverlapper->alignReadDuplicate(record, &endBlockList);

    // Search the block list for the exact match to the end read. This must exist
    OverlapBlockList::iterator matchIter = endBlockList.begin();
    while(matchIter != endBlockList.end())
    {
        if(matchIter->numDiff == 0 && !matchIter->flags.isQueryRev())
            break; // this block corresponds to the actual sequence of endRead
    }
    assert(matchIter != endBlockList.end());
    
    // Construct the canonical ID from the matching interval
    std::string endID = matchIter->toCanonicalID();

    Vertex* pVertex = m_pGraph->getVertex(endID);
    if(pVertex == NULL)
    {
      //pVertex = new(m_pGraph->getVertexAllocator()) Vertex(endID, record.seq.toString());
        pVertex = new Vertex(endID, record.seq.toString());
        m_pGraph->addVertex(pVertex);
    }
    return pVertex;
}

//
void StringGraphGenerator::resetContainmentFlags(Vertex* pVertex)
{
    if(!pVertex->isContained())
        return;
    pVertex->setContained(false);
    // Set the containment flag for all the vertices that have containment edges with this vertex
    EdgePtrVec edges = pVertex->getEdges();
    for(size_t i = 0; i < edges.size(); ++i)
    {
        Edge* pEdge = edges[i];
        if(pEdge->getOverlap().isContainment())
            pEdge->getEnd()->setContained(true);
    }
}

