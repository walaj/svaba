//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldWalk - a walk through a scaffold graph
//
#include "ScaffoldWalk.h"
#include "ScaffoldVertex.h"

ScaffoldWalk::ScaffoldWalk(ScaffoldVertex* pStartVertex) : m_pStartVertex(pStartVertex)
{
}

void ScaffoldWalk::addEdge(ScaffoldEdge* pEdge)
{
    m_edges.push_back(pEdge);
}

//
int64_t ScaffoldWalk::getGapSum() const
{
    int64_t sum = 0;
    ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
        sum += (*iter)->getDistance();

    return sum;
}

//
int64_t ScaffoldWalk::getContigLengthSum() const
{
    int64_t sum = m_pStartVertex->getSeqLen();
    ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
        sum += (*iter)->getEnd()->getSeqLen();
    return sum;
}

//
ScaffoldVertexPtrVector ScaffoldWalk::getVertices() const
{
    ScaffoldVertexPtrVector outVertices;
    outVertices.push_back(m_pStartVertex);
    ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
        outVertices.push_back((*iter)->getEnd());
    return outVertices;
}

//
ScaffoldEdgePtrVector ScaffoldWalk::getEdges() const
{
    return m_edges;
}

//
int ScaffoldWalk::findVertex(ScaffoldVertex* pVertex) const
{
    if(pVertex == m_pStartVertex)
        return 0;

    int idx = 1;
    for(ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
                                              iter != m_edges.end();
                                              ++iter) {
        if(pVertex == (*iter)->getEnd())
            return idx;
        ++idx;
    }
    return -1;
}

//
EdgeComp ScaffoldWalk::findOrientation(ScaffoldVertex* pVertex) const
{
    if(pVertex == m_pStartVertex)
        return EC_SAME;
    EdgeComp out = EC_SAME;
    for(ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
                                              iter != m_edges.end();
                                              ++iter) {
        if((*iter)->getComp() == EC_REVERSE)
            out = !out;

        if(pVertex == (*iter)->getEnd())
            return out;
    }

    std::cerr << "Error, ScaffoldWalk::findOrientation requires the vertex to be in the walk.\n";
    exit(EXIT_FAILURE);
    return EC_SAME;
}

// returns the sum of the gaps plus contig lengths from
// the starting vertex to pVertex.
// Precondition: pVertex is required to be in the walk
int ScaffoldWalk::getDistanceToVertex(ScaffoldVertex* pVertex) const
{
    int distance = 0;
    ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
    {
        distance += (*iter)->getDistance();
        if(pVertex == (*iter)->getEnd())
            return distance;
        distance += (*iter)->getEnd()->getSeqLen();
    }
    assert(false);
    return -1;
}


//
ScaffoldVertex* ScaffoldWalk::getStartVertex() const
{
    return m_pStartVertex;
}

//
ScaffoldVertex* ScaffoldWalk::getLastVertex() const
{
    if(m_edges.empty())
    {
        return m_pStartVertex;
    }
    else
    {
        return m_edges.back()->getEnd();
    }
}

// Print the edges in dot format
void ScaffoldWalk::printDot(std::ostream& out) const
{
    ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
    {
        ScaffoldEdge* pEdge = *iter;
        out << "\"" << pEdge->getStartID() << "\" -> \"" << pEdge->getEndID() << "\" [color=\"blue\"];\n";
    }
}

void ScaffoldWalk::print() const
{
    for(ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
                                              iter != m_edges.end();
                                              ++iter)
    {
        std::cout << *(*iter) << " ";
    }
    std::cout << "\n";
}
