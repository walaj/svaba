//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Vertex - Generic vertex class for bigraph
//
#include "Vertex.h"
#include "Edge.h"
#include <algorithm>

Vertex::~Vertex()
{
    EdgePtrVecIter iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
    {
        delete *iter;
        *iter = NULL;
    }
}

// Merging two string vertices has two parts
// First, the sequence of the vertex is extended
// by the the content of the edge label
// Then, all the edges that are pointing to this node
// must be updated to contain the extension of the vertex
void Vertex::merge(Edge* pEdge)
{
    Edge* pTwin = pEdge->getTwin();
    //std::cout << "Adding label to " << getID() << " str: " << pSE->getLabel() << "\n";

    // Merge the sequence
    DNAEncodedString label = pEdge->getLabel();
    size_t label_len = label.length();
    pEdge->updateSeqLen(m_seq.length() + label_len);
    bool prepend = false;

    if(pEdge->getDir() == ED_SENSE)
    {
        m_seq.append(label);
    }
    else
    {
        label.append(m_seq);
        std::swap(m_seq, label);
        prepend = true;
    }

    // Update the coverage value of the vertex
    m_coverage += pEdge->getEnd()->getCoverage();

    pEdge->extendMatch(label_len);
    pTwin->extendMatchFullLength();

    // All the SeqCoords for the edges must have their seqlen field updated
    // Also, if we prepended sequence to this edge, all the matches in the 
    // SENSE direction must have their coordinates offset
    size_t newLen = m_seq.length();
    for(EdgePtrVecIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        Edge* pUpdateEdge = *iter;
        pUpdateEdge->updateSeqLen(newLen);
        if(prepend && pUpdateEdge->getDir() == ED_SENSE && pEdge != pUpdateEdge)
            pUpdateEdge->offsetMatch(label_len);
    }

#ifdef VALIDATE
    VALIDATION_WARNING("Vertex::merge")
    validate();
#endif

}

void Vertex::validate() const
{
    for(EdgePtrVecConstIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        (*iter)->validate();
        /*
        std::string label = pSE->getLabel();
        StringVertex* pEnd = SV_CAST(pSE->getEnd());

        std::string vertSeq = pEnd->getSeq();
        EdgeDir suffixDir = (pSE->getComp() == EC_SAME) ? pSE->getDir() : !pSE->getDir();
        std::string vertSuffix = (suffixDir == ED_SENSE) ? vertSeq.substr(vertSeq.length() - label.length()) :
                                                                vertSeq.substr(0, label.length());

        if(pSE->getComp() == EC_REVERSE)
            vertSuffix = reverseComplement(vertSuffix);
        if(vertSuffix != label)
        {
            std::cerr << "Warning edge label " << label << " does not match vertex suffix " << vertSuffix << "\n";
        }
        */
    }
}

//
void Vertex::sortAdjListByID()
{
    EdgeIDComp comp;
    std::sort(m_edges.begin(), m_edges.end(), comp);
}

void Vertex::sortAdjListByLen()
{
    EdgeLenComp comp;
    std::sort(m_edges.begin(), m_edges.end(), comp);
}

// Mark duplicate edges with dupColor
// Returns true if a duplicate has been marked
bool Vertex::markDuplicateEdges(GraphColor dupColor)
{
    // Sort the edge lists by length
    sortAdjListByLen();
    bool hasDup = false;
    hasDup = markDuplicateEdges(ED_SENSE, dupColor) || hasDup;
    hasDup = markDuplicateEdges(ED_ANTISENSE, dupColor) || hasDup;
    return hasDup;
}

// Mark duplicate edges in the specified direction
bool Vertex::markDuplicateEdges(EdgeDir dir, GraphColor dupColor)
{
    for(EdgePtrVecIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        Edge* pEdge = *iter;
        if(pEdge->getDir() == dir)
        {
            Vertex* pY = pEdge->getEnd();
            if(pY->getColor() == GC_BLACK)
            {
                //std::cerr << getID() << " has a duplicate edge to " << pEdge->getEndID() << " in direction " << dir << "\n";

                // This vertex is the endpoint of some other (potentially longer) edge
                // Delete the edge
                Edge* pTwin = pEdge->getTwin();
                pTwin->setColor(dupColor);
                pEdge->setColor(dupColor);
            }
            else
            {
                assert(pY->getColor() == GC_WHITE);
                pY->setColor(GC_BLACK);
            }
        }
    }

    // Reset vertex colors
    for(EdgePtrVecIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
        (*iter)->getEnd()->setColor(GC_WHITE);
    return true;
}

// Get a multioverlap object representing the overlaps for this vertex
MultiOverlap Vertex::getMultiOverlap() const
{
    MultiOverlap mo(getID(), getSeq().toString());
    for(size_t i = 0; i < m_edges.size(); ++i)
    {
        Edge* pEdge = m_edges[i];
        mo.add(pEdge->getEnd()->getSeq().toString(), pEdge->getOverlap());
    }
    return mo;
}

// Add an edge
void Vertex::addEdge(Edge* ep)
{
    assert(ep->getStart() == this);

#ifdef VALIDATE
    for(EdgePtrVecConstIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        if((*iter)->getEndID() == ep->getEndID())
        {
            std::cout << "Attempted to add duplicate edge with ID: " << ep->getEndID() 
                        << " to vertex: " << ep->getStartID() << "\n";
            std::cout << "Added in desc: " << ep->getDesc() << " curr desc: " << (*iter)->getDesc() << "\n";
            //assert(false);
        }
    }
#endif
    m_edges.push_back(ep);

}

// Remove an edge from the edge list of the vertex
// and free the storage. This does not remove the twin
void Vertex::deleteEdge(Edge* pEdge)
{

    removeEdge(pEdge);
    delete pEdge;
    pEdge = NULL;
}

// Remove an edge but do not destroy it
void Vertex::removeEdge(Edge* pEdge)
{

    EdgePtrVecIter iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        if(*iter == pEdge)
            break;
        ++iter;
    }
    assert(iter != m_edges.end());
    m_edges.erase(iter);    
}

//
void Vertex::removeEdge(const EdgeDesc& ed)
{

    EdgePtrVecIter iter = findEdge(ed);
    if(iter == m_edges.end())
    {
        std::cout << "EDGE NOT FOUND: " << ed << "\n";
    }
    assert(iter != m_edges.end());
    m_edges.erase(iter);
}

// Delete all the edges, and their twins, from this vertex
void Vertex::deleteEdges()
{

    EdgePtrVecIter iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
    {
        Edge* pEdge = *iter;
        Edge* pTwin = pEdge->getTwin();
        Vertex* pPartner = pEdge->getEnd();
        pPartner->removeEdge(pTwin);
        delete pEdge;
        pEdge = NULL;
        delete pTwin;
        pTwin = NULL;
        *iter = NULL;
    }
    m_edges.clear();
}

// Delete edges that are marked
// This only deletes the edge and not its twin
int Vertex::sweepEdges(GraphColor c)
{

    int numRemoved = 0;
    EdgePtrVecIter iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        Edge* pEdge = *iter;
        if(pEdge->getColor() == c)
        {
            delete pEdge;
            pEdge = NULL;
            iter = m_edges.erase(iter);
            ++numRemoved;
        }
        else
            ++iter;
    }
    return numRemoved;
}

// Return the iterator to the edge matching edgedesc
EdgePtrVecIter Vertex::findEdge(const EdgeDesc& ed)
{
    for(EdgePtrVecIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        if((*iter)->getDesc() == ed)
            return iter;
    }
    return m_edges.end();
}

//
EdgePtrVecConstIter Vertex::findEdge(const EdgeDesc& ed) const
{
    for(EdgePtrVecConstIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        if((*iter)->getDesc() == ed)
            return iter;
    }
    return m_edges.end();
}

// Check for the presence of an edge
bool Vertex::hasEdge(Edge* pEdge) const
{
    return hasEdge(pEdge->getDesc());
}

//
bool Vertex::hasEdge(const EdgeDesc& ed) const
{
    return findEdge(ed) != m_edges.end();
}

//
bool Vertex::hasEdgeTo(const Vertex* pY) const
{
    assert(pY != NULL);
    EdgePtrVecConstIter iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
        if((*iter)->getEnd() == pY)
            return true;
    return false;
}

// Return the edge matching the descriptions
Edge* Vertex::getEdge(const EdgeDesc& ed)
{
     EdgePtrVecIter i = findEdge(ed);
     assert(i != m_edges.end());
     return *i;
}

// Find edges to the specified vertex
EdgePtrVec Vertex::findEdgesTo(VertexID id)
{
    EdgePtrVecConstIter iter = m_edges.begin();
    EdgePtrVec outEdges;
    for(; iter != m_edges.end(); ++iter)
    {
        if((*iter)->getEndID() == id)
            outEdges.push_back(*iter);
    }
    return outEdges;
}

// Returns the edge with the longest overlap length
// in direction dir
// Returns NULL if the vertex has no edges
Edge* Vertex::getLongestOverlapEdge(EdgeDir dir) const
{
    Edge* pOut = NULL;
    int maxOL = 0;
    EdgePtrVecConstIter iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
    {
        if((*iter)->getDir() != dir)
            continue;

        int currOL = (*iter)->getMatchLength();
        if(currOL > maxOL)
        {
            pOut = *iter;
            maxOL = currOL;
        }
    }
    return pOut;
}


//
// Get the edges in a particular direction
// This preserves the ordering of the edges
//
EdgePtrVec Vertex::getEdges(EdgeDir dir) const
{
    EdgePtrVecConstIter iter = m_edges.begin();
    EdgePtrVec outEdges;
    for(; iter != m_edges.end(); ++iter)
    {
        if((*iter)->getDir() == dir)
            outEdges.push_back(*iter);
    }
    return outEdges;
}


// Get the edges
EdgePtrVec Vertex::getEdges() const
{
    EdgePtrVec outEdges(m_edges.begin(), m_edges.end());
    return outEdges;    
}

void Vertex::setEdgeColors(GraphColor c) 
{ 
    for(EdgePtrVecIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
        (*iter)->setColor(c);
}

// Count the edges
// This function is not necessarily constant time
size_t Vertex::countEdges() const
{ 
    return m_edges.size(); 
}

//
size_t Vertex::countEdges(EdgeDir dir)
{
    EdgePtrVec ev = getEdges(dir);
    return ev.size();
}

// Calculate the difference in overlap lengths between
// the longest and second longest edge
int Vertex::getOverlapLengthDiff(EdgeDir dir) const
{
    int longest_len = 0;
    int second_longest_len = 0;
    EdgePtrVecConstIter iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
    {
        if((*iter)->getDir() != dir)
            continue;

        int currOL = (*iter)->getMatchLength();
        if(currOL > longest_len)
        {
            second_longest_len = longest_len;
            longest_len = currOL;
        }
        else if(currOL > second_longest_len)
        {
            second_longest_len = currOL;
        }
    }
    return longest_len - second_longest_len;
}


// Return the amount of memory this vertex is using, in bytes
size_t Vertex::getMemSize() const
{
    return sizeof(*this) + (m_edges.size() * sizeof(Edge*)) + m_seq.getMemSize();
}


// Output edges in graphviz format
void Vertex::writeEdges(std::ostream& out, int dotFlags) const
{
    EdgePtrVecConstIter iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
    {
        if(dotFlags & DF_UNDIRECTED)
        {
            if((*iter)->getStartID() < (*iter)->getEndID())
            {
                out << "\"" << (*iter)->getStartID() << "\" -- \"" << (*iter)->getEndID() << "\"";
            }
        }
        else
        {
            out << "\"" << (*iter)->getStartID() << "\" -> \"" << (*iter)->getEndID();
            std::string color = ((*iter)->getDir() == ED_SENSE) ? "black" : "red";
            std::string label = ((*iter)->getComp() == EC_SAME) ? "S" : "F";
            out << "\" [color=\"" << color << "\" ";
            out << "label=\"" << (*iter)->getMatchLength() << "\"";
            out << "];";
        }
        out << "\n";
    }
}

bool EdgeIDComp::operator()(const Edge* pA, const Edge* pB) 
{
       return pA->getEndID() < pB->getEndID();
}

// Compare string edge points by length
bool EdgeLenComp::operator()(const Edge* pA, const Edge* pB)
{
    return pA->getSeqLen() < pB->getSeqLen();
}

