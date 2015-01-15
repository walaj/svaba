//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGWalk - Data structure holding a walk through
// the graph
//
#include "SGWalk.h"

//
SGWalk::SGWalk(Vertex* pStartVertex, bool bIndexWalk) : m_pStartVertex(pStartVertex), m_extensionDistance(0), m_extensionFinished(false)
{
    if(bIndexWalk)
        m_pWalkIndex = new WalkIndex;
    else
        m_pWalkIndex = NULL;
}

//
SGWalk::~SGWalk()
{
    if(m_pWalkIndex != NULL)
    {
        delete m_pWalkIndex;
        m_pWalkIndex = NULL;
    }
}

//
SGWalk::SGWalk(const SGWalk& other)
{
    m_pStartVertex = other.m_pStartVertex;
    m_edges = other.m_edges;
    m_extensionDistance = other.m_extensionDistance;
    m_extensionFinished = other.m_extensionFinished;
    m_pWalkIndex = NULL;

    if(other.m_pWalkIndex != NULL)
        m_pWalkIndex = new WalkIndex(*other.m_pWalkIndex);
}

// Construct the walk structure from a vector of edges
SGWalk::SGWalk(const EdgePtrVec& edgeVec, bool bIndexWalk) : m_extensionDistance(0), m_extensionFinished(false)

{
    assert(!edgeVec.empty());

    if(bIndexWalk)
        m_pWalkIndex = new WalkIndex;
    else
        m_pWalkIndex = NULL;

    // The start vector is the start vertex of the first edge
    Edge* first = edgeVec.front();
    m_pStartVertex = first->getStart();

    for(EdgePtrVec::const_iterator iter = edgeVec.begin();
                                   iter != edgeVec.end();
                                   ++iter)
    {
        addEdge(*iter);
    }
}


//
SGWalk& SGWalk::operator=(const SGWalk& other)
{
    if(&other == this)
        return *this; // self-assign
   
    m_pStartVertex = other.m_pStartVertex;
    m_edges = other.m_edges;
    m_extensionDistance = other.m_extensionDistance;
    
    if(m_pWalkIndex != NULL)
        delete m_pWalkIndex;
    
    if(other.m_pWalkIndex != NULL)
        m_pWalkIndex = new WalkIndex(*other.m_pWalkIndex);
        
    return *this;
}

//
void SGWalk::addEdge(Edge* pEdge)
{
    m_edges.push_back(pEdge);
    m_extensionDistance += pEdge->getSeqLen();

    // Add the vertex ID to the index if necessary
    if(m_pWalkIndex != NULL)
        m_pWalkIndex->insert(m_edges.back()->getEndID());
}

//
void SGWalk::setFinished(bool b)
{
    m_extensionFinished = b;
}

bool SGWalk::isFinished() const
{
    return m_extensionFinished;
}

//
void SGWalk::popLast()
{
    m_edges.pop_back();
}

// 
Vertex* SGWalk::getStartVertex() const
{
    return m_pStartVertex;
}

//
size_t SGWalk::getNumVertices() const
{
    return getNumEdges() + 1;
}

//
size_t SGWalk::getNumEdges() const
{
    return m_edges.size();
}

bool SGWalk::isIndexed() const
{
    return m_pWalkIndex != NULL;
}

//
bool SGWalk::containsVertex(const VertexID& id) const
{
    if(m_pWalkIndex == NULL)
        assert(false);
    return m_pWalkIndex->count(id) > 0;
}

//
void SGWalk::truncate(const VertexID& id)
{
    EdgePtrVec::iterator iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        if((*iter)->getEndID() == id)
            break;
        ++iter;
    }
    assert(iter != m_edges.end());
    m_edges.erase(iter + 1, m_edges.end());
}

// Construct the extension string corresponding to the path
std::string SGWalk::getString(SGWalkType type, SGWalkVertexPlacementVector* pPlacementVector) const
{
    std::string out;

    // Append the full length of the starting vertex to the walk
    if(type == SGWT_START_TO_END || type == SGWT_INTERNAL)
    {
        out.append(m_pStartVertex->getSeq().toString());

        // Add the first record to the placement vector if required
        if(pPlacementVector != NULL)
        {
            SGWalkVertexPlacement firstPlace;
            firstPlace.pVertex = m_pStartVertex;
            firstPlace.position = 0; // 0-based coordinates
            firstPlace.isRC = false;
            pPlacementVector->push_back(firstPlace);
        }
    }

    // Determine if the string should go to the end of the last vertex
    // in the path
    size_t stop = m_edges.size();

    // The first edge is always in correct frame of reference 
    // so the comp is EC_SAME. This variable tracks where the 
    // string that is being added is different from the starting sequence
    // and needs to be flipped
    EdgeComp currComp = EC_SAME;

    // If the walk direction is antisense, we reverse every component and then
    // reverse the entire string to generate the final string
    bool reverseAll = !m_edges.empty() && m_edges[0]->getDir() == ED_ANTISENSE;
    if(reverseAll)
        out = reverse(out);
    
    for(size_t i = 0; i < stop; ++i)
    {
        Edge* pYZ = m_edges[i];
        // Append in the extension string
        std::string edge_str = pYZ->getLabel();
        assert(edge_str.size() != 0);

        // Determine whether this node is reverse complement wrt the string
        // we are building
        if(currComp == EC_REVERSE)
            edge_str = reverseComplement(edge_str);

        if(reverseAll)
            edge_str = reverse(edge_str);

        // Calculate the next comp, between X and Z
        EdgeComp ecYZ = pYZ->getComp();
        EdgeComp ecXZ;
        if(ecYZ == EC_SAME)
            ecXZ = currComp;
        else
            ecXZ = !currComp;

        out.append(edge_str);
        
        // Add this record into the placement vector
        if(pPlacementVector != NULL)
        {
            SGWalkVertexPlacement placement;
            placement.pVertex = pYZ->getEnd();
            placement.isRC = ecXZ == EC_REVERSE;
            placement.position = out.size() - pYZ->getEnd()->getSeqLen();
            pPlacementVector->push_back(placement);
        }

        currComp = ecXZ;
    }

    // If we want the internal portion of the string (which does not contain the endpoints
    // perform the truncation now. This needs to be done before the reversal.
    if(type == SGWT_INTERNAL)
    {
        if(pPlacementVector != NULL)
        {
            std::cerr << "Error: Vertex placement not supported for SGWT_INTERNAL walk types\n";
            exit(EXIT_FAILURE);
        }

        Edge* pFirstEdge = getFirstEdge();
        Edge* pLastEdge = getLastEdge();
        if(pFirstEdge == NULL || pLastEdge == NULL)
        {
            out.clear();
        }
        else
        {
            Vertex* pStart = m_pStartVertex;
            Vertex* pLast = getLastVertex();
            int start = pStart->getSeqLen() - pFirstEdge->getMatchLength();
            int end = out.size() - (pLast->getSeqLen() - pLastEdge->getMatchLength());

            if(end <= start)
                out.clear();
            else
            {
                std::string ss = out.substr(start, end - start);
                out = ss;
            }
        }
    }

    if(out.empty())
        std::cout << "No output for walk: " << pathSignature() << "\n";

    if(reverseAll)
    {
        out = reverse(out);

        // Reverse the placement vector too, including flipping the alignment coordinates
        if(pPlacementVector != NULL)
        {
            std::reverse(pPlacementVector->begin(), pPlacementVector->end());
            for(size_t i = 0; i < pPlacementVector->size(); ++i)
            {
                SGWalkVertexPlacement& item = pPlacementVector->at(i);
                item.position = out.size() - item.position - item.pVertex->getSeqLen();
            }
        }
    }

    return out;
}

// Returns a vector of EdgeComps of the orientation of each
// vertex in the path with respect to the start of the walk
std::vector<EdgeComp> SGWalk::getOrientationsToStart() const
{
    std::vector<EdgeComp> out;
    if(m_edges.empty())
        return out;
    
    // Tracking variable of the orientation between X (the start)
    // and the next vertex in the walk
    EdgeComp compXY = m_edges[0]->getComp();
    out.push_back(compXY);

    for(size_t i = 1; i < m_edges.size(); ++i)
    {
        // Calculate the orientation of XZ using YZ
        EdgeComp compYZ = m_edges[i]->getComp();

        // The direction flips if the XY/YZ are different
        EdgeComp compXZ = (compXY == compYZ) ? EC_SAME : EC_REVERSE;
        out.push_back(compXZ);
        compXY = compXZ;
    }
    return out;
}

// Get the substring of the full path string starting from position fromX
// to position toY on the first and last vertices, respectively.
// dirX is the direction along contig X towards vertex Y, vis-versa for dirY
std::string SGWalk::getFragmentString(const Vertex* pX, const Vertex* pY,
                                      int fromX, int toY, 
                                      EdgeDir dirX, EdgeDir dirY) const
{
    std::string out;

    // Calculate the portion of X that we should include in the string
    // If dirX is SENSE, we take the everything after position fromX
    // otherwise we take everything up to and including fromX
    SeqCoord xCoord(0,0,pX->getSeqLen());

    if(dirX == ED_SENSE)
    {
        xCoord.interval.start = fromX;
        xCoord.interval.end = pX->getSeqLen() - 1;
    }
    else
    {
        xCoord.interval.start = 0;
        xCoord.interval.end = fromX;
    }

    // Handle the trivial case where pX == pY and the walk is found immediately
    if(m_edges.empty() && pX == pY)
    {
        if(dirY == ED_SENSE)
        {
            xCoord.interval.start = toY;
        }
        else
        {
            xCoord.interval.end = toY;
        }
    }

    if(!xCoord.isValid())
        return "";

    //
    out.append(m_pStartVertex->getSeq().substr(xCoord.interval.start, xCoord.length()));

    // Determine if the string should go to the end of the last vertex
    // in the path
    size_t stop = m_edges.size();

    // The first edge is always in correct frame of reference 
    // so the comp is EC_SAME. This variable tracks where the 
    // string that is being added is different from the starting sequence
    // and needs to be flipped
    EdgeComp currComp = EC_SAME;

    // If the walk direction is antisense, we reverse every component and then
    // reverse the entire string to generate the final string
    bool reverseAll = !m_edges.empty() && m_edges[0]->getDir() == ED_ANTISENSE;
    if(reverseAll)
        out = reverse(out);

    for(size_t i = 0; i < stop; ++i)
    {
        Edge* pYZ = m_edges[i];
        bool isLast = i == (stop - 1);

        if(!isLast)
        {
            // Append the extension string without modification
            std::string edge_str = pYZ->getLabel();
            assert(edge_str.size() != 0);
            if(currComp == EC_REVERSE)
                edge_str = reverseComplement(edge_str);

            if(reverseAll)
                edge_str = reverse(edge_str);
            out.append(edge_str);
        }
        else
        {
            // 
            const Edge* pZY = pYZ->getTwin();
            
            // get the unmatched coordinates on pY
            SeqCoord unmatched = pZY->getMatchCoord().complement();

            // Now, we have to shrink the unmatched interval on Y to
            // only incude up to toY
            if(dirY == ED_SENSE)
                unmatched.interval.start = toY;
            else
                unmatched.interval.end = toY;

            if(!unmatched.isValid())
                return "";

            std::string seq = unmatched.getSubstring(pY->getStr());
            if(pYZ->getComp() != currComp)
                seq = reverseComplement(seq);
            
            if(reverseAll)
                seq = reverse(seq);
            out.append(seq);
        }

        // Calculate the next comp, between X and Z
        EdgeComp ecYZ = pYZ->getComp();
        EdgeComp ecXZ;
        if(ecYZ == EC_SAME)
            ecXZ = currComp;
        else
            ecXZ = !currComp;

        currComp = ecXZ;
    }

    if(reverseAll)
        out = reverse(out);
    return out;
}

//
int SGWalk::getExtensionDistance() const
{
    return m_extensionDistance;
}

// This is equivalent to the extension distance
int SGWalk::getEndToEndDistance() const
{
    return m_extensionDistance;
}

// 
int SGWalk::getStartToEndDistance() const
{
    return m_pStartVertex->getSeqLen() + getEndToEndDistance();
}

// Returns the distance from the end of pStart to the beginning of the last vertex in the path
// This can be negative if they overlap
int SGWalk::getEndToStartDistance() const
{
    if(m_edges.empty())
        return 0;

    int len_x = m_pStartVertex->getSeqLen();
    int len_y = getLastEdge()->getEnd()->getSeqLen();
    return getStartToEndDistance() - (len_x + len_y);
}

//
Edge* SGWalk::getFirstEdge() const
{
    if(m_edges.empty())
        return NULL;
    else
        return m_edges.front();
}

//
Edge* SGWalk::getLastEdge() const
{
    if(m_edges.empty())
        return NULL;
    else
        return m_edges.back();
}

//
Vertex* SGWalk::getLastVertex() const
{
    if(m_edges.empty())
        return NULL;
    else
        return getLastEdge()->getEnd();
}

//
Edge* SGWalk::getEdge(size_t idx) const
{
    return m_edges[idx];
}


//
Vertex* SGWalk::getVertex(size_t idx) const
{
    assert(idx < getNumVertices());
    if(idx == 0)
        return m_pStartVertex;
    else
        return m_edges[idx - 1]->getEnd();
}

//
VertexPtrVec SGWalk::getVertices() const
{
    VertexPtrVec out;
    out.push_back(m_pStartVertex);
    for(EdgePtrVec::const_iterator iter = m_edges.begin(); iter != m_edges.end(); ++iter)
        out.push_back((*iter)->getEnd());
    return out;
}

// Return a string describing this path
std::string SGWalk::pathSignature() const
{
    std::stringstream ss;
    ss << m_pStartVertex->getID() << ",";
    for(EdgePtrVec::const_iterator iter = m_edges.begin(); iter != m_edges.end(); ++iter)
        ss << (*iter)->getEndID() << ", ";
    return ss.str();
}

// 
void SGWalk::print() const
{
    std::cout << "Walk start: " << m_pStartVertex->getID() << "\nWalk: ";
    const Vertex* pLast = m_pStartVertex;
    for(EdgePtrVec::const_iterator iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        //std::cout << *(*iter) << " ";
        std::cout << (*iter)->getStartID() << " -- " << (*iter)->getEndID() << "," << (*iter)->getDir() << "," << (*iter)->getComp() << "\t";
        assert((*iter)->getStart() == pLast);
        pLast = (*iter)->getEnd();
    }
    std::cout << "\n";
}

// 
void SGWalk::printSimple() const
{
    std::cout << pathSignature() << "\n";
}

// 
bool SGWalk::compareByTotalLength(const SGWalk& a, const SGWalk& b)
{
    return a.getStartToEndDistance() > b.getStartToEndDistance();
}

