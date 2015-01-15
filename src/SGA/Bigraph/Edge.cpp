//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Base bidirectional edge class 
//
#include "Edge.h"
#include "Vertex.h"

// 
EdgeDesc Edge::getTwinDesc() const
{
    return EdgeDesc(getStart(), getTwinDir(), getComp());
}

// Get the edge's label
std::string Edge::getLabel() const
{
    const Edge* pTwin = getTwin();
    const Vertex* pEndpoint = m_pEnd;
    
    // get the unmatched coordinates in V2
    SeqCoord unmatched = pTwin->getMatchCoord().complement();
    std::string seq = unmatched.getSubstring(pEndpoint->getStr());

    if(getComp() == EC_REVERSE)
        seq = reverseComplement(seq);

    return seq;
}


// Join the edge pEdge into this edge, adding to the start
void Edge::join(const Edge* pEdge)
{
    // Update the match coordinate
    Match m12 = pEdge->getMatch();
    Match m23 = getMatch();
    m_matchCoord = m12.inverseTranslate(m23.coord[0]);

    if(pEdge->getComp() == EC_REVERSE)
        flip();

    // Now, update the twin of this edge to extend to the twin of pEdge
    m_pTwin->extend(pEdge->getTwin());
}

// Extend this edge by adding pEdge to the end
void Edge::extend(const Edge* pEdge)
{
    if(pEdge->getComp() == EC_REVERSE)
        flipComp();
    m_pEnd = pEdge->getEnd();
}

// return the mapping from V1 to V2 via this edge
// If necessary, flip the SC of V2 so that it is
// in the same coordinate system as V1
Match Edge::getMatch() const
{
    const SeqCoord& sc = getMatchCoord();
    const SeqCoord& tsc = m_pTwin->getMatchCoord();
    return Match(sc, tsc, getComp() == EC_REVERSE, -1);
}

// The overlap structure implied by this edge
Overlap Edge::getOverlap() const
{
    return Overlap(getStartID(), getEndID(), getMatch());
}

// Get the matching portion of V1 described by this edge
std::string Edge::getMatchStr() const
{
    return m_matchCoord.getSubstring(getStart()->getStr());
}

// Return the length of the sequence
size_t Edge::getSeqLen() const
{
    SeqCoord unmatched = m_pTwin->getMatchCoord().complement();
    return unmatched.length();
}

void Edge::offsetMatch(int offset)
{
    m_matchCoord.interval.start += offset;
    m_matchCoord.interval.end += offset;
}

void Edge::extendMatch(int ext_len)
{
    m_matchCoord.interval.end += ext_len;
}

// Bump the edges of the match outwards so it covers the entire 
// sequence
void Edge::extendMatchFullLength()
{
    if(m_matchCoord.isLeftExtreme())
        m_matchCoord.interval.end = m_matchCoord.seqlen - 1;
    else
        m_matchCoord.interval.start = 0;
}

void Edge::updateSeqLen(int newLen)
{
    m_matchCoord.seqlen = newLen;
}

// Validate that the edge members are sane
void Edge::validate() const
{
    const Edge* pTwin = getTwin();
    std::string m_v1 = getMatchStr();
    std::string m_v2 = pTwin->getMatchStr();

    if(getComp() == EC_REVERSE)
        m_v2 = reverseComplement(m_v2);

    bool error = false;
    if(m_v1.length() != m_v2.length())
    {
        std::cerr << "Error, matching strings are not the same length\n";
        error = true;
    }

    if(error)
    {
        std::cerr << "V1M: " << m_v1 << "\n";
        std::cerr << "V2M: " << m_v2 << "\n";
        std::cerr << "V1MC: " << getMatchCoord() << "\n";
        std::cerr << "V2MC: " << pTwin->getMatchCoord() << "\n";
        std::cerr << "V1: " << getStart()->getSeq() << "\n";
        std::cerr << "Validation failed for edge " << *this << "\n";
        assert(false);
    }
}

// Output
std::ostream& operator<<(std::ostream& out, const Edge& obj)
{
    out << obj.m_pEnd->getID() << "," << obj.getDir() << "," << obj.getComp();
    return out;
}

