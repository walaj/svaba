//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Base bidirectional edge class 
//

#ifndef EDGE_H
#define EDGE_H

#include <ostream>
#include "Match.h"
#include "Util.h"
#include "GraphCommon.h"
#include "EdgeDesc.h"
#include "Vertex.h"
#include "BitChar.h"
#include "SimpleAllocator.h"

// Packed structure holding the direction and comp of an edge
// The EdgeDir/EdgeComp enums (which only have values 0/1) are used 
// as the interface for this class so the set flags are cast to/from
// these types.
struct EdgeData
{
    public:
        EdgeData() {}

        // Setters 
        void setDir(EdgeDir dir)
        {
            if(dir == ED_ANTISENSE)
                m_data.set(DIR_BIT, true);
            else
                m_data.set(DIR_BIT, false);
        }

        void setComp(EdgeComp comp)
        {
            if(comp == EC_REVERSE)
                m_data.set(COMP_BIT, true);
            else
                m_data.set(COMP_BIT, false);
        }

        void flipDir() { m_data.flip(DIR_BIT); }
        void flipComp() { m_data.flip(COMP_BIT); }

        // Getters
        inline EdgeDir getDir() const
        {
            return m_data.test(DIR_BIT) ? ED_ANTISENSE : ED_SENSE;
        }

        inline EdgeComp getComp() const
        {
            return m_data.test(COMP_BIT) ? EC_REVERSE : EC_SAME;
        }

    private:
        static const size_t DIR_BIT = 0;
        static const size_t COMP_BIT = 1;
        BitChar m_data;
};

class Edge
{
    public:
        Edge(Vertex* end, EdgeDir dir, EdgeComp comp, SeqCoord m) : 
                 m_pEnd(end), m_pTwin(NULL), m_matchCoord(m), m_color(GC_WHITE), isTrusted(false)
        {
            m_edgeData.setDir(dir);
            m_edgeData.setComp(comp);
        }

        ~Edge() { }
        
        // High level modification functions
        
        // Join merges pEdge into this edge, with pEdge describing the starting point
        void join(const Edge* pEdge);

        // Extend merged pEdge into this edge, with pEdge describing the endpoint
        void extend(const Edge* pEdge);

        // Post merge update function
        void update() {}

        // Sequence Coordinate functions
        
        // update functions
        void extendMatchFullLength();
        void extendMatch(int ext_len);
        void offsetMatch(int offset);
        void updateSeqLen(int newLen);

        // access functions
        size_t getSeqLen() const;
        size_t getMatchLength() const { return m_matchCoord.length(); }
        const SeqCoord& getMatchCoord() const { return m_matchCoord; }
        std::string getMatchStr() const;
        Match getMatch() const;
        Overlap getOverlap() const;
        
        // setters
        void setTwin(Edge* pEdge) { m_pTwin = pEdge; }
        void setColor(GraphColor c) { m_color = c; }

        // getters
        VertexID getStartID() const { return getStart()->getID(); }
        VertexID getEndID() const { return m_pEnd->getID(); }
        inline Vertex* getStart() const { assert(m_pTwin != NULL); return m_pTwin->getEnd(); }
        inline Vertex* getEnd() const { return m_pEnd; }
        inline EdgeDir getDir() const { return m_edgeData.getDir(); }
        inline EdgeComp getComp() const { return m_edgeData.getComp(); }        
        inline Edge* getTwin() const { assert(m_pTwin != NULL); return m_pTwin; }
        EdgeDesc getTwinDesc() const;
        std::string getLabel() const;
        bool isSelf() const { return getStart() == getEnd(); }
        inline GraphColor getColor() const { return m_color; }
        size_t getMemSize() const { return sizeof(*this); }

        // Returns the direction of an edge that continues in the same direction
        // as this edge, corrected for complementary 
        inline EdgeDir getTransitiveDir() const { return (getComp() == EC_SAME) ? getDir() : !getDir(); }

        // Make the direction of the edge that its twin should point along 
        inline EdgeDir getTwinDir() const { return (getComp() == EC_SAME) ? !getDir() : getDir(); }
        inline EdgeDesc getDesc() const { return EdgeDesc(getEnd(), getDir(), getComp()); }
        
        // Flip the edge
        inline void flipComp() { m_edgeData.flipComp(); }
        inline void flipDir() { m_edgeData.flipDir(); }
        void flip() { flipComp(); flipDir(); }

        // Memory management
        //void* operator new(size_t /*size*/, SimpleAllocator<Edge>* pAllocator)
        //{
        //    return pAllocator->alloc();
        //}

        //void operator delete(void* /*target*/, size_t /*size*/)
        //{
            // Deletions are handled at the graph/pool level. The lifetime of an edge
            // is as long as the graph it belongs to, even if it is deleted before
            // the graph.
        //}

        // Validate that the edge is sane
        void validate() const;

        // Output
        friend std::ostream& operator<<(std::ostream& out, const Edge& obj);

    protected:
        
        // Global new is not allowed, allocation must go through the memory pool
        // belonging to the graph.
        //void* operator new(size_t size) { return malloc(size); } 
        
        Edge() {}; // Default constructor is not allowed

        Vertex* m_pEnd;
        Edge* m_pTwin;
        SeqCoord m_matchCoord;
        GraphColor m_color;
        EdgeData m_edgeData; // dir/comp member

    public:
        bool isTrusted;
};

#endif
