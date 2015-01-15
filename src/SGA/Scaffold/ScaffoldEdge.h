//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldEdge - An edge in a scaffold graph. 
// Wraps the ScaffoldLink data.
//
#ifndef SCAFFOLDEDGE_H
#define SCAFFOLDEDGE_H

#include "Edge.h"
#include "ScaffoldLink.h"

class ScaffoldVertex;

class ScaffoldEdge
{
    public:

        ScaffoldEdge(ScaffoldVertex* pEnd, ScaffoldLink link);

        //
        void setTwin(ScaffoldEdge* pEdge);
        void setLink(ScaffoldLink link);
        void setColor(GraphColor c);

        //
        VertexID getStartID() const;
        VertexID getEndID() const;
        ScaffoldVertex* getStart() const;
        ScaffoldVertex* getEnd() const;
        ScaffoldEdge* getTwin() const;
        const ScaffoldLink& getLink() const;

        //
        EdgeDir getDir() const;
        EdgeComp getComp() const;        
        ScaffoldLinkType getType() const;
        int getDistance() const;
        double getStdDev() const;
        char getTypeCode() const;
        GraphColor getColor() const;

        //
        std::string makeLinkString() const;

        //
        friend std::ostream& operator<<(std::ostream& out, const ScaffoldEdge& edge);

    private:
        ScaffoldVertex* m_pEnd;
        ScaffoldEdge* m_pTwin;
        ScaffoldLink m_link; // distance, direction information
        GraphColor m_color;
};

bool ScaffoldEdgePtrDistanceCompare(ScaffoldEdge* pXY, ScaffoldEdge* pXZ);


#endif
