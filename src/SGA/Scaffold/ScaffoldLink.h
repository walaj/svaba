//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldLink - Data structure representing
// an inferred relationship between two scaffold 
// components
//
#ifndef SCAFFOLDLINK_H
#define SCAFFOLDLINK_H

#include "GraphCommon.h"
#include "Edge.h"

enum ScaffoldLinkType
{
    SLT_DISTANCEEST,
    SLT_REFERENCE,
    SLT_INFERRED,
    SLT_NOTYPE
};

class ScaffoldLink
{
    public:
        ScaffoldLink() {}
        ScaffoldLink(const std::string& id, EdgeDir dir, EdgeComp comp, 
                     int dist, double sd, int np, int sl, ScaffoldLinkType slt);

        // Get the far endpoint of the link
        // Defined to be distance + seqLen
        int getEndpoint() const;

        EdgeDir getDir() const;
        EdgeDir getTwinDir() const;
        EdgeComp getComp() const;
        char getTypeCode() const;
        static ScaffoldLinkType getType(char tc);

        // IO
        void parse(const std::string& text);
        friend std::ostream& operator<<(std::ostream& out, const ScaffoldLink& link);

        // data 
        std::string endpointID;
        int distance;
        double stdDev;
        int numPairs;
        int seqLen; // length of the endpoint contig
        ScaffoldLinkType type;
        EdgeData edgeData; // dir/comp member
};

#endif
