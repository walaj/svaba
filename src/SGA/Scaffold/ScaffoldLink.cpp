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
#include "ScaffoldLink.h"

ScaffoldLink::ScaffoldLink(const std::string& id, 
                           EdgeDir dir, 
                           EdgeComp comp, 
                           int dist, 
                           double sd, 
                           int np, 
                           int sl, 
                           ScaffoldLinkType slt) : endpointID(id), 
                                                   distance(dist), 
                                                   stdDev(sd), 
                                                   numPairs(np), 
                                                   seqLen(sl),
                                                   type(slt)
{
    edgeData.setDir(dir);
    edgeData.setComp(comp);
}

int ScaffoldLink::getEndpoint() const
{
    assert(seqLen > 0);
    return distance + seqLen;
}

//
EdgeDir ScaffoldLink::getDir() const
{
    return edgeData.getDir();
}

// Return the direction of the twin edge
EdgeDir ScaffoldLink::getTwinDir() const
{
    return (getComp() == EC_SAME) ? !getDir() : getDir();     
}

//
EdgeComp ScaffoldLink::getComp() const
{
    return edgeData.getComp();
}

//
char ScaffoldLink::getTypeCode() const
{
    switch(type)
    {
        case SLT_DISTANCEEST:
            return 'D';
        case SLT_REFERENCE:
            return 'R';
        case SLT_INFERRED:
            return 'I';
        default:
            return 'N';
    }
}

//
ScaffoldLinkType ScaffoldLink::getType(char tc)
{
    switch(tc)
    {
        case 'D':
            return SLT_DISTANCEEST;
        case 'R':
            return SLT_REFERENCE;
        case 'I':
            return SLT_INFERRED;
        default:
            return SLT_NOTYPE;
    }
}

//
void ScaffoldLink::parse(const std::string& text)
{
    StringVector fields = split(text, ',');
    assert(fields.size() == 6);
    endpointID = fields[0];

    std::stringstream d_parser(fields[1]);
    d_parser >> distance;
    
    std::stringstream sd_parser(fields[2]);
    sd_parser >> stdDev;
    
    std::stringstream dir_parser(fields[3]);
    int dir;
    dir_parser >> dir;
    
    std::stringstream comp_parser(fields[4]);
    int comp;
    comp_parser >> comp;

    std::stringstream type_parser(fields[5]);
    char tc;
    type_parser >> tc;

    edgeData.setDir((EdgeDir)dir);
    edgeData.setComp((EdgeComp)comp);
    type = getType(tc);
}

//
std::ostream& operator<<(std::ostream& out, const ScaffoldLink& link)
{
    out << link.endpointID << "," << link.distance << "," << link.stdDev << ","
        << link.edgeData.getDir() << "," << link.edgeData.getComp() << "," << link.getTypeCode();
    return out;
}
