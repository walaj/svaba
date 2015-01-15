//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// GraphCommon - Definitions and functions common
// to the grpah implementation
//
#ifndef GRAPHCOMMON_H
#define GRAPHCOMMON_H

#include <vector>
#include "Util.h"

// The directions an edge can take.
// Do not change the values
enum EdgeDir
{
        ED_SENSE = 0,
        ED_ANTISENSE = 1,
        ED_COUNT = 2
};

// Flag indicating whether the sequences linked by an edge
// are from the same strand or not.
// Do not change the values
enum EdgeComp
{
        EC_SAME = 0,
        EC_REVERSE = 1
};

// Array holding the directions, for iteration
const EdgeDir EDGE_DIRECTIONS[ED_COUNT] = { ED_SENSE, ED_ANTISENSE };

// Flags specifying how the dot file should be drawn
enum DotFlags
{
    DF_UNDIRECTED = 0x01,
    DF_ANNOTATIONS = 0x02,
    DF_COLORED = 0x04,
    DF_NOID = 0x08
};

// GraphColors are generic flags that can be used to indicate state
typedef uint8_t GraphColor;

const GraphColor GC_WHITE = 0;
const GraphColor GC_GRAY = 1;
const GraphColor GC_BLACK = 2;
const GraphColor GC_BLUE = 3;
const GraphColor GC_RED = 4;

// Typedefs
typedef std::string VertexID;
typedef std::vector<VertexID> VertexIDVec;

//
// Edge Operations
//
inline EdgeDir operator!(const EdgeDir& dir)
{
        return (dir == ED_SENSE) ? ED_ANTISENSE : ED_SENSE;
}

inline EdgeComp operator!(const EdgeComp& comp)
{
        return (comp == EC_SAME) ? EC_REVERSE : EC_SAME;
}

// Correct an edges direction, given the relationship between the nodes
inline EdgeDir correctDir(EdgeDir dir, EdgeComp comp)
{
    return (comp == EC_SAME) ? dir : !dir;
}


#endif 
