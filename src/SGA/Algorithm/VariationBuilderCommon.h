///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// VariationBuilderCommon - Data structures
// commons to various algorithms for assembling
// variation haplotypes
//
#ifndef VARIATION_BUILDER_COMMON_H
#define VARIATION_BUILDER_COMMON_H
#include "BWT.h"
#include "BWTInterval.h"
#include "SGUtil.h"
#include "SGWalk.h"
#include <queue>

// Result of the bubble construction
// Used to track statistics
enum BubbleResultCode
{
    BRC_UNKNOWN,
    BRC_OK,
    BRC_SOURCE_BROKEN,
    BRC_SOURCE_BRANCH,
    BRC_TARGET_BROKEN,
    BRC_TARGET_BRANCH,
    BRC_WALK_FAILED,
    BRC_HB_FAILED,
    BRC_NO_SOLUTION
};

// The actual result structure
struct BubbleResult
{
    std::string targetString;
    std::string sourceString;
    double targetCoverage;
    double sourceCoverage;
    BubbleResultCode returnCode;
};

// A directed node in a graph that
// has not had it's neighbors visited
struct BuilderExtensionNode
{
    BuilderExtensionNode(Vertex* pX, EdgeDir d) : pVertex(pX), direction(d), distance(0) {}
    BuilderExtensionNode(Vertex* pX, EdgeDir dir, int dist) : pVertex(pX), direction(dir), distance(dist) {}

    Vertex* pVertex; // the vertex to extend
    EdgeDir direction; // the direction to extend to
    int distance; // the total number of nodes from the start to this node.
};
typedef std::queue<BuilderExtensionNode> BuilderExtensionQueue;
typedef std::map<std::string, int> StrIntMap;

// Common functions
namespace VariationBuilderCommon
{

// Count the number of extensions above the given threshold
size_t countValidExtensions(const AlphaCount64& ac, size_t threshold);

// Returns true if there are multiple characters above the threshold
bool hasMultipleBranch(const AlphaCount64& ac, size_t threshold);

// Filter out low counts in AlphaCount using a coverage threshold
// relative to the most frequent count. Returns the number of
// surviving counts
size_t filterLowFrequency(AlphaCount64& ac, double alpha);

// Make a de Bruijn graph string 
std::string makeDeBruijnVertex(const std::string& v, char edgeBase, EdgeDir direction);

// Add a de Bruijn graph edge to the given graph betwee pX and pY.
void addSameStrandDeBruijnEdges(StringGraph* pGraph, const Vertex* pX, const Vertex* pY, EdgeDir direction);

};

#endif
