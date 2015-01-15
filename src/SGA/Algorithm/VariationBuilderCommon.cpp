///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// VariantionVariationBuilderCommon -- Common functions and data
// structures for the abstract graph builders
//
#include "VariationBuilderCommon.h"
#include "SGAlgorithms.h"

// Count the number of extensions above the given threshold
size_t VariationBuilderCommon::countValidExtensions(const AlphaCount64& ac, size_t threshold)
{
    size_t n = 0;
    for(size_t i = 0; i < DNA_ALPHABET::size; ++i)
    {
        char b = DNA_ALPHABET::getBase(i);
        size_t count = ac.get(b);
        if(count >= threshold)
            n += 1;
    }
    return n;
}

// Filter out low counts in AlphaCount using a coverage threshold
// relative to the most frequent count. Returns the number of
// surviving counts
size_t VariationBuilderCommon::filterLowFrequency(AlphaCount64& ac, double alpha)
{
    size_t n = 0;
    size_t max = ac.getMaxCount();
    if(max == 0)
        return 0;

    for(size_t i = 0; i < DNA_ALPHABET::size; ++i)
    {
        char b = DNA_ALPHABET::getBase(i);
        size_t count = ac.get(b);
        double ratio = (double)count / max;
        if(ratio >= alpha)
            n += 1;
        else
            ac.set(b,0);
    }
    return n;
}


// Make a de Bruijn graph string 
std::string VariationBuilderCommon::makeDeBruijnVertex(const std::string& v, char edgeBase, EdgeDir direction)
{
    std::string w;
    size_t p = v.size() - 1;
    if(direction == ED_SENSE)
    {
        w = v.substr(1, p);
        w.append(1, edgeBase);
    }
    else
    {
        w.append(1, edgeBase);
        w.append(v.substr(0, p));
    }
    return w;
}

// Add a de Bruijn graph edge to the given graph betwee pX and pY.
// Assumes the sequences are from the same strand (not reverse complements)
void VariationBuilderCommon::addSameStrandDeBruijnEdges(StringGraph* pGraph, const Vertex* pX, const Vertex* pY, EdgeDir direction)
{
    assert(pX->getSeq().length() == pY->getSeq().length());
    
    // overlap length for a de bruijn edge
    size_t p = pX->getSeq().length() - 1;

    // Construct an overlap object for this relationship
    Overlap o;
    o.id[0] = pX->getID();
    o.id[1] = pY->getID();

    o.match.isReverse = false;
    o.match.numDiff = 0;

    if(direction == ED_SENSE)
    {
        // pX -> pY
        o.match.coord[0].interval.start = 1;
        o.match.coord[1].interval.start = 0;
    }
    else
    {
        // pY -> pX
        o.match.coord[0].interval.start = 0;
        o.match.coord[1].interval.start = 1;
    }

    o.match.coord[0].interval.end = o.match.coord[0].interval.start + p - 1; // inclusive coordinate
    o.match.coord[1].interval.end = o.match.coord[1].interval.start + p - 1;
    o.match.coord[0].seqlen = p + 1;
    o.match.coord[1].seqlen = p + 1;
    Edge* e = SGAlgorithms::createEdgesFromOverlap(pGraph, o, false);
    assert(e != NULL);
}

