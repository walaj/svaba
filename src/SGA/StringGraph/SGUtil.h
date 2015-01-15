//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SGUtils - Data structures/Functions related
// to building and manipulating string graphs
//
#ifndef SGUTIL_H
#define SGUTIL_H

#include "Bigraph.h"
#include "ASQG.h"

// typedefs
typedef Bigraph StringGraph;

namespace SGUtil
{
// Main string graph loading function
// The allowContainments flag forces the string graph to retain identical vertices
// Vertices that are substrings of other vertices (SS flag = 1) are never kept
StringGraph* loadASQG(const std::string& filename, const unsigned int minOverlap, bool allowContainments = false, size_t maxEdges = -1);
StringGraph* loadASQG(std::stringstream& pReader, const unsigned int minOverlap, bool allowContainments = false, size_t maxEdges = -1);

// Load a string graph from a fasta file.
// Returns a graph where each sequence in the fasta is a vertex but there are no edges in the graph.
StringGraph* loadFASTA(const std::string& filename);


};
#endif
