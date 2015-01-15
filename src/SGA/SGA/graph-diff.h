//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// graph-diff - Find strings that are only present
// in one of two input graphs
//
#ifndef GRAPH_DIFF_H
#define GRAPH_DIFF_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"
#include "Match.h"
#include "BWTAlgorithms.h"
#include "OverlapAlgorithm.h"

// functions

//
int graphDiffMain(int argc, char** argv);

// options
void parseGraphDiffOptions(int argc, char** argv);

#endif
