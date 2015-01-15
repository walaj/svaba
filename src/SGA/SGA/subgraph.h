//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// subgraph - extract a subgraph from an assembly graph
//
#ifndef SUBGRAPH_H
#define SUBGRAPH_H
#include <getopt.h>
#include "config.h"

// functions
int subgraphMain(int argc, char** argv);
void parseSubgraphOptions(int argc, char** argv);
void subgraph();

#endif
