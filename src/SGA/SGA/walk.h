//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// walk - traverse the graph from one vertex to another
// and write out the walk-strings to stdout
//
#ifndef WALK_H
#define WALK_H
#include <getopt.h>
#include "config.h"

// functions
int walkMain(int argc, char** argv);
void parseWalkOptions(int argc, char** argv);
void walk();

#endif
