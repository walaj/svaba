//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// gapfill - Fill intrascaffold gaps
//
#ifndef GAPFILL_H
#define GAPFILL_H
#include <getopt.h>
#include "config.h"

// functions

//
int gapfillMain(int argc, char** argv);

// options
void parseGapFillOptions(int argc, char** argv);

#endif
