//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// fm-merge - Merge reads/sequences using the FM-index
// without explicitly constructing the full graph
//
#ifndef FMMERGE_H
#define FMMERGE_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"
#include "Match.h"
#include "BWTAlgorithms.h"
#include "OverlapAlgorithm.h"

// functions

//
int FMMergeMain(int argc, char** argv);

// options
void parseFMMergeOptions(int argc, char** argv);

#endif
