//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// stats - Print statistics about the data set
//
#ifndef STATS_H
#define STATS_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"
#include "Match.h"
#include "BWTAlgorithms.h"
#include "OverlapAlgorithm.h"

// functions

//
int statsMain(int argc, char** argv);

// options
void parseStatsOptions(int argc, char** argv);

#endif
