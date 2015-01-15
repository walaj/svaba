//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// filter - remove reads from a data set based on various criteria
//
#ifndef FILTER_H
#define FILTER_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"
#include "Match.h"
#include "BWTAlgorithms.h"
#include "OverlapAlgorithm.h"

// functions

//
int filterMain(int argc, char** argv);

// options
void parseFilterOptions(int argc, char** argv);

#endif
