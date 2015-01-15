//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// merge - merge read files and their associated indices 
//
#ifndef MERGE_H
#define MERGE_H
#include <getopt.h>
#include "config.h"

int mergeMain(int argc, char** argv);
void parseMergeOptions(int argc, char** argv);

#endif
