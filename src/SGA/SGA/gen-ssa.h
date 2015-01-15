//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// gen-ssa - Build a sampled suffix array for a set of reads
// using its BWT
//
#ifndef GENSSA_H
#define GENSSA_H
#include <getopt.h>
#include "config.h"
#include "SuffixArray.h"

int genSSAMain(int argc, char** argv);
void parseGenSSAOptions(int argc, char** argv);

#endif
