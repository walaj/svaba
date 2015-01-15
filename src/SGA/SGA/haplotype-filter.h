//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// haplotype-filter - perform quality checks
// on haplotypes and their associated variants
//
#ifndef HAPLOTYPE_FILTER_H
#define HAPLOTYPE_FILTER_H
#include <getopt.h>
#include "config.h"

// functions

//
int haplotypeFilterMain(int argc, char** argv);

// options
void parseHaplotypeFilterOptions(int argc, char** argv);

#endif
