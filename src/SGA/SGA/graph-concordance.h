//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// vcf-read-count - count the number of reads
// supporting variant haplotypes
//
#ifndef VCF_READ_COUNT_H
#define VCF_READ_COUNT_H
#include <getopt.h>
#include "config.h"

// functions

//
int graphConcordanceMain(int argc, char** argv);

// options
void parseGraphConcordanceOptions(int argc, char** argv);

#endif
