//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// pre-QC - Perform pre-assembly quality checks on a set
//          of reads
//
#ifndef PREQC_H
#define PREQC_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"
#include "Match.h"
#include "BWTAlgorithms.h"
#include "OverlapAlgorithm.h"

// functions

//
int preQCMain(int argc, char** argv);

// options
void parsePreQCOptions(int argc, char** argv);

#endif
