//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// filterBAM - remove erroneous paired end connections
// from a bam file.
//
#ifndef FILTERBAM_H
#define FILTERBAM_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"
#include "Match.h"
#include "BWTAlgorithms.h"
#include "OverlapAlgorithm.h"

// functions

//
int filterBAMMain(int argc, char** argv);

// options
void parseFilterBAMOptions(int argc, char** argv);

#endif
