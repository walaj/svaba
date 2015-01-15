//-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Detect the power to detect variants at specific
// locations of a reference genome using k-mers
//
#ifndef VARIANT_DETECTABILITY_H
#define VARIANT_DETECTABILITY_H
#include <getopt.h>
#include "config.h"

int variantDetectabilityMain(int argc, char** argv);
void parseVarDetectOptions(int argc, char** argv);

#endif
