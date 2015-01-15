//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// correct - Correct sequencing errors in reads using the FM-index
//
#ifndef CORRECT_H
#define CORRECT_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"
#include "Match.h"
#include "BWTAlgorithms.h"
#include "OverlapAlgorithm.h"

// functions

//
int correctMain(int argc, char** argv);

// options
void parseCorrectOptions(int argc, char** argv);

#endif
