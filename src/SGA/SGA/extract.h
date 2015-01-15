//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// extract - Extract all sequences of a given length from a BWT
//
#ifndef EXTRACT_H
#define EXTRACT_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"
#include "Match.h"


// functions

int extractMain(int argc, char** argv);


// options
void parseExtractOptions(int argc, char** argv);

#endif
