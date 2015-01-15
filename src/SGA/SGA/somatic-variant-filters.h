//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// variant-filtering - apply various filters
// to a somatic VCF file
//
#ifndef SOMATIC_VARIANT_FILTERS_H
#define SOMATIC_VARIANT_FILTERS_H
#include <getopt.h>
#include "config.h"

// functions

//
int somaticVariantFiltersMain(int argc, char** argv);

// options
void parseSomaticVariantFiltersOptions(int argc, char** argv);

#endif
