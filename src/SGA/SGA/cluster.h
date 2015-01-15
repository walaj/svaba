//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// cluster - find connected components in an 
// overlap graph
//
#ifndef CLUSTER_H
#define CLUSTER_H
#include <getopt.h>
#include "config.h"

// functions
int clusterMain(int argc, char** argv);
void parseClusterOptions(int argc, char** argv);
void cluster();

#endif
