//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// assemble - Assemble reads into contigs
//
#ifndef ASSEMBLE_H
#define ASSEMBLE_H
#include <getopt.h>
#include "config.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "Util.h"
#include "contigs.h"

// functions
int assembleMain(int argc, char** argv);
void parseAssembleOptions(int argc, char** argv);
//void assemble();
void assemble(std::stringstream& asqg_stream, int minOverlap, int maxEdges, bool bExact, 
	      int trimLengthThreshold, bool bPerformTR, bool bValidate, int numTrimRounds, 
              int resolveSmallRepeatLen, int numBubbleRounds, double maxBubbleGapDivergence, 
              double maxBubbleDivergence, int maxIndelLength, int cutoff, std::string prefix, 
              ContigVector &contigs);


#endif
