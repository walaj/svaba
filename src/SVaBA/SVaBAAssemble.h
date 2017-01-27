//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// assemble - Assemble reads into contigs
//
#ifndef SVABA_ASSEMBLE_H
#define SVABA_ASSEMBLE_H
#include <getopt.h>
#include "config.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "Util.h"
#include "SGUtil.h"

#include "SeqLib/UnalignedSequence.h"


struct AssemblyOptions {
  
  unsigned int verbose = 1;
  std::string asqgFile;
  std::string outContigsFile;
  std::string outVariantsFile;
  std::string outGraphFile;
  
  unsigned int minOverlap;
  bool bEdgeStats = false;
  bool bSmoothGraph = false;
  int resolveSmallRepeatLen = -1;
  size_t maxEdges = 128;
  
  // Trim parameters
  int numTrimRounds = 10;
  size_t trimLengthThreshold = 300;
  
  // Bubble parameters
  int numBubbleRounds = 3;
  double maxBubbleDivergence = 0.05f;
  double maxBubbleGapDivergence = 0.01f;
  int maxIndelLength = 20;
  
  // 
  bool bValidate;
  bool bExact = true;
  bool bPerformTR = false;
  
};

StringGraph* assemble(std::stringstream& asqg_stream, int minOverlap, int maxEdges, bool bExact, 
	      int trimLengthThreshold, bool bPerformTR, bool bValidate, int numTrimRounds, 
              int resolveSmallRepeatLen, int numBubbleRounds, double maxBubbleGapDivergence, 
              double maxBubbleDivergence, int maxIndelLength, int cutoff, std::string prefix, 
		      SeqLib::UnalignedSequenceVector &contigs, bool walk_all, bool get_components);


#endif
