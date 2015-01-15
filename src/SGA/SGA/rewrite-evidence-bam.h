//-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// rewrite-evidence-bam - fill in read name and quality
// information in an SGA graph-diff evidence bam file
//
#ifndef REWRITE_EVIDENCE_BAM_H
#define REWRITE_EVIDENCE_BAM_H
#include <getopt.h>
#include "config.h"

// functions

//
int rewriteEvidenceBAMMain(int argc, char** argv);

// options
void parseRewriteEvidenceBAMOptions(int argc, char** argv);

#endif
