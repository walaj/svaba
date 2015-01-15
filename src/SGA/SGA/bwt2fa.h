//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// bwt2fa - Transform a bwt back into a set of sequences
//
#ifndef BWT2FA_H
#define BWT2FA_H
#include <getopt.h>
#include "config.h"

int bwt2faMain(int argc, char** argv);
void parseBWT2FAOptions(int argc, char** argv);

#endif
