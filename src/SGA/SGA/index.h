//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// index - Build a BWT/FM-index for a set of reads
//
#ifndef INDEX_H
#define INDEX_H
#include <getopt.h>
#include "config.h"
#include "SuffixArray.h"

int indexMain(int argc, char** argv);
void indexInMemorySAIS();
void indexInMemoryBCR();
void indexInMemoryRopebwt();
void indexOnDisk();
void buildIndexForTable(std::string outfile, const ReadTable* pRT, bool isReverse);
void parseIndexOptions(int argc, char** argv);

#endif
