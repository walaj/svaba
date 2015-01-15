//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// QuickBWT - Construct a simple BWT for a short input string
//
#include "BWT.h"

// Construct a BWT/suffix array of the given string
// Calling code is responsible for freeing the BWT/SA
void createQuickBWT(const std::string& str, BWT*& pBWT, SuffixArray*& pSA);
