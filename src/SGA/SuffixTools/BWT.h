//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWT - All functions that use a BWT include this file
// it simple typedefs the BWT name to the implementation
// of the BWT that we want, either the uncompressed version
// (SBWT) or the run-length encoded version (RLBWT). This could 
// be done using inheritence but the BWT is so used so much that 
// overhead of calling virtual functions is unwanted
//          
//
#ifndef BWT_H
#define BWT_H

#include "RLBWT.h"
#include "SBWT.h"

typedef RLBWT BWT;

#endif
