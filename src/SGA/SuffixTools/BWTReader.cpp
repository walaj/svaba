//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTReader - Abstract class for reading a BWT file
//
#include "BWTReader.h"
#include "BWTReaderBinary.h"
#include "BWTReaderAscii.h"

//
IBWTReader* BWTReader::createReader(const std::string& filename)
{
    return new BWTReaderBinary(filename);
}
