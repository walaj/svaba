//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SQG - Common functions for parsing 
// assembly and sequence graph files
//
#include "SQG.h"

namespace SQG
{

//
StringVector tokenizeRecord(const std::string& record)
{
    return split(record, FIELD_SEP);
}

//
StringVector tokenizeTagValue(const std::string& tagValue)
{
    return split(tagValue, TAG_SEP);
}

}
