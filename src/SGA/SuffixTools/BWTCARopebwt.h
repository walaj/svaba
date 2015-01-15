//-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// BWTCARopebwt - Construct the BWT for a set of reads
// using Heng Li's ropebwt implementation
//
#ifndef BWTCA_ROPEBWT_H
#define BWTCA_ROPEBWT_H

#include <string>

namespace BWTCA
{
    void runRopebwt(const std::string& input_filename, const std::string& bwt_out_name, 
                    bool use_threads, bool do_reverse);
};

#endif
