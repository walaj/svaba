//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SearchSeed.h - Data structure holding a partial 
// alignment to a FM-index. 
//
#include <stdio.h>
#include "SearchSeed.h"

//
void SearchSeed::print() const
{
    printf("li: %d ri: %d sl: %d dir: %d z: %d lrl: %d lru: %d rlr: %d rlu: %d\n", 
            left_index, right_index, seed_len, dir, z, 
            (int)ranges.interval[0].lower, (int)ranges.interval[0].upper, 
            (int)ranges.interval[1].lower, (int)ranges.interval[1].upper);
}

//
void SearchSeed::print(const std::string& w) const
{
    int range = (int)ranges.interval[0].upper - (int)ranges.interval[0].lower;
    printf("sub: %s li: %d ri: %d sl: %d dir: %d z: %d lrl: %d lru: %d range: %d rlr: %d rlu: %d\n", 
            w.substr(left_index, length()).c_str(), left_index, right_index,
            seed_len, dir, z, (int)ranges.interval[0].lower, (int)ranges.interval[0].upper, range,
            (int)ranges.interval[1].lower, (int)ranges.interval[1].upper);
}
