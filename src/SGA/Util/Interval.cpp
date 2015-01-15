//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Interval - A pair of integers denoting a range
//
#include "Interval.h"
#if 0
Interval Interval::intersect(const Interval& r1, const Interval& r2)
{
    Interval result;
    result.start = std::max(r1.start, r2.start);
    result.end = std::min(r1.end, r2.end);

    // Check for non-overlap
    if(result.end <= result.start)
    {
        result.start = 0;
        result.end = 0;
    }
    return result;
}
#endif

std::ostream& operator<<(std::ostream& out, const Interval& r)
{
    out << r.start << " " << r.end;
    return out;
}

std::istream& operator>>(std::istream& in, Interval& r)
{
    in >> r.start >> r.end;
    return in;
}
