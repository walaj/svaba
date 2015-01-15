//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Interval - A pair of integers denoting a closed interval
//

#ifndef INTERVAL_H
#define INTERVAL_H

#include "Util.h"

struct Interval
{
    // constructors
    Interval() : start(0), end(0) {}
    Interval(int s, int e) : start(s), end(e) {}

    // Find the intersection between [s1, e1] [s2, e2]
    // The intersection is placed in [s3, e3], if the two input
    // intervals do not intersect then e3 will be strictly less than s3
    template<class T>
    static void intersect(const T s1, const T e1, const T s2, const T e2, T& s3, T& e3)
    {
        assert(s1 <= e1 && s2 <= e2);
        s3 = std::max(s1, s2);
        e3 = std::min(e1, e2);
    }

    // Precondition: s1 >= e1 and s2 >= e2 
    // Return true if the coordinates intersect
    template<class T>
    static bool isIntersecting(const T s1, const T e1, const T s2, const T e2)
    {
        assert(s1 <= e1 && s2 <= e2);
        if(s2 > e1 || s1 > e2)
            return false;
        else
            return true;
    }

    // functions
    friend std::ostream& operator<<(std::ostream& out, const Interval& i);
    friend std::istream& operator>>(std::istream& in, Interval& i);    
    
    // data members
    int start;
    int end;
};

#endif
