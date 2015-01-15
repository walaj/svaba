//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// BWTInterval - Data structures for holding and manipulating
// the coordinates in a BWT/FM-index
// 
#ifndef BWTINTERVAL_H
#define BWTINTERVAL_H

#include <list>
#include <iostream>
#include <inttypes.h>

// A BWTInterval holds a pair of integers which delineate an alignment of some string
// to a BWT/Suffix Array
struct BWTInterval
{
    // Functions
    BWTInterval() : lower(0), upper(0) {}
    BWTInterval(int64_t l, int64_t u) : lower(l), upper(u) {}

    inline bool isValid() const { return lower <= upper; }
    inline int64_t size() const { return upper - lower + 1; }

    static inline bool compare(const BWTInterval& a, const BWTInterval& b)
    {
        if(a.lower == b.lower)
            return a.upper < b.upper;
        else
            return a.lower < b.lower;
    }

    static inline bool equal(const BWTInterval& a, const BWTInterval& b)
    {
        return a.lower == b.lower && a.upper == b.upper;
    }

    friend std::ostream& operator<<(std::ostream& out, const BWTInterval& a)
    {
        out << a.lower << " " << a.upper;
        return out;
    }

    friend std::istream& operator>>(std::istream& in, BWTInterval& a)
    {
        in >> a.lower >> a.upper;
        return in;
    }

    void write(std::ostream& out)
    {
        out.write((char*)&lower, sizeof(lower));
        out.write((char*)&upper, sizeof(upper));
    }

    void read(std::istream& in)
    {
        in.read((char*)&lower, sizeof(lower));
        in.read((char*)&upper, sizeof(upper));
    }

    // Data
    int64_t lower;
    int64_t upper;
};

// A pair of intervals, used for bidirectional searching a bwt/revbwt in lockstep
struct BWTIntervalPair
{
    // Functions
    BWTInterval& get(unsigned int idx) { return interval[idx]; }
    bool isValid() const { return interval[0].isValid() && interval[1].isValid(); }

    friend bool operator==(const BWTIntervalPair& a, const BWTIntervalPair& b)
    {
        return BWTInterval::equal(a.interval[0], b.interval[0]) && 
               BWTInterval::equal(a.interval[1], b.interval[1]);
    }

    // Sort an itnerval pair by interval[0]'s lower coordinate
    static bool sortFirstLower(const BWTIntervalPair& a, const BWTIntervalPair& b)
    {
        return a.interval[0].lower < b.interval[0].lower;
    }

    // Sort an itnerval pair by interval[1]'s lower coordinate
    static bool sortSecondLower(const BWTIntervalPair& a, const BWTIntervalPair& b)
    {
        return a.interval[1].lower < b.interval[1].lower;
    }    

    // I/O
    friend std::ostream& operator<<(std::ostream& out, const BWTIntervalPair& a)
    {
        out << a.interval[0] << " " << a.interval[1];
        return out;
    }

    friend std::istream& operator>>(std::istream& in, BWTIntervalPair& a)
    {
        in >> a.interval[0] >> a.interval[1];
        return in;
    }
    
    // Data
    BWTInterval interval[2];
};

#endif

