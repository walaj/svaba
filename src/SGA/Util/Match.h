//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Match - A pair of coordinates representing the overlapping
// regions of two sequences
//
#ifndef MATCH_H
#define MATCH_H

#include "Util.h"
#include "SeqCoord.h"

struct Match
{
    // functions
    Match() {}
    Match(const SeqCoord& sc1, const SeqCoord& sc2, bool isRC, int nd);
    Match(int s1, int e1, int l1, int s2, int e2, int l2, bool isRC, int nd);

    // Accessors
    inline bool isRC() const { return isReverse; }
    inline int getNumDiffs() const { return numDiff; }
    inline int getMinOverlapLength() const { return std::min(coord[0].length(), coord[1].length()); }
    inline int getMaxOverlapLength() const { return std::max(coord[0].length(), coord[1].length()); }

    void setNumDiffs(int n) { numDiff = n; }
    
    // Calculate the translation offset from coord[0] to coord[1]
    int calculateTranslation() const;
    int calculateInverseTranslation() const;

    // Translate the SeqCoord c from the frame of coord[0] to coord[1]
    SeqCoord translate(const SeqCoord& c) const;

    // Translate the SeqCoord c from the frame of coord[1] to coord[0]
    SeqCoord inverseTranslate(const SeqCoord& c) const;

    // Translate a single position from c[0] frame to c[1]
    int translate(int c) const;

    // Translate a single position from c[1] frame to c[0]
    int inverseTranslate(int c) const;

    // Swap the coords of this element
    void swap();

    // Infer an overlap between yz given matches xy and xz
    static Match infer(const Match& match_xy, const Match& match_xz);
    static bool doMatchesIntersect(const Match& match_xy, const Match& match_xz);

    // Expand a match outwards so each end is terminal for both coordinates
    void expand();

    // Count the number of differences between the strings assuming no indels
    int countDifferences(const std::string& s1, const std::string& s2) const;

    // Print the matched strings
    void printMatch(const std::string& s1, const std::string& s2) const;

    // Flip coord[1] if isReverse is true, effectively
    // bringing the matching strings into the same coordinate system
    void canonize();
    void decanonize();

    // Returns true if one overlap is contained within the other, ie they are not equal
    bool isSubstringContainment() const; 
    bool isContainment() const { return coord[0].isContained() || coord[1].isContained(); }

    // IO
    friend std::ostream& operator<<(std::ostream& out, const Match& m);
    friend std::istream& operator>>(std::istream& in, Match& m);

    // data
    SeqCoord coord[2];
    bool isReverse;
    int numDiff;
};

// Overlap
struct Overlap
{
    // constructors
    Overlap() {}
    Overlap(const std::string& i1, const std::string& i2, const Match& m);

    Overlap(const std::string& i1, const SeqCoord& sc1, 
            const std::string& i2, const SeqCoord& sc2, bool isRC, int nd); 

    Overlap(const std::string& i1, int s1, int e1, int l1, 
            const std::string& i2, int s2, int e2, int l2, bool isRC, int nd); 

    // functions
    
    // Swap the order of the elements
    void swap();

    //
    bool isContainment() const { return match.isContainment(); }
    bool isSubstringContainment() const { return match.isSubstringContainment(); }
    // Return the index (0 or 1) of the CONTAINED vertex (the discarded vertex of a containment)
    size_t getContainedIdx() const;
    std::string getContainedID() const;

    //
    int getOverlapLength(size_t idx) const;

    friend std::ostream& operator<<(std::ostream& out, const Overlap& o);
    friend std::istream& operator>>(std::istream& in, Overlap& o);

    // data
    std::string id[2];
    Match match;
};

// Collections
typedef std::vector<Overlap> OverlapVector;

#endif
