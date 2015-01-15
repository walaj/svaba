//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Match - A pair of coordinates representing the overlapping
// regions of two sequences
//
#include "Match.h"
#include <iostream>

// Constructor
Match::Match(const SeqCoord& sc1, const SeqCoord& sc2, bool isRC, int nd)
{
    coord[0] = sc1;
    coord[1] = sc2;
    isReverse = isRC;
    numDiff = nd;
}

// Constructor
Match::Match(int s1, int e1, int l1, int s2, int e2, int l2, bool isRC, int nd)
{
    coord[0] = SeqCoord(s1, e1, l1);
    coord[1] = SeqCoord(s2, e2, l2);
    isReverse = isRC;
    numDiff = nd;
}

// Swap the order of the elements
void Match::swap()
{
    SeqCoord tmp = coord[0];
    coord[0] = coord[1];
    coord[1] = tmp;
}

// Flip coord2 so that it is in the same frame as coord1
void Match::canonize()
{
    if(isRC())
        coord[1].flip();
    isReverse = false;
}

// Flip coord2 so it is out of frame of coord1
void Match::decanonize()
{
    if(!isRC())
        coord[1].flip();
    isReverse = true;
}

// Calculation the translation offset to shift
// a coord[0] position to a coord[1]. This must be calculated
// using canonical coordinates
int Match::calculateTranslation() const
{
    if(!isRC())
        return coord[1].interval.start - coord[0].interval.start;
    else
    {
        SeqCoord f = coord[1];
        f.flip();
        return f.interval.start - coord[0].interval.start;
    }
}

// Calculation the translation offset to shift
// a coord[1] position to a coord[0]. This must be calculated
// using canonical coordinates
int Match::calculateInverseTranslation() const
{
    if(!isRC())
        return coord[0].interval.start - coord[1].interval.start;
    else
    {
        SeqCoord f = coord[0];
        f.flip();
        return f.interval.start - coord[1].interval.start;
    }
}


// Translate the SeqCoord c from the frame of coord[0] to coord[1]
SeqCoord Match::translate(const SeqCoord& c) const
{
    assert(coord[0].length() == coord[1].length()); // ensure translation is valid
    int t = calculateTranslation();

    SeqCoord out;
    out.seqlen = coord[1].seqlen;
    out.interval.start = c.interval.start + t;
    out.interval.end = c.interval.end + t;
    if(isRC())
        out.flip();
    
    return out;
}

// Translate the SeqCoord c from the frame of coord[1] to coord[0]
SeqCoord Match::inverseTranslate(const SeqCoord& c) const
{
    assert(coord[0].length() == coord[1].length()); // ensure translation is valid
    int t = calculateInverseTranslation();
    
    SeqCoord out;
    out.seqlen = coord[0].seqlen;
    out.interval.start = c.interval.start + t;
    out.interval.end = c.interval.end + t;

    if(isRC())
        out.flip();

    return out;
}

// Translate a single position from c[0] frame to c[1]
int Match::translate(int c) const
{
    assert(!isRC());
    assert(coord[0].length() == coord[1].length()); // ensure translation is valid
    int t = calculateTranslation();
    return c + t;
}

// Translate a single position from c[1] frame to c[0]
int Match::inverseTranslate(int c) const
{
    assert(!isRC());
    assert(coord[0].length() == coord[1].length());
    int t = calculateInverseTranslation();
    return c + t;
}

// Given two matches, match_xy and match_xz infer the 
// match between yz. This requies coord[0] of the input 
// matches to be referring to the same sequence
// This returns the minimal matching region, it could possibly be extended
Match Match::infer(const Match& match_xy, const Match& match_xz)
{
    assert(match_xy.coord[0].seqlen == match_xz.coord[0].seqlen);

    // Calculate the max/min start/end coordinates of coord[0]
    int s = std::max(match_xy.coord[0].interval.start, match_xz.coord[0].interval.start);
    int e = std::min(match_xy.coord[0].interval.end, match_xz.coord[0].interval.end);

    // These are the coordinates in frame X
    SeqCoord r_y(s, e, match_xy.coord[0].seqlen);
    SeqCoord r_z(s, e, match_xz.coord[0].seqlen);

    // Translate into the desired frame
    SeqCoord t_y = match_xy.translate(r_y);
    SeqCoord t_z = match_xz.translate(r_z);

    // Set the new lengths of the seqcoords
    t_y.seqlen = match_xy.coord[1].seqlen;
    t_z.seqlen = match_xz.coord[1].seqlen;

    return Match(t_y, t_z, match_xy.isRC() != match_xz.isRC(), -1);
}

// Given two matches, match_xy and match_xz, return
// true if the match coordinates intersect. Assumes
// that the matches have a common basis which is coord[0]
bool Match::doMatchesIntersect(const Match& match_xy, const Match& match_xz)
{
    int s = std::max(match_xy.coord[0].interval.start, match_xz.coord[0].interval.start);
    int e = std::min(match_xy.coord[0].interval.end, match_xz.coord[0].interval.end);
    if(s > e)
        return false;
    else
        return true;    
}

// Expand the match outwards so one sequence is left terminal and one sequence
// is right terminal. This makes it a "proper" overlap
// This function assumes the coordinates are valid (see SeqCoord)
void Match::expand()
{
    assert(coord[0].isValid() && coord[1].isValid());

    // This is simpler if the coordinates are in canonical form, so here were canonize
    // them and decanonize after
    bool flipped = false;
    if(isRC())
    {
        flipped = true;
        canonize();
    }

    // left expansion
    if(!coord[0].isLeftExtreme() || !coord[1].isLeftExtreme())
    {
        int dist = std::min(coord[0].getLeftDist(), coord[1].getLeftDist());
        coord[0].interval.start -= dist;
        coord[1].interval.start -= dist;
    }

    // right expansion
    if(!coord[0].isRightExtreme() || !coord[1].isRightExtreme())
    {
        int dist = std::min(coord[0].getRightDist(), coord[1].getRightDist());
        coord[0].interval.end += dist;
        coord[1].interval.end += dist;
    }

    if(flipped)
        decanonize();
    assert(coord[0].isValid() && coord[1].isValid());
}

// Count the number of differences between the two strings based on the match
int Match::countDifferences(const std::string& s1, const std::string& s2) const
{
    std::string matched1 = coord[0].getSubstring(s1);
    std::string matched2 = coord[1].getSubstring(s2);

    if(isReverse)
        matched2 = reverseComplement(matched2);

    return ::countDifferences(matched1, matched2, matched1.length());
}

bool Match::isSubstringContainment() const
{ 
    return isContainment() && (!coord[0].isContained() || !coord[1].isContained()); 
}


// Print the matched strings
void Match::printMatch(const std::string& s1, const std::string& s2) const
{
    std::string matched1 = coord[0].getSubstring(s1);
    std::string matched2 = coord[1].getSubstring(s2);

    if(isReverse)
        matched2 = reverseComplement(matched2);

    std::cout << "M1: " << matched1 << "\n";
    std::cout << "M2: " << matched2 << "\n";
}

// Output
std::ostream& operator<<(std::ostream& out, const Match& m)
{
    out << m.coord[0] << " " << m.coord[1] << " " << m.isReverse << " " << m.numDiff;
    return out;
}

// Input
std::istream& operator>>(std::istream& in, Match& m)
{
    in >> m.coord[0] >> m.coord[1] >> m.isReverse >> m.numDiff;
    return in;
}

//
// Overlap
//
Overlap::Overlap(const std::string& i1, const std::string& i2, const Match& m) : match(m)
{
    id[0] = i1;
    id[1] = i2;
}

Overlap::Overlap(const std::string& i1, const SeqCoord& sc1, 
                 const std::string& i2, const SeqCoord& sc2, bool isRC, int nd)
{
    id[0] = i1;
    id[1] = i2;
    match = Match(sc1, sc2, isRC, nd);
}

Overlap::Overlap(const std::string& i1, int s1, int e1, int l1,
                 const std::string& i2, int s2, int e2, int l2, bool isRC, int nd)
{
    id[0] = i1;
    id[1] = i2;
    match = Match(s1, e1, l1, s2, e2, l2, isRC, nd);
}

// Return the index of the contained vertex
// Precondition: the overlap is a containment
size_t Overlap::getContainedIdx() const
{
    // The verts are mutually contained, return the lexographically lower id
    if(match.coord[0].isContained() && match.coord[1].isContained())
    {
        if(id[0] < id[1])
            return 1;
        else
            return 0;
    }
    else if(match.coord[0].isContained())
    {
        return 0;
    }
    else
    {
        assert(match.coord[1].isContained());
        return 1;
    }
}

// Return the ID of the contained vertex
// Precondition: The overlap is a containment
std::string Overlap::getContainedID() const
{
    return id[getContainedIdx()];
}

// Swap the elements
void Overlap::swap()
{
    id[0].swap(id[1]);
    match.swap();
}

// 
int Overlap::getOverlapLength(size_t idx) const
{
    assert(idx <= 1);
    return match.coord[idx].length();
}

// Output
std::ostream& operator<<(std::ostream& out, const Overlap& o)
{
    out << o.id[0] << " " << o.id[1] << " " << o.match;
    return out;
}

// Input
std::istream& operator>>(std::istream& in, Overlap& o)
{
    in >> o.id[0] >> o.id[1] >> o.match;
    return in;
}
