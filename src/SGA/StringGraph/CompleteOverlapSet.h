//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// CompleteOverlapSet - A collection of all
// valid overlaps for a given vertex. This is
// inferred from the structure of a transitive reduced
// graph
//
#ifndef COMPLETEOVERLAPSET_H
#define COMPLETEOVERLAPSET_H

#include "Bigraph.h"
#include "SGAlgorithms.h"

// Structure used for iterative exploration of the graph
// The exploration starts at some vertex X, each element
// holds a possible overlap between X and another vertex Y
struct ExploreElement
{
    ExploreElement(const EdgeDesc& e, const Overlap& o) : ed(e), ovr(o) {}
    EdgeDesc ed;
    Overlap ovr;
};

// Comparison operator used to compare ExploreElements
// by the length of the overlap on vertex X
struct CompareExploreElemOverlapLength
{
    bool operator()(const ExploreElement& elemXY, const ExploreElement& elemXZ)
    {
        return elemXY.ovr.match.coord[0].length() < elemXZ.ovr.match.coord[0].length();
    }
};

//
typedef std::priority_queue<ExploreElement, 
                            std::vector<ExploreElement>, 
                            CompareExploreElemOverlapLength> ExplorePriorityQueue;
typedef std::list<EdgeDesc> EdgeDescList;
typedef std::queue<ExploreElement> ExploreQueue;

class CompleteOverlapSet
{
    public:

        CompleteOverlapSet(const Vertex* pVertex, double maxER, int minLength);

        void getDiffMap(SGAlgorithms::EdgeDescOverlapMap& missingMap, SGAlgorithms::EdgeDescOverlapMap& extraMap);
        void removeOverlapsTo(Vertex* pRemove);

        // Remove all the transitive and containment relationships, leaving only the irreducible
        // If the input pointers are not NULL, the transitive and containment overlaps will be placed
        // in the appropriate map
        void computeIrreducible(SGAlgorithms::EdgeDescOverlapMap* pTransitive, 
                                SGAlgorithms::EdgeDescOverlapMap* pContainments);
        
        void constructMap();
        void recursiveConstruct(const EdgeDesc& edXY, const Overlap& ovrXY, int depth, int distance);

        // Partition the OverlapMap into edges that are containments, irreducible and transitive
        // If the pointer for an output map is NULL, simply discard the edges
        void partitionOverlaps(SGAlgorithms::EdgeDescOverlapMap* pIrreducible, 
                               SGAlgorithms::EdgeDescOverlapMap* pTransitive, 
                               SGAlgorithms::EdgeDescOverlapMap* pContainment) const;

        void resetParameters(double maxER, int minLength);
        SGAlgorithms::EdgeDescOverlapMap getOverlapMap() const { return m_overlapMap; }
        size_t size() const { return m_overlapMap.size(); }
        size_t getCost() const { return m_cost; }
    private:

        // functions
        void iterativeConstruct();
        void constructBFS();

        // data
        SGAlgorithms::EdgeDescOverlapMap m_overlapMap;

        // the vertex this set is centered on
        const Vertex* m_pX;

        double m_maxER;
        int m_minLength;
        mutable size_t m_cost;
};

#endif
