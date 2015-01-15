//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldGroup - A group of links that
// can be potentially ordered into a scaffold 
//
#ifndef SCAFFOLDGROUP_H
#define SCAFFOLDGROUP_H

#include "ScaffoldVertex.h"
#include "ScaffoldLink.h"


struct LinkVertexPair
{
    ScaffoldLink link;
    ScaffoldVertex* pEndpoint;
};

typedef std::list<LinkVertexPair> LinkList;
typedef LinkList::iterator LinkListIterator;
typedef std::vector<ScaffoldLink> LinkVector;
typedef std::vector<LinkVertexPair> LinkPairVector;
typedef LinkPairVector::iterator LinkVectorPairIterator;
typedef LinkPairVector::const_iterator LinkVectorPairConstIterator;

class ScaffoldGroup
{
    public:
        ScaffoldGroup(const ScaffoldVertex* pRootVertex, int maxOverlap);

        void addLink(const ScaffoldLink& link, ScaffoldVertex* pVertex);

        bool isOrderAmbiguous();
        bool markPolymorphic(double p_cutoff, double cn_cutoff);

        bool hasConsistentLayout();

        // Is the ordering between linkA and linkB ambiguous?
        // Returns true if the probabilty that the most likely
        // ordering is incorrect is greater than p
        bool areLinksAmbiguous(const ScaffoldLink& linkA,
                               const ScaffoldLink& linkB,
                               double p);

        void getSecondaryLinks();

        void computeBestOrdering();
        std::string getBestOrderingString() const;

        int calculateLongestOverlap();
        double calculateProbACloserThanB(const ScaffoldLink& linkA,
                                         const ScaffoldLink& linkB);

        
        // Construct a set of links between successive elements
        // of the ordered scaffold
        void getLinearLinks(LinkVector& outLinks);

    private:

        double normCDF(double x, double mean, double sd);
        int scoreLinkPlacement(const ScaffoldLink& link, const LinkList& unplacedList);

        const ScaffoldVertex* m_pRootVertex;
        int m_maxOverlap;
        LinkPairVector m_links;

        bool m_isOrdered;

};

#endif
