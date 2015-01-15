//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldGroup - A set of links that
// can be potentially ordered into a scaffold 
//
#include "ScaffoldGroup.h"
#include "Interval.h"

//
ScaffoldGroup::ScaffoldGroup(const ScaffoldVertex* pRootVertex, 
                             int maxOverlap) : m_pRootVertex(pRootVertex),
                                               m_maxOverlap(maxOverlap),
                                               m_isOrdered(false)
{

}

//
void ScaffoldGroup::addLink(const ScaffoldLink& link, ScaffoldVertex* pVertex)
{
    LinkVertexPair pair = {link, pVertex};
    m_links.push_back(pair);
}

// 
bool ScaffoldGroup::isOrderAmbiguous()
{
    double ambiguity_p = 0.01f;
    LinkVectorPairIterator i = m_links.begin();
    for(; i != m_links.end(); ++i)
    {
        LinkVectorPairIterator j = i + 1;
        for(; j != m_links.end(); ++j)
        {
            bool isAmbiguous = areLinksAmbiguous(i->link, j->link, ambiguity_p);
            if(isAmbiguous)
            {
                return true;
            }
        }
    }
    return false;
}

// 
bool ScaffoldGroup::markPolymorphic(double p_cutoff, double cn_cutoff)
{
    //std::cout << "Checking " << m_pRootVertex->getID() << " for polymorphic links\n";
    LinkVectorPairIterator i = m_links.begin();
    for(; i != m_links.end(); ++i)
    {
        LinkVectorPairIterator j = i + 1;
        for(; j != m_links.end(); ++j)
        {
            bool isAmbiguous = areLinksAmbiguous(i->link, j->link, p_cutoff);
            if(isAmbiguous)
            {
                double sum = i->pEndpoint->getEstCopyNumber() + j->pEndpoint->getEstCopyNumber();
                /*
                std::cout << "\tLinks " << i->link << " and " << j->link << " are ambiguous\n";
                std::cout << "\tECN I: "<< i->pEndpoint->getEstCopyNumber() << " ECN J: " << j->pEndpoint->getEstCopyNumber() << "\n";
                std::cout << "\tsum: " << sum << "\n";
                */
                if(sum < cn_cutoff)
                {
                    if(i->pEndpoint->getEstCopyNumber() < j->pEndpoint->getEstCopyNumber())
                        i->pEndpoint->setClassification(SVC_POLYMORPHIC);
                    else
                        j->pEndpoint->setClassification(SVC_POLYMORPHIC);
                    return true;
                }
            }
        }
    }
    
    // No nodes marked as polymorphic
    return false;
}

//
bool ScaffoldGroup::areLinksAmbiguous(const ScaffoldLink& linkA,
                                      const ScaffoldLink& linkB,
                                      double p)
{
    // Calculate the probabilty A comes before B and vice-versa
    double p_AB = calculateProbACloserThanB(linkA, linkB);
    double p_BA = 1.0f - p_AB;
    double best = std::max(p_AB, p_BA);
    double p_wrong = 1.0f - best;
    return p_wrong > p;
}

// Return true if the scaffoldgroup has a consistent layout
bool ScaffoldGroup::hasConsistentLayout()
{
    return true;
}

// Calculate the longest overlap between any
// pair of nodes in the scaffold
// If two nodes are significantly overlapped, it
// may indicate that the root node is a repeat or connected
// to polymorphic nodes
int ScaffoldGroup::calculateLongestOverlap()
{
    int longestOverlap = 0;
    LinkVectorPairIterator i = m_links.begin();
    for(; i != m_links.end(); ++i)
    {
        // Compute interval for this link using closed coordionates [start, end]
        Interval interval_i(i->link.distance, i->link.getEndpoint() - 1);

        LinkVectorPairIterator j = i + 1;
        for(; j != m_links.end(); ++j)
        {
            Interval interval_j(j->link.distance, j->link.getEndpoint() - 1);
            if(Interval::isIntersecting(interval_i.start, interval_i.end,
                                        interval_j.start, interval_j.end)) {
                Interval intersection;
                Interval::intersect(interval_i.start, interval_i.end,
                                    interval_j.start, interval_j.end,
                                    intersection.start, intersection.end);

                int overlap = intersection.end - intersection.start + 1;
                if(overlap > longestOverlap)
                    longestOverlap = overlap;
            }
        }
    }

    return longestOverlap;
}

//
void ScaffoldGroup::computeBestOrdering()
{
    // We compute the best ordering of the links for a given vertex
    // with a greedy algorithm as follows. Initially, all links
    // are unplaced. At each step, we test and score each unplaced
    // link as a candidate for being the next-closest contig
    // in the scaffold to the root. The lowest scoring contig
    // is removed from the unplaced list. This continues
    // until all contigs have been placed in the putative scaffold.
    LinkList unplacedLinks(m_links.begin(), m_links.end());
    LinkPairVector placedLinks;

    int totalScore = 0;
    while(!unplacedLinks.empty())
    {
        int bestScore = std::numeric_limits<int>::max();
        LinkListIterator bestLink = unplacedLinks.end();

        for(LinkListIterator iter = unplacedLinks.begin();
                             iter != unplacedLinks.end();
                             ++iter)
        {
            int score = scoreLinkPlacement(iter->link, unplacedLinks);
            if(score < bestScore)
            {
                bestScore = score;
                bestLink = iter;
            }
        }

        assert(bestLink != unplacedLinks.end());

        placedLinks.push_back(*bestLink);
        unplacedLinks.erase(bestLink);
        totalScore += bestScore;
    }

    m_links = placedLinks;
    m_isOrdered = true;
}

// Get a simple representation of the best ordering as a string
std::string ScaffoldGroup::getBestOrderingString() const
{
    assert(m_isOrdered);
    std::stringstream ss;
    for(LinkVectorPairConstIterator iter = m_links.begin(); iter != m_links.end(); ++iter)
    {
        ss << iter->link.endpointID << ":" << iter->link.distance << "-" << iter->link.getEndpoint() << " ";
    }
    return ss.str();
}

// Construct a set of links between successive elements of the
// ordered group
void ScaffoldGroup::getLinearLinks(LinkVector& outLinks)
{
    assert(m_isOrdered);

    if(m_links.empty())
        return;


    // The first link is added without modification
    LinkPairVector::const_iterator iter = m_links.begin();
    ScaffoldLink previousLink = iter->link;
    outLinks.push_back(previousLink);
    ++iter;
    for(; iter != m_links.end(); ++iter)
    {
        // Infer the distance and orientation of the current contig in the scaffold
        // with the previous
        const ScaffoldLink& currentLink = iter->link;
        EdgeComp orientation = previousLink.getComp() == currentLink.getComp() ? EC_SAME : EC_REVERSE;

        EdgeDir dir;
        if(previousLink.getComp() == EC_SAME)
            dir = previousLink.getDir();
        else
            dir = !previousLink.getDir();

        int distance = currentLink.distance - previousLink.getEndpoint();
        double sd = sqrt(pow(currentLink.stdDev, 2.0) + pow(previousLink.stdDev, 2.0));
        ScaffoldLink inferred(currentLink.endpointID, dir, orientation, distance, sd, 0, currentLink.seqLen, SLT_INFERRED);
        outLinks.push_back(inferred);
        previousLink = currentLink;
    }
}

// Attempt to place contigs in the gaps of the current scaffold using secondary links
void ScaffoldGroup::getSecondaryLinks()
{
    assert(m_isOrdered);

    if(m_links.empty())
        return;

    LinkPairVector::const_iterator iter = m_links.begin();
    for(; iter != m_links.end(); ++iter)
    {
        ScaffoldLink xyLink = iter->link;
        ScaffoldVertex* pY = iter->pEndpoint;
        ScaffoldEdgePtrVector yEdges = pY->getEdges();

        // Infer a link to this contig

        // Get the direction from contig y back to contig x (the root)
        int xyDistance = xyLink.distance;
        EdgeDir yxDir = xyLink.getTwinDir();
        EdgeComp yxComp = xyLink.getComp();
        for(size_t i = 0; i < yEdges.size(); ++i)
        {
            ScaffoldEdge* pYZ = yEdges[i];
            ScaffoldVertex* pZ = pYZ->getEnd();
            EdgeDir yzDir = pYZ->getDir();
            EdgeComp yzComp = pYZ->getComp();
            int yzDistance = pYZ->getDistance();

            // Calculate the distance between x and z
            // Two possible orientations of the three contigs:
            //  ---x---  ----z---- ----y---- (case 1)
            //  ---x---  ----y---- ----z---- (case 2)
            // If the direction y->z is the same as y -> x, its case 1
            // otherwise, case 2
            int xzDistance;
            if(yzDir == yxDir)
            {
                // In same direction as x
                xzDistance = xyDistance - (pZ->getSeqLen() + yzDistance);
                printf("SAME xy %d - (%d + %d) = %d\n", xyDistance, (int)pZ->getSeqLen(), yzDistance, xzDistance);
            }
            else
            {
                // opposite direction of x
                xzDistance = xyDistance + pY->getSeqLen() + yzDistance;
                printf("DIFF xy %d + %d + %d = %d\n", xyDistance, (int)pY->getSeqLen(), yzDistance, xzDistance);
            }

            std::cout << "yzLink: " << pYZ->getLink() << "\n";

            WARN_ONCE("Recalculate xzDir");
            EdgeDir xzDir = xyLink.getDir();
            double xzSD = sqrt(pow(xyLink.stdDev, 2.0) + pow(pYZ->getStdDev(), 2.0));
            EdgeComp xzComp = (yxComp == yzComp) ? EC_SAME : EC_REVERSE;
            ScaffoldLink xzLink(pZ->getID(), xzDir, xzComp, xzDistance, xzSD, 0, pZ->getSeqLen(), SLT_INFERRED);
            std::cout << "Inferred secondary edge: " << xzLink << "\n";
        }
    }
}

// Calculate the score of placing the given link before
// all the other links in the unplaced list
int ScaffoldGroup::scoreLinkPlacement(const ScaffoldLink& link,
                                      const LinkList& unplacedList)
{
    // Set the range of coordinates
    int end = link.getEndpoint();

    int score = 0;
    for(LinkList::const_iterator iter = unplacedList.begin();
                                 iter != unplacedList.end();
                                 ++iter)
    {
        if(link.endpointID != iter->link.endpointID)
        {
            // Compute the offset of this link based on input link preceding it
            // If the interval [start,end] includes the distance estimate
            // of the unplaced link, shift the estimate to the first valid
            // position after the placed contig (given by end - max overlap)
            int next_valid_position = end - m_maxOverlap - 1;
            if(iter->link.distance < next_valid_position)
            {
                score += next_valid_position - iter->link.distance;
            }
        }
    }
    return score;
}

// Given the links A and B, calculate the probabilty that A
// comes before B in the scaffold.
// In other words, if D_A and D_B are normal random variables, 
// what is the probability that D_B < D_A. The difference
// D_B - D_A is a random variable with mean u_A - u_B and variance
// v_A + v_B
double ScaffoldGroup::calculateProbACloserThanB(const ScaffoldLink& linkA,
                                                const ScaffoldLink& linkB)
{
    // calculate the mean and variance of the normal rv B - A
    double mean = linkA.distance - linkB.distance;
    double variance = pow(linkA.stdDev, 2.0f) + pow(linkB.stdDev, 2.0f);

    // The probabilty that A is closer than B is given by probabilty that
    // the random variable has a value less than zero.
    return normCDF(0.0f, mean, sqrt(variance));
}

double ScaffoldGroup::normCDF(double x, double mean, double sd)
{
    double t = (x - mean) / (sd * sqrt(2.0f));
    return 0.5f * (1 + erf(t)); 
}
