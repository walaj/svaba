//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldVisitors - Functors that perform
// some operation on a ScaffoldVertex/Graph
//
#include "ScaffoldVisitors.h"
#include "ScaffoldRecord.h"
#include "ScaffoldGroup.h"
#include "SGSearch.h"
#include "ScaffoldSearch.h"
#include <limits.h>

//
// Returns true if the distance estimates for the two edges
// are close enough that the ordering cannot be resolved. Assumes
// that pXY is closer than pXZ
bool ScaffoldAlgorithms::areEdgesAmbiguous(ScaffoldEdge* pXY, ScaffoldEdge* pXZ)
{
    const double AMBIGUOUS_TOLERANCE = 2.0f;
    assert(pXY->getDistance() <= pXZ->getDistance());
    int translated = pXZ->getDistance() - static_cast<int>(AMBIGUOUS_TOLERANCE*pXZ->getStdDev());
    if(translated < pXY->getDistance())
        return true;
    else
        return false;
}

//
// Returns the edge properties inferred from the edges X->Y and X->Z
void ScaffoldAlgorithms::inferScaffoldEdgeYZ(ScaffoldEdge* pXY, ScaffoldEdge* pXZ,
                                             int& dist, double& sd, EdgeDir& dir_yz, 
                                             EdgeDir& dir_zy, EdgeComp& comp)
{
    comp = (pXY->getComp() == pXZ->getComp() ? EC_SAME : EC_REVERSE);
    dir_yz = !pXY->getTwin()->getDir(); // Opposite direction of X <- Y
    dir_zy = pXZ->getTwin()->getDir();  // Same direction as X <- Z since Z is after Y
    dist = pXZ->getDistance() - (pXY->getDistance() + pXY->getEnd()->getSeqLen());
    sd = sqrt(pow(pXY->getStdDev(), 2.0) + pow(pXZ->getStdDev(), 2.0));
}

// ScaffoldStatsVisitor
void ScaffoldStatsVisitor::previsit(ScaffoldGraph* /*pGraph*/)
{
    m_numVertices = 0;
    m_numEdges = 0;
}

bool ScaffoldStatsVisitor::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    ++m_numVertices;
    m_numEdges += pVertex->getNumEdges();

    return false;
}

void ScaffoldStatsVisitor::postvisit(ScaffoldGraph* /*pGraph*/)
{
    printf("Scaffold Stats -- Num vertices: %zu Num edges: %zu\n", m_numVertices, m_numEdges);
}

// ScaffoldAStatistic
ScaffoldAStatisticVisitor::ScaffoldAStatisticVisitor(double uniqueThreshold, 
                                                     double minCopyNumber) : m_uniqueThreshold(uniqueThreshold),
                                                                             m_minCopyNumber(minCopyNumber)

{
                                                                               
}

//
void ScaffoldAStatisticVisitor::previsit(ScaffoldGraph* /*pGraph*/)
{
    m_sumUnique = 0;
    m_sumRepeat = 0;
    m_sumLowCN = 0;
    m_numUnique = 0;
    m_numRepeat = 0;
    m_numLowCN = 0;
}

//
bool ScaffoldAStatisticVisitor::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    // Never re-classify repeats
    if(pVertex->getClassification() != SVC_REPEAT)
    {
        if(pVertex->getAStatistic() > m_uniqueThreshold)
        {
            pVertex->setClassification(SVC_UNIQUE);
            ++m_numUnique;
            m_sumUnique += pVertex->getSeqLen();
        }
        else
        {
            pVertex->setClassification(SVC_REPEAT);
            ++m_numRepeat;
            m_sumRepeat += pVertex->getSeqLen();
        }

        if(pVertex->getEstCopyNumber() < m_minCopyNumber)
        {
            pVertex->setClassification(SVC_REPEAT);
            ++m_numLowCN;
            m_sumLowCN += pVertex->getSeqLen();
        }
    }
    return false;   
}

//
void ScaffoldAStatisticVisitor::postvisit(ScaffoldGraph* /*pGraph*/)
{
    printf("A-statistic: classified %zu vertices as unique (%.2lf Mbp)\n", m_numUnique, (double)m_sumUnique / 1000000);
    printf("A-statistic: classified %zu vertices as repeat (%.2lf Mbp)\n", m_numRepeat, (double)m_sumRepeat / 1000000);
    printf("A-statistic: classified %zu vertices as spurious (%.2lf Mbp)\n", m_numLowCN, (double)m_sumLowCN / 1000000);
}

//
ScaffoldLinkValidator::ScaffoldLinkValidator(int maxOverlap, 
                                             double threshold,
                                             int verbose) : m_maxOverlap(maxOverlap),
                                                            m_threshold(threshold),
                                                            m_verbose(verbose)
{

}

//
void ScaffoldLinkValidator::previsit(ScaffoldGraph* pGraph)
{
    m_numUnique = 0;
    m_numRepeat = 0;
    m_numCut = 0;
    pGraph->setEdgeColors(GC_WHITE);
}

//
bool ScaffoldLinkValidator::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    // Never re-classify repeats
    if(pVertex->getClassification() == SVC_REPEAT)
        return false;

    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() > 1)
        {
            ScaffoldGroup group(pVertex, m_maxOverlap);
            for(size_t i = 0; i < edgeVec.size(); ++i)
                group.addLink(edgeVec[i]->getLink(), edgeVec[i]->getEnd());
        
            bool isAmbiguous = group.isOrderAmbiguous();
            int longestOverlap = group.calculateLongestOverlap();
            group.computeBestOrdering();
            bool overlapCheck = longestOverlap < 400;
            bool result = overlapCheck;

            if(m_verbose >= 1)
            {
                std::string ambiStr = isAmbiguous ? "ambiguous" : "unambiguous";
                std::string overStr = overlapCheck ? "good-ordering" : "no-ordering";
                std::string resultStr = result ? "PASS" : "FAIL";
                std::string orderStr = group.getBestOrderingString();
                printf("LV %s %d CL:%d EC:%2.2lf AS:%2.2lf %s %s %s LO:%d BO:%s\n", pVertex->getID().c_str(), 
                                                                                  (int)idx, 
                                                                                  (int)pVertex->getSeqLen(), 
                                                                                  pVertex->getEstCopyNumber(),
                                                                                  pVertex->getAStatistic(),
                                                                                  ambiStr.c_str(), 
                                                                                  overStr.c_str(), 
                                                                                  resultStr.c_str(), 
                                                                                  longestOverlap,
                                                                                  orderStr.c_str());
            }

            // If the link did not pass validation, cut the scaffold at this point for all involved vertices
            if(!result)
            {
                for(size_t i = 0; i < edgeVec.size(); ++i)
                {
                    edgeVec[i]->setColor(GC_BLACK);
                    edgeVec[i]->getEnd()->markEdgesInDir(edgeVec[i]->getTwin()->getDir(), GC_BLACK);
                }
                m_numCut += 1;
            }
        }
    }
    return false;   
}

//
void ScaffoldLinkValidator::postvisit(ScaffoldGraph* pGraph)
{
    printf("ScaffoldLinkValidator removed edges for %zu vertices\n", m_numCut);
    pGraph->deleteEdgesByColor(GC_BLACK);
}

//
ScaffoldChainVisitor::ScaffoldChainVisitor(int maxOverlap) : m_maxOverlap(maxOverlap)
{
}

//
void ScaffoldChainVisitor::previsit(ScaffoldGraph* /*pGraph*/)
{
}

// 
bool ScaffoldChainVisitor::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    // Never try to make chains from a repeat
    if(pVertex->getClassification() == SVC_REPEAT)
        return false;

    bool changed_graph = false;
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() > 1)
        {
            // Try to resolve this vertex edge set
            // by creating edges between the vertices that it
            // is linked to
            
            // Process the edges in order of length
            std::sort(edgeVec.begin(), edgeVec.end(), ScaffoldEdgePtrDistanceCompare);
            ScaffoldEdge* pXY = edgeVec[0];
            ScaffoldEdge* pXZ = edgeVec[1];
            assert(pXY->getDistance() <= pXZ->getDistance());
            ScaffoldVertex* pY = pXY->getEnd();
            ScaffoldVertex* pZ = pXZ->getEnd();

            if(pY->isRepeat() || pZ->isRepeat())
            {
                std::cerr << "Warning, skipping repeat\n";
                continue;
            }


            std::cout << "CV " << pVertex->getID() << "\t XY: " << *pXY << " XZ: " << *pXZ << "\n";

            // Check if the ordering of the contigs is ambiguous
            bool isAmbiguous = ScaffoldAlgorithms::areEdgesAmbiguous(pXY, pXZ);
            if(isAmbiguous)
            {
                std::cout << "\tEdges " << *pXY << " and " << *pXZ << 
                              " are ambiguously ordered\n";
            }

            // Infer the edge data for the edge Y->Z
            int dist;
            double sd;
            EdgeDir dir_yz;
            EdgeDir dir_zy; 
            EdgeComp comp;
            ScaffoldAlgorithms::inferScaffoldEdgeYZ(pXY, pXZ, dist, sd, dir_yz, dir_zy, comp);

            // Sanity checks
            bool isConsistent = (dist > -1*m_maxOverlap);
            if(!isConsistent)
            {
                std::cout << "\tEdge is not consistent with max overlap: " << dist << "\n";
            }

            // Check if an edge between Y->Z already exists
            ScaffoldEdge* pCheckEdge = pY->findEdgeTo(pZ->getID(), dir_yz, comp);
            if(pCheckEdge == NULL)
            {
                // Create the new edges
                ScaffoldLink linkYZ(pZ->getID(), dir_yz, comp, dist, sd, 0, pZ->getSeqLen(), SLT_INFERRED);
                ScaffoldLink linkZY(pY->getID(), dir_zy, comp, dist, sd, 0, pY->getSeqLen(), SLT_INFERRED);
                ScaffoldEdge* pYZ = new ScaffoldEdge(pZ, linkYZ);
                ScaffoldEdge* pZY = new ScaffoldEdge(pY, linkZY);
                pYZ->setTwin(pZY);
                pZY->setTwin(pYZ);

                std::cout << "\tCreating edge: " << *pYZ << " " << *pZY << "\n";
                pY->addEdge(pYZ);
                pZ->addEdge(pZY);
            }
            else
            {
                std::cout << "\tEdge already exists: " << *pCheckEdge << " inferred dist: " << dist << "\n";
            }

            // Remove the edges pXZ and pZX
            ScaffoldEdge* pZX = pXZ->getTwin();

            std::cout << "\tDeleting edge: " << *pXZ << "\n";
            pVertex->deleteEdge(pXZ); // deallocates
            pZ->deleteEdge(pZX); // deallocates
            
            pXZ = 0;
            pZX = 0;
            changed_graph = true;
        }
    }

    return changed_graph;
}

//
void ScaffoldChainVisitor::postvisit(ScaffoldGraph* /*pGraph*/)
{

}

//
bool ScaffoldMultiEdgeRemoveVisitor::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() > 1)
        {
            pVertex->deleteEdgesAndTwins(dir);
        }
    }
    return false;
}

// 
// ScaffoldTransitiveReductionVisitor - Remove transitive edges from the graph
//
ScaffoldTransitiveReductionVisitor::ScaffoldTransitiveReductionVisitor()
{

}

//
void ScaffoldTransitiveReductionVisitor::previsit(ScaffoldGraph* pGraph)
{
    pGraph->setEdgeColors(GC_WHITE);
}

//
bool ScaffoldTransitiveReductionVisitor::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    // Never try to make chains from a repeat
    if(pVertex->getClassification() == SVC_REPEAT)
        return false;

    bool changed_graph = false;
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() <= 1)
            continue;

        // Try to make a transitive set out of these edges
        // Search the graph for a variation walk that contains a common
        // endpoint for all links
        ScaffoldWalkVector walkVector;
        ScaffoldSearch::findVariantWalks(pVertex, dir, 1000000, 100, walkVector);

        // Search the walks to see if one contains all the links.
        int walkIdx = -1;
        int lowestIdxInWalk = std::numeric_limits<int>::max();
        int lowestIdxInVec = -1;
        for(size_t i = 0; i < walkVector.size(); ++i)
        {
            bool allContained = true;
            lowestIdxInWalk = std::numeric_limits<int>::max();
            lowestIdxInVec = -1;

            for(size_t j = 0; j < edgeVec.size(); ++j)
            {
                EdgeComp expected_orientation = edgeVec[j]->getComp();
                int idxInWalk = walkVector[i].findVertex(edgeVec[j]->getEnd());

                if(idxInWalk == -1)
                {
                    // vertex not found
                    allContained = false;
                    break;
                }

                // The vertex is in the walk, check that the orientation in
                // the walk matches the expected orientation
                EdgeComp walk_orientation = walkVector[i].findOrientation(edgeVec[j]->getEnd());
                if(walk_orientation != expected_orientation)
                {
                    // Invalid solution
                    allContained = false;
                    break;
                }
                
                // This is a valid solution
                if(idxInWalk < lowestIdxInWalk)
                {
                    lowestIdxInWalk = idxInWalk;
                    lowestIdxInVec = j;
                }
            }

            if(allContained)
            {
                walkIdx = i;
                break;
            }
        }

        if(walkIdx == -1)
        {
            // no path found
            continue;
        }

        /*
        std::cout << "Walk " << walkIdx << " contains all links ";
        walkVector[walkIdx].print();
        std::cout << "keeping link: " << lowestIdxInVec << " " << edgeVec[lowestIdxInVec]->getLink() << "\n";
        */

        // Keep the first link in the walk and discard all others
        // Mark the links that should not be deleted
        for(int i = 0; i < (int)edgeVec.size(); ++i)
        {
            if(i != lowestIdxInVec)
            {
                edgeVec[i]->setColor(GC_BLACK);
                edgeVec[i]->getTwin()->setColor(GC_BLACK);
                changed_graph = true;
            }
        }
    }
    return changed_graph;
}

//
void ScaffoldTransitiveReductionVisitor::postvisit(ScaffoldGraph* pGraph)
{
    pGraph->deleteEdgesByColor(GC_BLACK);
}

//
// ScaffoldPolymorphismVisitor - Detect and remove polymorpism in the scaffold
// graph
//
ScaffoldPolymorphismVisitor::ScaffoldPolymorphismVisitor(int maxOverlap) : m_maxOverlap(maxOverlap)
{

}

//
void ScaffoldPolymorphismVisitor::previsit(ScaffoldGraph* /*pGraph*/)
{
    m_numMarked = 0;
}

//
bool ScaffoldPolymorphismVisitor::visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex)
{
    (void)pGraph;
    (void)pVertex;
    if(pVertex->getClassification() == SVC_REPEAT)
        return false;

    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() > 1)
        {
            ScaffoldGroup group(pVertex, m_maxOverlap);

            for(size_t i = 0; i < edgeVec.size(); ++i)
            {
                group.addLink(edgeVec[i]->getLink(), edgeVec[i]->getEnd());
            }

            // Check if any links to this scaffold are polymorphic
            // This is done be seeing if there are a pair of primary links
            // that have an ambiguous ordering (the probability they are 
            // incorrectly ordered by the distance estimates is greater than
            // p_cutoff). If they are ambiguous and the sum of their 
            // estimated copy number is less than cn_cutoff, we mark
            // the lower-covered node as polymorphic.
            double p_cutoff = 0.01;
            double cn_cutoff = 1.5f;
            bool nodeMarked = group.markPolymorphic(p_cutoff, cn_cutoff);
            if(nodeMarked)
                m_numMarked++;
            return nodeMarked;
        }
    }    
    return false;
}

//
void ScaffoldPolymorphismVisitor::postvisit(ScaffoldGraph* pGraph)
{
    printf("Marked %d nodes as polymorphic\n", m_numMarked);
    pGraph->deleteVertices(SVC_POLYMORPHIC);
}


// 
// ScaffoldTransitiveReductionVisitor - Detect and resolve structural variation in the scaffold graph
//
ScaffoldSVVisitor::ScaffoldSVVisitor(int maxSize) : m_maxSVSize(maxSize)
{

}

//
void ScaffoldSVVisitor::previsit(ScaffoldGraph* pGraph)
{
    pGraph->setEdgeColors(GC_WHITE);
    m_numMarked = 0;
}

//
bool ScaffoldSVVisitor::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    // Never try to make chains from a repeat
    if(pVertex->getClassification() == SVC_REPEAT)
        return false;

    bool changed_graph = false;
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() <= 1)
            continue;

        // Try to make a transitive set out of these edges
        // Search the graph for a variation walk that contains a common
        // endpoint for all links
        ScaffoldWalkVector walkVector;
        ScaffoldSearch::findVariantWalks(pVertex, dir, 50000, 100, walkVector);

        // Search the walks to see if one contains all the links.
        ScaffoldVertexPtrVector linkedVertices;
        for(size_t i = 0; i < edgeVec.size(); ++i)
            linkedVertices.push_back(edgeVec[i]->getEnd());

        int coveringWalkIdx = ScaffoldSearch::findCoveringWalk(walkVector, linkedVertices);
        if(coveringWalkIdx == -1)
        {
            // no valid found
            continue;
        }

        // Keep the first link in the walk and discard all others
        //std::cout << "Walk " << coveringWalkIdx << " contains all links ";
        //walkVector[coveringWalkIdx].print();
        
        // Calculate the difference between the distance in the walk and the distance
        // estimate for each linked vertex
        int maxDiff = 0;
        int closestIdx = 0;
        int closestDistance = std::numeric_limits<int>::max();
        for(int i = 0; i < (int)edgeVec.size(); ++i)
        {
            int distance = walkVector[coveringWalkIdx].getDistanceToVertex(edgeVec[i]->getEnd());
            int diff = abs(distance - edgeVec[i]->getDistance());
          
            /*
            printf("SVW: %s -> %s DE: %d WD: %d DF: %d TL: %.2lf\n", edgeVec[i]->getStartID().c_str(), 
                                                                     edgeVec[i]->getEndID().c_str(),
                                                                     edgeVec[i]->getDistance(), distance, 
                                                                     diff, 
                                                                     3*edgeVec[i]->getStdDev());
            */

            if(diff > maxDiff)
                maxDiff = diff;

            if(distance < closestDistance)
            {
                closestDistance = distance;
                closestIdx = i;
            }
        }

        //printf("SVW: %s MD: %d OR: %d\n", pVertex->getID().c_str(), maxDiff, hasOOR);
        if(maxDiff > m_maxSVSize)
            continue;

        // Mark the edges for removal that are not the closest
        for(int i = 0; i < (int)edgeVec.size(); ++i)
        {
            if(i != closestIdx)
            {
                edgeVec[i]->setColor(GC_BLACK);
                edgeVec[i]->getTwin()->setColor(GC_BLACK);
                m_numMarked += 1;
            }
        }
    }
    return changed_graph;
}

//
void ScaffoldSVVisitor::postvisit(ScaffoldGraph* pGraph)
{
    pGraph->deleteEdgesByColor(GC_BLACK);
    std::cout << "Structural variation resolver marked " << m_numMarked << " edges\n";
}

//
//
//
void ScaffoldConflictingVisitor::previsit(ScaffoldGraph*)
{
    m_numMarked = 0;
}

//
bool ScaffoldConflictingVisitor::visit(ScaffoldGraph*, ScaffoldVertex* pVertex)
{
    if(pVertex->hasConflictingLink())
    {
        pVertex->setClassification(SVC_REPEAT);
        m_numMarked += 1;
    }
    return false;
}

//
void ScaffoldConflictingVisitor::postvisit(ScaffoldGraph*)
{
    std::cout << "[conflict] marked: " << m_numMarked << "\n";
}

ScaffoldLayoutVisitor::ScaffoldLayoutVisitor()
{

}
  
//      
void ScaffoldLayoutVisitor::previsit(ScaffoldGraph* pGraph)
{
    pGraph->setEdgeColors(GC_WHITE);
}

//
bool ScaffoldLayoutVisitor::visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex)
{
    // Never re-classify repeats
    if(pVertex->getClassification() == SVC_REPEAT)
        return false;

    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() <= 1)
            continue;

        ScaffoldGroup group(pVertex, 100);
        for(size_t i = 0; i < edgeVec.size(); ++i)
            group.addLink(edgeVec[i]->getLink(), edgeVec[i]->getEnd());

        group.computeBestOrdering();
        int longestOverlap = group.calculateLongestOverlap();

        if(longestOverlap > 150)
        {
            std::cout << "cannot scaffold from " << pVertex->getID() << " OL: " << longestOverlap << "\n";
            continue;
        }

        // Calculate the ordered set of links that would make up this scaffold
        std::cout << "Computing layout for " << pVertex->getID() << "\n";
        LinkVector linearLinks;
        group.getLinearLinks(linearLinks);

        std::cout << "Linearized links:\n";

        // Iterate through the linearized links and see if they are sane
        bool allResolved = true;
        ScaffoldVertex* pCurrStartVertex = pVertex;
        for(size_t i = 0; i < linearLinks.size(); ++i)
        {
            std::cout << "\t" << linearLinks[i] << "\n";

            // Does a real link exist between the two vertices that are predicted to be linked?
            ScaffoldEdge* pRealEdge = pCurrStartVertex->findEdgeTo(linearLinks[i].endpointID, linearLinks[i].getDir(), linearLinks[i].getComp());
            if(pRealEdge != NULL)
            {
                std::cout << "\t\tresolved:" << pRealEdge->getLink() << "\n";
            }
            else
            {
                pCurrStartVertex->setClassification(SVC_REPEAT);
                allResolved = false;
                std::cout << "\t\tnotfound\n";
            }

            ScaffoldVertex* pCurrEndVertex = pGraph->getVertex(linearLinks[i].endpointID);
            assert(pCurrEndVertex != NULL);
            pCurrStartVertex = pCurrEndVertex;
        }

        group.getSecondaryLinks();

        if(allResolved)
        {
            // Removed all the edges for current vertex except for the first
            std::string firstID = linearLinks[0].endpointID;
            for(size_t i = 0; i < edgeVec.size(); ++i)
            {
                if(edgeVec[i]->getEndID() == firstID)
                {
                    edgeVec[i]->setColor(GC_BLACK);
                    edgeVec[i]->getTwin()->setColor(GC_BLACK);
                }
            }
        }
    }
    return false;
}

//
void ScaffoldLayoutVisitor::postvisit(ScaffoldGraph* pGraph)
{
    pGraph->deleteEdgesByColor(GC_BLACK);
}


//
// ScaffoldDistanceRefinementVisitor - refine distance estimates using a string graph
// of contig connectivity.
//
ScaffoldDistanceRefinementVisitor::ScaffoldDistanceRefinementVisitor(const StringGraph* pStringGraph) : m_pStringGraph(pStringGraph)
{

}

//
void ScaffoldDistanceRefinementVisitor::previsit(ScaffoldGraph* /*pGraph*/)
{
}

//
bool ScaffoldDistanceRefinementVisitor::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);

        for(size_t i = 0; i < edgeVec.size(); ++i)
        {
            Vertex* pX = m_pStringGraph->getVertex(edgeVec[i]->getStartID());
            Vertex* pY = m_pStringGraph->getVertex(edgeVec[i]->getEndID());
            assert(pX != NULL && pY != NULL);

            SGWalkVector walks;
            SGSearch::findWalks(pX, pY, edgeVec[i]->getDir(), edgeVec[i]->getDistance() + pY->getSeqLen() + 1000, true, 100000, walks); 

            // Select the walk closest to the distance estimate
            if(walks.size() > 0)
            {
                int closest = std::numeric_limits<int>::max();
                int est = edgeVec[i]->getDistance();

                size_t idx = -1;
                for(size_t j = 0; j < walks.size(); ++j)
                {
                    int diff = abs(walks[j].getEndToStartDistance() - est);
                    if(diff < closest)
                    {
                        closest = diff;
                        idx = j;
                    }
                }

                printf("%s -> %s\t%d\t%d\t%s\n", pX->getID().c_str(), 
                                                 pY->getID().c_str(), 
                                                 edgeVec[i]->getDistance(), 
                                                 walks[idx].getEndToStartDistance(),
                                                 walks[idx].pathSignature().c_str());
            }
            else
            {
                printf("%s -> %s\t%d\tN/A\n", pX->getID().c_str(), pY->getID().c_str(), edgeVec[i]->getDistance());
            }

        }
    }    
    return false;
}

//
void ScaffoldDistanceRefinementVisitor::postvisit(ScaffoldGraph* /*pGraph*/)
{
}

//
ScaffoldWriterVisitor::ScaffoldWriterVisitor(const std::string& filename)
{
    m_pWriter = createWriter(filename);
}

//
ScaffoldWriterVisitor::~ScaffoldWriterVisitor()
{
    delete m_pWriter;
}

//
void ScaffoldWriterVisitor::previsit(ScaffoldGraph* pGraph)
{
    pGraph->setVertexColors(GC_WHITE);
    m_statsVector.clear();
}

//
bool ScaffoldWriterVisitor::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    if(pVertex->getColor() == GC_RED)
        return false; //already output

    ScaffoldEdgePtrVector edges = pVertex->getEdges();
    
    if(edges.size() <= 1)
    {
        // Start of a chain found, traverse it to the other end
        size_t num_contigs = 0;
        size_t bases = 0; // number of bases in contigs
        size_t span = 0; // number of bases plus gaps

        pVertex->setColor(GC_RED);
    
        ScaffoldRecord record;
        record.setRoot(pVertex->getID());
        
        bases += pVertex->getSeqLen();
        span += pVertex->getSeqLen();
        num_contigs += 1;
        
        if(edges.size() == 1)
        {
            // Write the linked contigs
            ScaffoldEdge* pStartEdge = edges[0];
            ScaffoldEdge* pXY = pStartEdge;
            while(1)
            {
                record.addLink(pXY->getLink());
                ScaffoldVertex* pY = pXY->getEnd();
                if(pY->getColor() == GC_RED)
                    break; // loop found

                pY->setColor(GC_RED);
                bases += pY->getSeqLen();
                span += pY->getSeqLen() + pXY->getDistance();
                num_contigs += 1;

                // get the next direction
                EdgeDir nextDir = !pXY->getTwin()->getDir();
                ScaffoldEdgePtrVector nextEdges = pY->getEdges(nextDir);
                
                if(nextEdges.size() == 1)
                    pXY = nextEdges[0];
                else
                    break;
            }
        }
        record.writeScaf(m_pWriter);

        ScaffoldStats stats;
        stats.numContigs = num_contigs;
        stats.bases = bases;
        stats.span = span;
        m_statsVector.push_back(stats);
    }
    return false;
}

void ScaffoldWriterVisitor::postvisit(ScaffoldGraph*)
{
    size_t totalScaffolds = m_statsVector.size();
    size_t singleton = 0;
    size_t sumSpan = 0;
    size_t sumBases = 0;
    size_t totalContigs = 0;
    std::vector<ScaffoldStats>::iterator iter = m_statsVector.begin();
    for(; iter != m_statsVector.end(); ++iter)
    {
        sumSpan += iter->span;
        sumBases += iter->bases;       
        totalContigs += iter->numContigs;
        if(iter->numContigs == 1)
            ++singleton;
    }


    std::sort(m_statsVector.begin(), m_statsVector.end(), ScaffoldStats::sortSpanDesc);

    size_t targetSpan = sumSpan / 2;
    size_t spanN50 = -1;
    size_t maxSpan = 0;
    size_t runningSum = 0;
    iter = m_statsVector.begin();
    for(; iter != m_statsVector.end(); ++iter)
    {
        if((size_t)iter->span > maxSpan)
            maxSpan = iter->span;
        runningSum += iter->span;
        if(runningSum > targetSpan)
        {
            spanN50 = iter->span;
            break;
        }
    }

    printf("\n======\n");
    printf("Constructed %zu scaffolds from %zu contigs (%zu singleton scaffolds)\n", totalScaffolds, totalContigs, singleton);
    printf("Total bases: %.2lfMbp\n", (double)sumBases / 1000000);
    printf("Total span: %.2lfMbp\n", (double)sumSpan / 1000000);
    
    printf("Max scaffold span: %zu bp\n", maxSpan);
    printf("N50 scaffold span: %zu bp\n", spanN50);
    printf("Mean scaffold span: %.0lf bp\n", (double)sumSpan / totalScaffolds);
}
