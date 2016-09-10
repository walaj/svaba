//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGVisitors - Algorithms that visit
// each vertex in the graph and perform some
// operation
//
#include "SGVisitors.h"
#include "CompleteOverlapSet.h"
#include "SGSearch.h"
#include "stdaln.h"

//
// SGFastaVisitor - output the vertices in the graph in 
// fasta format
//
bool SGFastaVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    m_fileHandle << ">" << pVertex->getID() << " " <<  pVertex->getSeq().length() 
                 << " " << 0 << "\n";
    m_fileHandle << pVertex->getSeq() << "\n";
    return false;
}

// jwala added
bool SGVisitorContig::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
  m_ct.push_back(SeqLib::UnalignedSequence(pVertex->getID(), pVertex->getSeq().toString(),std::string()));
  return false;
}


//
// SGTransRedVisitor - Perform a transitive reduction about this vertex
// This uses Myers' algorithm (2005, The fragment assembly string graph)
// Precondition: the edge list is sorted by length (ascending)
void SGTransitiveReductionVisitor::previsit(StringGraph* pGraph)
{
    // The graph must not have containments
    assert(!pGraph->hasContainment());

    // Set all the vertices in the graph to "vacant"
    pGraph->setColors(GC_WHITE);
    pGraph->sortVertexAdjListsByLen();

    marked_verts = 0;
    marked_edges = 0;
}

//
bool SGTransitiveReductionVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    size_t trans_count = 0;
    static const size_t FUZZ = 10; // see myers

    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        EdgePtrVec edges = pVertex->getEdges(dir); // These edges are already sorted
        if(edges.size() == 0)
            continue;

        for(size_t i = 0; i < edges.size(); ++i)
            (edges[i])->getEnd()->setColor(GC_GRAY);

        Edge* pLongestEdge = edges.back();
        size_t longestLen = pLongestEdge->getSeqLen() + FUZZ;
        
        // Stage 1
        for(size_t i = 0; i < edges.size(); ++i)
        {
            Edge* pVWEdge = edges[i];
            Vertex* pWVert = pVWEdge->getEnd();

            EdgeDir transDir = !pVWEdge->getTwinDir();
            if(pWVert->getColor() == GC_GRAY)
            {
                EdgePtrVec w_edges = pWVert->getEdges(transDir);
                for(size_t j = 0; j < w_edges.size(); ++j)
                {
                    Edge* pWXEdge = w_edges[j];
                    size_t trans_len = pVWEdge->getSeqLen() + pWXEdge->getSeqLen();
                    if(trans_len <= longestLen)
                    {
                        if(pWXEdge->getEnd()->getColor() == GC_GRAY)
                        {
                            // X is the endpoint of an edge of V, therefore it is transitive
                            pWXEdge->getEnd()->setColor(GC_BLACK);
                        }
                    }
                    else
                        break;
                }
            }
        }
        
        // Stage 2
        for(size_t i = 0; i < edges.size(); ++i)
        {
            Edge* pVWEdge = edges[i];
            Vertex* pWVert = pVWEdge->getEnd();

            EdgeDir transDir = !pVWEdge->getTwinDir();
            EdgePtrVec w_edges = pWVert->getEdges(transDir);
            for(size_t j = 0; j < w_edges.size(); ++j)
            {
                Edge* pWXEdge = w_edges[j];
                size_t len = pWXEdge->getSeqLen();

                if(len < FUZZ || j == 0)
                {
                    if(pWXEdge->getEnd()->getColor() == GC_GRAY)
                    {
                        // X is the endpoint of an edge of V, therefore it is transitive
                        pWXEdge->getEnd()->setColor(GC_BLACK);
                    }
                }
                else
                {
                    break;
                }
            }
        }

        for(size_t i = 0; i < edges.size(); ++i)
        {
            if(edges[i]->getEnd()->getColor() == GC_BLACK)
            {
                // Mark the edge and its twin for removal
                if(edges[i]->getColor() != GC_BLACK || edges[i]->getTwin()->getColor() != GC_BLACK)
                {
                    edges[i]->setColor(GC_BLACK);
                    edges[i]->getTwin()->setColor(GC_BLACK);
                    marked_edges += 2;
                    trans_count++;
                }
            }
            edges[i]->getEnd()->setColor(GC_WHITE);
        }
    }

    if(trans_count > 0)
        ++marked_verts;

    return false;
}

// Remove all the marked edges
void SGTransitiveReductionVisitor::postvisit(StringGraph* pGraph)
{
    //printf("TR marked %d verts and %d edges\n", marked_verts, marked_edges);
    pGraph->sweepEdges(GC_BLACK);
    pGraph->setTransitiveFlag(false);
    assert(pGraph->checkColors(GC_WHITE));
}

//
// SGIdenticalRemoveVisitor - Removes identical vertices
// from the graph. This algorithm is less complex 
// than SGContainRemoveVisitor because we do not have to
// remodel the graph in this case because no irreducible
// edges need to be moved.
//
void SGIdenticalRemoveVisitor::previsit(StringGraph* pGraph)
{
    pGraph->setColors(GC_WHITE);
    count = 0;
}

// 
bool SGIdenticalRemoveVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    if(!pVertex->isContained())
        return false;

    // Check if this vertex is identical to any other vertex
    EdgePtrVec neighborEdges = pVertex->getEdges();
    for(size_t i = 0; i < neighborEdges.size(); ++i)
    {
        Edge* pEdge = neighborEdges[i];
        Vertex* pOther = pEdge->getEnd();
        if(pVertex->getSeqLen() != pOther->getSeqLen())
            continue;
        
        Overlap ovr = pEdge->getOverlap();
        if(!ovr.isContainment() || ovr.getContainedIdx() != 0)
            continue;

        if(pVertex->getSeq() == pOther->getSeq())
        {
            pVertex->setColor(GC_BLACK);
            ++count;
            break;
        }
    }
            
    return false;
}

void SGIdenticalRemoveVisitor::postvisit(StringGraph* pGraph)
{
    pGraph->sweepVertices(GC_BLACK);
}

//
// SGContainRemoveVisitor - Removes contained
// vertices from the graph
//
void SGContainRemoveVisitor::previsit(StringGraph* pGraph)
{
    pGraph->setColors(GC_WHITE);

    // Clear the containment flag, if any containments are added
    // during this algorithm the flag will be reset and another
    // round must be re-run
    pGraph->setContainmentFlag(false);    
}

//
bool SGContainRemoveVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{

    if(!pVertex->isContained())
        return false;

    //cout << pVertex->getID() << endl; // debug
    // Add any new irreducible edges that exist when pToRemove is deleted
    // from the graph
    EdgePtrVec neighborEdges = pVertex->getEdges();
    
    // If the graph has been transitively reduced, we have to check all
    // the neighbors to see if any new edges need to be added. If the graph is a
    // complete overlap graph we can just remove the edges to the deletion vertex
    if(!pGraph->hasTransitive() && !pGraph->isExactMode())
    {
        // This must be done in order of edge length or some transitive edges
        // may be created
        EdgeLenComp comp;
        std::sort(neighborEdges.begin(), neighborEdges.end(), comp);

        for(size_t j = 0; j < neighborEdges.size(); ++j)
        {
            Vertex* pRemodelVert = neighborEdges[j]->getEnd();
            Edge* pRemodelEdge = neighborEdges[j]->getTwin();
            SGAlgorithms::remodelVertexForExcision(pGraph, 
                                                   pRemodelVert, 
                                                   pRemodelEdge);
        }
    }
            
    // Delete the edges from the graph
    for(size_t j = 0; j < neighborEdges.size(); ++j)
    {
        Vertex* pRemodelVert = neighborEdges[j]->getEnd();
        Edge* pRemodelEdge = neighborEdges[j]->getTwin();
        pRemodelVert->deleteEdge(pRemodelEdge);
        pVertex->deleteEdge(neighborEdges[j]);
    }
    pVertex->setColor(GC_BLACK);
    return false;
}

void SGContainRemoveVisitor::postvisit(StringGraph* pGraph)
{
    pGraph->sweepVertices(GC_BLACK);
}

//
// Validate the structure of the graph by detecting missing
// or erroneous edges
//

//
bool SGValidateStructureVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
    SGAlgorithms::EdgeDescOverlapMap irreducibleMap;
    SGAlgorithms::EdgeDescOverlapMap transitiveMap;
    
    // Construct the set of overlaps reachable within the current parameters
    CompleteOverlapSet vertexOverlapSet(pVertex, pGraph->getErrorRate(), pGraph->getMinOverlap());
    vertexOverlapSet.computeIrreducible(NULL, NULL);

    SGAlgorithms::EdgeDescOverlapMap missingMap;
    SGAlgorithms::EdgeDescOverlapMap extraMap;
    
    vertexOverlapSet.getDiffMap(missingMap, extraMap);  
    
    if(!missingMap.empty())
    {
        std::cout << "Missing irreducible for " << pVertex->getID() << ":\n";
        SGAlgorithms::printOverlapMap(missingMap);
    }

    if(!extraMap.empty())
    {
        std::cout << "Extra irreducible for " << pVertex->getID() << ":\n";
        SGAlgorithms::printOverlapMap(extraMap);
    }
    return false;
}

//
// SGTrimVisitor - Remove "dead-end" vertices from the graph
//
void SGTrimVisitor::previsit(StringGraph* pGraph)
{
    num_island = 0;
    num_terminal = 0;
    pGraph->setColors(GC_WHITE);
}

// Mark any nodes that either dont have edges or edges in only one direction for removal
bool SGTrimVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    if(pVertex->countEdges() == 0)
    {
        // Is an island, remove if the sequence length is less than the threshold
        if(pVertex->getSeqLen() < m_minLength)
        {
            pVertex->setColor(GC_BLACK);
            ++num_island;
        }
    }
    else
    {
        // Check if this node is a dead-end
        for(size_t idx = 0; idx < ED_COUNT; idx++)
        {
            EdgeDir dir = EDGE_DIRECTIONS[idx];
            if(pVertex->countEdges(dir) == 0 && pVertex->getSeqLen() < m_minLength)
            {
                pVertex->setColor(GC_BLACK);
                ++num_terminal;
            }
        }
    }

    return false;
}

// Remove all the marked edges
void SGTrimVisitor::postvisit(StringGraph* pGraph)
{
    pGraph->sweepVertices(GC_BLACK);
    //printf("StringGraphTrim: Removed %d island and %d dead-end short vertices\n", num_island, num_terminal);
}

//
// SGDuplicateVisitor - Detect and remove duplicate edges
//
void SGDuplicateVisitor::previsit(StringGraph* pGraph)
{
    assert(pGraph->checkColors(GC_WHITE));
    (void)pGraph;
    m_hasDuplicate = false;
}

bool SGDuplicateVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    m_hasDuplicate = pVertex->markDuplicateEdges(GC_RED) || m_hasDuplicate;
    return false;
}


void SGDuplicateVisitor::postvisit(StringGraph* pGraph)
{
    assert(pGraph->checkColors(GC_WHITE));
    //if(m_hasDuplicate)
    // {
    //int numRemoved = pGraph->sweepEdges(GC_RED);
        //if(!m_bSilent)
        //    std::cerr << "Warning: removed " << numRemoved << " duplicate edges\n";
	//}
}

//
// Small repeat resolver - Remove edges induced from small (sub-read length)
// repeats
// 
void SGSmallRepeatResolveVisitor::previsit(StringGraph*)
{

}

//
bool SGSmallRepeatResolveVisitor::visit(StringGraph* /*pGraph*/, Vertex* pX)
{
    bool changed = false;

    // If the vertex has more than MAX_EDGES, do not
    // attempt to resolve
    size_t MAX_EDGES = 10;

    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        EdgePtrVec x_edges = pX->getEdges(dir); // These edges are already sorted

        if(x_edges.size() < 2 || x_edges.size() > MAX_EDGES)
            continue;

        // Try to eliminate the shortest edge from this vertex (let this be X->Y)
        // If Y has a longer edge than Y->X in the same direction, we remove X->Y

        // Edges are sorted by length so the last edge is the shortest
        Edge* pXY = x_edges.back();
        size_t xy_len = pXY->getOverlap().getOverlapLength(0);
        size_t x_longest_len = x_edges.front()->getOverlap().getOverlapLength(0);
        if(xy_len == x_longest_len)
            continue;

        Edge* pYX = pXY->getTwin();
        Vertex* pY = pXY->getEnd();

        EdgePtrVec y_edges = pY->getEdges(pYX->getDir());
        if(y_edges.size() > MAX_EDGES)
            continue;

        size_t yx_len = pYX->getOverlap().getOverlapLength(0);

        size_t y_longest_len = 0;
        for(size_t i = 0; i < y_edges.size(); ++i)
        {
            Edge* pYZ = y_edges[i];
            if(pYZ == pYX)
                continue; // skip Y->X

            size_t yz_len = pYZ->getOverlap().getOverlapLength(0);
            if(yz_len > y_longest_len)
                y_longest_len = yz_len;
        }


        if(y_longest_len > yx_len)
        {
            // Delete the edge if the difference between the shortest and longest is greater than minDiff
            int x_diff = x_longest_len - xy_len;
            int y_diff = y_longest_len - yx_len;

            if(x_diff > m_minDiff && y_diff > m_minDiff)
            {
                pX->deleteEdge(pXY);
                pY->deleteEdge(pYX);
                changed = true;
            }
        }
    }

    return changed;
}

//
void SGSmallRepeatResolveVisitor::postvisit(StringGraph* pGraph)
{
    pGraph->sweepEdges(GC_RED);
}

//
// OverlapRatio Visitor - Only keep edges in the graph when the ratio
// between their length and the length of the longest overlap meets a minimum cutoff
// 
void SGOverlapRatioVisitor::previsit(StringGraph*)
{

}

//
bool SGOverlapRatioVisitor::visit(StringGraph* /*pGraph*/, Vertex* pX)
{
    bool changed = false;

    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        EdgePtrVec x_edges = pX->getEdges(dir); // These edges are already sorted

        if(x_edges.size() < 2)
            continue;

        size_t x_longest_len = x_edges.front()->getOverlap().getOverlapLength(0);
        for(size_t i = 1; i < x_edges.size(); ++i)
        {
            size_t curr_len = x_edges[i]->getOverlap().getOverlapLength(0);
            double ratio = (double)curr_len / x_longest_len;
            if(ratio < m_minRatio)
            {
                x_edges[i]->setColor(GC_RED);
                x_edges[i]->getTwin()->setColor(GC_RED);
                changed = true;
            }
        }
    }

    return changed;
}

//
void SGOverlapRatioVisitor::postvisit(StringGraph* pGraph)
{
    pGraph->sweepEdges(GC_RED);
}

//
// SGSuperRepeatVisitor
//

// Remove all edges of nodes that have been marked as super repeats
void SGSuperRepeatVisitor::previsit(StringGraph*)
{
    m_num_superrepeats = 0;
}

//
bool SGSuperRepeatVisitor::visit(StringGraph*, Vertex* pVertex)
{
    if(pVertex->isSuperRepeat())
    {
        pVertex->deleteEdges();
        m_num_superrepeats += 1;
        return true;
    }
    return false;
}

//
void SGSuperRepeatVisitor::postvisit(StringGraph*)
{
  //printf("Deleted edges for %zu super repetitive vertices\n", m_num_superrepeats); 
}

//
// SGSmoothingVisitor - Find branches in the graph
// which arise from variation and remove them
//
void SGSmoothingVisitor::previsit(StringGraph* pGraph)
{
    pGraph->setColors(GC_WHITE);
    m_simpleBubblesRemoved = 0;
    m_complexBubblesRemoved = 0;
}

//
bool SGSmoothingVisitor::visit(StringGraph* pGraph, Vertex* pVertex)
{
    (void)pGraph;
    if(pVertex->getColor() == GC_RED)
        return false;

    bool found = false;
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        EdgePtrVec edges = pVertex->getEdges(dir);
        if(edges.size() <= 1)
            continue;

        for(size_t i = 0; i < edges.size(); ++i)
        {
            if(edges[i]->getEnd()->getColor() == GC_RED)
                return false;
        }

        //std::cout << "Smoothing " << pVertex->getID() << "\n";

        const int MAX_WALKS = 10;
        const int MAX_DISTANCE = 5000;
        bool bIsDegenerate = false;
        bool bFailGapCheck = false;
        bool bFailDivergenceCheck = false;
        bool bFailIndelSizeCheck = false;

        SGWalkVector variantWalks;
        SGSearch::findVariantWalks(pVertex, dir, MAX_DISTANCE, MAX_WALKS, variantWalks);

        if(variantWalks.size() > 0)
        {
            found = true;
            size_t selectedIdx = -1;
            size_t selectedCoverage = 0;

            // Calculate the minimum amount overlapped on the start/end vertex.
            // This is used to properly extract the sequences from walks that represent the variation.
            int minOverlapX = std::numeric_limits<int>::max();
            int minOverlapY = std::numeric_limits<int>::max();

            for(size_t i = 0; i < variantWalks.size(); ++i)
            {
                if(variantWalks[i].getNumEdges() <= 1)
                    bIsDegenerate = true;

                // Calculate the walk coverage using the internal vertices of the walk. 
                // The walk with the highest coverage will be retained
                size_t walkCoverage = 0;
                for(size_t j = 1; j < variantWalks[i].getNumVertices() - 1; ++j)
                    walkCoverage += variantWalks[i].getVertex(j)->getCoverage();

                if(walkCoverage > selectedCoverage || selectedCoverage == 0)
                {
                    selectedIdx = i;
                    selectedCoverage = walkCoverage;
                }
                
                Edge* pFirstEdge = variantWalks[i].getFirstEdge();
                Edge* pLastEdge = variantWalks[i].getLastEdge();

                if((int)pFirstEdge->getMatchLength() < minOverlapX)
                    minOverlapX = pFirstEdge->getMatchLength();

                if((int)pLastEdge->getTwin()->getMatchLength() < minOverlapY)
                    minOverlapY = pLastEdge->getTwin()->getMatchLength();
            }

            // Calculate the strings for each walk that represent the region of variation
            StringVector walkStrings;
            for(size_t i = 0; i < variantWalks.size(); ++i)
            {
                Vertex* pStartVertex = variantWalks[i].getStartVertex();
                Vertex* pLastVertex = variantWalks[i].getLastVertex();
                assert(pStartVertex != NULL && pLastVertex != NULL);
                
                std::string full = variantWalks[i].getString(SGWT_START_TO_END);
                int posStart = 0;
                int posEnd = 0;

                if(dir == ED_ANTISENSE)
                {
                    // pLast   -----------
                    // pStart          ------------
                    // full    --------------------
                    // out             ----
                    posStart = pLastVertex->getSeqLen() - minOverlapY;
                    posEnd = full.size() - (pStartVertex->getSeqLen() - minOverlapX);
                }
                else
                {
                    // pStart         --------------
                    // pLast   -----------
                    // full    ---------------------
                    // out            ----
                    posStart = pStartVertex->getSeqLen() - minOverlapX; // match start position
                    posEnd = full.size() - (pLastVertex->getSeqLen() - minOverlapY); // match end position
                }
                
                std::string out;
                if(posEnd > posStart)
                    out = full.substr(posStart, posEnd - posStart);
                walkStrings.push_back(out);
            }

            assert(selectedIdx != (size_t)-1);
            SGWalk& selectedWalk = variantWalks[selectedIdx];
            assert(selectedWalk.isIndexed());

            // Check the divergence of the other walks to this walk
            StringVector cigarStrings;
            std::vector<int> maxIndel;
            std::vector<double> gapPercent; // percentage of matching that is gaps
            std::vector<double> totalPercent; // percent of total alignment that is mismatch or gap

            cigarStrings.resize(variantWalks.size());
            gapPercent.resize(variantWalks.size());
            totalPercent.resize(variantWalks.size());
            maxIndel.resize(variantWalks.size());

            for(size_t i = 0; i < variantWalks.size(); ++i)
            {
                if(i == selectedIdx)
                    continue;

                // We want to compute the total gap length, total mismatches and percent
                // divergence between the two paths.
                int matchLen = 0;
                int totalDiff = 0;
                int gapLength = 0;
                int maxGapLength = 0;
                // We have to handle the degenerate case where one internal string has zero length
                // this can happen when there is an isolated insertion/deletion and the walks are like:
                // x -> y -> z
                // x -> z
                if(walkStrings[selectedIdx].empty() || walkStrings[i].empty())
                {
                    matchLen = std::max(walkStrings[selectedIdx].size(), walkStrings[i].size());
                    totalDiff = matchLen;
                    gapLength = matchLen;
                }
                else
                {
                    AlnAln *aln_global;
                    aln_global = aln_stdaln(walkStrings[selectedIdx].c_str(), walkStrings[i].c_str(), &aln_param_blast, 1, 1);

                    // Calculate the alignment parameters
                    while(aln_global->outm[matchLen] != '\0')
                    {
                        if(aln_global->outm[matchLen] == ' ')
                            totalDiff += 1;
                        matchLen += 1;
                    }

                    std::stringstream cigarSS;
                    for (int j = 0; j != aln_global->n_cigar; ++j)
                    {
                        char cigarOp = "MID"[aln_global->cigar32[j]&0xf];
                        int cigarLen = aln_global->cigar32[j]>>4;
                        if(cigarOp == 'I' || cigarOp == 'D')
                        {
                            gapLength += cigarLen;
                            if(gapLength > maxGapLength)
                                maxGapLength = gapLength;
                        }

                        cigarSS << cigarLen;
                        cigarSS << cigarOp;
                    }
                    cigarStrings[i] = cigarSS.str();
                    aln_free_AlnAln(aln_global);
                }

                double percentDiff = (double)totalDiff / matchLen;
                double percentGap = (double)gapLength / matchLen;

                if(percentDiff > m_maxTotalDivergence)
                    bFailDivergenceCheck = true;
                
                if(percentGap > m_maxGapDivergence)
                    bFailGapCheck = true;

                if(maxGapLength > m_maxIndelLength)
                    bFailIndelSizeCheck = true;

                gapPercent[i] = percentGap;
                totalPercent[i] = percentDiff;
                maxIndel[i] = maxGapLength;
            }

            if(bIsDegenerate || bFailGapCheck || bFailDivergenceCheck || bFailIndelSizeCheck)
                continue;

            // Write the selected path to the variants file as variant 0
            int variantIdx = 0;
            std::string selectedSequence = selectedWalk.getString(SGWT_START_TO_END);
            std::stringstream ss;
            ss << "variant-" << m_numRemovedTotal << "/" << variantIdx++;
            writeFastaRecord(&m_outFile, ss.str(), selectedSequence);


            // The vertex set for each walk is not necessarily disjoint,
            // the selected walk may contain vertices that are part
            // of other paths. We handle this be initially marking all
            // vertices of the 
            for(size_t i = 0; i < variantWalks.size(); ++i)
            {
                if(i == selectedIdx)
                    continue;

                SGWalk& currWalk = variantWalks[i];
                for(size_t j = 0; j < currWalk.getNumEdges() - 1; ++j)
                {
                    Edge* currEdge = currWalk.getEdge(j);
                    
                    // If the vertex is also on the selected path, do not mark it
                    Vertex* currVertex = currEdge->getEnd();
                    if(!selectedWalk.containsVertex(currVertex->getID()))
                    {
                        currEdge->getEnd()->setColor(GC_RED);
                    }
                }

                // Write the variant to a file
                std::string variantSequence = currWalk.getString(SGWT_START_TO_END);
                std::stringstream variantID;
                std::stringstream ss;
                ss << "variant-" << m_numRemovedTotal << "/" << variantIdx++;
                ss << " IGD:" << (double)gapPercent[i] << " ITD:" << totalPercent[i] << " MID: " << maxIndel[i] << " InternalCigar:" << cigarStrings[i];
                writeFastaRecord(&m_outFile, ss.str(), variantSequence);
            }

            if(variantWalks.size() == 2)
                m_simpleBubblesRemoved += 1;
            else
                m_complexBubblesRemoved += 1;
            ++m_numRemovedTotal;
        }
    }
    return found;
}

// Remove all the marked edges
void SGSmoothingVisitor::postvisit(StringGraph* pGraph)
{
    pGraph->sweepVertices(GC_RED);
    assert(pGraph->checkColors(GC_WHITE));

    //printf("VariationSmoother: Removed %d simple and %d complex bubbles\n", m_simpleBubblesRemoved, m_complexBubblesRemoved);
}

//
// SGGraphStatsVisitor - Collect summary stasitics
// about the graph
//
void SGGraphStatsVisitor::previsit(StringGraph* /*pGraph*/)
{
    num_terminal = 0;
    num_island = 0;
    num_monobranch = 0;
    num_dibranch = 0;
    num_simple = 0;
    num_edges = 0;
    num_vertex = 0;
    sum_edgeLen = 0;
}

// Find bubbles (nodes where there is a split and then immediate rejoin) and mark them for removal
bool SGGraphStatsVisitor::visit(StringGraph* /*pGraph*/, Vertex* pVertex)
{
    int s_count = pVertex->countEdges(ED_SENSE);
    int as_count = pVertex->countEdges(ED_ANTISENSE);
    if(s_count == 0 && as_count == 0)
    {
        ++num_island;
    }
    else if(s_count == 0 || as_count == 0)
    {
        ++num_terminal;
    }

    if(s_count > 1 && as_count > 1)
        ++num_dibranch;
    else if(s_count > 1 || as_count > 1)
        ++num_monobranch;

    if(s_count == 1 || as_count == 1)
        ++num_simple;

    num_edges += (s_count + as_count);
    ++num_vertex;

    EdgePtrVec edges = pVertex->getEdges();
    for(size_t i = 0; i < edges.size(); ++i)
        sum_edgeLen += edges[i]->getSeqLen();

    return false;
}

//
void SGGraphStatsVisitor::postvisit(StringGraph* /*pGraph*/)
{
  printf("Vertices: %d Edges: %d Islands: %d Tips: %d Monobranch: %d Dibranch: %d Simple: %d\n", num_vertex, num_edges, 
  num_island, num_terminal,
                                                                                                   num_monobranch, num_dibranch, num_simple);
  
}
