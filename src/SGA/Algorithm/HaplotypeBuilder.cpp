///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// HaplotypeBuilder - Construct candidate
// haplotypes from a pair of k-mer seeds.
//
#include "HaplotypeBuilder.h"
#include "BWTAlgorithms.h"
#include "SGSearch.h"
#include "SGAlgorithms.h"
#include "Profiler.h"

//
//
//
HaplotypeBuilder::HaplotypeBuilder() : m_pStartVertex(NULL), m_pJoinVertex(NULL), m_kmerThreshold(1), m_kmerSize(51)
{
    m_pGraph = new StringGraph;
}

//
HaplotypeBuilder::~HaplotypeBuilder()
{
    delete m_pGraph;
}

//
void HaplotypeBuilder::setKmerParameters(size_t k, size_t t)
{
    m_kmerSize = k;
    m_kmerThreshold = t;
}

// The source string is the string the bubble starts from
void HaplotypeBuilder::setTerminals(const AnchorSequence& leftAnchor, const AnchorSequence& rightAnchor)
{
    // Create a new vertex for the source sequence
    // As we are creating a de Bruijn graph, we use the sequence
    // of the vertex as its ID
    assert(leftAnchor.sequence != rightAnchor.sequence);
    Vertex* pLeftVertex = new(m_pGraph->getVertexAllocator()) Vertex(leftAnchor.sequence, leftAnchor.sequence);
    addVertex(pLeftVertex, leftAnchor.count);

    Vertex* pRightVertex = new(m_pGraph->getVertexAllocator()) Vertex(rightAnchor.sequence, rightAnchor.sequence);
    addVertex(pRightVertex, rightAnchor.count);

    // Add the vertex to the extension queue
    m_queue.push(BuilderExtensionNode(pLeftVertex, ED_SENSE));
    
    m_pStartVertex = pLeftVertex;
    m_pJoinVertex = pRightVertex;
}

// The source index is the index that the contains the source string
void HaplotypeBuilder::setIndex(const BWT* pBWT, const BWT* pRBWT)
{
    m_pBWT = pBWT;
    m_pRevBWT = pRBWT;
}

// Run the bubble construction process
HaplotypeBuilderReturnCode HaplotypeBuilder::run()
{
    PROFILE_FUNC("HaplotypeBuilder::run")
    assert(m_queue.size() == 1);
    assert(m_pJoinVertex != NULL);
    assert(m_pBWT != NULL);

    size_t MAX_ITERATIONS = 2000;
    size_t MAX_SIMULTANEOUS_BRANCHES = 20;
    size_t MAX_TOTAL_BRANCHES = 50;
   
    // Tracking stats
    size_t total_branches = 0;
    size_t iterations = 0;

    while(!m_queue.empty())
    {
        if(iterations > MAX_ITERATIONS || m_queue.size() > MAX_SIMULTANEOUS_BRANCHES || total_branches > MAX_TOTAL_BRANCHES)
            return HBRC_TOO_MANY_VERTICES;

        iterations += 1;
        BuilderExtensionNode curr = m_queue.front();
        m_queue.pop();

        // Calculate de Bruijn extensions for this node
        std::string vertStr = curr.pVertex->getSeq().toString();
        AlphaCount64 extensionCounts;
        if(m_pRevBWT != NULL)
            extensionCounts = BWTAlgorithms::calculateDeBruijnExtensions(vertStr, m_pBWT, m_pRevBWT, curr.direction);
        else
            extensionCounts = BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(vertStr, m_pBWT, curr.direction);
        
        size_t num_added = 0;
        for(size_t i = 0; i < DNA_ALPHABET::size; ++i)
        {
            char b = DNA_ALPHABET::getBase(i);
            size_t count = extensionCounts.get(b);
            if(count < m_kmerThreshold)
                continue;

            std::string newStr = makeDeBruijnVertex(vertStr, b, curr.direction);
            Vertex* pVertex = m_pGraph->getVertex(newStr);
            
            // Check if we have found the vertex we are assembling to
            bool joinFound = pVertex != NULL && pVertex == m_pJoinVertex;
            if(!joinFound && pVertex == NULL)
            {
                pVertex = new(m_pGraph->getVertexAllocator()) Vertex(newStr, newStr);
                addVertex(pVertex, count);
                m_queue.push(BuilderExtensionNode(pVertex, curr.direction));
                num_added += 1;
            }
            
            // Create the new edge in the graph
            addDeBruijnEdges(curr.pVertex, pVertex, curr.direction);

            // If we've found the join vertex, we have completed the target half of the bubble
            if(joinFound)
                return HBRC_OK;
        }

        if(num_added > 0)
            total_branches += (num_added - 1);
    }

    // no path between the nodes found
    return HBRC_NO_PATH;
}

// After the bubble has been built into the graph, this function
// finds and compares the two sequences
HaplotypeBuilderReturnCode HaplotypeBuilder::parseWalks(HaplotypeBuilderResult& results) const
{
    // Parse walks from the graph that go through the bubbles
    SGWalkVector outWalks;
    bool success = SGSearch::findWalks(m_pStartVertex,
                                       m_pJoinVertex,
                                       ED_SENSE,
                                       10000, // max distance to search
                                       10000, // max nodes to search
                                       true, // exhaustive search
                                       outWalks);
    if(!success)
        return HBRC_WALK_FAILED;

    // Convert the walks into strings
    for(size_t i = 0; i < outWalks.size(); ++i)
    {
        std::string walkStr = outWalks[i].getString(SGWT_START_TO_END);
        results.haplotypes.push_back(walkStr);
    }
    
    return HBRC_OK;
}

// Add a vertex to the graph
void HaplotypeBuilder::addVertex(Vertex* pVertex, int coverage)
{
    m_pGraph->addVertex(pVertex);
    m_vertexCoverageMap[pVertex->getSeq().toString()] = coverage;
}

// Add a new edge to the graph denoting the relationship between pX and pY.
// Assumes pX and pY are already present in the m_pGraph
void HaplotypeBuilder::addDeBruijnEdges(const Vertex* pX, const Vertex* pY, EdgeDir direction)
{
    assert(pX->getSeq().length() == pY->getSeq().length());
    
    // overlap length for a de bruijn edge
    size_t p = pX->getSeq().length() - 1;

    // Construct an overlap object for this relationship
    Overlap o;
    o.id[0] = pX->getID();
    o.id[1] = pY->getID();

    o.match.isReverse = false;
    o.match.numDiff = 0;

    if(direction == ED_SENSE)
    {
        // pX -> pY
        o.match.coord[0].interval.start = 1;
        o.match.coord[1].interval.start = 0;
    }
    else
    {
        // pY -> pX
        o.match.coord[0].interval.start = 0;
        o.match.coord[1].interval.start = 1;
    }

    o.match.coord[0].interval.end = o.match.coord[0].interval.start + p - 1; // inclusive coordinate
    o.match.coord[1].interval.end = o.match.coord[1].interval.start + p - 1;
    o.match.coord[0].seqlen = p + 1;
    o.match.coord[1].seqlen = p + 1;
    Edge* e = SGAlgorithms::createEdgesFromOverlap(m_pGraph, o, false);
    assert(e != NULL);
}


// Make the sequence of a new deBruijn vertex using the edge details
std::string HaplotypeBuilder::makeDeBruijnVertex(const std::string& v, char edgeBase, EdgeDir direction)
{
    std::string w;
    size_t p = v.size() - 1;
    if(direction == ED_SENSE)
    {
        w = v.substr(1, p);
        w.append(1, edgeBase);
    }
    else
    {
        w.append(1, edgeBase);
        w.append(v.substr(0, p));
    }
    return w;
}

// Count the number of extensions of a de Bruijn node that are above
// the required k-mer coverage
size_t HaplotypeBuilder::countValidExtensions(const AlphaCount64& ac) const
{
    size_t n = 0;
    for(size_t i = 0; i < DNA_ALPHABET::size; ++i)
    {
        char b = DNA_ALPHABET::getBase(i);
        size_t count = ac.get(b);
        if(count >= m_kmerThreshold)
            n += 1;
    }
    return n;
}
