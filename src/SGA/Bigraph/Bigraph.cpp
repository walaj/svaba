//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Bidirectional graph 
//
#include <assert.h>
#include <ostream>
#include <fstream>
#include <iostream>
#include "Bigraph.h"
#include "Timer.h"
#include "ASQG.h"

//
//
//
Bigraph::Bigraph() : m_hasContainment(false), m_hasTransitive(false), m_isExactMode(false), m_minOverlap(0), m_errorRate(0.0f)
{
  //JEREMIAH
    // Set up the memory pools for the graph
    //m_pEdgeAllocator = new SimpleAllocator<Edge>();
    //m_pVertexAllocator = new SimpleAllocator<Vertex>();

    //m_vertices.set_deleted_key(""); // JEREMIAH
    //WARN_ONCE("HARDCODED HASH TABLE MAX SIZE");
    //m_vertices.resize(600000000);
}

//
//
//
Bigraph::~Bigraph()
{
    VertexPtrMap::iterator iter = m_vertices.begin(); 
    for(; iter != m_vertices.end(); ++iter)
    {
        delete iter->second;
        iter->second = NULL;
    }

    //JEREMIAH
    // Clean up the memory pools
    //delete m_pEdgeAllocator;
    //delete m_pVertexAllocator;
}

//
// Add a vertex
//
void Bigraph::addVertex(Vertex* pVert)
{
  //    std::pair<VertexPtrMapIter, bool> result = 
  //    m_vertices_tokeep.insert(std::make_pair(pVert->getID(), pVert));//JEREMIAH
  std::pair<VertexPtrMapIter, bool> result = m_vertices.insert(std::make_pair(pVert->getID(), pVert));

    if(!result.second)
    {
        std::cerr << "Error: Attempted to insert vertex into graph with a duplicate id: " <<
                     pVert->getID() << "\n";
        std::cerr << "All reads must have a unique identifier\n";
        exit(1);
    }
}

//
// Remove a vertex that is not connected to any other
//
void Bigraph::removeIslandVertex(Vertex* pVertex)
{
    assert(pVertex->countEdges() == 0);

    // Remove the vertex from the collection
    VertexID id = pVertex->getID();
    delete pVertex;
    m_vertices.erase(id);
    //m_vertices_tokeep.erase(id); //JEREMIAH
}

//
// Remove a vertex that is not connected to any other
//
void Bigraph::removeIslandVertexNoErase(Vertex* pVertex)
{
    assert(pVertex->countEdges() == 0);

    // Remove the vertex from the collection
    VertexID id = pVertex->getID();
    delete pVertex;
    m_vertices.erase(id);
}


//
// Remove a vertex
//
void Bigraph::removeConnectedVertex(Vertex* pVertex)
{
    // Remove the edges pointing to this Vertex
    pVertex->deleteEdges();

    // Remove the vertex from the collection
    VertexID id = pVertex->getID();
    delete pVertex;
    m_vertices.erase(id);
}


//
// Check for the existance of a vertex
//
bool Bigraph::hasVertex(VertexID id)
{
    VertexPtrMapIter iter = m_vertices.find(id);
    return iter != m_vertices.end();
}

//
// Get a vertex
//
Vertex* Bigraph::getVertex(VertexID id) const
{
    VertexPtrMapConstIter iter = m_vertices.find(id);
    if(iter == m_vertices.end())
        return NULL;
    return iter->second;
}

//
// Add an edge
//
void Bigraph::addEdge(Vertex* pVertex, Edge* pEdge)
{
    assert(pEdge->getStart() == pVertex);
    pVertex->addEdge(pEdge);
}

//
// Remove an edge
//
void Bigraph::removeEdge(const EdgeDesc& ed)
{
    ed.pVertex->removeEdge(ed);
}

//
// High level merge function that does not specify an edge
//
void Bigraph::mergeVertices(VertexID id1, VertexID id2)
{
    Vertex* pVert1 = getVertex(id1);

    // Get the edges from vertex1 to vertex2
    EdgePtrVec edgesTo = pVert1->findEdgesTo(id2);

    if(edgesTo.empty())
    {
        std::cerr << "mergeVertices: vertices are not connected\n";
        return;
    }

    if(edgesTo.size() > 1)
    {
        std::cerr << "mergeVertces: cannot merge because of ambigious edges\n";
        return;
    }

    // There is a single unique edge between the vertices
    Edge* mergeEdge = *edgesTo.begin();

    // Call the real merging function
    merge(pVert1, mergeEdge);
}

//
// Merge two vertices along the specified edge
//
void Bigraph::merge(Vertex* pV1, Edge* pEdge)
{
    Vertex* pV2 = pEdge->getEnd();
    //std::cout << "Merging " << pV1->getID() << " with " << pV2->getID() << "\n";

    // jeremiah connected component tracking
    if (m_get_components) {
      int lab1 = m_connected_components[pV1->getID()];
      int lab2 = m_connected_components[pV2->getID()];
      ++m_label;
      m_connected_components[pV1->getID()] = m_label;
      m_connected_components[pV2->getID()] = m_label;
      for (auto& i : m_connected_components)
	if (i.second == lab1 || i.second == lab2)
	  i.second = m_label;
    }

    //JEREMIAH
    //pV1->addParent(pV2->getID());
    //std::cerr << "LEN: " << pV1->getParentID().size() << std::endl;

    // Merge the data
    pV1->merge(pEdge);

    // Get the twin edge (the edge in v2 that points to v1)
    Edge* pTwin = pEdge->getTwin();

    // Ensure v2 has the twin edge
    assert(pV2->hasEdge(pTwin));

    // Get the edge set opposite of the twin edge (which will be the new edges in this direction for V1)
    EdgePtrVec transEdges = pV2->getEdges(!pTwin->getDir());

    // Move the edges from pV2 to pV1
    for(EdgePtrVecIter iter = transEdges.begin(); iter != transEdges.end(); ++iter)
    {
        Edge* pTransEdge = *iter;

        // Remove the edge from V2, this does not destroy the edge
        pV2->removeEdge(pTransEdge);

        // Join pEdge to the start of transEdge
        // This updates the starting point of pTransEdge to be V1
        // This calls Edge::extend on the twin edge
        pTransEdge->join(pEdge);
        assert(pTransEdge->getDir() == pEdge->getDir());
        pV1->addEdge(pTransEdge); // add to V1

        // Notify the edges they have been updated
        pTransEdge->update();
        pTransEdge->getTwin()->update();
    }

    // Remove the edge from pV1 to pV2
    pV1->removeEdge(pEdge);
    delete pEdge;
    pEdge = 0;

    // Remove the edge from pV2 to pV1
    pV2->removeEdge(pTwin);
    delete pTwin;
    pEdge = 0;

    // Remove V2
    // It is guarenteed to not be connected
    removeIslandVertexNoErase(pV2);
    //validate();
}

//
int Bigraph::sweepVertices(GraphColor c)
{
    int numRemoved = 0;
    VertexPtrMapIter iter = m_vertices.begin();
    while(iter != m_vertices.end())
    {
        VertexPtrMapIter next = iter;
        ++next;
        if(iter->second->getColor() == c)
        {
            removeConnectedVertex(iter->second); 
            ++numRemoved;
        }
        iter = next;
    }
    return numRemoved;
}


//
int Bigraph::sweepEdges(GraphColor c)
{
    int numRemoved = 0;
    for(VertexPtrMapIter iter = m_vertices.begin(); iter != m_vertices.end(); ++iter)
        numRemoved += iter->second->sweepEdges(c);
    return numRemoved;
}

//    Simplify the graph by compacting singular edges
void Bigraph::simplify()
{
    assert(!hasContainment());
    simplify(ED_SENSE);
    simplify(ED_ANTISENSE);
}

// Simplify the graph by compacting edges in the given direction
void Bigraph::simplify(EdgeDir dir)
{
    bool graph_changed = true;
    while(graph_changed)
    {
        graph_changed = false;
        VertexPtrMapIter iter = m_vertices.begin(); 
        while(iter != m_vertices.end())
        {



            // Get the edges for this direction
            EdgePtrVec edges = iter->second->getEdges(dir);

	    //std::cout << "--VERTEX: dir " << dir << " " << iter->first << " seq len " << iter->second->getSeqLen() << " eges size " << edges.size() << std::endl;

	    // for (auto& i : edges)
	    // std::cerr << "    EDGE: " << i->getStartID() << " " << i->getEndID()<< std::endl;

            // If there is a single edge in this direction, merge the vertices
            // Don't merge singular self edges though
            if(edges.size() == 1 && !edges.front()->isSelf())
            {
                // Check that the edge back is singular as well
                Edge* pSingle = edges.front();
                Edge* pTwin = pSingle->getTwin();
                Vertex* pV2 = pSingle->getEnd();
		//std::cerr << "...edge back edge counts " << pV2->countEdges(pTwin->getDir()) << std::endl;
                if(pV2->countEdges(pTwin->getDir()) == 1)
                {
		  //std::cerr << "merging " << std::endl;
                    merge(iter->second, pSingle);
                    graph_changed = true;
                }
            }

            ++iter;
        }
    } 
}

//
// Rename the vertices to have a sequential idx
// starting with prefix. This requires extra memory
// as we need to keep a vector of vertex pointers to reconstruct
// the graph after the vertices are renamed. This should
// only be done after the string graph has been simplified or
// else it could require a lot of memory.
//
void Bigraph::renameVertices(const std::string& prefix)
{
    size_t currIdx = 0;
    size_t numVertices = m_vertices.size();
    std::vector<Vertex*> vertexPtrVec(numVertices, 0);

    VertexPtrMapIter iter = m_vertices.begin();
    while(iter != m_vertices.end())
    {
        std::stringstream ss;
        ss << prefix << currIdx;
	std::string newlab = ss.str();

	// jeremiah
	if (m_get_components) {
	  int label = m_connected_components[iter->second->getID()];
	  if (label) { // is a connected component
	    for (auto& i : m_connected_components) {
	      if (i.second == label)
		m_reads_on_contigs[newlab].push_back(i.first); // i.first is read name of merged read
	    }
	  }
	}

        iter->second->setID(ss.str());
        vertexPtrVec[currIdx] = iter->second;
        ++iter;
        ++currIdx;
    }

    // Clear the old graph
    m_vertices.clear();
    
    // Re-add the vertices
    for(size_t i = 0; i < numVertices; ++i)
        addVertex(vertexPtrVec[i]);
}

//
// Sort the adjacency list for each vertex by length
//
void Bigraph::sortVertexAdjListsByLen()
{
    VertexPtrMapIter iter = m_vertices.begin();
    for(; iter != m_vertices.end(); ++iter)
        iter->second->sortAdjListByLen();
}


//
// Sort the adjacency list for each vertex by ID
//
void Bigraph::sortVertexAdjListsByID()
{
    VertexPtrMapIter iter = m_vertices.begin();
    for(; iter != m_vertices.end(); ++iter)
        iter->second->sortAdjListByID();
}

//
// Validate that the graph is sane
//
void Bigraph::validate()
{
    VertexPtrMapIter iter = m_vertices.begin();
    for(; iter != m_vertices.end(); ++iter)
    {
        iter->second->validate();
    }
}

//
// Flip a vertex
//
void Bigraph::flip(VertexID /*id*/)
{
    assert(false);
#if 0
    // TODO: update this code
    Vertex* pVertex = getVertex(id);
    EdgePtrVec edges = pVertex->getEdges();

    for(EdgePtrVecIter iter = edges.begin(); iter != edges.end(); ++iter)
    {
        // Get the old twin
        GraphEdgeType twin = iter->getTwin();
        
        GraphEdgeType flipped = *iter; 
        flipped.flip();

        // Remove the edge from the source ver
        pVertex->removeEdge(*iter);
        pVertex->addEdge(flipped);

        // Update the partner by deleting the old twin and 
        Vertex* pV2 = getVertex(twin.getStart());
        pV2->removeEdge(twin);
        pV2->addEdge(flipped.getTwin());
    }
#endif
}

//
// Get the IDs of the vertices that do not branch (both sense/antisense degree <= 1)
//
VertexIDVec Bigraph::getNonBranchingVertices() const
{
    VertexIDVec out;
    VertexPtrMapConstIter iter = m_vertices.begin(); 
    for(; iter != m_vertices.end(); ++iter)
    {
        int senseEdges = iter->second->countEdges(ED_SENSE);
        int antisenseEdges = iter->second->countEdges(ED_ANTISENSE);
        if(antisenseEdges <= 1 && senseEdges <= 1)
        {
            out.push_back(iter->second->getID());
        }
    }
    return out;
}


//
// Get all the paths corresponding to the linear components of the graph
// Precondition: all vertices have in/out degree at most 1 (no branching)
//
PathVector Bigraph::getLinearComponents()
{
    PathVector outPaths;
    setColors(GC_WHITE);
    VertexPtrMapIter iter = m_vertices.begin(); 
    for(; iter != m_vertices.end(); ++iter)
    {
        // Output the linear path containing this vertex if it hasnt been visited already
        if(iter->second->getColor() != GC_BLACK)
        {
            outPaths.push_back(constructLinearPath(iter->second->getID()));
        }
    }
    assert(checkColors(GC_BLACK));
    return outPaths;
}

//
// Return all the path of nodes that can be linearally reached from this node
// The path expands in both directions so the first node in the path is not necessarily the source
//
Path Bigraph::constructLinearPath(VertexID id)
{
    Path sensePath;
    Path antisensePath;
    followLinear(id, ED_SENSE, sensePath);
    followLinear(id, ED_ANTISENSE, antisensePath);

    // Construct the final path 
    Path final = reversePath(antisensePath);
    final.insert(final.end(), sensePath.begin(), sensePath.end());
    return final;
}

//
// Recursively follow the graph in the specified direction without branching
// outPath is an out-parameter of the edges that were followed
//
void Bigraph::followLinear(VertexID id, EdgeDir dir, Path& outPath)
{
    Vertex* pVertex = getVertex(id);
    EdgePtrVec edges = pVertex->getEdges(dir);

    // Color the vertex
    pVertex->setColor(GC_BLACK);
    
    if(edges.size() == 1)
    {
        Edge* pSingle = edges.front();
        outPath.push_back(pSingle);
        // Correct the direction for the comp
        assert(pSingle->getDir() == dir);
        EdgeDir corrected_dir = correctDir(pSingle->getDir(), pSingle->getComp());

        // Recurse
        followLinear(pSingle->getEndID(), corrected_dir, outPath);
    }
}

//
Path Bigraph::reversePath(const Path& path)
{
    Path out;
    for(Path::const_reverse_iterator iter = path.rbegin(); iter != path.rend(); ++iter)
        out.push_back((*iter)->getTwin());
    return out;
}

//
Vertex* Bigraph::getFirstVertex() const
{
    if(m_vertices.empty())
        return NULL;
    else
        return m_vertices.begin()->second;
}

// Returns a vector of pointers to the vertices
VertexPtrVec Bigraph::getAllVertices() const
{
    VertexPtrVec out;
    VertexPtrMapConstIter iter = m_vertices.begin();
    for(; iter != m_vertices.end(); ++iter)
        out.push_back(iter->second);
    return out;
}

// Returns a vector of pointers to the vertices
/*VertexPtrVec Bigraph::getAllKeepVertices() const
{
    VertexPtrVec out;
    VertexPtrMapConstIter iter = m_vertices_tokeep.begin();
    for(; iter != m_vertices_tokeep.end(); ++iter)
        out.push_back(iter->second);
    return out;
    }*/

//
// Append vertex sequences to the vector
void Bigraph::getVertexSequences(std::vector<std::string>& outSequences) const
{
    VertexPtrMapConstIter iter = m_vertices.begin();
    for(; iter != m_vertices.end(); ++iter)
        outSequences.push_back(iter->second->getSeq().toString());
}


//
// Visit each vertex and call the function pointer
// The function returns TRUE if it modifies the graph
//
bool Bigraph::visit(VertexVisitFunction f)
{
    bool modified = false;
    VertexPtrMapConstIter iter = m_vertices.begin(); 
    for(; iter != m_vertices.end(); ++iter)
    {
        modified = f(this, iter->second) || modified;
    }
    return modified;
}


//
// Set all the vertices and edges in the graph to the given color
//
void Bigraph::setColors(GraphColor c)
{
    VertexPtrMapIter iter = m_vertices.begin(); 
    for(; iter != m_vertices.end(); ++iter)
    {
        iter->second->setColor(c);
        iter->second->setEdgeColors(c);
    }
}

//
// Get/Set the containment flag
//
void Bigraph::setContainmentFlag(bool b)
{
    m_hasContainment = b;
}

//
bool Bigraph::hasContainment() const
{
    return m_hasContainment;
}

//
// Get/Set the exact mode flag
//
void Bigraph::setExactMode(bool b)
{
    m_isExactMode = b;
}

//
bool Bigraph::isExactMode() const
{
    return m_isExactMode;
}

//
// Get/Set the transitive flag
//
void Bigraph::setTransitiveFlag(bool b)
{
    m_hasTransitive = b;
}

//
bool Bigraph::hasTransitive() const
{
    return m_hasTransitive;
}

//
//
//
void Bigraph::setMinOverlap(int mo)
{
    m_minOverlap = mo;
}

//
int Bigraph::getMinOverlap() const
{
    return m_minOverlap;
}

//
//
//
void Bigraph::setErrorRate(double er)
{
    m_errorRate = er;
}

//
double Bigraph::getErrorRate() const
{
    return m_errorRate;
}

//
// Check if all the vertices in the graph are the given color
//
bool Bigraph::checkColors(GraphColor c)
{
    VertexPtrMapIter iter = m_vertices.begin(); 
    for(; iter != m_vertices.end(); ++iter)
    {
        if(iter->second->getColor() != c)
        {
            std::cerr << "Warning vertex " << iter->second->getID() << " is color " << iter->second->getColor() << " expected " << c << "\n";
            return false;
        }
    }
    return true;
}

//
// Output simple stats
//
void Bigraph::stats() const
{
    int numVerts = 0;
    int numEdges = 0;

    VertexPtrMapConstIter iter = m_vertices.begin(); 
    for(; iter != m_vertices.end(); ++iter)
    {
        numEdges += iter->second->countEdges();
        ++numVerts;
    }

    std::cout << "Graph has " << numVerts << " vertices and " << numEdges << " edges\n";
}

//
//
//
size_t Bigraph::getNumVertices() const
{
    return m_vertices.size();
}

//
// Output mem stats
//
void Bigraph::printMemSize() const
{
    size_t numVerts = 0;
    size_t vertMem = 0;

    size_t numEdges = 0;
    size_t edgeMem = 0;

    VertexPtrMapConstIter iter = m_vertices.begin(); 
    for(; iter != m_vertices.end(); ++iter)
    {
        ++numVerts;
        vertMem += iter->second->getMemSize();

        EdgePtrVec edges = iter->second->getEdges();
        for(EdgePtrVecIter edgeIter = edges.begin(); edgeIter != edges.end(); ++edgeIter)
        {
            ++numEdges;
            edgeMem += (*edgeIter)->getMemSize();
        }
    }
    //printf("num verts: %zu using %zu bytes (%.2lf per vert)\n", numVerts, vertMem, double(vertMem) / numVerts);
    //printf("num edges: %zu using %zu bytes (%.2lf per edge)\n", numEdges, edgeMem, double(edgeMem) / numEdges);
    //printf("total: %zu\n", edgeMem + vertMem);
}

//
// Write the graph to a dot file
//
void Bigraph::writeDot(const std::string& filename, int dotFlags) const
{
    std::ofstream out(filename.c_str());
    
    std::string graphType = (dotFlags & DF_UNDIRECTED) ? "graph" : "digraph";

    out << graphType << " G\n{\n";
    VertexPtrMapConstIter iter = m_vertices.begin(); 
    for(; iter != m_vertices.end(); ++iter)
    {
        VertexID id = iter->second->getID();
        std::string label = (dotFlags & DF_NOID) ? "" : id;
        
        out << "\"" << id << "\" [ label=\"" << label << "\" ";
        if(dotFlags & DF_COLORED)
            out << " style=\"filled\" fillcolor=\"" << getColorString(iter->second->getColor()) << "\" ";
        out << "];\n";
        iter->second->writeEdges(out, dotFlags);
    }
    out << "}\n";
    out.close();
}

//
// Write the graph to an ASQG file
//
void Bigraph::writeASQG(const std::string& filename) const
{
    std::ostream* pWriter = createWriter(filename);
    
    // Header
    ASQG::HeaderRecord headerRecord;
    headerRecord.setOverlapTag(m_minOverlap);
    headerRecord.setErrorRateTag(m_errorRate);
    headerRecord.setTransitiveTag(m_hasTransitive);
    headerRecord.setContainmentTag(m_hasContainment);
    headerRecord.write(*pWriter);


    VertexPtrMapConstIter iter; 

    // Vertices
    for(iter = m_vertices.begin(); iter != m_vertices.end(); ++iter)
    {
        ASQG::VertexRecord vertexRecord(iter->second->getID(), iter->second->getSeq().toString());
        vertexRecord.write(*pWriter);
    }

    // Edges
    for(iter = m_vertices.begin(); iter != m_vertices.end(); ++iter)
    {
        EdgePtrVec edges = iter->second->getEdges();
        for(EdgePtrVecIter edgeIter = edges.begin(); edgeIter != edges.end(); ++edgeIter)
        {
            // We write one record for every bidirectional edge so only write edges
            // that are in canonical form (where id1 < id2)
            Overlap ovr = (*edgeIter)->getOverlap();
            if(ovr.id[0] <= ovr.id[1])
            {
                // Containment edges are in both directions so only output one
                // record if it is a containment
                if(!ovr.isContainment() || ((*edgeIter)->getDir() == ED_SENSE))
                {
                    ASQG::EdgeRecord edgeRecord(ovr);
                    edgeRecord.write(*pWriter);
                }
            }
        }
    }
    delete pWriter;
}

//
std::string Bigraph::getColorString(GraphColor c)
{
    switch(c)
    {
        case GC_WHITE:
            return "white";
        case GC_GRAY:
            return "gray";
        case GC_BLACK:
            return "black";
        case GC_BLUE:
            return "blue";
        case GC_RED:
            return "red";
        default:
            return "black";
    }
}
