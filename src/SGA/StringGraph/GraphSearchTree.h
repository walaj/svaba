//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// GraphSearchTree - Generic data structure used for implementing
// a breadth-first search of a bidirectional graph. It is designed
// to return all possible walks between the given start
// and end vertices, up to a given distance. Used to search a
// string graph or scaffold graph.
//
#ifndef GRAPHSEARCHTREE_H
#define GRAPHSEARCHTREE_H

#include "Bigraph.h"
#include "SGWalk.h"
#include <deque>
#include <queue>

template<typename VERTEX, typename EDGE, typename DISTANCE>
class GraphSearchNode
{
    // Typedefs
    public:
        typedef std::deque<GraphSearchNode<VERTEX,EDGE,DISTANCE>* > GraphSearchNodePtrDeque;
        typedef std::vector<EDGE*> _EDGEPtrVector;

    public:
        GraphSearchNode(VERTEX* pVertex, EdgeDir expandDir, GraphSearchNode* pParent, EDGE* pEdgeFromParent, int distance);
        ~GraphSearchNode();

        // Reduce the child count by 1
        void decrementChildren();

        // Create the children of this node and place pointers to their nodes
        // on the queue. Returns the number of children created;
        int createChildren(GraphSearchNodePtrDeque& outQueue, const DISTANCE& distanceFunc);

        GraphSearchNode* getParent() const { return m_pParent; }
        VERTEX* getVertex() const { return m_pVertex; }
        int64_t getDistance() const { return m_distance; }
        EDGE* getEdgeFromParent() const { return m_pEdgeFromParent; }
        int getNumChildren() const { return m_numChildren; }

    private:

        // data
        VERTEX* m_pVertex;
        EdgeDir m_expandDir;
        GraphSearchNode* m_pParent;
        EDGE* m_pEdgeFromParent;

        int m_numChildren;
        int64_t m_distance;
};

template<typename VERTEX, typename EDGE, typename DISTANCE>
class GraphSearchTree
{
    // typedefs
    typedef GraphSearchNode<VERTEX,EDGE,DISTANCE> _SearchNode;
    typedef typename _SearchNode::GraphSearchNodePtrDeque _SearchNodePtrDeque;
    typedef typename std::set<_SearchNode*> _SearchNodePtrSet;
    typedef std::vector<EDGE*> WALK; // list of edges defines a walk through the graph
    typedef std::vector<WALK> WALKVector; // vector of walks
    
    typedef std::vector<VERTEX*> VertexPtrVector;
    typedef std::vector<VertexPtrVector> VertexPtrVectorVector;

    public:

        GraphSearchTree(VERTEX* pStartVertex, 
                     VERTEX* pEndVertex,
                     EdgeDir searchDir,
                     int64_t distanceLimit,
                     size_t nodeLimit);

        ~GraphSearchTree();

        // Find connected components in the graph
        // Takes in a vector of all the vertices in the graph
        static void connectedComponents(VertexPtrVector allVertices, VertexPtrVectorVector& connectedComponents);

        // Returns true if the search has converged on a single vertex. In
        // other words, all walks from the start node share a common vertex,
        // which is represented by one of the nodes waiting expansion.
        // This function is the key to finding collapsed walks that represent
        // complex variation bubbles. It is a heavy operation however, as the
        // entire tree is search for each node in the expand queue. It should
        // only be used on small trees.
        // Returns true if the search converged and the pointer to the vertex
        // is return in pConvergedVertex.
        bool hasSearchConverged(VERTEX*& pConvergedVertex);

        // build walks to a given vertex
        template<typename BUILDER>
        void buildWalksContainingVertex(VERTEX* pTarget, BUILDER& walkBuilder);

        // Construct walks representing every path from the start node
        template<typename BUILDER>
        void buildWalksToAllLeaves(BUILDER& walkBuilder);

        // Construct walks representing every path from the start vertex to the goal vertex
        template<typename BUILDER>
        void buildWalksToGoal(BUILDER& walkBuilder);
                
        // Expand all nodes in the queue a single time
        // Returns false if the search has stopped
        bool stepOnce();

        // Check whether the search was aborted or not
        bool wasSearchAborted() const { return m_searchAborted; }

    private:

        // Search the branch from pNode to the root for pX.  
        bool searchBranchForVertex(_SearchNode* pNode, VERTEX* pX, _SearchNode*& pFoundNode) const;

        // Build the walks from the root to the leaves in the queue
        template<typename BUILDER>
        void _buildWalksToLeaves(const _SearchNodePtrDeque& queue, BUILDER& walkBuilder);

        //
        void addEdgesFromBranch(_SearchNode* pNode, WALK& outEdges);

        // Build a queue with all the leaves in it
        void _makeFullLeafQueue(_SearchNodePtrDeque& completeQueue) const;

        // print the branch sequence
        void printBranch(_SearchNode* pNode) const;

        // We keep the pointers to the search nodes
        // in one of three queues. 
        // The goal queue contains the nodes representing the vertex we are searching for.
        // The expand queue contains nodes that have not yet been explored.
        // The done queue contains non-goal nodes that will not be expanded further.
        // Together, they represent all leaves of the tree
        _SearchNodePtrDeque m_goalQueue;
        _SearchNodePtrDeque m_expandQueue;
        _SearchNodePtrDeque m_doneQueue;
        size_t m_totalNodes; // The total number of nodes in the search tree
    
        _SearchNode* m_pRootNode;
        VERTEX* m_pGoalVertex;
        EdgeDir m_initialDir;

        int64_t m_distanceLimit;
        size_t m_nodeLimit;

        // Flag indicating the search was aborted
        bool m_searchAborted;

        // Distance functor
        DISTANCE m_distanceFunc;
};

//
template<typename VERTEX, typename EDGE, typename DISTANCE>
GraphSearchNode<VERTEX,EDGE,DISTANCE>::GraphSearchNode(VERTEX* pVertex,
                           EdgeDir expandDir,
                           GraphSearchNode<VERTEX,EDGE,DISTANCE>* pParent,
                           EDGE* pEdgeFromParent,
                           int distance) : m_pVertex(pVertex),
                                                    m_expandDir(expandDir),
                                                    m_pParent(pParent), 
                                                    m_pEdgeFromParent(pEdgeFromParent),
                                                    m_numChildren(0)
{
    // Set the extension distance
    if(m_pParent == NULL)
    {
        // this is the root node with distance 0
        m_distance = 0;
    }
    else
    {
        assert(m_pEdgeFromParent != NULL);
        m_distance = m_pParent->m_distance + distance;
    }
}

// Delete this node and decrement the number of children
// in the parent node. All children of a node must
// be deleted before the parent
template<typename VERTEX, typename EDGE, typename DISTANCE>
GraphSearchNode<VERTEX,EDGE,DISTANCE>::~GraphSearchNode()
{
    assert(m_numChildren == 0);
    if(m_pParent != NULL)
        m_pParent->decrementChildren();
}

//
template<typename VERTEX, typename EDGE, typename DISTANCE>
void GraphSearchNode<VERTEX,EDGE,DISTANCE>::decrementChildren()
{
    assert(m_numChildren != 0);
    m_numChildren -= 1;
}

// creates nodes for the children of this node
// and place pointers to them in the queue.
// Returns the number of nodes created
template<typename VERTEX, typename EDGE, typename DISTANCE>
int GraphSearchNode<VERTEX,EDGE,DISTANCE>::createChildren(GraphSearchNodePtrDeque& outDeque, const DISTANCE& distanceFunc)
{
    assert(m_numChildren == 0);

    _EDGEPtrVector edges = m_pVertex->getEdges(m_expandDir);

    for(size_t i = 0; i < edges.size(); ++i)
    {
        EdgeDir childExpandDir = !edges[i]->getTwin()->getDir();
        GraphSearchNode* pNode = new GraphSearchNode(edges[i]->getEnd(), childExpandDir, this, edges[i], distanceFunc(edges[i]));
        outDeque.push_back(pNode);
        m_numChildren += 1;
    }
    return edges.size();
}

//
// GraphSearchTree
//
template<typename VERTEX, typename EDGE, typename DISTANCE>
GraphSearchTree<VERTEX,EDGE,DISTANCE>::GraphSearchTree(VERTEX* pStartVertex, 
                                                       VERTEX* pEndVertex, 
                                                       EdgeDir searchDir,
                                                       int64_t distanceLimit,
                                                       size_t nodeLimit) : m_pGoalVertex(pEndVertex),
                                                                           m_distanceLimit(distanceLimit),
                                                                           m_nodeLimit(nodeLimit),
                                                                           m_searchAborted(false)
{
    // Create the root node of the search tree
    m_pRootNode = new GraphSearchNode<VERTEX, EDGE, DISTANCE>(pStartVertex, searchDir, NULL, NULL, 0);

    // add the root to the expand queue
    m_expandQueue.push_back(m_pRootNode);

    m_totalNodes = 1;
}

template<typename VERTEX, typename EDGE, typename DISTANCE>
GraphSearchTree<VERTEX,EDGE,DISTANCE>::~GraphSearchTree()
{
    // Delete the tree
    // We delete each leaf and recurse up the tree iteratively deleting
    // parents with a single child node. This ensure that each parent is
    // deleted after all its children
    _SearchNodePtrDeque completeLeafNodes;
    _makeFullLeafQueue(completeLeafNodes);

    size_t totalDeleted = 0;
    for(typename _SearchNodePtrDeque::iterator iter = completeLeafNodes.begin(); 
                                               iter != completeLeafNodes.end();
                                               ++iter)    
    {
        _SearchNode* pCurr = *iter;

        // loop invariant: pCurr is a deletable node
        // the loop stops when the parent is NULL or has 
        // a child other than pCurr
        do
        {
            assert(pCurr->getNumChildren() == 0);
            _SearchNode* pNext = pCurr->getParent();
            
            delete pCurr; // decrements pNext's child count
            totalDeleted += 1;

            pCurr = pNext;
        } while(pCurr && pCurr->getNumChildren() == 0);
    }

    assert(totalDeleted == m_totalNodes);
}

// Perform one step of the BFS
template<typename VERTEX, typename EDGE, typename DISTANCE>
bool GraphSearchTree<VERTEX,EDGE,DISTANCE>::stepOnce()
{
    if(m_expandQueue.empty())
        return false;

    if(m_totalNodes > m_nodeLimit)
    {
        // Move all nodes in the expand queue to the done queue
        m_doneQueue.insert(m_doneQueue.end(), m_expandQueue.begin(), m_expandQueue.end());
        m_expandQueue.clear();

        // Set a flag indicating the search was aborted
        m_searchAborted = true;
        return false;
    }

    // Iterate over the expand queue. If the path to the node
    // is outside the depth limit, move that node to the done queue. It cannot
    // yield a valid path to the goal. Otherwise, add the children of the node 
    // to the incoming queue
    _SearchNodePtrDeque incomingQueue;
    while(!m_expandQueue.empty())
    {
        _SearchNode* pNode = m_expandQueue.front();
        m_expandQueue.pop_front();
        
        if(pNode->getVertex() == m_pGoalVertex)
        {
            // This node represents the goal, add it to the goal queue
            m_goalQueue.push_back(pNode);
            continue;
        }

        if(pNode->getDistance() > m_distanceLimit)
        {
            // Path to this node is too long, expand it no further
            m_doneQueue.push_back(pNode);
        }
        else
        {
            // Add the children of this node to the queue
            int numCreated = pNode->createChildren(incomingQueue, m_distanceFunc);
            m_totalNodes += numCreated;

            if(numCreated == 0)
            {
                // No children created, add this node to the done queue
                m_doneQueue.push_back(pNode);
            }
        }
    }

    m_expandQueue = incomingQueue;
    return true;
}

// Return true if all the walks from the root converge
// to one vertex (ie if the search from pX converged
// to pY, then ALL paths from pX must go through pY).
template<typename VERTEX, typename EDGE, typename DISTANCE>
bool GraphSearchTree<VERTEX,EDGE,DISTANCE>::hasSearchConverged(VERTEX*& pConvergedVertex)
{
    // Construct a set of all the leaf nodes
    _SearchNodePtrDeque completeLeafNodes;
    _makeFullLeafQueue(completeLeafNodes);

    // Search all the tree for all the nodes in the expand queue
    for(typename _SearchNodePtrDeque::iterator iter = m_expandQueue.begin(); 
                                               iter != m_expandQueue.end();
                                               ++iter)
    {
        _SearchNode* pNode = *iter;
        // If this node has the same vertex as the root skip it
        // We do not want to collapse at the root
        if(pNode->getVertex() == m_pRootNode->getVertex())
            continue;

        bool isInAllBranches = true;
        for(typename _SearchNodePtrDeque::iterator leafIter = completeLeafNodes.begin();
                                                   leafIter != completeLeafNodes.end();
                                                   ++leafIter)
        {
            // Search the current branch from this leaf node to the root
            _SearchNode* pFoundNode = NULL;
            bool isInBranch = searchBranchForVertex(*leafIter, pNode->getVertex(), pFoundNode);
            if(!isInBranch)
            {
                isInAllBranches = false;
                break;
            }
        }
    
        // search has converted
        if(isInAllBranches)
        {
            pConvergedVertex = pNode->getVertex();
            return true;
        }
    }

    // search has not converted
    pConvergedVertex = NULL;
    return false;
}

// Construct walks representing every path from the start node
template<typename VERTEX, typename EDGE, typename DISTANCE>
template<typename BUILDER>
void GraphSearchTree<VERTEX,EDGE,DISTANCE>::buildWalksToAllLeaves(BUILDER& walkBuilder)
{
    // Construct a queue with all leaf nodes in it
    _SearchNodePtrDeque completeLeafNodes;
    _makeFullLeafQueue(completeLeafNodes);

    _buildWalksToLeaves(completeLeafNodes, walkBuilder);
}

// Construct walks representing every path from the start vertex to the goal vertex
template<typename VERTEX, typename EDGE, typename DISTANCE>
template<typename BUILDER>
void GraphSearchTree<VERTEX,EDGE,DISTANCE>::buildWalksToGoal(BUILDER& walkBuilder)
{
    _buildWalksToLeaves(m_goalQueue, walkBuilder);
}

// Build all the walks that contain pTarget.
template<typename VERTEX, typename EDGE, typename DISTANCE>
template<typename BUILDER>
void GraphSearchTree<VERTEX,EDGE,DISTANCE>::buildWalksContainingVertex(VERTEX* pTarget, BUILDER& walkBuilder)
{
    _SearchNodePtrDeque completeLeafNodes;
    _makeFullLeafQueue(completeLeafNodes);

    // Search upwards from each leaf until pTarget is found.
    // When it is found, insert the pointer to the search node
    // in the set
    _SearchNodePtrSet leafSet;

    // Find pTarget in each branch of the graph
    for(typename _SearchNodePtrDeque::const_iterator iter = completeLeafNodes.begin();
                                                     iter != completeLeafNodes.end();
                                                     ++iter)
    {
        _SearchNode* pFoundNode = NULL;
        searchBranchForVertex(*iter, pTarget, pFoundNode);
        assert(pFoundNode != NULL);
        leafSet.insert(pFoundNode);
    }

    // Construct all the walks to the found leaves
    completeLeafNodes.clear();
    completeLeafNodes.insert(completeLeafNodes.end(), leafSet.begin(), leafSet.end());
    _buildWalksToLeaves(completeLeafNodes, walkBuilder);
}

// Main function for constructing a vector of walks from a set of leaves
template<typename VERTEX, typename EDGE, typename DISTANCE>
template<typename BUILDER>
void GraphSearchTree<VERTEX,EDGE,DISTANCE>::_buildWalksToLeaves(const _SearchNodePtrDeque& queue, BUILDER& walkBuilder)
{
    for(typename _SearchNodePtrDeque::const_iterator iter = queue.begin();
                                                     iter != queue.end();
                                                     ++iter)
    {

        // Recursively travel the tree from the leaf to the root collecting the edges in the vector
        WALK currWalk;
        addEdgesFromBranch(*iter, currWalk);

        // Reverse the walk and write it to the output structure
        walkBuilder.startNewWalk(m_pRootNode->getVertex());
        for(typename WALK::reverse_iterator iter = currWalk.rbegin(); iter != currWalk.rend(); ++iter)
            walkBuilder.addEdge(*iter);
        walkBuilder.finishCurrentWalk();
    }
}

// Return true if the vertex pX is found somewhere in the branch 
// from pNode to the root. If it is found, pFoundNode is set
// to the furtherest instance of pX from the root.
template<typename VERTEX, typename EDGE, typename DISTANCE>
bool GraphSearchTree<VERTEX,EDGE,DISTANCE>::searchBranchForVertex(_SearchNode* pNode, VERTEX* pX, _SearchNode*& pFoundNode) const
{
    if(pNode == NULL)
    {
        pFoundNode = NULL;
        return false;
    }

    if(pNode->getVertex() == pX && pNode != m_pRootNode)
    {
        pFoundNode = pNode;
        return true;
    }
    return searchBranchForVertex(pNode->getParent(), pX, pFoundNode);
}

//
template<typename VERTEX, typename EDGE, typename DISTANCE>
void GraphSearchTree<VERTEX,EDGE,DISTANCE>::addEdgesFromBranch(_SearchNode* pNode, WALK& outEdges)
{
    // Terminate the recursion at the root node and dont add an edge
    if(pNode->getParent() != NULL)
    {
        outEdges.push_back(pNode->getEdgeFromParent());
        return addEdgesFromBranch(pNode->getParent(), outEdges);
    }
}

//
template<typename VERTEX, typename EDGE, typename DISTANCE>
void GraphSearchTree<VERTEX,EDGE,DISTANCE>::_makeFullLeafQueue(_SearchNodePtrDeque& completeQueue) const
{
    completeQueue.insert(completeQueue.end(), m_expandQueue.begin(), m_expandQueue.end());
    completeQueue.insert(completeQueue.end(), m_goalQueue.begin(), m_goalQueue.end());
    completeQueue.insert(completeQueue.end(), m_doneQueue.begin(), m_doneQueue.end());
}

//
template<typename VERTEX, typename EDGE, typename DISTANCE>
void GraphSearchTree<VERTEX,EDGE,DISTANCE>::printBranch(_SearchNode* pNode) const
{
    if(pNode != NULL)
    {
        std::cout << pNode->getVertex()->getID() << ",";
        printBranch(pNode->getParent());
    }
}

template<typename VERTEX, typename EDGE, typename DISTANCE>
void GraphSearchTree<VERTEX,EDGE,DISTANCE>::connectedComponents(VertexPtrVector allVertices, 
                                                                VertexPtrVectorVector& connectedComponents)
{
    // Set the color of each vertex to be white signalling its not visited
    typename VertexPtrVector::iterator iter = allVertices.begin();
    for(; iter != allVertices.end(); ++iter)
    {
        assert((*iter)->getColor() == GC_WHITE);
        (*iter)->setColor(GC_WHITE);
    }

    // 
    iter = allVertices.begin();
    for(; iter != allVertices.end(); ++iter)
    {
        // Do nothing if this vertex is already part of a CC
        if((*iter)->getColor() == GC_BLACK)
            continue;

       // Start a new CC
       VertexPtrVector currComponent;

       std::queue<VERTEX*> exploreQueue;
       (*iter)->setColor(GC_GRAY); // queued color
       exploreQueue.push(*iter);

       while(!exploreQueue.empty())
       {
            VERTEX* pCurr = exploreQueue.front();
            exploreQueue.pop();

            assert(pCurr->getColor() != GC_BLACK);
            currComponent.push_back(pCurr);
            pCurr->setColor(GC_BLACK); //done with this vertex

            // Enqueue edges if they havent been visited already
            std::vector<EDGE*> edges = pCurr->getEdges();
            for(size_t i = 0; i < edges.size(); ++i)
            {
                EDGE* pEdge = edges[i];
                VERTEX* pNext = pEdge->getEnd();
                if(pNext->getColor() == GC_WHITE)
                {
                    pNext->setColor(GC_GRAY); // queued
                    exploreQueue.push(pNext);
                }
            }
       }

        connectedComponents.push_back(currComponent);
    }

    iter = allVertices.begin();
    for(; iter != allVertices.end(); ++iter)
    {
        (*iter)->setColor(GC_WHITE);
    }

    // Sanity check
    size_t totalVertices = allVertices.size();
    size_t totalInComponents = 0;

    for(size_t i = 0; i < connectedComponents.size(); ++i)
    {
        totalInComponents += connectedComponents[i].size();
    }
    assert(totalVertices == totalInComponents);
    //std::cout << "[CC] total: " << totalVertices << " num components: " << connectedComponents.size() << "\n";
    //std::cout << "[CC] total vertices in components: " << totalInComponents << "\n";
}

#endif
