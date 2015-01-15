///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// StringThreader - Iteratively construct a
// string representing a walk through an assembly graph
// matching a query sequence. 
//
// The assembly graph is abstractly represented as
// an FM-index.
//
#include "StringThreader.h"
#include "BWTAlgorithms.h"
#include "LRAlignment.h"
#include "StdAlnTools.h"

//
// StringThreaderNode
//
StringThreaderNode::StringThreaderNode(const std::string* pQuery,
                                       StringThreaderNode* parent) : m_pQuery(pQuery),
                                                                     m_pParent(parent)
{

}

// Destructor, recurisvely delete the children of the node
StringThreaderNode::~StringThreaderNode()
{
    // Delete children
    for(STNodePtrList::iterator iter = m_children.begin(); iter != m_children.end(); ++iter)
        delete *iter;

    // Delete alignment columns
    for(size_t i = 0; i < m_alignmentColumns.size(); ++i)
        delete m_alignmentColumns[i];
}

// Return a suffix of length l of the path from the root to this node
std::string StringThreaderNode::getSuffix(size_t l) const
{
    size_t n = m_label.size();
    if(l <= n)
    {
        return m_label.substr(n - l, l);
    }
    else 
    {
        assert(m_pParent != NULL);
        return m_pParent->getSuffix(l - n) + m_label;
    }
}

// Return the full string of the path from the root to this node
std::string StringThreaderNode::getFullString() const
{
    if(m_pParent == NULL)
        return m_label;
    else
        return m_pParent->getFullString() + m_label;
}

// Create a new child node with the given label. Returns a pointer to the new node.
StringThreaderNode* StringThreaderNode::createChild(const std::string& label)
{
    StringThreaderNode* pAdded = new StringThreaderNode(m_pQuery, this);
    m_children.push_back(pAdded);

    assert(!m_alignmentColumns.empty());
    pAdded->computeExtendedAlignment(label, m_alignmentColumns.back());
    return pAdded;
}

// Extend the label of this node
void StringThreaderNode::extend(const std::string& ext)
{
    assert(!ext.empty());
    assert(!m_alignmentColumns.empty());
    computeExtendedAlignment(ext, m_alignmentColumns.back());
}
   
// Update the alignment columns for this node.
void StringThreaderNode::computeExtendedAlignment(const std::string& ext, const BandedDPColumn* pPrevColumn)
{
    // Calculate a new alignment column for each base of the label
    for(size_t i = 0; i < ext.size(); ++i)
    {
        BandedDPColumn* pNewColumn = ExtensionDP::createNewColumn(ext[i], *m_pQuery, pPrevColumn);
        m_alignmentColumns.push_back(pNewColumn);
        pPrevColumn = pNewColumn;
    }
    m_label.append(ext);
}

// Initialize the alignment columns. 
void StringThreaderNode::computeInitialAlignment(const std::string& initialLabel, int queryAlignmentEnd, int bandwidth)
{
    // Create the initial alignment columnds between label and query
    m_label = initialLabel;
    assert(!m_label.empty());
    ExtensionDP::createInitialAlignment(m_label, m_pQuery->substr(0, queryAlignmentEnd), bandwidth, m_alignmentColumns);
}

// Calculate error rate over last context bases of the alignment
double StringThreaderNode::getLocalErrorRate(int context) const
{
    return ExtensionDP::calculateLocalEditPercentage(m_alignmentColumns.back(), context);
}

// Calculate error rate over the entire alignment
double StringThreaderNode::getGlobalErrorRate() const
{
    return ExtensionDP::calculateGlobalEditPercentage(m_alignmentColumns.back());
}

// Calculate the edit distance between the thread and query
int StringThreaderNode::getEditDistance() const
{
    int edits, alignLength;
    ExtensionDP::countEditsAndAlignLength(m_alignmentColumns.back(), edits, alignLength);
    return edits;
}

// Returns true if the extension has terminated
bool StringThreaderNode::hasExtensionTerminated() const
{
    return ExtensionDP::isExtensionTerminated(m_alignmentColumns.back(), 2);
}

// Return the best alignment between the string represented by this node and the query
StringThreaderResult StringThreaderNode::getAlignment() const
{
    ExtensionDPAlignment alignment = ExtensionDP::findGlocalAlignment(m_alignmentColumns.back());
    StringThreaderResult result;
    result.query_align_length = alignment.query_align_length;
    result.thread =  getFullString().substr(0, alignment.target_align_length);
    return result;
}

// Print the alignment
void StringThreaderNode::printFullAlignment() const
{
    std::string fullString = getFullString();
    ExtensionDP::printAlignment(fullString, *m_pQuery, m_alignmentColumns.back());
}

// Print the string(s) represented by this node and its children
void StringThreaderNode::printAllStrings(const std::string& parent) const
{
    if(m_children.empty())
    {
        std::cout << "S: " << parent + m_label << "\n";
    }
    else
    {
        for(STNodePtrList::const_iterator iter = m_children.begin(); iter != m_children.end(); ++iter)
            (*iter)->printAllStrings(parent + m_label);
    }
}

//
// StringTheader
//
StringThreader::StringThreader(const std::string& seed, 
                               const std::string* pQuery,
                               int queryAlignmentEnd,
                               int kmer, 
                               const BWT* pBWT) : m_pBWT(pBWT), m_kmer(kmer), m_pQuery(pQuery)
{
    // Create the root node containing the seed string
    m_pRootNode = new StringThreaderNode(pQuery, NULL);
    m_pRootNode->computeInitialAlignment(seed, queryAlignmentEnd, 50);
    m_leaves.push_back(m_pRootNode);
}

//
StringThreader::~StringThreader()
{
    // Recursively destroy the tree
    delete m_pRootNode;
}

// Run the threading algorithm
void StringThreader::run(StringThreaderResultVector& results)
{
    // Extend the leaf nodes
    while(!m_leaves.empty())
    {
        extendLeaves();
        cullLeavesByEdits();
        checkTerminated(results);
    }
}

// Print the string represented by every node
void StringThreader::printAll()
{
    std::cout << "Print all: \n";
    m_pRootNode->printAllStrings("");
}

// Extend each leaf node
void StringThreader::extendLeaves()
{
    STNodePtrList newLeaves;
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        StringVector extensions = getDeBruijnExtensions(*iter);

        // Either extend the current node or branch it
        // If no extension, do nothing and this node
        // is no longer considered a leaf
        if(extensions.size() == 1)
        {
            // Single extension, do not branch
            (*iter)->extend(extensions.front());
            newLeaves.push_back(*iter);
        }
        else if(extensions.size() > 1)
        {
            // Branch
            for(size_t i = 0; i < extensions.size(); ++i)
            {
                StringThreaderNode* pAdded = (*iter)->createChild(extensions[i]);
                newLeaves.push_back(pAdded);
            }
        }
    }

    m_leaves.clear();
    m_leaves = newLeaves;
}

// Remove leaves that are have a high local error rate
// These leaves probably represent branches that are not
// the query sequence.
void StringThreader::cullLeavesByLocalError()
{
    STNodePtrList newLeaves;

    int context = 20;
    double threshold = 0.3f;

    // Calculate the local error rate of the alignments to each new leaf
    // If it is less than threshold, add the leaf to the node
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        double ler = (*iter)->getLocalErrorRate(context);

        if(ler < threshold)
            newLeaves.push_back(*iter);
    }
    m_leaves = newLeaves;
}

// Remove leaves that have an edit distance much higher than
// the best leaf
void StringThreader::cullLeavesByEdits()
{
    STNodePtrList newLeaves;

    // Calculate the local error rate of the alignments to each new leaf
    // If it is less than threshold, add the leaf to the node
    int bestEdits = std::numeric_limits<int>::max();
    IntVector editsVector;
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        int edits = (*iter)->getEditDistance();
        if(edits < bestEdits)
            bestEdits = edits;
         editsVector.push_back(edits);
    }

    int leafID = 0;
    int threshold = 1;
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        if(editsVector[leafID] <= bestEdits + threshold)
            newLeaves.push_back(*iter);
        leafID += 1;
    }

    m_leaves = newLeaves;
}

// Check for leaves whose extension has terminated. If the leaf has
// terminated, its alignment result is pushed to the result vector
void StringThreader::checkTerminated(StringThreaderResultVector& results)
{
    STNodePtrList newLeaves;
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        if((*iter)->hasExtensionTerminated())
            results.push_back((*iter)->getAlignment());
        else
            newLeaves.push_back(*iter);
    }
    m_leaves = newLeaves;
}

// Calculate the successors of this node in the implicit deBruijn graph
StringVector StringThreader::getDeBruijnExtensions(StringThreaderNode* pNode)
{
    // Get the last k-1 bases of the node
    std::string pmer = pNode->getSuffix(m_kmer - 1);
    AlphaCount64 extensions = BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(pmer, m_pBWT, ED_SENSE, NULL);

    // Loop over the DNA symbols, if there is are more than two characters create a branch
    // otherwise just perform an extension.
    bool hasExtension = extensions.hasDNAChar();

    StringVector out;
    if(hasExtension)
    {
        for(int i = 0; i < DNA_ALPHABET::size; ++i)
        {
            char b = DNA_ALPHABET::getBase(i);
            if(extensions.get(b) >= 3)
            {
                // extend to b
                std::string tmp;
                tmp.append(1,b);
                out.push_back(tmp);
            }
        }
    }

    return out;
}
