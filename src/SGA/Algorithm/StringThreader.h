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
#ifndef STRING_THREADER_H
#define STRING_THREADER_H

#include <list>
#include "BWT.h"
#include "ExtensionDP.h"

// Typedefs
class StringThreaderNode;
typedef std::list<StringThreaderNode*> STNodePtrList;

// Object to hold the result of the threading process
struct StringThreaderResult
{
    std::string thread;
    int query_align_length;
};
typedef std::vector<StringThreaderResult> StringThreaderResultVector;

// A node in the threading tree
class StringThreaderNode
{
    public:

        //
        // Functions
        // 
        StringThreaderNode(const std::string* pQuery, StringThreaderNode* parent);
        ~StringThreaderNode();
      
        // Add a child node to this node with the given label
        // Returns a pointer to the created node
        StringThreaderNode* createChild(const std::string& label);

        // Extend the label of this node by l
        void extend(const std::string& ext);
        
        // Return a suffix of length l of the string represented by this node
        std::string getSuffix(size_t l) const;

        // Return the complete sequence of the string represented by the branch
        std::string getFullString() const;

        // Get alignment mismatch/error rate metrics
        double getLocalErrorRate(int context) const;
        double getGlobalErrorRate() const;
        int getEditDistance() const;
        
        // Initialize or update the alignment data
        void computeInitialAlignment(const std::string& initialLabel, int queryAlignmentEnd, int bandwidth);
        void computeExtendedAlignment(const std::string& ext, const BandedDPColumn* pPrevColumn);

        // Return the best alignment for the string this node represents
        StringThreaderResult getAlignment() const;

        // Check if this node can be extended any further
        bool hasExtensionTerminated() const;

        // Print the alignment between the string represented by the
        // path from the root of the tree to this node
        void printFullAlignment() const;
        
        // Recursive function to print all the strings represented
        // by this node and all its children.
        void printAllStrings(const std::string& parent) const;


    private:
        
        //
        // Data
        //
        
        // The extension string from the parent
        std::string m_label;

        // The query string being threaded through the graph
        const std::string* m_pQuery;

        // The parent node, can be NULL
        StringThreaderNode* m_pParent;
        STNodePtrList m_children;

        // Alignment information between the label of this node and the query sequence
        // One column per label base
        BandedDPColumnPtrVector m_alignmentColumns;
};

class StringThreader
{
    public:

        //
        // Functions
        //
        StringThreader(const std::string& seed, 
                       const std::string* pQuery,
                       int queryAlignmentEnd,
                       int kmer,
                       const BWT* pBWT);

        ~StringThreader();

        // Run the threading process. Valid alignments are pushed to the results
        // vector
        void run(StringThreaderResultVector& results);

        // Print all the strings represented by the tree
        void printAll();

    private:

        //
        // Functions
        //
        void extendLeaves();

        // Leaf removal heuristics
        void cullLeavesByLocalError();
        void cullLeavesByEdits();
        
        // Check if the leaves can be extended no further
        // If so, the best alignment is pushed to results
        void checkTerminated(StringThreaderResultVector& results);

        // Calculate the successors of this node in the implicit deBruijn graph
        StringVector getDeBruijnExtensions(StringThreaderNode* pNode);
        
        //
        // Data
        //
        const BWT* m_pBWT; 
        int m_kmer;
        const std::string* m_pQuery;
        StringThreaderNode* m_pRootNode;
        STNodePtrList m_leaves;
};

#endif
