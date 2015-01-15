//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//
// SearchHistory - data structures to track
// the history of a FM-index search
//-----------------------------------------------
#ifndef SEARCHHISTORY_H
#define SEARCHHISTORY_H

#include "Util.h"

// Base, Position pair indicating a divergence during the search
struct SearchHistoryItem
{
    SearchHistoryItem(int p, char b) : pos(p), base(b) {}
    int pos;
    char base;

    static bool sortPos(const SearchHistoryItem& a, const SearchHistoryItem& b)
    {
        return a.pos < b.pos;
    }

    friend std::ostream& operator<<(std::ostream& out, const SearchHistoryItem& hi)
    {
        out << hi.pos << "," << hi.base;
        return out;
    }
};
typedef std::vector<SearchHistoryItem> HistoryItemVector;

// A vector of history items that can be compared with other histories
class SearchHistoryVector
{
    public:
        
        //
        SearchHistoryVector() {}

        //
        void add(int pos, char base);
        void add(SearchHistoryItem& item);
        void normalize(bool doComplement);
        size_t size() const { return m_history.size(); }
        std::string getBaseString() const;

        // Transform the original string using the history into a string representing
        // the sequence of bases searched to create the history
        std::string transform(const std::string& original, bool queryReversed) const;

        //
        static int countDifferences(const SearchHistoryVector& a, const SearchHistoryVector& b, int maxPos);

        //
        friend std::ostream& operator<<(std::ostream& out, const SearchHistoryVector& hist);

    private:

        HistoryItemVector m_history;
};


class SearchHistoryNode;

// A SearchHistoryLink is a reference-counted wrapper of a 
// search node. This is the external interface to the SearchHistoryNodes
// This allows the SearchHistoryNodes to be automatically cleaned up when 
// they are no longer referred to
class SearchHistoryLink
{
    public:
        
        //
        SearchHistoryLink();
        SearchHistoryLink(SearchHistoryNode* ptr);

        // We need to handle the copy constructor and the assignment operator
        // for the reference counting to be correct
        SearchHistoryLink(const SearchHistoryLink& link);        
        SearchHistoryLink& operator=(SearchHistoryLink const& link);
        ~SearchHistoryLink();


        SearchHistoryNode* operator->() const { assert(pNode != NULL); return pNode; }
        SearchHistoryNode& operator*() const { assert(pNode != NULL); return *pNode; }

    private:
        SearchHistoryNode* pNode;
};

// A search history node is one link in a chain of history items. createChild makes a 
// new element in the chain, indicating a divergence from the parent history
class SearchHistoryNode
{
    public:

        SearchHistoryLink createChild(int var_pos, char var_base);
        static SearchHistoryLink createRoot(); // Create the root node of the history tree
        SearchHistoryVector getHistoryVector();

    private:

        friend class SearchHistoryLink;
        friend class SearchTree;

        // The nodes should only be constructed/destructed through the links
        SearchHistoryNode(SearchHistoryNode* pParent, 
                          int var_pos, char var_base) : m_parentLink(pParent), 
                                                        m_variant(var_pos, var_base),
                                                        m_refCount(0) {}
        
        ~SearchHistoryNode() { assert(m_refCount == 0); }

        inline void increment() { ++m_refCount; }
        inline void decrement() { --m_refCount; }
        inline int getCount() const { return m_refCount; }

        SearchHistoryLink m_parentLink;
        SearchHistoryItem m_variant;
        int m_refCount;

        static const char ROOT_CHAR = '0';
};

#endif
