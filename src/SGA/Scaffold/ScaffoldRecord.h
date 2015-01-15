//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldRecord - A scaffold consisting of a 
// starting component and a vector of ordered
// links
//
#ifndef SCAFFOLDRECORD_H
#define SCAFFOLDRECORD_H

#include "ScaffoldLink.h"
#include "ScaffoldSequenceCollection.h"
#include "SGUtil.h"

// Gap resolution statistics
struct ResolveStats
{
    ResolveStats()
    {
        numGapsResolved = 0;
        numGapsAttempted = 0;
        numScaffolds = 0;
        
        graphWalkFound = 0;
        graphWalkTooMany = 0;
        graphWalkNoPath = 0;

        overlapFound = 0;
        overlapFailed = 0;
    }

    ~ResolveStats() { print(); }
    
    void print() const
    {
        printf("Num scaffolds: %d\n", numScaffolds);
        printf("Num gaps attempted: %d\n", numGapsAttempted);
        printf("Num gaps resolved: %d\n", numGapsResolved);

        printf("Num gaps resolved by graph walk: %d\n", graphWalkFound);
        printf("Num graph walks failed because of too many solutions: %d\n", graphWalkTooMany);
        printf("Num graph walks failed because of no path: %d\n", graphWalkNoPath);

        printf("Num gaps resolved by overlap: %d\n", overlapFound);
        printf("Num overlaps failed: %d\n", overlapFailed);
    }

    int numGapsResolved;
    int numGapsAttempted;
    int numScaffolds;

    // Graph resolve stats
    int graphWalkFound;
    int graphWalkTooMany;
    int graphWalkNoPath;

    // Overlap resolve stats
    int overlapFound;
    int overlapFailed;
};

// Parameter object for the scaffold record generateString function
struct ResolveParams
{
    ScaffoldSequenceCollection* pSequenceCollection;
    StringGraph* pGraph;

    int minOverlap;
    int maxOverlap;
    double maxErrorRate;
    int resolveMask;
    int minGapLength;
    double distanceFactor;
    ResolveStats* pStats;
};

// Flags indicating what level of gap resolution should be performed
const int RESOLVE_GRAPH_UNIQUE = 1;
const int RESOLVE_GRAPH_BEST = 2;
const int RESOLVE_OVERLAP = 4;

class ScaffoldRecord
{
    public:
        ScaffoldRecord();

        void setRoot(const std::string& root);
        void addLink(const ScaffoldLink& link);
        size_t getNumComponents() const;

        // Generate a sequence string representing the constructed scaffold
        std::string generateString(const ResolveParams& params, StringVector& ids) const;

        // Resolve a link by find walks through the graph
        bool graphResolve(const ResolveParams& params, const std::string& startID, 
                          const ScaffoldLink& link, std::string& extensionString) const;

        // Resolve a predicted overlap between s1 and s2 by aligning the ends of the sequences
        bool overlapResolve(const ResolveParams& params, const std::string& s1, const std::string& s2, 
                            const ScaffoldLink& link, std::string& outString) const;

        // Resolve a link by introducing a gap
        bool introduceGap(int minGapLength, const std::string& contigString, const ScaffoldLink& link, std::string& outString) const;

        void parse(const std::string& text);
        void writeScaf(std::ostream* pWriter);

    private:
        
        typedef std::vector<ScaffoldLink> LinkVector;
        
        std::string m_rootID;
        LinkVector m_links;
        ResolveStats* m_pStats;
};

#endif
