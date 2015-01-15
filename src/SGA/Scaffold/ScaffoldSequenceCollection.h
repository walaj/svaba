//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldSequenceCollection - A set of input
// sequences that are being scaffolded. Can
// either be held in a string graph or just a map
// from ID to sequence. Supports setting particular
// sequences as being placed in a scaffold
//
#ifndef SCAFFOLDSEQUENCECOLLECTION_H
#define SCAFFOLDSEQUENCECOLLECTION_H

#include "SGUtil.h"
#include <string>
#include <map>

//
class ScaffoldSequenceCollection
{
    public:
        ScaffoldSequenceCollection() {}
        
        virtual ~ScaffoldSequenceCollection() {}

        // Returns the sequence with the given ID
        virtual std::string getSequence(const std::string& id) const = 0;

        // Mark the sequence with id as being placed in a scaffold
        virtual void setPlaced(const std::string& id) = 0;
        
        // write the unplaced sequences of length at least minLength using pWriter
        virtual void writeUnplaced(std::ostream* pWriter, int minLength) = 0;
};

// ScaffoldSequenceCollection implemented using a string graph
// as the base storage
class GraphSequenceCollection : public ScaffoldSequenceCollection
{
    public:

        //
        GraphSequenceCollection(StringGraph* pGraph);
        ~GraphSequenceCollection() {}
        
        // Returns the sequence with the given ID
        std::string getSequence(const std::string& id) const;

        // 
        void setPlaced(const std::string& id);

        // write the unplaced sequences of length at least minLength using pWriter
        void writeUnplaced(std::ostream* pWriter, int minLength);

    private:
        
        StringGraph* m_pGraph;
};

// ScaffoldSequenceCollection implemented using a simple map
class MapSequenceCollection : public ScaffoldSequenceCollection
{
    public:

        // Read the sequences from a fasta file
        MapSequenceCollection(std::string filename);
        ~MapSequenceCollection() {}
        
        // Returns the sequence with the given ID
        std::string getSequence(const std::string& id) const;

        // 
        void setPlaced(const std::string& id);

        // write the unplaced sequences of length at least minLength using pWriter
        void writeUnplaced(std::ostream* pWriter, int minLength);

    private:
        
        struct SequenceMapData
        {
            std::string sequence;
            bool isPlaced;
        };
        typedef std::map<std::string, SequenceMapData> SMPMap;
        SMPMap m_map;
};

#endif
