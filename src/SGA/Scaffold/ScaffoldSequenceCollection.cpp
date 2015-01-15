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
// sequences as being visited
//
#include "ScaffoldSequenceCollection.h"
#include "SeqReader.h"

//
GraphSequenceCollection::GraphSequenceCollection(StringGraph* pGraph) : m_pGraph(pGraph)
{

}
        
// Returns the sequence with the given ID
std::string GraphSequenceCollection::getSequence(const std::string& id) const
{
    Vertex* pVertex = m_pGraph->getVertex(id);
    if(pVertex == NULL)
    {
        std::cerr << "Error: could not find the contig record with id " << id << "\n";
        std::cerr << "Check to make sure you have loaded the correct set of contigs\n";
        exit(EXIT_FAILURE);
    }

    return pVertex->getSeq().toString();
}

// 
void GraphSequenceCollection::setPlaced(const std::string& id)
{
    Vertex* pVertex = m_pGraph->getVertex(id);
    assert(pVertex != NULL);
    pVertex->setColor(GC_BLACK);
}

// Small visitor function to write nodes of the graph that are not colored black
struct UnplacedVisitor
{
    UnplacedVisitor(std::ostream* pWriter, int minLength) : m_pWriter(pWriter), m_numUnplaced(0), m_minLength(minLength) {}

    void previsit(StringGraph*) {}
    
    bool visit(StringGraph*, Vertex* pVertex)
    {
        if(pVertex->getColor() == GC_BLACK)
            return false;
        
        if((int)pVertex->getSeqLen() >= m_minLength)
        {
            std::stringstream idss;
            idss << "unplaced-" << m_numUnplaced++;
            writeFastaRecord(m_pWriter, idss.str(), pVertex->getSeq().toString());
        }
        return false;
    }

    void postvisit(StringGraph*) {}

    std::ostream* m_pWriter;
    int m_numUnplaced;
    int m_minLength;
};

// Write the unplaced sequences using a graph visitor
void GraphSequenceCollection::writeUnplaced(std::ostream* pWriter, int minLength)
{
    UnplacedVisitor uvisit(pWriter, minLength);
    m_pGraph->visit(uvisit);
}

// Read the sequences from the file
MapSequenceCollection::MapSequenceCollection(std::string filename)
{
    SeqReader reader(filename, SRF_NO_VALIDATION | SRF_KEEP_CASE);
    SeqRecord record;
    while(reader.get(record))
    {
        SequenceMapData& md = m_map[record.id];
        md.sequence = record.seq.toString();
        md.isPlaced = false;
    }
}

// Returns the sequence with the given ID
std::string MapSequenceCollection::getSequence(const std::string& id) const
{
    SMPMap::const_iterator iter = m_map.find(id);
    if(iter == m_map.end())
    {
        std::cerr << "Error: Sequence with id " << id << " not found in input sequence collection\n";
        exit(EXIT_FAILURE);
    }
    assert(iter != m_map.end());
    return iter->second.sequence;
}

// 
void MapSequenceCollection::setPlaced(const std::string& id)
{
    SMPMap::iterator iter = m_map.find(id);
    assert(iter != m_map.end());
    iter->second.isPlaced = true;
}

// write the unplaced sequences of length at least minLength using pWriter
void MapSequenceCollection::writeUnplaced(std::ostream* pWriter, int minLength)
{
    int numUnplaced = 0;
    for(SMPMap::iterator iter = m_map.begin(); iter != m_map.end(); ++iter)
    {
        if(!iter->second.isPlaced && (int)iter->second.sequence.size() >= minLength)
        {
            std::stringstream idss;
            idss << "unplaced-" << numUnplaced++;
            writeFastaRecord(pWriter, idss.str(), iter->second.sequence);
        }
    }
}

