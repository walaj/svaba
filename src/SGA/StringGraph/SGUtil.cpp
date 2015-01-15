//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGUtils - Data structures/Functions related
// to building and manipulating string graphs
//
#include "SGUtil.h"
#include "SeqReader.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"

StringGraph* SGUtil::loadASQG(const std::string& filename, const unsigned int minOverlap, 
                              bool allowContainments, size_t maxEdges)
{
    // Initialize graph
    StringGraph* pGraph = new StringGraph;

    std::istream* pReader = createReader(filename);

    int stage = 0;
    int line = 0;
    std::string recordLine;
    while(getline(*pReader, recordLine))
    {
        ASQG::RecordType rt = ASQG::getRecordType(recordLine);
        switch(rt)
        {
            case ASQG::RT_HEADER:
            {
                if(stage != 0)
                {
                    std::cerr << "Error: Unexpected header record found at line " << line << "\n";
                    exit(EXIT_FAILURE);
                }

                ASQG::HeaderRecord headerRecord(recordLine);
                const SQG::IntTag& overlapTag = headerRecord.getOverlapTag();
                if(overlapTag.isInitialized())
                    pGraph->setMinOverlap(overlapTag.get());
                else
                    pGraph->setMinOverlap(0);

                const SQG::FloatTag& errorRateTag = headerRecord.getErrorRateTag();
                if(errorRateTag.isInitialized())
                    pGraph->setErrorRate(errorRateTag.get());
                
                const SQG::IntTag& containmentTag = headerRecord.getContainmentTag();
                if(containmentTag.isInitialized())
                    pGraph->setContainmentFlag(containmentTag.get());
                else
                    pGraph->setContainmentFlag(true); // conservatively assume containments are present

                const SQG::IntTag& transitiveTag = headerRecord.getTransitiveTag();
                if(!transitiveTag.isInitialized())
                {
                    std::cerr << "Warning: ASQG does not have transitive tag\n";
                    pGraph->setTransitiveFlag(true);
                }
                else
                {
                    pGraph->setTransitiveFlag(transitiveTag.get());
                }

                break;
            }
            case ASQG::RT_VERTEX:
            {
                // progress the stage if we are done the header
                if(stage == 0)
                    stage = 1;

                if(stage != 1)
                {
                    std::cerr << "Error: Unexpected vertex record found at line " << line << "\n";
                    exit(EXIT_FAILURE);
                }

                ASQG::VertexRecord vertexRecord(recordLine);
                const SQG::IntTag& ssTag = vertexRecord.getSubstringTag();

		//                Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(vertexRecord.getID(), vertexRecord.getSeq());
		Vertex* pVertex = new Vertex(vertexRecord.getID(), vertexRecord.getSeq());
                if(ssTag.isInitialized() && ssTag.get() == 1)
                {
                    // Vertex is a substring of some other vertex, mark it as contained
                    pVertex->setContained(true);
                    pGraph->setContainmentFlag(true);
                }
                pGraph->addVertex(pVertex);
                break;
            }
            case ASQG::RT_EDGE:
            {
                if(stage == 1)
                    stage = 2;
                
                if(stage != 2)
                {
                    std::cerr << "Error: Unexpected edge record found at line " << line << "\n";
                    exit(EXIT_FAILURE);
                }

                ASQG::EdgeRecord edgeRecord(recordLine);
                const Overlap& ovr = edgeRecord.getOverlap();

                // Add the edge to the graph
                if(ovr.match.getMinOverlapLength() >= (int)minOverlap)
                    SGAlgorithms::createEdgesFromOverlap(pGraph, ovr, allowContainments, maxEdges);
                break;
            }
        }
        ++line;
    }

    // Completely delete the edges for all nodes that were marked as super-repetitive in the graph
    SGSuperRepeatVisitor superRepeatVisitor;
    pGraph->visit(superRepeatVisitor);

    // Remove any duplicate edges
    SGDuplicateVisitor dupVisit;
    pGraph->visit(dupVisit);

    SGGraphStatsVisitor statsVisit;
    pGraph->visit(statsVisit);
    // Remove identical vertices
    // This is much cheaper to do than remove via
    // SGContainRemove as no remodelling needs to occur
   /*
    SGIdenticalRemoveVisitor irv;
    pGraph->visit(irv);

    // Remove substring vertices
    while(pGraph->hasContainment())
    {
        SGContainRemoveVisitor crv;
        pGraph->visit(crv);
    }
*/
    delete pReader;
    return pGraph;
}

//JEREMIAH
StringGraph* SGUtil::loadASQG(std::stringstream& pReader, const unsigned int minOverlap, 
                              bool allowContainments, size_t maxEdges)
{
    // Initialize graph
    StringGraph* pGraph = new StringGraph;

    //std::istream* pReader = createReader(filename);

    int stage = 0;
    int line = 0;
    std::string recordLine;
    while(getline(pReader, recordLine))
    {
        ASQG::RecordType rt = ASQG::getRecordType(recordLine);
        switch(rt)
        {
            case ASQG::RT_HEADER:
            {
                if(stage != 0)
                {
                    std::cerr << "Error: Unexpected header record found at line " << line << "\n";
                    exit(EXIT_FAILURE);
                }

                ASQG::HeaderRecord headerRecord(recordLine);
                const SQG::IntTag& overlapTag = headerRecord.getOverlapTag();
                if(overlapTag.isInitialized())
                    pGraph->setMinOverlap(overlapTag.get());
                else
                    pGraph->setMinOverlap(0);

                const SQG::FloatTag& errorRateTag = headerRecord.getErrorRateTag();
                if(errorRateTag.isInitialized())
                    pGraph->setErrorRate(errorRateTag.get());
                
                const SQG::IntTag& containmentTag = headerRecord.getContainmentTag();
                if(containmentTag.isInitialized())
                    pGraph->setContainmentFlag(containmentTag.get());
                else
                    pGraph->setContainmentFlag(true); // conservatively assume containments are present

                const SQG::IntTag& transitiveTag = headerRecord.getTransitiveTag();
                if(!transitiveTag.isInitialized())
                {
                    std::cerr << "Warning: ASQG does not have transitive tag\n";
                    pGraph->setTransitiveFlag(true);
                }
                else
                {
                    pGraph->setTransitiveFlag(transitiveTag.get());
                }

                break;
            }
            case ASQG::RT_VERTEX:
            {
                // progress the stage if we are done the header
                if(stage == 0)
                    stage = 1;

                if(stage != 1)
                {
                    std::cerr << "Error: Unexpected vertex record found at line " << line << "\n";
                    exit(EXIT_FAILURE);
                }

                ASQG::VertexRecord vertexRecord(recordLine);
                const SQG::IntTag& ssTag = vertexRecord.getSubstringTag();

		//                Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(vertexRecord.getID(), vertexRecord.getSeq());
                Vertex* pVertex = new Vertex(vertexRecord.getID(), vertexRecord.getSeq());
                if(ssTag.isInitialized() && ssTag.get() == 1)
                {
                    // Vertex is a substring of some other vertex, mark it as contained
                    pVertex->setContained(true);
                    pGraph->setContainmentFlag(true);
                }
                pGraph->addVertex(pVertex);
                break;
            }
            case ASQG::RT_EDGE:
            {
                if(stage == 1)
                    stage = 2;
                
                if(stage != 2)
                {
                    std::cerr << "Error: Unexpected edge record found at line " << line << "\n";
                    exit(EXIT_FAILURE);
                }

                ASQG::EdgeRecord edgeRecord(recordLine);
                const Overlap& ovr = edgeRecord.getOverlap();

                // Add the edge to the graph
                if(ovr.match.getMinOverlapLength() >= (int)minOverlap)
                    SGAlgorithms::createEdgesFromOverlap(pGraph, ovr, allowContainments, maxEdges);
                break;
            }
        }
        ++line;
    }

    // Completely delete the edges for all nodes that were marked as super-repetitive in the graph
    SGSuperRepeatVisitor superRepeatVisitor;
    pGraph->visit(superRepeatVisitor);

    // Remove any duplicate edges
    SGDuplicateVisitor dupVisit;
    pGraph->visit(dupVisit);

    //SGGraphStatsVisitor statsVisit;
    //pGraph->visit(statsVisit);
    // Remove identical vertices
    // This is much cheaper to do than remove via
    // SGContainRemove as no remodelling needs to occur
   /*
    SGIdenticalRemoveVisitor irv;
    pGraph->visit(irv);

    // Remove substring vertices
    while(pGraph->hasContainment())
    {
        SGContainRemoveVisitor crv;
        pGraph->visit(crv);
    }
*/
    //delete pReader;
    return pGraph;
}

// Load a graph (with no edges) from a fasta file
StringGraph* SGUtil::loadFASTA(const std::string& filename)
{
    StringGraph* pGraph = new StringGraph;
    SeqReader reader(filename);
    SeqRecord record;

    while(reader.get(record))
    {
      //Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(record.id, record.seq.toString());
        Vertex* pVertex = new Vertex(record.id, record.seq.toString());
        pGraph->addVertex(pVertex);
    }
    return pGraph;
}
