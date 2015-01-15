//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldGraph - A graph representing long-distance
// relationships between contigs.
//
#include "ScaffoldGraph.h"
#include "SeqReader.h"

ScaffoldGraph::ScaffoldGraph()
{
    m_vertices.set_deleted_key("");
}

//
ScaffoldGraph::~ScaffoldGraph()
{
    for(ScaffoldVertexMap::iterator iter = m_vertices.begin();
         iter != m_vertices.end(); ++iter)
    {
        iter->second->deleteEdges();
        delete iter->second;
        iter->second = NULL;
    }
}

//
void ScaffoldGraph::addVertex(ScaffoldVertex* pVertex)
{
    m_vertices.insert(std::make_pair(pVertex->getID(), pVertex));
}

//
void ScaffoldGraph::addEdge(ScaffoldVertex* pVertex, ScaffoldEdge* pEdge)
{
    assert(pVertex != NULL);
    pVertex->addEdge(pEdge);
}

ScaffoldVertexPtrVector ScaffoldGraph::getAllVertices() const
{
    ScaffoldVertexPtrVector outVertices;
    ScaffoldVertexMap::const_iterator iter = m_vertices.begin();
    for(; iter != m_vertices.end(); ++iter)
        outVertices.push_back(iter->second);
    return outVertices;
}

//
ScaffoldVertex* ScaffoldGraph::getVertex(VertexID id) const
{
    ScaffoldVertexMap::const_iterator iter = m_vertices.find(id);
    if(iter == m_vertices.end())
        return NULL;
    return iter->second;
}

//
void ScaffoldGraph::deleteVertices(ScaffoldVertexClassification classification)
{
    ScaffoldVertexMap::iterator iter = m_vertices.begin(); 
    while(iter != m_vertices.end())
    {
        if(iter->second->getClassification() == classification)
        {
            iter->second->deleteEdgesAndTwins();
            delete iter->second;
            iter->second = NULL;
            m_vertices.erase(iter++);
        }
        else
        {
            ++iter;
        }
    }
}

//
void ScaffoldGraph::deleteEdgesByColor(GraphColor c)
{
    ScaffoldVertexMap::iterator iter = m_vertices.begin(); 
    while(iter != m_vertices.end())
    {
        iter->second->deleteEdgesAndTwinsByColor(c);
        ++iter;
    }
}



//
void ScaffoldGraph::setVertexColors(GraphColor c)
{
    ScaffoldVertexMap::iterator iter = m_vertices.begin(); 
    while(iter != m_vertices.end())
    {
        iter->second->setColor(c);
        ++iter;
    }
}

//
void ScaffoldGraph::setEdgeColors(GraphColor c)
{
    ScaffoldVertexMap::iterator iter = m_vertices.begin(); 
    while(iter != m_vertices.end())
    {
        iter->second->setEdgeColors(c);
        ++iter;
    }
}

// 
void ScaffoldGraph::loadVertices(const std::string& filename, int minLength)
{
    SeqReader reader(filename, SRF_NO_VALIDATION | SRF_KEEP_CASE);
    SeqRecord sr;
    while(reader.get(sr))
    {
        int contigLength = sr.seq.length();
        if(contigLength >= minLength)
        {
            ScaffoldVertex* pVertex = new ScaffoldVertex(sr.id, sr.seq.length());
            addVertex(pVertex);
        }
    }    
}

// 
void ScaffoldGraph::loadDistanceEstimateEdges(const std::string& filename, bool isMatePair, int verbose)
{
    std::cout << "Reading distance estimates from " << filename << "\n";
    std::istream* pReader = createReader(filename);
    std::string line;

    while(getline(*pReader, line))
    {
        assert(line.substr(0,4) != "Mate");
        StringVector fields = split(line, ' ');
        assert(fields.size() >= 1);

        std::string rootID = fields[0];
        EdgeDir currDir = ED_SENSE; // abyss distance estimate outputs the sense contigs first

        for(size_t i = 1; i < fields.size(); ++i)
        {
            std::string record = fields[i];
            if(record == ";")
            {
                currDir = !currDir;
                continue;
            }

            std::string id;
            EdgeComp comp;
            int distance;
            int numPairs;
            double stdDev;
            parseDERecord(record, id, comp, distance, numPairs, stdDev);

            // Get the vertices that are linked
            ScaffoldVertex* pVertex1 = getVertex(rootID);
            ScaffoldVertex* pVertex2 = getVertex(id);

            if(pVertex1 != NULL && pVertex2 != NULL)
            {
                if(pVertex1 == pVertex2)
                {
                    std::cout << "Self-edges not allowed\n";
                    assert(false);
                    continue;
                }

                ScaffoldLink link1(id, currDir, comp, distance, stdDev, numPairs, pVertex2->getSeqLen(), SLT_DISTANCEEST);
                ScaffoldLink link2(rootID, !correctDir(currDir, comp), comp, distance, stdDev, numPairs, pVertex1->getSeqLen(), SLT_DISTANCEEST);

                // Check if there already exists a DistanceEstimate edge between these vertices
                ScaffoldEdge* pEdge = pVertex1->findEdgeTo(id, SLT_DISTANCEEST);
                if(pEdge != NULL)
                {
                    // An edge to this vertex already exists
                    // If the current estimate is mate-pair link, never overwrite
                    // the current estimate
                    if(!isMatePair && pEdge->getLink().stdDev < stdDev)
                    {
                        pEdge->setLink(link1);
                        pEdge->getTwin()->setLink(link2);
                    }
                    else
                    {
                        if(abs(pEdge->getDistance() - link1.distance) > 100)
                        {
                            if(verbose >= 1)
                            {
                                printf("LL skipped from %s to %s. Distance1: %d Distance2: %d\n", pVertex1->getID().c_str(), 
                                                                                                  link1.endpointID.c_str(), 
                                                                                                  link1.distance, 
                                                                                                  pEdge->getDistance());
                            }
                            pVertex1->setConflictingFlag(true);
                            pVertex2->setConflictingFlag(true);
                        }
                    }
                }
                else
                {
                    ScaffoldEdge* pEdge1 = new ScaffoldEdge(pVertex2, link1);
                    ScaffoldEdge* pEdge2 = new ScaffoldEdge(pVertex1, link2);

                    pEdge1->setTwin(pEdge2);
                    pEdge2->setTwin(pEdge1);

                    addEdge(pVertex1, pEdge1);
                    addEdge(pVertex2, pEdge2);
                }
            }
        }
    }

    delete pReader;
}

void ScaffoldGraph::loadAStatistic(const std::string& filename)
{
    std::istream* pReader = createReader(filename);
    std::string line;

    while(getline(*pReader, line))
    {
        StringVector fields = split(line, '\t');
        assert(fields.size() == 6);

        VertexID id = fields[0];

        std::stringstream cn_parser(fields[4]);
        double cn;
        cn_parser >> cn;

        std::stringstream as_parser(fields[5]);
        double as;
        as_parser >> as;

        ScaffoldVertex* pVertex = getVertex(id);
        if(pVertex != NULL)
        {
            pVertex->setAStatistic(as);
            pVertex->setEstCopyNumber(cn);
        }
    }
    delete pReader;
}


//
void ScaffoldGraph::parseDERecord(const std::string& record, std::string& id, 
                                  EdgeComp& comp, int& distance, int& numPairs, double& stdDev)
{
    StringVector fields = split(record, ',');
    if(fields.size() != 4)
    {
        std::cerr << "Distance Estimate record is not formatted correctly: " << record << "\n";
        exit(1);
    }

    // Parse the ID and its orientation
    id = fields[0].substr(0, fields[0].size() - 1);
    comp = (fields[0][fields[0].size() - 1] == '+' ? EC_SAME : EC_REVERSE);

    std::stringstream d_parser(fields[1]);
    d_parser >> distance;
    
    std::stringstream np_parser(fields[2]);
    np_parser >> numPairs;

    std::stringstream sd_parser(fields[3]);
    sd_parser >> stdDev;
}

//
void ScaffoldGraph::writeDot(const std::string& outFile) const
{
    std::ostream* pWriter = createWriter(outFile);
    
    std::string graphType = "digraph";

    *pWriter << graphType << " G\n{\n";
    ScaffoldVertexMap::const_iterator iter = m_vertices.begin(); 
    for(; iter != m_vertices.end(); ++iter)
    {
        iter->second->writeDot(pWriter);
    }
    *pWriter << "}\n";
    delete pWriter;
}
