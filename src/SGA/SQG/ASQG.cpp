//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ASQG - Definitions and functions
// for handling ASGQ files
//
#include "ASQG.h"
#include <iterator>

namespace ASQG
{

const int HEADER_VERSION = 1;

// Record ID tags
static int RECORD_TAG_SIZE = 2; // do not include null terminator 
static char HEADER_TAG[] = "HT";
static char VERTEX_TAG[] = "VT";
static char EDGE_TAG[] = "ED";

// Header tags
const int FIELD_TAG_SIZE = 2;
static char VERSION_TAG[] = "VN";
static char OVERLAP_TAG[] = "OL";
static char INPUT_FILE_TAG[] = "IN";
static char ERROR_RATE_TAG[] = "ER";
static char CONTAINMENT_TAG[] = "CN"; // 1 if the graph has containment edges/vertices
static char TRANSITIVE_TAG[] = "TE"; // 1 if the graph has transitive edges

// Vertex tags
static char SUBSTRING_TAG[] = "SS";

//
// Header Record
//
HeaderRecord::HeaderRecord()
{
    setVersionTag(HEADER_VERSION);
}

//
HeaderRecord::HeaderRecord(const std::string& recordLine)
{
    parse(recordLine);
}

//
void HeaderRecord::setVersionTag(int version)
{
    m_versionTag.set(version);
}

//
void HeaderRecord::setOverlapTag(int overlapLen)
{
    m_overlapTag.set(overlapLen);
}

//
void HeaderRecord::setInputFileTag(const std::string& name)
{
    m_infileTag.set(name);
}

//
void HeaderRecord::setErrorRateTag(float errorRate)
{
    m_errorRateTag.set(errorRate);
}

//
void HeaderRecord::setContainmentTag(int v)
{
    m_containmentTag.set(v);
}

//
void HeaderRecord::setTransitiveTag(int v)
{
    m_transitiveTag.set(v);
}

//
void HeaderRecord::write(std::ostream& out)
{
    StringVector fields;
    fields.push_back(HEADER_TAG);

    // Check for mandatory tags
    if(!m_versionTag.isInitialized())
    {
        std::cerr << "Error: Header version tag not set, aborting." << std::endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        fields.push_back(m_versionTag.toTagString(VERSION_TAG));
    }

    if(m_errorRateTag.isInitialized())
        fields.push_back(m_errorRateTag.toTagString(ERROR_RATE_TAG));

    if(m_overlapTag.isInitialized())
        fields.push_back(m_overlapTag.toTagString(OVERLAP_TAG));
    
    if(m_infileTag.isInitialized())
        fields.push_back(m_infileTag.toTagString(INPUT_FILE_TAG));
    
    if(m_containmentTag.isInitialized())
        fields.push_back(m_containmentTag.toTagString(CONTAINMENT_TAG));
    
    if(m_transitiveTag.isInitialized())
        fields.push_back(m_transitiveTag.toTagString(TRANSITIVE_TAG));

    writeFields(out, fields);
}

//JEREMIAH
void HeaderRecord::write(std::stringstream& out)
{
    StringVector fields;
    fields.push_back(HEADER_TAG);

    // Check for mandatory tags
    if(!m_versionTag.isInitialized())
    {
        std::cerr << "Error: Header version tag not set, aborting." << std::endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        fields.push_back(m_versionTag.toTagString(VERSION_TAG));
    }

    if(m_errorRateTag.isInitialized())
        fields.push_back(m_errorRateTag.toTagString(ERROR_RATE_TAG));

    if(m_overlapTag.isInitialized())
        fields.push_back(m_overlapTag.toTagString(OVERLAP_TAG));
    
    if(m_infileTag.isInitialized())
        fields.push_back(m_infileTag.toTagString(INPUT_FILE_TAG));
    
    if(m_containmentTag.isInitialized())
        fields.push_back(m_containmentTag.toTagString(CONTAINMENT_TAG));
    
    if(m_transitiveTag.isInitialized())
        fields.push_back(m_transitiveTag.toTagString(TRANSITIVE_TAG));

    writeFields(out, fields);
}


//
void HeaderRecord::parse(const std::string& record)
{
    if(record.size() < 2)
    {
        std::cerr << "Error: Record is not valid\n";
        exit(EXIT_FAILURE);
    }

    // Tokenize record
    StringVector tokens = SQG::tokenizeRecord(record);

    // Ensure the first token indicates this is a valid header record
    if(tokens[0].compare(0, RECORD_TAG_SIZE, HEADER_TAG) != 0)
    {
        std::cerr << "Error: Record does not have a header tag" << std::endl;
        std::cerr << "Record: " << record << std::endl;
        exit(EXIT_FAILURE);
    }

    for(size_t i = 1; i < tokens.size(); ++i)
    {
        static char VERSION_TAG[] = "VN";
        static char OVERLAP_TAG[] = "OL";
        static char INPUT_FILE_TAG[] = "IN";
        static char ERROR_RATE_TAG[] = "ER";
        
        if(tokens[i].compare(0, FIELD_TAG_SIZE, VERSION_TAG) == 0)
            m_versionTag.fromString(tokens[i]);

        if(tokens[i].compare(0, FIELD_TAG_SIZE, OVERLAP_TAG) == 0)
            m_overlapTag.fromString(tokens[i]);

        if(tokens[i].compare(0, FIELD_TAG_SIZE, INPUT_FILE_TAG) == 0)
            m_infileTag.fromString(tokens[i]);

        if(tokens[i].compare(0, FIELD_TAG_SIZE, ERROR_RATE_TAG) == 0)
            m_errorRateTag.fromString(tokens[i]);

        if(tokens[i].compare(0, FIELD_TAG_SIZE, CONTAINMENT_TAG) == 0)
            m_containmentTag.fromString(tokens[i]);
        
        if(tokens[i].compare(0, FIELD_TAG_SIZE, TRANSITIVE_TAG) == 0)
            m_transitiveTag.fromString(tokens[i]);

    }
}

//
// Vertex Record
//
//
VertexRecord::VertexRecord(const std::string& recordLine)
{
    parse(recordLine);
}

//
void VertexRecord::setSubstringTag(bool b)
{
    m_substringTag.set(b);
}

//
void VertexRecord::write(std::ostream& out)
{
    StringVector fields;
    fields.push_back(VERTEX_TAG);
    fields.push_back(m_id);
    fields.push_back(m_seq);

    if(m_substringTag.isInitialized())
        fields.push_back(m_substringTag.toTagString(SUBSTRING_TAG));

    writeFields(out, fields);
}

//JEREMIAH
void VertexRecord::write(std::stringstream& out)
{
    StringVector fields;
    fields.push_back(VERTEX_TAG);
    fields.push_back(m_id);
    fields.push_back(m_seq);

    if(m_substringTag.isInitialized())
        fields.push_back(m_substringTag.toTagString(SUBSTRING_TAG));

    writeFields(out, fields);
}



//
void VertexRecord::parse(const std::string& record)
{
    if(record.size() < 2)
    {
        std::cerr << "Error: Record is not valid\n";
        exit(EXIT_FAILURE);
    }

    // Tokenize record
    StringVector tokens = SQG::tokenizeRecord(record);

    if(tokens.size() < 3)
    {
        std::cerr << "Error: Vertex record is incomplete.\n";
        std::cerr << "Record: " << record << std::endl;
        exit(EXIT_FAILURE);
    }

    // Ensure the first token indicates this is a valid vertex record
    if(tokens[0].compare(0, RECORD_TAG_SIZE, VERTEX_TAG) != 0)
    {
        std::cerr << "Error: Record does not have a vertex tag" << std::endl;
        std::cerr << "Record: " << record << std::endl;
        exit(EXIT_FAILURE);
    }

    m_id = tokens[1];
    m_seq = tokens[2];
    
    for(size_t i = 3; i < tokens.size(); ++i)
    {
        if(tokens[i].compare(0, FIELD_TAG_SIZE, SUBSTRING_TAG) == 0)
            m_substringTag.fromString(tokens[i]);
    }    
}



//
// EdgeRecord
//
EdgeRecord::EdgeRecord(const std::string& recordLine)
{
    parse(recordLine);
}

//
void EdgeRecord::write(std::ostream& out)
{
    StringVector fields;
    fields.push_back(EDGE_TAG);
    std::stringstream ss;
    ss << m_overlap;
    fields.push_back(ss.str());
    writeFields(out, fields);
}

//
void EdgeRecord::parse(const std::string& record)
{
    if(record.size() < 2)
    {
        std::cerr << "Error: Record is not valid\n";
        exit(EXIT_FAILURE);
    }

    // Tokenize record
    StringVector tokens = SQG::tokenizeRecord(record);

    if(tokens.size() < 2)
    {
        std::cerr << "Error: Edge record is incomplete.\n";
        std::cerr << "Record: " << record << std::endl;
        exit(EXIT_FAILURE);
    }

    // Ensure the first token indicates this is a valid edge record
    if(tokens[0].compare(EDGE_TAG) != 0)
    {
        std::cerr << "Error: Record does not have an edge tag" << std::endl;
        std::cerr << "Record: " << record << std::endl;
        exit(EXIT_FAILURE);
    }

    // Parse the overlap
    std::stringstream ssparser(tokens[1]);
    ssparser >> m_overlap;
}

//
//
//
RecordType getRecordType(const std::string& record)
{
    // the record type is the first two characters of the record
    if(record.size() < 2)
    {
        std::cerr << "Error: record does not have a valid record-type tag" << std::endl;
        std::cerr << "Record: " << record << std::endl;
        exit(EXIT_FAILURE);
    }

    char recordTag[RECORD_TAG_SIZE];
    record.copy(recordTag, RECORD_TAG_SIZE);

    if(strncmp(recordTag, HEADER_TAG, RECORD_TAG_SIZE) == 0)
    {
        return RT_HEADER;
    }

    if(strncmp(recordTag, VERTEX_TAG, RECORD_TAG_SIZE) == 0)
    {
        return RT_VERTEX;
    }

    if(strncmp(recordTag, EDGE_TAG, RECORD_TAG_SIZE) == 0)
    {
        return RT_EDGE;
    }

    // Unknown tag
    std::cerr << "Error: Unknown ASQG file record tag: " << recordTag << std::endl;
    exit(EXIT_FAILURE);
}

void writeFields(std::ostream& out, const StringVector& fields)
{
    for(size_t i = 0; i < fields.size(); ++i)
    {
        out << fields[i];
        if(i != fields.size() - 1)
            out << SQG::FIELD_SEP;
    }
    out << "\n";
}

}

