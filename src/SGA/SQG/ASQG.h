//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ASQG - Definitions and functions
// for handling ASGQ files
//
#ifndef ASQG_H
#define ASQG_H

#include "SQG.h"
#include "Match.h"

namespace ASQG
{
    enum RecordType
    {
        RT_HEADER = 0,
        RT_VERTEX,
        RT_EDGE
    };

    // A header record is just a tag:value pairs
    struct HeaderRecord
    {
        public:
            HeaderRecord();
            HeaderRecord(const std::string& recordLine);
            
            void setOverlapTag(int overlapLen);
            void setInputFileTag(const std::string& name);
            void setErrorRateTag(float errorRate);
            void setContainmentTag(int v);
            void setTransitiveTag(int v);

            const SQG::IntTag& getVersionTag() const { return m_versionTag; }
            const SQG::FloatTag& getErrorRateTag() const { return m_errorRateTag; }
            const SQG::StringTag& getInfileTag() const { return m_infileTag; }
            const SQG::IntTag& getOverlapTag() const { return m_overlapTag; }
            const SQG::IntTag& getContainmentTag() const { return m_containmentTag; };
            const SQG::IntTag& getTransitiveTag() const { return m_transitiveTag; };

            void write(std::ostream& out);
            void write(std::stringstream& out); // JEREMIAH
            void parse(const std::string& record);

        private:

            void setVersionTag(int version);
            
            SQG::IntTag m_versionTag;
            SQG::FloatTag m_errorRateTag;
            SQG::StringTag m_infileTag;
            SQG::IntTag m_overlapTag;
            SQG::IntTag m_containmentTag;
            SQG::IntTag m_transitiveTag;
    };

    // A vertex record is an id, sequence and an array of
    // tag:value 
    struct VertexRecord
    {
        public:
            VertexRecord() {}
            VertexRecord(const std::string& recordLine);
            VertexRecord(const std::string& i, const std::string& s) : m_id(i), m_seq(s) {}

            void setSubstringTag(bool b);
            
            const std::string& getID() const { return m_id; }
            const std::string& getSeq() const { return m_seq; }
            const SQG::IntTag& getSubstringTag() const { return m_substringTag; }

            void write(std::ostream& out);
            void write(std::stringstream& out);
            void parse(const std::string& record);

        private:

            std::string m_id;
            std::string m_seq;
            SQG::IntTag m_substringTag;
    };

    // An edge record is just an overlap object and tag:values
    struct EdgeRecord
    {
        public:
            EdgeRecord() {}
            EdgeRecord(const std::string& recordLine);
            EdgeRecord(const Overlap& o) : m_overlap(o) {}

            const Overlap& getOverlap() const { return m_overlap; }

            void write(std::ostream& out);
            void parse(const std::string& record);

        private:
        
            Overlap m_overlap;
    };

    // Parsing functions
    RecordType getRecordType(const std::string& record);
    HeaderRecord parseHeaderRecord(const std::string& record);
    VertexRecord parseVertexRecord(const std::string& record);
    EdgeRecord parseEdgeRecord(const std::string& record);

    // Writing functions
    void writeFields(std::ostream& out, const StringVector& fields);
};

#endif
