#include "Contig.h"

//
// Statics
//
static const char* IDFIELD = "ID";
static const char* LENFIELD = "LN";
static const char* COVFIELD = "CV";
static const char* SEQFIELD = "SQ";
static const char* UNQFIELD = "UN";
static const char* ANTFIELD = "AN";

Contig::Contig() :    m_length(0), 
                    m_coverage(0.0f), 
                    m_uniqueFlag(UF_UNKNOWN) {}

//
//
//
std::istream& readFasta(std::istream& in, Contig& c)
{
    in.ignore(1); // skip ">"
    in >> c.m_id >> c.m_length >> c.m_coverage; // read header
    in.ignore(1); // skip newline
    in >> c.m_seq; // read seq
    in.ignore(1); // place the ifstream at the next line
    return in;
}

std::istream& readCAF(std::istream& in, Contig& c)
{
    std::string line;
    getline(in, line);
    std::stringstream tokenizer(line);

    std::string token;
    while(tokenizer >> token)
    {
        std::string key;
        std::string value;
        splitKeyValue(token, key, value);
        c.setFromKeyValue(key, value);
    }
    assert(c.getSequence().size() == c.getLength());
    return in;
}

//
//
//
SequenceVector Contig::getVariants() const
{
    SequenceVector out;
    out.push_back(getSequence());
    out.insert(out.end(), m_variants.begin(), m_variants.end());
    return out;
}

//
//
//
void Contig::addVariants(const SequenceVector& seqs)
{
    m_variants.insert(m_variants.end(), seqs.begin(), seqs.end());
}


//
//
//
std::ostream& writeCAF(std::ostream& out, Contig& c)
{
    std::vector<std::string> fields;
        
    // ID
    fields.push_back(makeKeyValue(IDFIELD, c.m_id));
    fields.push_back(makeKeyValue(LENFIELD, c.m_length));
    fields.push_back(makeKeyValue(COVFIELD, c.m_coverage));
    fields.push_back(makeKeyValue(SEQFIELD, c.m_seq));
    fields.push_back(makeKeyValue(UNQFIELD, c.m_uniqueFlag));
    std::copy(fields.begin(), fields.end(), std::ostream_iterator<std::string>(out, "\t"));
    return out;
}

//
// Set a field based on a key/value pair
//
void Contig::setFromKeyValue(std::string& key, std::string& value)
{
    std::stringstream parser(value);
    if(key == IDFIELD)
        parser >> m_id;
    else if(key == LENFIELD)
        parser >> m_length;
    else if(key == COVFIELD)
        parser >> m_coverage;
    else if(key == SEQFIELD)
        parser >> m_seq;
    else if(key == ANTFIELD)
        parser >> m_annotation;
    else if(key == UNQFIELD)
    {
        int v;
        parser >> v;
        m_uniqueFlag = (UniqueFlag)v;
    }
    else
        assert(false && "Unknown token found while parsing contig");
}
