//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Util - Common data structures and functions
//
#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <list>
#include <string>
#include <istream>
#include <fstream>
#include <cassert>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <stdint.h>
#include <iostream>
#include "DNAString.h"
#include "gzstream.h"
#include "Quality.h"

#define CAF_SEP ':'
#define FUNCTION_TIMER Timer functionTimer(__PRETTY_FUNCTION__);
#define VALIDATION_WARNING(x) static bool validation_warn = true; if(validation_warn) \
                               printf("[%s] Warning validation is on\n", (x)); validation_warn = false;

#define WARN_ONCE(x) static bool _warn_once = true; if(_warn_once) \
                     printf("WARNING: [%s]\n", (x)); _warn_once = false;

#define GZIP_EXT ".gz"

//
// Typedef
//
typedef std::string Sequence;
typedef std::string ContigID;
typedef std::vector<int> IntVec;
typedef std::vector<double> DoubleVec;
typedef std::vector<std::string> StringVector;
typedef std::list<std::string> StringList;
typedef std::vector<Sequence> SequenceVector;

// SeqItem is just an id, sequence pair
struct SeqItem
{
    // data
    std::string id;
    DNAString seq;

    void write(std::ostream& out, const std::string& meta = "") const
    {
        out << ">" << id << (meta.empty() ? "" : " ") << meta << "\n";
        out << seq.toString() << "\n";
    }    
};
typedef std::vector<SeqItem> SeqItemVector;

// SeqRecord is a id,sequence pair with an associated quality value
struct SeqRecord
{
    SeqItem toSeqItem() const
    {
        SeqItem r;
        r.id = id;
        r.seq = seq;
        return r;
    }

    // Write the sequence to the file handle. If the 
    void write(std::ostream& out, const std::string& meta = "") const
    {
        // If there is a quality string write the record as fastq, otherwise fasta
        if(!qual.empty())
        {
            out << "@" << id << (meta.empty() ? "" : " ") << meta << "\n";
            out << seq.toString() << "\n";
            out << "+" << "\n";
            out << qual << "\n";
        }
        else
        {
            toSeqItem().write(out, meta);
        }
    }

    bool hasQuality() const
    {
        return !qual.empty();
    }

    // Get the phred score for base i. If there are no quality string, returns DEFAULT_QUAL_SCORE
    int getPhredScore(size_t i) const
    {
        if(!qual.empty())
        {
            assert(i < qual.size());
            return Quality::char2phred(qual[i]);
        }
        else
        {
            return DEFAULT_QUAL_SCORE;
        }
    }

    // data
    std::string id;
    DNAString seq;
    std::string qual;
};
typedef std::vector<SeqRecord> SeqRecordVector;

//
// Functions
//
std::string stripFilename(const std::string& filename);
std::string stripExtension(const std::string& filename);
std::string stripGzippedExtension(const std::string& filename);
std::string stripDirectories(const std::string& filename);
std::string getFileExtension(const std::string& filename);
bool isGzip(const std::string& filename);
bool isFastq(const std::string& filename);
std::ifstream::pos_type getFilesize(const std::string& filename);

// Write out a fasta record
void writeFastaRecord(std::ostream* pWriter, const std::string& id, const std::string& seq, size_t maxLength = 80);

// Wrapper function for opening a reader of compressed or uncompressed file
std::istream* createReader(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in);
std::ostream* createWriter(const std::string& filename, 
                           std::ios_base::openmode mode = std::ios_base::out);

void assertFileOpen(std::ifstream& fh, const std::string& fn);
void assertFileOpen(std::ofstream& fh, const std::string& fn);
void assertGZOpen(gzstreambase& gh, const std::string& fn);

char randomBase();

// Key-value operations
template <class C>
std::string makeKeyValue(std::string key, C value)
{
    std::stringstream ss;
    ss << key << CAF_SEP << value;
    return ss.str();
}

StringVector split(std::string in, char delimiter);
void splitKeyValue(std::string in, std::string& key, std::string& value);

std::string getPairBasename(const std::string& id);
std::string getPairID(const std::string& id);

// Returns 0 if the id indicates the first read in a pair, 1 otherwise
int getPairIndex(const std::string& id);

// Debug function to get the distance between two reads based on their names, which 
// encodes the positions
size_t debug_getReadDistFromNames(const std::string& name1, const std::string& name2);

//
// Return the lexographic value for the given base
//
static const uint8_t s_lexoRankLUT[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,1,0,2,0,0,0,3,0,0,0,0,0,0,0,0,
    0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

inline static uint8_t getBaseRank(char b)
{
    return s_lexoRankLUT[static_cast<uint8_t>(b)];
}

//
// Sequence operations
//
std::string reverseComplement(const std::string& seq);
std::string complement(const std::string& seq);
std::string reverse(const std::string& seq);

// Reverse/complement functions, allowing a full IUPAC alphabet
std::string reverseComplementIUPAC(const std::string& seq);
std::string complementIUPAC(const std::string& seq);

// Return the prefix/suffix of seq of length len 
std::string prefix(const std::string& seq, const unsigned int len);
std::string suffix(const std::string& seq, const unsigned int len);

// Calculate the dust score for the given sequence
double maxDustWindow(const std::string& seq, size_t windowSize = 64, size_t minWindow = 64);
double calculateDustScore(const std::string& seq);

// Count the number of differences between s1 and s2 over the first n chars
int countDifferences(const std::string& s1, const std::string& s2, size_t n);

// Construct a string encoding the differences between s1 and s2
std::string getDiffString(const std::string& s1, const std::string& s2);

// Complement a base
inline char complement(char base)
{
    switch(base)
    {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        case 'N':
            return 'N';
        default:
            assert(false && "Unknown base!");
            return 'N';
    }
}

// Complement a base using the full IUPAC alphabet
// Also allows for lowercase bases
inline char complementIUPAC(char c)
{
    char cmp = '\0';
    bool is_lc = std::islower(c);

    switch(std::toupper(c)) {
        case 'A': cmp = 'T'; break;
        case 'C': cmp = 'G'; break;
        case 'G': cmp = 'C'; break;
        case 'T': cmp = 'A'; break;
        case 'M': cmp = 'K'; break;
        case 'R': cmp = 'Y'; break;
        case 'W': cmp = 'W'; break;
        case 'S': cmp = 'S'; break;
        case 'Y': cmp = 'R'; break;
        case 'K': cmp = 'M'; break;
        case 'V': cmp = 'B'; break;
        case 'H': cmp = 'D'; break;
        case 'D': cmp = 'H'; break;
        case 'B': cmp = 'V'; break;
        case 'N': cmp = 'N'; break;
        default:
            assert(false);
    }

    if(is_lc)
        cmp = std::tolower(cmp);
    return cmp;
}

// Wrapper function to determine whether a calculated error rate is within tolerance of a
// threshold. A fixed tolerance is used as this does not need to be super precise, it only
// needs to account for small differences in fp calculations
inline bool isErrorRateAcceptable(double er, double threshold)
{
    static double tolerance = 0.000001f;
    if(er - tolerance <= threshold)
        return true;
    else
        return false;
}


#endif
