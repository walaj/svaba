//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Util - Common data structures and functions
//
#include <iostream>
#include <math.h>
#include <map>
#include "Util.h"

//
// Sequence operations
//

// Reverse complement a sequence
std::string reverseComplement(const std::string& seq)
{
    std::string out(seq.length(), 'A');
    size_t last_pos = seq.length() - 1;
    for(int i = last_pos; i >= 0; --i)
    {
        out[last_pos - i] = complement(seq[i]);
    }
    return out;
}

// Reverse complement a sequence using the full iupac alphabet
std::string reverseComplementIUPAC(const std::string& seq)
{
    std::string out(seq.length(), 'A');
    size_t last_pos = seq.length() - 1;
    for(int i = last_pos; i >= 0; --i)
    {
        out[last_pos - i] = complementIUPAC(seq[i]);
    }
    return out;
}

// Reverse a sequence
std::string reverse(const std::string& seq)
{
    return std::string(seq.rbegin(), seq.rend());
}

// Complement a sequence
std::string complement(const std::string& seq)
{
    std::string out(seq.length(), 'A');
    size_t l = seq.length();
    for(size_t i = 0; i < l; ++i)
        out[i] = complement(seq[i]);
    return out;
}

// Complement a sequence over the IUPAC alphabet
std::string complementIUPAC(const std::string& seq)
{
    std::string out(seq.length(), 'A');
    size_t l = seq.length();
    for(size_t i = 0; i < l; ++i)
        out[i] = complementIUPAC(seq[i]);
    return out;
}

// 
std::string prefix(const std::string& seq, const unsigned int len)
{
    assert(seq.length() >= len);
    return seq.substr(0, len);
}

// 
std::string suffix(const std::string& seq, const unsigned int len)
{
    assert(seq.length() >= len);
    return seq.substr(seq.length() - len);
}

//
// Dust scoring scheme as given by:
// Morgulis A. "A fast and symmetric DUST implementation to Mask
// Low-Complexity DNA Sequences". J Comp Bio.
double calculateDustScore(const std::string& seq)
{
    std::map<std::string, int> scoreMap;
    
    // Cannot calculate dust scores on very short reads
    if(seq.size() < 3)
        return 0.0f;

    // Slide a 3-mer window over the sequence and insert the sequences into the map
    for(size_t i = 0; i < seq.size() - 3; ++i)
    {
        std::string triMer = seq.substr(i, 3);
        scoreMap[triMer]++;
    }

    // Calculate the score by summing the square of every element in the map
    double sum = 0;
    std::map<std::string, int>::iterator iter = scoreMap.begin();
    for(; iter != scoreMap.end(); ++iter)
    {
        int tc = iter->second;
        double score = (double)(tc * (tc - 1)) / 2.0f;
        sum += score;
    }
    return sum / (seq.size() - 2);
}

// Returns the window over seq with the highest dust score
double maxDustWindow(const std::string& seq, size_t windowSize, size_t minWindow)
{
    double maxScore = 0.0f;
    for(size_t i = 0; i < seq.size(); i += 1)
    {
        size_t r = seq.size() - i;
        size_t w = r < windowSize ? r : windowSize;
        if(w >= minWindow) // don't calculate score for small windows
        {
            double s = calculateDustScore(seq.substr(i, w));
            if(s > maxScore)
                maxScore = s;
        }
    }
    return maxScore;
}

// count the differences between s1 and s2 over the first n bases
int countDifferences(const std::string& s1, const std::string& s2, size_t n)
{
    int numDiff = 0;
    for(size_t i = 0; i < n; ++i)
    {
        if(s1[i] != s2[i])
            numDiff++;
    }
    return numDiff;
}

// count the differences between s1 and s2 over the first n bases
std::string getDiffString(const std::string& s1, const std::string& s2)
{
    std::string out;
    size_t stop = std::min(s1.size(), s2.size());
    for(size_t i = 0; i < stop; ++i)
        out.push_back(s1[i] == s2[i] ? '.' : s1[i]);
    return out;
}

// 
char randomBase()
{
    int i = rand() % 4;
    switch(i)
    {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        default:
            assert(false);
    }
    return 'A';        
}

// Strip the leading directories and
// the last trailling suffix from a filename
std::string stripFilename(const std::string& filename)
{
    std::string out = stripDirectories(filename);
    // Remove the gzip extension if necessary
    if(isGzip(out))
        out = stripExtension(out);
    return stripExtension(out);
}

// Remove a single file extension from the filename
std::string stripExtension(const std::string& filename)
{
    size_t suffixPos = filename.find_last_of('.');
    if(suffixPos == std::string::npos)
        return filename; // no suffix
    else
        return filename.substr(0, suffixPos);
}

// Remove an extension from a filename, including the .gz extension
// file.fastq will return file
// file.fastq.gz will return file
std::string stripGzippedExtension(const std::string& filename)
{
    if(isGzip(filename))
        return stripExtension(stripExtension(filename));
    else
        return stripExtension(filename);
}


// Strip the leadering directories from a filename
std::string stripDirectories(const std::string& filename)
{
    size_t lastDirPos = filename.find_last_of('/');
    
    if(lastDirPos == std::string::npos)
        return filename; // no directories
    else
        return filename.substr(lastDirPos + 1);
}

// Return the file extension
std::string getFileExtension(const std::string& filename)
{
    size_t suffixPos = filename.find_last_of('.');
    
    if(suffixPos == std::string::npos)
        return "";
    else
        return filename.substr(suffixPos + 1);
}

// Write out a fasta record
void writeFastaRecord(std::ostream* pWriter, const std::string& id, const std::string& seq, size_t maxLineLength)
{
    *pWriter << ">" << id << " " << seq.length() << "\n";
    for(size_t i = 0; i < seq.size(); i += maxLineLength)
    {
        size_t span = std::min(seq.size() - i, maxLineLength);
        *pWriter << seq.substr(i, span) << "\n";
    }
}

// Returns true if the filename has an extension indicating it is compressed
bool isGzip(const std::string& filename)
{
    size_t suffix_length = sizeof(GZIP_EXT) - 1;

    // Assume files without an extension are not compressed
    if(filename.length() < suffix_length)
        return false;

    std::string extension = suffix(filename, suffix_length);
    return extension == GZIP_EXT;
}

// Returns true if the filename has an extension indicating it is fastq
bool isFastq(const std::string& filename)
{
    return filename.find(".fastq") != std::string::npos || filename.find(".fq") != std::string::npos;
}


// Returns the size of the file. Code from stackoverflow.
std::ifstream::pos_type getFilesize(const std::string& filename)
{
    std::ifstream in(filename.c_str(), std::ifstream::in | std::ifstream::binary);
    in.seekg(0, std::ifstream::end);
    return in.tellg();
}

// Open a file that may or may not be gzipped for reading
// The caller is responsible for freeing the handle
std::istream* createReader(const std::string& filename, std::ios_base::openmode mode)
{
    if(isGzip(filename))
    {
        igzstream* pGZ = new igzstream(filename.c_str(), mode);
        assertGZOpen(*pGZ, filename);
        return pGZ;
    }
    else
    {
        std::ifstream* pReader = new std::ifstream(filename.c_str(), mode);
        assertFileOpen(*pReader, filename);
        return pReader;
    }
}

// Open a file that may or may not be gzipped for writing
// The caller is responsible for freeing the handle
std::ostream* createWriter(const std::string& filename,
                           std::ios_base::openmode mode)
{
    if(isGzip(filename))
    {
        ogzstream* pGZ = new ogzstream(filename.c_str(), mode);
        assertGZOpen(*pGZ, filename);
        return pGZ;
    }
    else
    {
        std::ofstream* pWriter = new std::ofstream(filename.c_str(), mode);
        assertFileOpen(*pWriter, filename);
        return pWriter;
    }
}

// Ensure a filehandle is open
void assertFileOpen(std::ifstream& fh, const std::string& fn)
{
    if(!fh.is_open())
    {
        std::cerr << "Error: could not open " << fn << " for read\n";
        exit(EXIT_FAILURE);
    }    
}

// Ensure a filehandle is open
void assertFileOpen(std::ofstream& fh, const std::string& fn)
{
    if(!fh.is_open())
    {
        std::cerr << "Error: could not open " << fn << " for write\n";
        exit(EXIT_FAILURE);
    }    
}

//
void assertGZOpen(gzstreambase& gh, const std::string& fn)
{
    if(!gh.good())
    {
        std::cerr << "Error: could not open " << fn << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Split a string into parts based on the delimiter
StringVector split(std::string in, char delimiter)
{
    StringVector out;
    size_t lastPos = 0;
    size_t pos = in.find_first_of(delimiter);

    while(pos != std::string::npos)
    {
        out.push_back(in.substr(lastPos, pos - lastPos));
        lastPos = pos + 1;
        pos = in.find_first_of(delimiter, lastPos);
    }
    out.push_back(in.substr(lastPos));
    return out;
}

// Split a key-value pair
void splitKeyValue(std::string in, std::string& key, std::string& value)
{
    StringVector parts = split(in, CAF_SEP);
    if(parts.size() != 2)
    {
        std::cerr << "Invalid key-value pair " << in << std::endl;
        assert(false);
    }

    key = parts[0];
    value = parts[1];

    assert(key.size() > 0 && value.size() > 0 && "Invalid key-value pair");
}

// Get the shared component of a read pair name
// This is the part of the name preceding the "/"
std::string getPairBasename(const std::string& id)
{
    assert(!id.empty());

    size_t pos = id.find_last_of('/');
    if(pos == std::string::npos)
        return id;
    else
        return id.substr(0, pos);
}

// Get the ID of the pair of a given read
std::string getPairID(const std::string& id)
{
    assert(!id.empty());
    std::string pid(id);

    size_t li = id.length() - 1;
    char last = id[li];

    if(last == 'A')
        pid[li] = 'B';
    else if(last == 'B')
        pid[li] = 'A';
    else if(last == '1')
        pid[li] = '2';
    else if(last == '2')
        pid[li] = '1';
    else if(last == 'f')
        pid[li] = 'r';
    else if(last == 'r')
        pid[li] = 'f';
    else
        pid = "";
    return pid;
}

// Return 0 if the id indicates the first read of a pair, 1 otherwise
int getPairIndex(const std::string& id)
{
    size_t li = id.length() - 1;
    char last = id[li];
    if(last == 'A' || last == '1' || last == 'f')
        return 0;
    else if(last == 'B' || last == '2' || last == 'r')
        return 1;
    else
    {
        std::cerr << "Unrecognized pair format: " << id << "\n";
        exit(EXIT_FAILURE);
    }
}

// Debug function to parse the distance between two reads
// based on their names
// This assumes is that the read names are just the positions the reads
// were sampled from
size_t debug_getReadDistFromNames(const std::string& name1, const std::string& name2)
{
    std::string id;
    std::string pos;
    StringVector parts = split(name1, ':');
    if(parts.size() != 2)
        return 0;

    int p1 = atoi(parts[1].c_str());
    
    parts = split(name2, ':');
    if(parts.size() != 2)
        return 0;

    int p2 = atoi(parts[1].c_str());
    int dist = int(abs(p1 - p2));
    return dist;
}

