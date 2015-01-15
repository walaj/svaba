///-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// StdAlnTools - Collection of wrappers around the
// stdaln dynamic programming alignment functions
#include <assert.h>
#include "StdAlnTools.h"
#include "stdaln.h"
#include "Alphabet.h"

// Perform a global alignment between the given strings
int StdAlnTools::globalAlignment(const std::string& target, const std::string& query, bool bPrint)
{
    path_t* path;
    int path_len = 0;
    int score = 0;
    createGlobalAlignmentPath(target, query, &path, &path_len, &score);

    if(bPrint)
    {
        std::string paddedTarget, paddedQuery, paddedMatch;
        makePaddedStringsFromPath(target, query, path, path_len, paddedTarget, paddedQuery, paddedMatch);
        printPaddedStrings(paddedTarget, paddedQuery, paddedMatch);

        std::cout << "CIGAR: " << makeCigar(path, path_len) << "\n";
        std::cout << "Global alignment score: " << score << "\n";
    }

    free(path);

    return score;
}

// Perform a global alignment between the given strings and return the CIGAR string
std::string StdAlnTools::globalAlignmentCigar(const std::string& target, const std::string& query)
{
    std::string cigar;
    int score = 0;
    globalAlignment(target, query, cigar, score);
    return cigar;
}

// Perform a global alignment between the given strings and return the CIGAR string and score
void StdAlnTools::globalAlignment(const std::string& target, const std::string& query, std::string& cigar, int& score)
{
    path_t* path;
    int path_len = 0;
    createGlobalAlignmentPath(target, query, &path, &path_len, &score);
    cigar = makeCigar(path, path_len);
    free(path);
}


// Perform a local alignment
LocalAlignmentResult StdAlnTools::localAlignment(const std::string& target, const std::string& query)
{
    // Set up global alignment parameters and data structures
    GlobalAlnParams params;
    int max_path_length = target.size() + query.size();
    path_t *path = (path_t*)calloc(max_path_length, sizeof(path_t));
    int path_len;

    AlnParam par;
    int matrix[25];
    StdAlnTools::setAlnParam(par, matrix, params);

    // Make a packed version of the query and target
    uint8_t* pQueryT = createPacked(query);
    uint8_t* pTargetT = createPacked(target);
    
    LocalAlignmentResult result;
    result.score = aln_local_core(pTargetT, target.size(), pQueryT, query.size(), &par, path, &path_len, 1, 0);
    assert(path_len <= max_path_length);

    result.cigar = makeCigar(path, path_len);

    // Calculate the aligned coordinates
    // This returns inclusive coordinates of the substrings aligned
    path_t* p = path + path_len - 1;
    result.targetStartIndex = (p->i ? p->i : 1) - 1;
    result.targetEndIndex = path->i - 1;
    result.queryStartIndex = (p->j ? p->j : 1) - 1;
    result.queryEndIndex = path->j - 1;
    
    // Update cigar
    if(result.queryStartIndex > 0)
    {
        std::stringstream csa;
        csa << result.queryStartIndex << "S";
        result.cigar.insert(0, csa.str());
    }

    int query_left_overhang = query.length() - (result.queryEndIndex + 1);

    if(query_left_overhang > 0)
    {
        std::stringstream csa;
        csa << query_left_overhang << "S";
        result.cigar.append(csa.str());
    }

    // Clean up
    delete [] pQueryT;
    delete [] pTargetT;
    free(path);

    return result;
}

// Expand a CIGAR string into a character code for each symbol of the alignment
std::string StdAlnTools::expandCigar(const std::string& cigar)
{
    std::string expanded;
    std::stringstream parser(cigar);
    int length;
    char code;
    while(parser >> length)
    {
        bool success = parser >> code;
        expanded.append(length, code);
        assert(success);
        (void)success;
    }
    return expanded;
}

// Compact an expanded CIGAR string into a regular cigar string
std::string StdAlnTools::compactCigar(const std::string& expanded_cigar)
{
    if(expanded_cigar.empty())
        return "";

    std::stringstream parser(expanded_cigar);
    std::stringstream writer;

    char prev_code = '\0';
    char curr_code;
    int curr_length = 0;
    while(parser >> curr_code)
    {
        if(curr_code == prev_code)
        {
            curr_length += 1;
        }
        else
        {
            // write the previous cigar character
            if(prev_code != '\0')
                writer << curr_length << prev_code;
            prev_code = curr_code;
            curr_length = 1;
        }
    }

    // Write the last symbol
    writer << curr_length << prev_code;

    return writer.str();
}


// Remove padding characters from str
std::string StdAlnTools::unpad(const std::string& str)
{
    std::string out;
    for(size_t i = 0; i < str.size(); ++i)
    {
        if(str[i] != '-')
            out.push_back(str[i]);   
    }
    return out;
}

// Create a path array representing the global alignment between target and query
// Caller must free the path
void StdAlnTools::createGlobalAlignmentPath(const std::string& target, const std::string& query,
                                            path_t** path, int* path_len, int* score)
{
    // Set up global alignment parameters and data structures
    GlobalAlnParams params;
    int max_path_length = target.size() + query.size();
    *path = (path_t*)calloc(max_path_length, sizeof(path_t));

    AlnParam par;
    int matrix[25];
    StdAlnTools::setAlnParam(par, matrix, params);

    // Make a packed version of the query and target
    uint8_t* pQueryT = createPacked(query);
    uint8_t* pTargetT = createPacked(target);
    *score = aln_global_core(pTargetT, target.size(), pQueryT, query.size(), &par, *path, path_len);
    assert(*path_len <= max_path_length);

    // Clean up
    delete [] pQueryT;
    delete [] pTargetT;
}

// Convert a std::string into the stdaln required packed format.
// This function allocates memory which the caller must free
uint8_t* StdAlnTools::createPacked(const std::string& s, size_t start, size_t length)
{
    if(length == std::string::npos)
        length = s.size();
    assert(length <= s.size());

    uint8_t* pBuffer = new uint8_t[length];
    for(size_t j = 0; j < length; ++j)
        pBuffer[j] = DNA_ALPHABET::getBaseRank(s[start + j]);
    return pBuffer;
}


// Calculate the maximum target length for a query of length ql
size_t StdAlnTools::calculateMaxTargetLength(int ql, const GlobalAlnParams& params)
{
    size_t mt = ((ql + 1) / 2 * params.match + params.gap_ext) / params.gap_ext + ql;
    return mt;
}

// Fill in the stdaln AlnParam data, using GlobalAlnParams
void StdAlnTools::setAlnParam(AlnParam& par, int matrix[25], const GlobalAlnParams& params)
{
    par.matrix = matrix;

    // Set matrix
	for (size_t i = 0; i < 25; ++i) par.matrix[i] = -params.mismatch;
	for (size_t i = 0; i < 4; ++i) par.matrix[i*5+i] = params.match; 
	par.gap_open = params.gap_open;
    par.gap_ext = params.gap_ext;
    par.gap_end = params.gap_ext;
	par.row = 5; 
    par.band_width = params.bandwidth;
}

// Convert a path to a CIGAR string
std::string StdAlnTools::makeCigar(path_t* path, int path_len)
{
    int cigarLen = 0;
	uint32_t* pCigar = aln_path2cigar32(path, path_len, &cigarLen); 
    std::stringstream cigarSS;
    for (int j = 0; j != cigarLen; ++j)
    {
        cigarSS << (pCigar[j]>>4);
        cigarSS << "MID"[pCigar[j]&0xf];
    }
    free(pCigar);
    return cigarSS.str();
}

// Make a pair of padded alignment strings from the pair of sequences
void StdAlnTools::makePaddedStrings(const std::string& s1, const std::string& s2, 
                                    std::string& out1, std::string& out2)
{
    path_t* path;
    int path_len = 0;
    int score = 0;
    createGlobalAlignmentPath(s1, s2, &path, &path_len, &score);
    std::string outm;
    makePaddedStringsFromPath(s1, s2, path, path_len, out1, out2, outm);
    free(path);
}

// Convert a dynamic programming path to a set of padded strings
// Algorithm ported from stdaln.c:aln_stdaln_aux
void StdAlnTools::makePaddedStringsFromPath(const std::string& s1, const std::string& s2, path_t* path, int path_len, 
                                            std::string& out1, std::string& out2, std::string& outm)
{
    out1.resize(path_len, 'A');
    out2.resize(path_len, 'A');
    outm.resize(path_len, 'A');

    path_t* p = path + path_len - 1;

    for (int l = 0; p >= path; --p, ++l) 
    {
        assert(l < (int)out1.size());
        int idx_i = p->i - 1;
        int idx_j = p->j - 1;

        switch (p->ctype) 
        {
            // p->i and p->j are indices into the DP matrix, which is
            // the string coordinate - 1
            case FROM_M:
                out1[l] = s1[idx_i]; 
                out2[l] = s2[idx_j];
                outm[l] = (s1[idx_i] == s2[idx_j] /*&& s1[p->i] != ap->row*/)? '|' : ' ';
                break;
            case FROM_I: 
                out1[l] = '-'; 
                out2[l] = s2[idx_j]; 
                outm[l] = ' '; 
                break;
            case FROM_D: 
                out1[l] = s1[idx_i]; 
                out2[l] = '-'; 
                outm[l] = ' '; 
                break;
        }
    }    
}

// 
void StdAlnTools::printPaddedStrings(const std::string& s1, const std::string& s2, const std::string& m, int colSize)
{
    assert(s1.size() == s2.size() && s1.size() == m.size());
    size_t len = s1.size();
    for(size_t l = 0; l < len; l += colSize)
    {
        int diff = len - l;
        int stop = diff < colSize ? diff : colSize;
        printf("1\t%s\n", s1.substr(l, stop).c_str());   
        printf("2\t%s\n", s2.substr(l, stop).c_str());   
        printf("m\t%s\n\n", m.substr(l, stop).c_str());
    }
}
