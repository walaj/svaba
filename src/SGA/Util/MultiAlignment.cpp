//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// MultiAlignment - Class constructing a multiple
// alignment between a root sequence and a set of sequences
//
#include <iostream>
#include <algorithm>
#include <map>
#include <limits>
#include "MultiAlignment.h"
#include "StdAlnTools.h"

struct CigarIter
{
    const MAlignData* pData;
    int symIdx; // the index of the cigar string we are parsing
    int baseIdx; // the index of the current base 

    char updateAndEmit(char updateMode)
    {
        char cigSym = getCigarSymbol();
        char outBase;

        if(updateMode == 'I')
        {
            if(cigSym == 'I')
            {
                symIdx += 1;
                outBase = '-';
            }
            else if(cigSym == 'D')
            {
                assert(false);
            }
            else
            {
                symIdx += 1;
                outBase = getOutputBase(cigSym);
                if(outBase != '.' && outBase != '-')
                    baseIdx += 1;
            }
        }
        else if(updateMode == 'D')
        {
            if(cigSym == 'D')
            {
                symIdx += 1;
                outBase = getOutputBase(cigSym);
                if(outBase != '.')
                    baseIdx += 1;
            }
            else if(cigSym == 'I')
            {
                outBase = '-';
            }
            else if(cigSym == 'S')
            {
                outBase = '.';
            }
            else
            {
                outBase = '-';
                //symIdx += 1;
            }
        }
        else
        {
            symIdx += 1;
            outBase = getOutputBase(cigSym);
            if(outBase != '.')
                baseIdx += 1;
        }

        return outBase;
    }

    char getCigarSymbol()
    {
        if(symIdx >= (int)pData->expandedCigar.size())
            return 'S';
        return pData->expandedCigar[symIdx];
    }

    char getOutputBase(char cigSym)
    {
        if(cigSym == 'S' || cigSym == 'N')
            return '.';
        else
            return pData->str[baseIdx];
    }

    static bool sortPosition(const CigarIter& a, const CigarIter& b)
    {
        return a.pData->position < b.pData->position;
    }
};
typedef std::vector<CigarIter> CigarIterVector;

MultiAlignment::MultiAlignment(std::string rootStr, const MAlignDataVector& inData, std::string rootName)
{
    m_verbose = 1;
    // Build a padded multiple alignment from the pairwise alignments to the root
    CigarIterVector iterVec;

    // Create entry for the root
    MAlignData rootData;
    rootData.str = rootStr;
    rootData.expandedCigar = std::string(rootStr.size(), 'M');
    rootData.position = 0;
    rootData.name = rootName;

    m_alignData.push_back(rootData);
    m_alignData.insert(m_alignData.end(), inData.begin(), inData.end());

    CigarIter tmp = {&rootData, 0, 0};
    iterVec.push_back(tmp);

    if(m_verbose > 1)
        printf("%d\t%s\n", 0, tmp.pData->expandedCigar.c_str());

    for(size_t i = 0; i < inData.size(); ++i)
    {
        CigarIter iter = {&inData[i], 0, 0};
        iterVec.push_back(iter);
        if(m_verbose > 1)
            printf("%zu\t%s\n", i+1, iter.pData->expandedCigar.c_str());
    }

    // Build the padded strings by inserting '-' as necessary
    bool done = false;
    while(!done)
    {
        // Check if any strings have a deletion or insertion at this base
        bool hasDel = false;
        bool hasInsert = false;
        bool hasMatch = false;

        for(size_t i = 0; i < iterVec.size(); ++i)
        {
            char sym = iterVec[i].getCigarSymbol();
            if(sym == 'D')
                hasDel = true;
            else if(sym == 'I')
                hasInsert = true;
            else if(sym == 'M')
                hasMatch = true;
        }

        done = !hasDel && !hasInsert && !hasMatch;
        if(done)
            break;

        char updateMode = hasDel ? 'D' : (hasInsert ? 'I' : 'M');
        for(size_t i = 0; i < iterVec.size(); ++i)
        {
            char outSym = iterVec[i].updateAndEmit(updateMode);
            m_alignData[i].padded.append(1,outSym);
        }
    }
}

// 
size_t MultiAlignment::getIdxByName(const std::string& name) const
{
    size_t max = std::numeric_limits<size_t>::max();
    size_t idx = max;
    for(size_t i = 0; i < m_alignData.size(); ++i)
    {
        if(m_alignData[i].name == name)
        {
            if(idx == max)
            {
                idx = i;
            }
            else
            {
                std::cerr << "Error in MultiAlignment::getIdxByName: duplicate rows found for " << name << "\n";
                exit(EXIT_FAILURE);
            }
        }
    }

    if(idx == max)
    {
        std::cerr << "Error in MultiAlignment::getIdxByName: row not found for " << name << "\n";
        exit(EXIT_FAILURE);
    }
    return idx;
}

//
size_t MultiAlignment::getNumColumns() const
{
    assert(!m_alignData.empty());
    return m_alignData.front().padded.size();
}

std::string MultiAlignment::getCigar(size_t rowIdx) const
{
    assert(rowIdx < m_alignData.size());
    std::cout << "MA has " << m_alignData.size() << " elements\n";
    return StdAlnTools::compactCigar(m_alignData[rowIdx].expandedCigar);
}

//
char MultiAlignment::getSymbol(size_t rowIdx, size_t colIdx) const
{
    assert(rowIdx < m_alignData.size());
    assert(colIdx < m_alignData[rowIdx].padded.size());
    return m_alignData[rowIdx].padded[colIdx];
}

size_t MultiAlignment::getBaseIdx(size_t rowIdx, size_t colIdx) const
{
    assert(rowIdx < m_alignData.size());
    assert(colIdx < m_alignData[rowIdx].padded.size());

    // Convert the column index to the baseIndex in the sequence at

    // Ensure this is a real base on the target string.
    assert(getSymbol(rowIdx, colIdx) != '-');

    // Substract the number of padding charcters from the column index to get the 
    // base index
    size_t padSyms = 0;
    for(size_t i = 0; i < colIdx; ++i)
    {
        if(getSymbol(rowIdx, i) == '-')
            padSyms += 1;
    }
    assert(padSyms <= colIdx);
    return colIdx - padSyms;
}


//
std::string MultiAlignment::getPaddedSubstr(size_t rowIdx, size_t start, size_t length) const
{
    assert(rowIdx < m_alignData.size());
    assert(start < m_alignData[rowIdx].padded.size());
    assert(start + length <= m_alignData[rowIdx].padded.size());
    return m_alignData[rowIdx].padded.substr(start, length);
}

// Count the length of the homopolymer starting at from
size_t MultiAlignment::countHomopolymer(size_t rowIdx, int from, int to) const
{
    assert(rowIdx < m_alignData.size());
    assert(from < (int)m_alignData[rowIdx].padded.size());
    
    //
    if(to == -1)
        to = m_alignData[rowIdx].padded.size() - 1; // inclusive

    // Determine iteration direction
    int step = from <= to ? 1 : -1;
    
    // The first or last column of the multiple alignment was requested
    // This can only be a homopolymer of length 1
    if(from == (int)m_alignData[rowIdx].padded.size() - 1 || (from == 0 && step < 0))
        return 1;

    // Get the base of the homopolymer
    // If it is a padding symbol we use the next non-padded base in the sequence
    char b = '\0';
    int max_position = (int)m_alignData[rowIdx].padded.size() - 1;
    while(from >= 0 && from <= max_position)
    {
        b = getSymbol(rowIdx, from);
        if(b == '-')
            from += step;
        else
            break;
    }

    size_t length = 1;
    do
    {
        from += step;
        if(from < 0 || from > max_position)
            return length;

        char s = getSymbol(rowIdx, from);
        if(s == '-')
            continue;
        if(s == b)
            length += 1;
        else
            break;
    }
    while(from != to);
    return length;
}

// Generate simple consensus string from the multiple alignment
std::string MultiAlignment::generateConsensus()
{
    assert(!m_alignData.empty() && !m_alignData.front().padded.empty());

    std::string consensus;
    std::string paddedConsensus;

    std::map<char, int> baseMap;

    size_t num_rows = m_alignData.size() - 1; // do not include root
    size_t num_cols = m_alignData.front().padded.size();
    for(size_t i = 0; i < num_cols; ++i)
    {
        baseMap.clear();
        for(size_t j = 1; j < num_rows; ++j)
        {
            char b = m_alignData[j].padded[i];
            if(b != '.')
                baseMap[b]++;
        }

        int max = 0;
        char maxBase = '.';
        for(std::map<char,int>::iterator iter = baseMap.begin(); iter != baseMap.end(); ++iter)
        {
            if(iter->second > max)
            {
                max = iter->second;
                maxBase = iter->first;
            }
        }

        paddedConsensus.append(1,maxBase);

        if(maxBase == '.' && consensus.empty())
            continue; // skip no-call at beginning
        else if(maxBase == '.')
            break; // non-called position in middle of read, stop consensus generation
        else if(maxBase != '-') //
            consensus.append(1, maxBase);
    }

    if(m_verbose > 0)
        print(80, &paddedConsensus);

    return consensus;
}

// Generate a string representing the columns that differ between the strings
std::string MultiAlignment::generateMatchString() const
{
    assert(!m_alignData.empty() && !m_alignData.front().padded.empty());

    std::string matchString;

    size_t num_rows = m_alignData.size();
    size_t num_cols = m_alignData.front().padded.size();
    for(size_t i = 0; i < num_cols; ++i)
    {
        char rootChar = m_alignData[0].padded[i];
        bool allMatch = true;
        for(size_t j = 1; j < num_rows; ++j)
        {
            char b = m_alignData[j].padded[i];
            if(b != rootChar)
                allMatch = false;
        }
        matchString.append(1, allMatch ? '*' : ' ');
    }

    return matchString;
}

//
void MultiAlignment::filterByEditDistance(int max_edit_distance)
{
    size_t num_rows = m_alignData.size();
    assert(num_rows > 0);
    size_t num_cols = m_alignData.front().padded.size();
    std::vector<bool> keep_vector(num_rows, true);
    for(size_t i = 1; i < num_rows; ++i)
    {
        int current_edit_distance = 0;
        for(size_t j = 0; j < num_cols; ++j)
        {
            char root_symbol = getSymbol(0, j);
            char seq_symbol = getSymbol(i, j);
            if(root_symbol != '.' && seq_symbol != '.' && seq_symbol != '-' && root_symbol != seq_symbol)
                current_edit_distance++;
        }

        if(current_edit_distance > max_edit_distance)
            keep_vector[i] = false;
    }

    MAlignDataVector tmp_ma;
    for(size_t i = 0; i < num_rows; ++i)
    {
        if(keep_vector[i])
            tmp_ma.push_back(m_alignData[i]);
    }

    printf("Removed %zu\n", m_alignData.size() - tmp_ma.size());
    m_alignData.swap(tmp_ma);
}

//
void MultiAlignment::print(int col_size, const std::string* pConsensus, bool sorted, bool masked) const
{
    assert(!m_alignData.empty() && !m_alignData.front().padded.empty());

    std::string matchString = generateMatchString();

    // Create a copy of the m_alignData and sort it by position
    MAlignDataVector sortedAlignments = m_alignData;
    if(sorted)
        std::stable_sort(sortedAlignments.begin(), sortedAlignments.end(), MAlignData::sortPosition);
    col_size = 100000;
    size_t len = sortedAlignments[0].padded.size();
    for(size_t l = 0; l < len; l += col_size)
    {
        // Print the consensus if requested
        if(pConsensus != NULL)
        {
            int diff = pConsensus->size() - l;
            int stop = diff < col_size ? diff : col_size;
            if(stop > 0)
                printf("C\t%s\n", pConsensus->substr(l,stop).c_str());
            else
                printf("C\n");
        }
        
        // Print each row
        for(size_t i = 0; i < sortedAlignments.size(); ++i)
        {
            const MAlignData& mad = sortedAlignments[i];
            int diff = mad.padded.size() - l;
            int stop = diff < col_size ? diff : col_size;

            std::string s = mad.padded.substr(l, stop).c_str();
            if(masked)
            {
                for(size_t j = 0; j < s.size(); ++j)
                {
                    // get base row symbol
                    char rb = getSymbol(0, l + j);
                    if(s[j] == rb && s[j] != '-')
                        s[j] = '=';
                }
            }
            
            printf("%zu\t%s\t%s\n", i, s.c_str(), mad.name.c_str());
        }
    
        // Print the matched columns
        int diff = matchString.size() - l;
        int stop = diff < col_size ? diff : col_size;
        
        printf("M\t%s\n", matchString.substr(l, stop).c_str());
        std::cout << "\n";
    }
}

// Construct a multiple alignment of the input strings from global alignments
// to the first element in the vector
MultiAlignment MultiAlignmentTools::alignSequencesGlobal(const SeqItemVector& sequences)
{
    assert(!sequences.empty());
    const std::string& base = sequences[0].seq.toString();

    MAlignDataVector madVector;
    for(size_t i = 1; i < sequences.size(); ++i)
    {
        std::string seq = sequences[i].seq.toString();
        std::string cigar = StdAlnTools::globalAlignmentCigar(seq, base);
    
        // Initialize the multiple alignment data
        MAlignData maData;
        maData.position = 0;
        maData.str = seq;
        maData.expandedCigar = StdAlnTools::expandCigar(cigar);

        maData.name = sequences[i].id;
        madVector.push_back(maData);
    }

    return MultiAlignment(base, madVector, sequences[0].id);
}

// Construct a multiple alignment of the input strings from local alignments
// to the first element of the vector
MultiAlignment MultiAlignmentTools::alignSequencesLocal(const SeqItemVector& sequences)
{
    assert(!sequences.empty());
    const std::string& base = sequences[0].seq.toString();

    MAlignDataVector madVector;
    for(size_t i = 1; i < sequences.size(); ++i)
    {
        // The alignment is (slightly counterintuitively) with respect to the 
        // base string as the query so that all the cigar strings used
        // are based on edit operations to the base.
        std::string seq = sequences[i].seq.toString();
        LocalAlignmentResult result = StdAlnTools::localAlignment(seq,base);
        
        // Initialize the multiple alignment data
        MAlignData maData;
        maData.position = result.queryStartIndex;
        // If the non-base sequence overhangs the base, clip it
        if(result.targetStartIndex > 0)
            maData.str = seq.substr(result.targetStartIndex);
        else
            maData.str = seq;

        // Pad out the cigar with reference skip characters
        maData.expandedCigar = StdAlnTools::expandCigar(result.cigar);

        maData.name = sequences[i].id;
        madVector.push_back(maData);
    }

    return MultiAlignment(base, madVector, sequences[0].id);
}
