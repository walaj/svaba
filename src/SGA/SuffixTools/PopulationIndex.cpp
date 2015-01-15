//-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// PopulationIndex - A structure mapping a read
// index to a member of a population
//
#include "PopulationIndex.h"
#include "Util.h"
#include <algorithm>

PopulationIndex::PopulationIndex(const std::string& filename)
{
    std::istream* reader = createReader(filename);
    std::string line;
    while(getline(*reader, line))
    {
        PopulationMember member = str2member(line);
        m_population.push_back(member);
    }
    delete reader;
    reader = NULL;

    // Make sure the index is properly formatted
    assert(m_population.size() != 0);
    assert(m_population[0].start == 0);
    for(size_t i = 1; i < m_population.size(); ++i)
    {
        assert(m_population[i].start > m_population[i-1].start);
        assert(m_population[i].start == m_population[i-1].end + 1);
    }
    
    /* Test the index
    size_t start = 0;
    size_t end = m_population.back().end + 2;
    for(size_t i = start; i <= end; ++i)
    {
        std::string name = getName(i);
        printf("%zu = %s\n", i, name.c_str());
    }
    */
}

//
std::string PopulationIndex::getName(size_t read_index) const
{
    std::vector<PopulationMember>::const_iterator iter = getIterByReadIndex(read_index);
    return iter->name;
}

//
size_t PopulationIndex::getSampleIndex(size_t read_index) const
{
    std::vector<PopulationMember>::const_iterator iter = getIterByReadIndex(read_index);
    return iter - m_population.begin();
}

//
std::vector<PopulationMember>::const_iterator PopulationIndex::getIterByReadIndex(size_t read_index) const
{
    assert(read_index <= m_population.back().end);
    // Perform a binary search to identify the first element of the population
    // that is strictly greater than the index. This returns an iterator that is 
    // guarenteed to be one past the individual that index belongs to.
    PopulationMember search = { read_index, read_index, "" };
    std::vector<PopulationMember>::const_iterator iter = std::upper_bound(m_population.begin(), 
                                                                          m_population.end(), 
                                                                          search, 
                                                                          PopulationMember::compareByStart);
    assert(iter != m_population.begin());
    iter--;
    assert(read_index >= iter->start && read_index <= iter->end);
    return iter;
}

//
StringVector PopulationIndex::getSamples() const
{
    StringVector out;
    for(size_t i = 0; i < m_population.size(); ++i)
        out.push_back(m_population[i].name);
    return out;
}

//
void PopulationIndex::mergeIndexFiles(const std::string& file1, const std::string& file2, const std::string& outfile)
{
    std::ostream* writer = createWriter(outfile);
    
    // Copy the first index to the output unmodified but track the number of elements read
    size_t num_file_1 = 0;
    std::istream* reader = createReader(file1);
    std::string line;
    while(getline(*reader, line))
    {
        // Copy
        *writer << line << "\n";
        
        // Parse
        PopulationMember member = str2member(line);
        num_file_1 += (member.end - member.start + 1);
    }    
    delete reader;

    // Copy the second index, offsetting by the number of reads in file1
    reader = createReader(file2);
    while(getline(*reader, line))
    {
        PopulationMember member = str2member(line);
        member.start += num_file_1;
        member.end += num_file_1;

        // Copy
        *writer << member.start << "\t" << member.end << "\t" << member.name << "\n";
    } 
    delete reader;
    delete writer;
}

//
PopulationMember PopulationIndex::str2member(const std::string& line)
{
    std::stringstream parser(line);
    PopulationMember member;
    parser >> member.start;
    parser >> member.end;
    parser >> member.name;
    return member;
}
