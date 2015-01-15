//-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// PopulationIndex - A structure mapping a read
// index to a member of a population
//
#ifndef POPULATION_INDEX_H
#define POPULATION_INDEX_H

#include <vector>
#include <string>

typedef std::vector<std::string> StringVector;

struct PopulationMember
{
    size_t start;
    size_t end;
    std::string name;

    //
    static bool compareByStart(const PopulationMember& a,
                               const PopulationMember& b)
    {
        return a.start < b.start;
    }
};

class PopulationIndex
{
    public:
        PopulationIndex(const std::string& filename);

        // Return the name of the individual for the read at index
        std::string getName(size_t read_index) const;
        
        // Return the index of the sample containing the given read index
        size_t getSampleIndex(size_t read_index) const;
        
        // Get the names of all samples in the collection
        StringVector getSamples() const;
        
        // Returns the number of samples in the population
        size_t getNumSamples() const { return m_population.size(); }

        // Merge two population index files into a new one
        static void mergeIndexFiles(const std::string& file1, const std::string& file2, const std::string& outfile);

    private:
        
        // Returns an iterator pointing to the sample that the contains the given index
        std::vector<PopulationMember>::const_iterator getIterByReadIndex(size_t index) const;

        // Parse a string from a .popidx file into a population member
        static PopulationMember str2member(const std::string& line);
        
        // Data
        std::vector<PopulationMember> m_population;
};

#endif
