///----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// DBGPathGuide - Determine whether edges in a de Bruijn
// graph are supported by a subset of sequence reads
//
#ifndef DBGPATHGUIDE_H
#define DBGPATHGUIDE_H

#include <string>
#include <set>

class DBGPathGuide
{
    public:
        DBGPathGuide(size_t k);

        // Add the (k+1)-mers in seq to the set
        void addSequence(const std::string& seq);
        
        // Check whether the given p-mer is in the sest
        bool hasPmer(const std::string& str);

        // Write some information to stdout
        void printStats() const;

    private:

        size_t m_k;
        size_t m_pmers_checked;
        size_t m_pmers_passed;
        std::set<std::string> m_pmer_set;
};

#endif
