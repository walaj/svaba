///----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// DBGPathGuide - Determine whether edges in a de Bruijn
// graph are supported by a subset of sequence reads
//
#include "DBGPathGuide.h"
#include <stdio.h>

DBGPathGuide::DBGPathGuide(size_t k) : m_k(k), m_pmers_checked(0), m_pmers_passed(0)
{

}

//
void DBGPathGuide::addSequence(const std::string& str)
{
    size_t p = m_k + 1;
    if(str.size() < p)
        return;

    for(size_t i = 0; i < str.size() - p + 1; i++)
        m_pmer_set.insert(str.substr(i, p));
}

//
bool DBGPathGuide::hasPmer(const std::string& str)
{
    m_pmers_checked += 1;
    bool passed = m_pmer_set.find(str) != m_pmer_set.end();
    m_pmers_passed += passed;
    return passed;
}

//
void DBGPathGuide::printStats() const
{
    printf("DBGPathGuide has %zu pmers. %zu out of %zu have passed the check\n", m_pmer_set.size(), m_pmers_passed, m_pmers_checked);
}
