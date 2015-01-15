//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// PrimerScreen - Singleton class to filter sequences 
// that match a database of primer sequences
//
#include "PrimerScreen.h"

// Hardcoded primer sequences that we wish to filter against

// Primers used in the Sanger's pcr-free library prep
// See Kozarewa et al. (http://www.nature.com/nmeth/journal/v6/n4/full/nmeth.1311.html)
// supplemental info. We do not use the whole primer sequences.
#define ILLUMINA_SANGER_PCR_FREE_A "AATGATACGGCGACCACCGAGATCTACA"
#define ILLUMINA_SANGER_PCR_FREE_B "GATCGGAAGAGCGGTTCAGCAGGAATGC"

PrimerScreen::PrimerScreen()
{
    m_db.push_back(ILLUMINA_SANGER_PCR_FREE_A);
    m_db.push_back(ILLUMINA_SANGER_PCR_FREE_B);
}

// Check seq against the primer database
bool PrimerScreen::containsPrimer(const std::string& seq)
{
    static PrimerScreen screener; // initializes singleton object if necessary
    
    // For now we only check if the first 14
    // bases of seq is a perfect match to any sequence in the db
    // This is sufficient to get rid of the vast majority of the primer
    // contamination
    const size_t check_size = 14;
    std::string check = seq.substr(0, check_size);
    for(size_t i = 0; i < screener.m_db.size(); ++i)
    {
        if(screener.m_db[i].find(check) != std::string::npos)
            return true;
    }
    return false;
}
