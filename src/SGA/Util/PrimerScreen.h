//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// PrimerScreen - Singleton class to filter sequences 
// that match a database of primer sequences
//
#ifndef PRIMERSCREEN_H
#define PRIMERSCREEN_H

#include "Util.h"

class PrimerScreen
{

    public:
        
        // Return true if the sequence fails the primer check
        static bool containsPrimer(const std::string& seq);

    private:
        
        PrimerScreen();

        StringVector m_db;
};

#endif
