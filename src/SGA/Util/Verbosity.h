///----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Verbosity - Singleton object to only
// hold the verbosity level for the entire
// program
//
#ifndef VERBOSITY_H
#define VERBOSITY_H

class Verbosity
{
    public:
        Verbosity() : m_maxPrintLevel(0) {}

        // Get singleton object
        static Verbosity& Instance()
        {
            static Verbosity instance;
            return instance;
        }
     
        // Get/set the print level
        int getPrintLevel() const { return m_maxPrintLevel; }
        void setPrintLevel(int l) { m_maxPrintLevel = l; }

    private:
        int m_maxPrintLevel;
};


#endif
