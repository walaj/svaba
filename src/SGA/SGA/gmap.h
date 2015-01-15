//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// gmap - Map sequences to the vertices of a graph
//
#ifndef GMAP_H
#define GMAP_H
#include <getopt.h>
#include "config.h"

// struct
struct GmapRecord
{
    std::string readID;
    std::string readSeq;
    std::string mappedID;
    bool isRC;

    bool isMapped() const
    {
        if(mappedID != "-" && mappedID != "MM")
            return true;
        else
            return false;
    }

    friend std::ostream& operator<<(std::ostream& out, const GmapRecord& record)
    {
        out << record.readID << "\t" << record.readSeq << "\t" << record.mappedID << "\t" << record.isRC;
        return out;
    }

    friend std::istream& operator>>(std::istream& in, GmapRecord& record)
    {
        in >> record.readID >> record.readSeq >> record.mappedID >> record.isRC;
        return in;
    }
};

// functions
int gmapMain(int argc, char** argv);
void parseGmapOptions(int argc, char** argv);
void gmap();
void parseGmapHits(const StringVector& hitsFilenames);

#endif
