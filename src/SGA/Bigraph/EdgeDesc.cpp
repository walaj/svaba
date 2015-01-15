//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// EdgeDesc - A unique description of an edge 
//
#include "EdgeDesc.h"
#include "Vertex.h"

// Operators
bool EdgeDesc::operator<(const EdgeDesc& obj) const
{
    assert(pVertex != NULL && obj.pVertex != NULL);
    if(pVertex->getID() < obj.pVertex->getID())
        return true;
    else if(pVertex->getID() > obj.pVertex->getID())
        return false;
    else if(dir < obj.dir)
        return true;
    else if(dir > obj.dir)
        return false;
    else if(comp < obj.comp)
        return true;
    else if(comp > obj.comp)
        return false;
    return false;
}

bool EdgeDesc::operator==(const EdgeDesc& obj) const
{
    assert(pVertex != NULL && obj.pVertex != NULL);
    return pVertex->getID() == obj.pVertex->getID() && dir == obj.dir && comp == obj.comp;
}

std::ostream& operator<<(std::ostream& out, const EdgeDesc& ed)
{
    assert(ed.pVertex != NULL);
    out << ed.pVertex->getID() << "," << ed.dir << "," << ed.comp;
    return out;
}

