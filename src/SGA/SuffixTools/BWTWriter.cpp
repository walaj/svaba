//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTWriter - Abstract class for writing a BWT file to disk
//
#include "BWTWriter.h"
#include "BWTWriterBinary.h"
#include "BWTWriterAscii.h"

//
void IBWTWriter::write(const SuffixArray* pSA, const ReadTable* pRT)
{
    size_t num_symbols = pSA->getSize();
    size_t num_strings = pSA->getNumStrings();
    writeHeader(num_strings, num_symbols, BWF_NOFMI);

    for(size_t i = 0; i < num_symbols; ++i)
    {
        SAElem saElem = pSA->get(i);
        const SeqItem& si = pRT->getRead(saElem.getID());

        // Get the position of the start of the suffix
        uint64_t f_pos = saElem.getPos();
        uint64_t l_pos = (f_pos == 0) ? si.seq.length() : f_pos - 1;
        char b = (l_pos == si.seq.length()) ? '$' : si.seq.get(l_pos);
        writeBWChar(b);
    }
    finalize();
}

//
IBWTWriter* BWTWriter::createWriter(const std::string& filename)
{
    return new BWTWriterBinary(filename);
}
