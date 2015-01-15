//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// QuickBWT - Construct a simple BWT for a short input string
//
#include "QuickBWT.h"

// Construct a BWT for a small string.
// Calling function is responsible for freeing the memory
void createQuickBWT(const std::string& str, BWT*& pBWT, SuffixArray*& pSA)
{
    ReadTable rt;
    SeqItem si = { "a", str };
    rt.addRead(si);  

    pSA = new SuffixArray(&rt, 1, true);
    pBWT = new BWT(pSA, &rt);

    /*
    std::cout << "QBWT -- str len: " << str.length() << "\n";
    std::cout << "QBWT -- bwt len: " << pBWT->getBWLen() << "\n";
    pBWT->printInfo();
    */
}
