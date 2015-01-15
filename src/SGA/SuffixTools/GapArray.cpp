//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// GapArray - Data structure and functions used to count 
// the number of times a suffix of a given rank occurs in a data set
//
#include "GapArray.h"
#include "SparseGapArray.h"
#if 0
// SimpleGapArray
SimpleGapArray::SimpleGapArray()
{

}

//
void SimpleGapArray::resize(size_t n)
{
    m_data.resize(n);
}

//
void SimpleGapArray::increment(size_t i)
{
    static size_t max_gap_count = std::numeric_limits<GAP_TYPE>::max();
    assert(i < m_data.size());    
    assert(m_data[i] < max_gap_count);
    ++m_data[i];
}

//
size_t SimpleGapArray::get(size_t i) const
{
    return m_data[i];
}

// 
size_t SimpleGapArray::size() const
{
    return m_data.size();
}
#endif

// Construct a gap array for the given underlying storage storage
GapArray* createGapArray(int storage)
{
    switch(storage)
    {
        case 1:
            return new SparseGapArray1;
        case 4:
            return new SparseGapArray4;
        case 8:
            return new SparseGapArray8;
        case 16:
            return new SparseGapArray16;
        case 32:
            return new SparseGapArray32;
        default:
        {
            std::cerr << "Invalid gap array storage parameter: " << storage << "\n";
            exit(1);
        }
    }
}

// Increment the gap array for each suffix of seq. Not thread safe.
void updateGapArray(const DNAString& w, const BWT* pBWTInternal, GapArray* pGapArray)
{
    (void)w;
    (void)pBWTInternal;
    (void)pGapArray;

    assert(false);

#if 0
    size_t l = w.length();
    int i = l - 1;

    // Compute the rank of the last character of seq. We consider $ to be implicitly
    // terminated by a $ character. The first rank calculated is for this and it is given
    // by the C(a) array in BWTInternal
    int64_t rank = pBWTInternal->getPC('$'); // always zero
    pGapArray->incrementSerial(rank);

    // Compute the starting rank for the last symbol of w
    char c = w.get(i);
    rank = pBWTInternal->getPC(c);
    pGapArray->incrementSerial(rank);
    --i;

    // Iteratively compute the remaining ranks
    while(i >= 0)
    {
        char c = w.get(i);
        rank = pBWTInternal->getPC(c) + pBWTInternal->getOcc(c, rank - 1);
        pGapArray->incrementSerial(rank);
        --i;
    }
#endif
}

//
void analyzeGapArray(GapArray* pGapArray)
{
    size_t thresh_8 = 255;
    size_t thresh_16 = 65535;
    size_t count_8 = 0;
    size_t count_16 = 0;
    size_t count_0 = 0;
    size_t count_1 = 0;
    size_t count_2 = 0;
    size_t count_4 = 0;

    for(size_t i = 0; i < pGapArray->size(); ++i)
    {
        size_t c = pGapArray->get(i);
        printf("GA\t%zu\t%zu\n", i, c);
        if(c >= thresh_8)
            ++count_8;
        if(c >= thresh_16)
            ++count_16;
        if(c == 0)
            ++count_0;
        if(c == 1)
            ++count_1;
        if(c == 2)
            ++count_2;
        if(c >= 16)
            ++count_4;
    }

    printf("Num >= %zu: %zu\n", thresh_8, count_8);
    printf("Num >= %zu: %zu\n", thresh_16, count_16);
    printf("Num >= 16: %zu\n", count_4);
    printf("Num == 0: %zu\n", count_0);
    printf("Num == 1: %zu\n", count_1);
    printf("Num == 2: %zu\n", count_2);
}
