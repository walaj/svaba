//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SuffixCompareRadix 
//
#include "SuffixCompare.h"

//
SuffixCompareRadix::SuffixCompareRadix(const ReadTable* pRT) : m_pRT(pRT), m_pNumSuffixLUT(0)
{
    m_bucketOffset = 0;
    m_bucketLen = 6;
    initializeNumSuffixLUT();
}

//
SuffixCompareRadix::SuffixCompareRadix(const ReadTable* pRT, int bucket_len) : m_pRT(pRT),
                                                                                 m_bucketOffset(0),
                                                                                 m_bucketLen(bucket_len),
                                                                                 m_pNumSuffixLUT(0)

{
    if(m_bucketOffset % m_bucketLen != 0)
    {
        std::cerr << "Error, bucket offset must be a multiple of bucket length\n";
        assert(false);
    }
    initializeNumSuffixLUT();
}

// 
SuffixCompareRadix::~SuffixCompareRadix()
{
    delete [] m_pNumSuffixLUT;
}

void SuffixCompareRadix::initializeNumSuffixLUT()
{
    m_pNumSuffixLUT = new int[m_bucketLen];
    for(int i = 0; i < (int)m_bucketLen; ++i)
    {
        m_pNumSuffixLUT[i] = calcNumSuffixes(i);
    }
}

// Get the bucket for a particular SAElem
int SuffixCompareRadix::getBucket(SAElem x) const
{
    //std::cout << "Finding bucket for " << x << "\n";
    const DNAString& read = m_pRT->getRead(x.getID()).seq;

    size_t suffix_start = x.getPos() + m_bucketOffset;
    const char* suffix = read.getSuffix(suffix_start);
    size_t suffix_len = read.getSuffixLength(suffix_start);

    size_t stop = std::min(m_bucketLen, suffix_len);

    int rank = 0;
    for(size_t i = 0; i < stop; ++i)
    {
        char b = suffix[i];
        rank += numPredSuffixes(b, m_bucketLen - i);
    }

    return rank;
}

//
void SuffixCompareRadix::printElem(SAElem& x) const
{
    const DNAString& read = m_pRT->getRead(x.getID()).seq;
    std::cout << x << " = " << read.getSuffixString(x.getPos()) << "\n";
}

//
void SuffixCompareRadix::setBucketDepth(int depth)
{
    m_bucketOffset = m_bucketLen * depth;
}

//
int SuffixCompareRadix::calcNumSuffixes(int maxLen) const
{
    int r = 0;
    for(int i = 0; i <= maxLen; ++i)
        r += (1 << 2*i);
    return r;
}

bool SuffixCompareRadix::isBucketDegenerate(int index) const
{
    if(index == 0)
        return true;
    else
    {
        for(int i = m_bucketLen; i >= 2; --i)
        {
            --index;
            int block_size = calcNumSuffixes(i - 1);
            if(index % block_size == 0)
                return true;
            else
            {
                int block = index / block_size;
                index -= block * block_size;
            }
        }
    }
    return false;
}

//
int SuffixCompareRadix::getNumBuckets() const
{
    return calcNumSuffixes(m_bucketLen);
}

// Compare two suffixes by ID
// Precondition: the suffixes are sorted by sequence already
bool SuffixCompareID::operator()(SAElem x, SAElem y) const
{ 
    const SeqItem& rx = m_pRT->getRead(x.getID());
    const SeqItem& ry = m_pRT->getRead(y.getID());
    return rx.id < ry.id;
}


