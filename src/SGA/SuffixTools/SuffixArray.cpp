//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SuffixArray - Generalized suffix array
//
#include "SuffixArray.h"
#include "InverseSuffixArray.h"
#include "bucketSort.h"
#include "SuffixCompare.h"
#include "mkqs.h"
#include "SACAInducedCopying.h"
#include "Timer.h"
#include "SAReader.h"
#include "SAWriter.h"
#include "BWTWriter.h"

// Read a suffix array from a file
SuffixArray::SuffixArray(const std::string& filename)
{
    SAReader reader(filename);
    reader.read(this);
}

// Construct the suffix array for a table of reads
SuffixArray::SuffixArray(const ReadTable* pRT, int numThreads, bool silent)
{
    //Timer timer("SuffixArray Construction", silent);
    saca_induced_copying(this, pRT, numThreads, silent);
}

// Initialize a suffix array for the strings in RT
void SuffixArray::initialize(const ReadTable& rt)
{
    size_t n = rt.countSumLengths() + rt.getCount(); 
    initialize(n, rt.getCount());

    // Fill the data table with the linear ordering of the suffixes
    size_t count = 0;
    for(size_t i = 0; i < rt.getCount(); ++i)
    {
        // + 1 below is for the empty suffix (is it actually needed?)
        for(size_t j = 0; j < rt.getRead(i).seq.length() + 1; ++j)
        {
            m_data[count++] = SAElem(i, j);
        }
    }
}

void SuffixArray::initialize(size_t num_suffixes, size_t num_strings)
{
    m_data.resize(num_suffixes);
    m_numStrings = num_strings;
}

// Validate the suffix array using the read table
void SuffixArray::validate(const ReadTable* pRT) const
{
    size_t maxIdx = pRT->getCount();
    size_t n = m_data.size();

    // Exit if there is nothing to do
    if(n == 0)
        return;

    // Compute the ISA
    InverseSuffixArray isa(*this);

    // Validate the ISA is a permutation of 1..n, this implies that the id,pos pairs of the SA are valid
    isa.validate();

    size_t empty_count = 0;
    // Ensure that the suffix at pos i is lexographically lower than the suffix at i + 1 using the full string
    for(size_t i = 0; i < n - 1; ++i)
    {
        SAElem id1 = m_data[i];
        SAElem id2 = m_data[i+1];

        assert(id1.getID() < maxIdx);
        std::string suffix1 = pRT->getRead(id1.getID()).seq.getSuffixString(id1.getPos());
        std::string suffix2 = pRT->getRead(id2.getID()).seq.getSuffixString(id2.getPos());

        if(suffix1.length() == 1)
            ++empty_count;

        bool suffixValidated = true;
        if(suffix1 == suffix2)
        {
            suffixValidated = id1.getID() < id2.getID();
        }
        else
        {
            suffixValidated = suffix1 < suffix2;
        }

        if(!suffixValidated)
        {
            std::cerr << "Validation failure: " << suffix1 << " is not less than " << suffix2
                        << " ids: " << id1.getID() << "," << id2.getID() << "\n";
            assert(suffix1 < suffix2);
        }
    }

    assert(m_numStrings == empty_count);
}

// 
void SuffixArray::removeReads(const NumericIDSet& idSet)
{
    // Early exit if the idset is empty
    if(idSet.empty())
        return;

    SAElemVector newData;
    newData.reserve(m_data.size());

    for(size_t idx = 0; idx < m_data.size(); ++idx)
    {
        SAElem id = m_data[idx];
        if(idSet.find(id.getID()) == idSet.end())
        {
            // not on the delete list
            newData.push_back(id);
        }
    }

    m_data.swap(newData);
    m_numStrings -= idSet.size();
}


// Get the suffix cooresponding to idx using the read table
std::string SuffixArray::getSuffix(size_t idx, const ReadTable* pRT) const
{
    SAElem id = m_data[idx];
    return pRT->getRead(id.getID()).seq.getSuffixString(id.getPos());
}

// Return the length of the suffix corresponding to elem 
size_t SuffixArray::getSuffixLength(const ReadTable* pRT, const SAElem elem) const
{
    size_t readLength = pRT->getReadLength(elem.getID());
    return readLength - elem.getPos();
}


// Print the suffix array
void SuffixArray::print(const ReadTable* pRT) const
{
    std::cout << "i\tSA(i)\n";
    for(size_t i = 0; i < m_data.size(); ++i)
    {
        SAElem id1 = m_data[i];
        std::string suffix = !id1.isEmpty() ? getSuffix(i, pRT) : "";
        int pos = id1.getPos();
        char b;
        if(pos == 0)
            b = '$';
        else
            b =  pRT->getRead(id1.getID()).seq.get(pos - 1);
        std::cout << i << "\t" << id1 << "\t" << b << "\t" << suffix << "\n";
    }
}

// Print the suffix array
void SuffixArray::print() const
{
    std::cout << "i\tSA(i)\n";
    for(size_t i = 0; i < m_data.size(); ++i)
    {
        std::cout << i << "\t" << m_data[i] << "\n";
    }
}

// write the suffix array to a file
void SuffixArray::write(const std::string& filename)
{
    SAWriter writer(filename);
    writer.write(this);
}

// write the suffix array to a file
void SuffixArray::writeBWT(const std::string& filename, const ReadTable* pRT)
{
    //RLBWTWriter writer(filename);
    IBWTWriter* pWriter = BWTWriter::createWriter(filename);
    pWriter->write(this, pRT);
    delete pWriter;
}

// write the index of the suffix array to a file
void SuffixArray::writeIndex(std::string& filename)
{
    SAWriter writer(filename);
    writer.writeHeader(m_numStrings, m_numStrings);
    for(size_t i = 0; i < m_data.size(); ++i)
    {
      if(m_data[i].isFull()) 
            writer.writeElem(m_data[i]);
    }
}

/*template<typename T>
inline std::vector<T> erase_indices(const std::vector<T>& data, std::vector<size_t>& indicesToDelete)
{
  if(indicesToDelete.empty())
    return data;

  std::vector<T> ret;
  ret.reserve(data.size() - indicesToDelete.size());

  std::sort(indicesToDelete.begin(), indicesToDelete.end());

  // new we can assume there is at least 1 element to delete. copy blocks at a time.
  std::vector<T>::const_iterator itBlockBegin = data.begin();
  for(std::vector<size_t>::const_iterator it = indicesToDelete.begin(); it != indicesToDelete.end(); ++ it)
    {
      std::vector<T>::const_iterator itBlockEnd = data.begin() + *it;
      if(itBlockBegin != itBlockEnd)
        {
	  std::copy(itBlockBegin, itBlockEnd, std::back_inserter(ret));
        }
      itBlockBegin = itBlockEnd + 1;
    }

  // copy last block.
  if(itBlockBegin != data.end())
    {
      std::copy(itBlockBegin, data.end(), std::back_inserter(ret));
    }

  return ret;
}*/

//JEREMIAH
void SuffixArray::writeIndex()
{
    SAElemVector sae;
    for(int i = 1; i < m_data.size(); ++i)
    {
      if(m_data[i].isFull()) 
        sae.push_back(m_data[i]);
    }

    m_data = sae;
}


// Output operator
std::ostream& operator<<(std::ostream& out, const SuffixArray& sa)
{
    // Write the size and number of strings
    out << sa.m_data.size() << "\n";
    out << sa.m_numStrings << "\n";

    for(size_t i = 0; i < sa.m_data.size(); ++i)
    {
        out << sa.m_data[i] << "\n";
    }
    return out;
}


// Input operator
std::istream& operator>>(std::istream& in, SuffixArray& sa)
{
    // Read the size and number of strings
    size_t n;
    in >> n;
    in >> sa.m_numStrings;

    sa.m_data.resize(n);
    size_t i = 0;
    SAElem id;
    while(in >> id)
        sa.m_data[i++] = id;
    assert(i == n);
    return in;
}
