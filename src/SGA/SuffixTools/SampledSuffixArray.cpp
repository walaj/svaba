//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SampledSuffixArray - Data structure holding
// a subset of suffix array positions. From the
// sampled positions, the other entries of the
// suffix array can be calculated.
//
#include "SampledSuffixArray.h"
#include "SAReader.h"
#include "SAWriter.h"
#include "config.h"

#if HAVE_OPENMP
#include <omp.h>
#endif

static const uint32_t SSA_MAGIC_NUMBER = 12412;
#define SSA_READ(x) pReader->read(reinterpret_cast<char*>(&(x)), sizeof((x)));
#define SSA_READ_N(x,n) pReader->read(reinterpret_cast<char*>(&(x)), (n));

#define SSA_WRITE(x) pWriter->write(reinterpret_cast<const char*>(&(x)), sizeof((x)));
#define SSA_WRITE_N(x,n) pWriter->write(reinterpret_cast<const char*>(&(x)), (n));

//
SampledSuffixArray::SampledSuffixArray()
{

}

SampledSuffixArray::SampledSuffixArray(const std::string& filename, SSAFileType filetype)
{
    // Read the sampled suffix array from a file - either from a .ssa or .sai file
    if(filetype == SSA_FT_SSA)
        readSSA(filename);
    else
        readSAI(filename);
}

// 
SAElem SampledSuffixArray::calcSA(int64_t idx, const BWT* pBWT) const
{
    size_t offset = 0;
    SAElem elem;

    while(1)
    {
        // Check if this position is sampled. If the sample rate is zero we are using the lexo. index only
        if(m_sampleRate > 0 && idx % m_sampleRate == 0 && !m_saSamples[idx / m_sampleRate].isEmpty())
        {
            // A valid sample is stored for this idx
            elem = m_saSamples[idx / m_sampleRate];
            break;
        }

        // A sample does not exist for this position, perform a backtracking step
        char b = pBWT->getChar(idx);
        idx = pBWT->getPC(b) + pBWT->getOcc(b, idx - 1);

        if(b == '$')
        {
            // idx (before the update) corresponds to the start of a read.
            // We can directly look up the saElem for idx from the lexicographic index
            assert(idx < (int64_t)m_saLexoIndex.size());
            elem.setID(m_saLexoIndex[idx]);
            elem.setPos(0);
            break;
        }
        else
        {
            // A backtracking step is performed, increment offset
            offset += 1;
        }
    }

    elem.setPos(elem.getPos() + offset);
    return elem;
}

// Returns the ID of the read with lexicographic rank r
size_t SampledSuffixArray::lookupLexoRank(size_t r) const
{
    return m_saLexoIndex[r];
}

// 
void SampledSuffixArray::build(const BWT* pBWT, const ReadInfoTable* pRIT, int sampleRate)
{
    m_sampleRate = sampleRate;

    size_t numStrings = pRIT->getCount();
    m_saLexoIndex.resize(numStrings);

    size_t MAX_ELEMS = std::numeric_limits<SSA_INT_TYPE>::max();
    if(numStrings > MAX_ELEMS)
    {
        std::cerr << "Error: Only " << MAX_ELEMS << " reads are allowed in the sampled suffix array\n";
        std::cerr << "Number of reads in your index: " << numStrings << "\n";
        std::cerr << "Contact sga-users@googlegroups.com for help\n";
        exit(EXIT_FAILURE);
    }

    // Set the size of the sampled vector
    size_t numElems = (pBWT->getBWLen() / m_sampleRate) + 1;
    m_saSamples.resize(numElems);

    // For each read, start from the end of the read and backtrack through the suffix array/BWT.
    // For every idx that is divisible by the sample rate, store the calculate SAElem
    for(size_t i = 0; i < numStrings; ++i)
    {
        // The suffix array positions for the ends of reads are ordered
        // by their position in the read information table, therefore
        // the starting suffix array index is i
        size_t idx = i;

        // The ID of the read is i. The position coordinate is inclusive but 
        // since the read information table does not store the '$' symbol
        // the starting position equals the read length
        SAElem elem(i, pRIT->getReadLength(i));

        while(1)
        {
            if(idx % m_sampleRate == 0)
            {
                // store this SAElem
                m_saSamples[idx / m_sampleRate] = elem;
            }

            char b = pBWT->getChar(idx);
            idx = pBWT->getPC(b) + pBWT->getOcc(b, idx - 1);
            if(b == '$')
            {
                // we have hit the beginning of this string
                // store the SAElem for the beginning of the read
                // in the lexicographic index
                if(elem.getPos() != 0)
                    std::cout << "elem: " << elem << " i: " << i << "\n";

                assert(elem.getPos() == 0);
                size_t id = elem.getID();
                assert(id < MAX_ELEMS);
                m_saLexoIndex[idx] = id;
                break; // done;
            }
            else
            {
                // Decrease the position of the elem
                elem.setPos(elem.getPos() - 1);
            }
        }
    }
}

// A streamlined version of the above function
void SampledSuffixArray::buildLexicoIndex(const BWT* pBWT, int num_threads)
{
    int64_t numStrings = pBWT->getNumStrings();
    m_saLexoIndex.resize(numStrings);
    int64_t MAX_ELEMS = std::numeric_limits<SSA_INT_TYPE>::max();
    assert(numStrings < MAX_ELEMS);

    (void)num_threads;
    // Parallelize this computaiton using openmp, if the compiler supports it
#if HAVE_OPENMP
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
#endif
    for(int64_t read_idx = 0; read_idx < numStrings; ++read_idx)
    {
        // For each read, start from the end of the read and backtrack through the suffix array/BWT
        // to calculate its lexicographic rank in the collection
        size_t idx = read_idx;
        while(1)
        {
            char b = pBWT->getChar(idx);
            idx = pBWT->getPC(b) + pBWT->getOcc(b, idx - 1);
            if(b == '$')
            {
                // There is a one-to-one mapping between read_index and the element
                // of the array that is set - therefore we can perform this operation
                // without a lock.
                m_saLexoIndex[idx] = read_idx;
                break; // done;
            }
        }
    }
}

// Validate the sampled suffix array values are correct
void SampledSuffixArray::validate(const std::string filename, const BWT* pBWT)
{
    ReadTable* pRT = new ReadTable(filename);
    SuffixArray* pSA = new SuffixArray(pRT, 1);
    
    std::cout << "Validating sampled suffix array entries\n";

    for(size_t i = 0; i < pSA->getSize(); ++i)
    {
        SAElem calc = calcSA(i, pBWT);
        SAElem real = pSA->get(i);
        if(calc.getID() != real.getID() || calc.getPos() != real.getPos())
        {
            std::cout << "Error SA elements do not match for " << i << "\n";
            std::cout << "Calc: " << calc << "\n";
            std::cout << "Real: " << real << "\n";
            exit(1);
        }
    }
    
    std::cout << "All calculate SA values are correct\n";

    delete pRT;
    delete pSA;
}

// Save the SA to disc
void SampledSuffixArray::writeSSA(std::string filename)
{
    std::ostream* pWriter = createWriter(filename, std::ios::out | std::ios::binary);
    
    // Write a magic number
    SSA_WRITE(SSA_MAGIC_NUMBER)

    // Write sample rate
    SSA_WRITE(m_sampleRate)

    // Write number of lexicographic index entries
    size_t n = m_saLexoIndex.size();
    SSA_WRITE(n)

    // Write lexo index
    SSA_WRITE_N(m_saLexoIndex.front(), sizeof(SSA_INT_TYPE) * n)
    
    // Write number of samples
    n = m_saSamples.size();
    SSA_WRITE(n)

    // Write samples
    SSA_WRITE_N(m_saSamples.front(), sizeof(SAElem) * n)

    delete pWriter;
}

// Save just the lexicographic index portion of the SSA to disk as plaintext
void SampledSuffixArray::writeLexicoIndex(const std::string& filename)
{
    SAWriter writer(filename);
    size_t num_strings = m_saLexoIndex.size();
    writer.writeHeader(num_strings, num_strings);
    for(size_t i = 0; i < m_saLexoIndex.size(); ++i) 
    {
        SAElem elem(m_saLexoIndex[i], 0);
        writer.writeElem(elem);
    }
}


void SampledSuffixArray::readSSA(std::string filename)
{
    std::istream* pReader = createReader(filename, std::ios::binary);
    
    // Write a magic number
    uint32_t magic = 0;
    SSA_READ(magic)
    assert(magic == SSA_MAGIC_NUMBER);

    // Read sample rate
    SSA_READ(m_sampleRate)

    // Read number of lexicographic index entries
    size_t n = 0;
    SSA_READ(n)
    m_saLexoIndex.resize(n);

    // Read lexo index
    SSA_READ_N(m_saLexoIndex.front(), sizeof(SSA_INT_TYPE) * n)
    
    // Read number of samples
    n = 0;
    SSA_READ(n)
    m_saSamples.resize(n);

    // Read samples
    SSA_READ_N(m_saSamples.front(), sizeof(SAElem) * n)

    delete pReader;
}

void SampledSuffixArray::readSAI(std::string filename)
{
    SAReader reader(filename);
    size_t num_strings, num_elems;
    reader.readHeader(num_strings, num_elems);
    assert(num_strings == num_elems);
    m_saLexoIndex.reserve(num_strings);
    reader.readElems(m_saLexoIndex);

    // Set the sample rate to zero to signify there are no samples
    m_sampleRate = 0;
}

// Print memory usage information
void SampledSuffixArray::printInfo() const
{
    double mb = (double)(1024*1024);
    double lexoSize = (double)(sizeof(SSA_INT_TYPE) * m_saLexoIndex.capacity()) / mb;
    double sampleSize = (double)(sizeof(SAElem) * m_saSamples.capacity()) / mb;
    
    printf("SampledSuffixArray info:\n");
    printf("Sample rate: %d\n", m_sampleRate);
    printf("Contains %zu entries in lexicographic array (%.1lf MB)\n", m_saLexoIndex.size(), lexoSize);
    printf("Contains %zu entries in sample array (%.1lf MB)\n", m_saSamples.size(), sampleSize);
    printf("Total size: %.1lf\n", lexoSize + sampleSize);
}
