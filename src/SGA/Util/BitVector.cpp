//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// BitVector - Vector of bits 
//
#include "BitVector.h"
#include <assert.h>
#include <cstdlib>

//
BitVector::BitVector()
{
    initializeMutex();
}

//
BitVector::BitVector(size_t n)
{
    initializeMutex();
    resize(n);
}

void BitVector::initializeMutex()
{
    int ret = pthread_mutex_init(&m_mutex, NULL);
    if(ret != 0)
    {
        std::cerr << "Mutex initialization in BitVector failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }
}

//
BitVector::~BitVector()
{
    int ret = pthread_mutex_destroy(&m_mutex);
    if(ret != 0)
    {
        std::cerr << "Mutex destruction in BitVector failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }
}

//
void BitVector::resize(size_t n)
{
    size_t num_bytes = (n % 8 == 0) ? n / 8 : (n / 8) + 1;
    m_data.resize(num_bytes);
}

//
void BitVector::lock()
{
    pthread_mutex_lock(&m_mutex);
}

//
void BitVector::unlock()
{
    pthread_mutex_unlock(&m_mutex);
}

//
bool BitVector::updateCAS(size_t i, bool oldValue, bool newValue)
{
    size_t byte = i / 8;
    assert(byte < m_data.size());
    size_t offset = i - byte * 8;
    return m_data[byte].updateCAS(offset, oldValue, newValue);
}


// Set bit at position i to value v
void BitVector::set(size_t i, bool v)
{
    size_t byte = i / 8;
    assert(byte < m_data.size());
    size_t offset = i - byte * 8;
    m_data[byte].set(offset, v);
}

// Test bit i
bool BitVector::test(size_t i) const
{
    size_t byte = i / 8;
    size_t offset = i - byte * 8;
    return m_data[byte].test(offset);
}

