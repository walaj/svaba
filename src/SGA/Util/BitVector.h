//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// BitVector - Vector of bits. The structure
// can be locked by a mutex to guarentee atomic access.
//
#ifndef BITVECTOR_H
#define BITVECTOR_H

#include "BitChar.h"
#include <vector>

class BitVector
{
    public:
    
        BitVector();
        BitVector(size_t n);
        ~BitVector();

        // Functions to acquire/release the mutex
        // The client code is responsible for acquiring
        // the lock before calling set(). Reading a bit with 
        // test() is ok however.
        void lock();
        void unlock();

        // Update the bit at position i from oldValue to newValue using a
        // compare and swap operation. Returns true if the update is successful.
        bool updateCAS(size_t i, bool oldValue, bool newValue);

        void resize(size_t n);
        void set(size_t i, bool v);
        bool test(size_t i) const;

        size_t capacity() const { return m_data.size() * 8; }

    private:

        void initializeMutex();

        std::vector<BitChar> m_data;
        pthread_mutex_t m_mutex;
};

#endif
