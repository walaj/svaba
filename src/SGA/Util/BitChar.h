//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// BitChar - Fixed-size bitset of 8 bits
//
#ifndef BITCHAR_H
#define BITCHAR_H

#include <iostream>

struct BitChar
{
    public:
        BitChar() : m_data(0) {}
        
        // Update the value of bit idx from oldValue to newValue using an atomic compare and swap operation
        // Returns true if the update was successfully performed. This is the only thread-safe way to update the BitChar
        bool updateCAS(unsigned char idx, bool oldValue, bool newValue);

        // set the bit at idx to the value u
        void set(unsigned char idx, bool v);

        // returns true if bit at idx is set
        bool test(unsigned char idx) const;

        // flips the bit at idx
        void flip(unsigned idx);

        // I/O
        friend std::ostream& operator<<(std::ostream& out, const BitChar& bc);
        friend std::istream& operator>>(std::istream& in, BitChar& bc);

        void write(std::ostream& out);
        void read(std::istream& in);

    private:
        volatile unsigned char m_data;
};

#endif
