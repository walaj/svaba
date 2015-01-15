//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// BitChar - Fixed-size bitset of 8 bits
//

#include <iostream>
#include <assert.h>
#include "BitChar.h"

static unsigned char bc_mask[]={0x80,0x40,0x20,0x10,0x08,0x04,0x02,0x01};

//
void BitChar::set(unsigned char idx, bool v)
{
    if(v)
        m_data = m_data | bc_mask[idx];
    else
        m_data = m_data & ~bc_mask[idx];
}

//
bool BitChar::test(unsigned char idx) const
{
    return m_data & bc_mask[idx];
}

//
void BitChar::flip(unsigned idx)
{
    m_data = m_data ^ bc_mask[idx];
}

//
void printBinary(std::ostream& out, const BitChar& bc)
{
    for(int i = 7; i >= 0; i--)
    {
        if(bc.test(i))
            out << "1";
        else
            out << "0";
    }
}

// Update the value of bit idx from oldValue to newValue using an atomic compare and swap operation
// Returns true if the update was successfully performed.
bool BitChar::updateCAS(unsigned char idx, bool oldValue, bool newValue)
{
    assert(oldValue != newValue);
    
    // Iterate attempts of the CAS operation until the desired bit has the correct value.
    while(1)
    {
        // Get the current value of the data. 
        // If any thread updates d after this call the CAS will fail.
        unsigned char oldData = m_data;
        
        // Check if the bit is already set to newValue in our copy of the data. If so,
        // some other thread must have updated the value already so we return false.
        bool currValue = (oldData & bc_mask[idx]);
        if(currValue == newValue)
            return false;

        // Calculate the value of the data with the bit set to the correct value
        unsigned char newData;
        if(newValue)
            newData = oldData | bc_mask[idx]; //set the bit
        else
            newData = oldData & ~bc_mask[idx]; //clear the bit

        // perform the atomic CAS. this returns false if the real value of d does not match oldData
        bool success = __sync_bool_compare_and_swap(&m_data, oldData, newData);

        // If the update was successfully performed, return true.
        // Otherwise, some other thread must have updated a bit in this char so we iterate again.
        if(success)
            return true;
    }
}

//
std::ostream& operator<<(std::ostream& out, const BitChar& bc)
{
    int temp = bc.m_data;
    out << temp;
    return out;
}

//
std::istream& operator>>(std::istream& in, BitChar& bc)
{
    int temp;
    in >> temp;
    bc.m_data = (unsigned char)temp;
    return in;
}

// 
void BitChar::write(std::ostream& out)
{
    out.write((char*)&m_data, sizeof(m_data));
}

//
void BitChar::read(std::istream& in)
{
    in.read((char*)&m_data, sizeof(m_data));
}
