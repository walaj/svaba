//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// DNACodec - Encoder/decoder for a DNA alphabet
// of A,C,G,T using 2 bits per symbol
//
#ifndef DNACODEC_H
#define DNACODEC_H

#include "Alphabet.h"

static const unsigned char dna_offset_mask[]={0xC0,0x30,0x0C,0x03};

class DNACodec
{
    public:
        typedef uint8_t UNIT_TYPE;
        static const int SYMBOLS_PER_UNIT = 4;

        // Encoded the character b into a value
        inline uint8_t encode(char b)
        {
            return DNA_ALPHABET::getBaseRank(b);
        }

        // Decode the value c into a character
        inline char decode(uint8_t c)
        {
            return DNA_ALPHABET::getBase(c);
        }

        // Returns the number of units required to encode
        // a string of length n
        inline size_t getRequiredUnits(size_t n)
        {
            return (n + SYMBOLS_PER_UNIT - 1) / SYMBOLS_PER_UNIT;
        }

        // Returns the number of symbols that can be stored
        // in n units
        inline size_t getCapacity(size_t n)
        {
            return n * SYMBOLS_PER_UNIT;
        }

        // Returns the index of the unit to store the
        // i-th symbol of the string
        inline size_t getUnitIndex(size_t i)
        {
            return i / SYMBOLS_PER_UNIT;
        }

        // Returns the position within a unit that the i-th symbol 
        // should be stored in
        inline size_t getUnitOffset(size_t i)
        {
            // this position is the k-th symbol of the unit
            return i % SYMBOLS_PER_UNIT;
        }

        // Return the amount that a value must be shifted to
        // store a code at a given offset
        inline uint8_t getOffsetShift(size_t offset)
        {
            return 2*(3 - offset);
        }

        // Store the value v at the i-th encoded position in the data array
        inline void store(UNIT_TYPE* pData, size_t i, char b)
        {
            UNIT_TYPE& unit = pData[getUnitIndex(i)];
            size_t offset = getUnitOffset(i);
            uint8_t shift = getOffsetShift(offset);

            // Clear the 
            unsigned char mask = dna_offset_mask[offset];

            // Clear position
            unit &= ~mask;

            // Set position
            uint8_t code = encode(b);
            code <<= shift;
            unit |= code; 
        }

        // get the character stored at position i
        inline char get(const UNIT_TYPE* pData, size_t i)
        {
            const UNIT_TYPE& unit = pData[getUnitIndex(i)];
            size_t offset = getUnitOffset(i);
            unsigned char mask = dna_offset_mask[offset];
            UNIT_TYPE code = unit & mask;
            uint8_t shift = getOffsetShift(offset);
            code >>= shift;
            return decode(code);
        }
};

#endif
