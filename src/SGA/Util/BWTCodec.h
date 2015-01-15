//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// BWTCodec - Encoder/decoder for a BWT
// string over the alphabet $ACGT using
// 3 bits per symbol
//
#ifndef BWTCODEC_H
#define BWTCODEC_H

#include "Alphabet.h"

static const uint16_t bwt_offset_mask[]={0xE000,0x1C00,0x0380,0x0070,0x000E};
static const uint8_t  bwt_offset_shift[]={13, 10, 7, 4, 1};

class BWTCodec
{
    public:
        // The data is stored 5 characters per uint16_t
        // The data pattern is:
        // 11100000 00000000 first symbol
        // 00011100 00000000 second symbol, etc
        typedef uint16_t UNIT_TYPE;
        static const int SYMBOLS_PER_UNIT = 5;

        // Encoded the character b into a value
        inline uint8_t encode(char b) const
        {
            return BWT_ALPHABET::getRank(b);
        }

        // Decode the value c into a character
        inline char decode(uint8_t c) const
        {
            return BWT_ALPHABET::getChar(c);
        }

        // Returns the number of units required to encode
        // a string of length n
        inline size_t getRequiredUnits(size_t n) const
        {
            return (n + SYMBOLS_PER_UNIT - 1) / SYMBOLS_PER_UNIT;
        }

        // Returns the number of symbols that can be stored
        // in n units
        inline size_t getCapacity(size_t n) const
        {
            return n * SYMBOLS_PER_UNIT;
        }

        // Returns the index of the unit to store the
        // i-th symbol of the string
        inline size_t getUnitIndex(size_t i) const
        {
            return i / SYMBOLS_PER_UNIT;
        }

        // Returns the position within a unit that the i-th symbol 
        // should be stored in
        inline size_t getUnitOffset(size_t i) const
        {
            // this position is the k-th symbol of the unit
            return i % SYMBOLS_PER_UNIT;
        }

        // Return the amount that a value must be shifted to
        // store a code at a given offset
        inline uint8_t getOffsetShift(size_t offset) const
        {
            return 1 + 3*(4 - offset);
        }

        // Store the value v at the i-th encoded position in the data array
        inline void store(UNIT_TYPE* pData, size_t i, char b) const
        {
            UNIT_TYPE& unit = pData[getUnitIndex(i)];
            size_t offset = getUnitOffset(i);
            uint8_t shift = getOffsetShift(offset);
            
            // Clear the currrent position
            uint16_t mask = bwt_offset_mask[offset];
            unit &= ~mask;

            // Set position
            UNIT_TYPE code = encode(b);
            code <<= shift;
            unit |= code; 
        }

        // get the character stored at position i
        inline char get(const UNIT_TYPE* pData, size_t i) const
        {
            const UNIT_TYPE& unit = pData[getUnitIndex(i)];
            size_t offset = getUnitOffset(i);
            uint16_t mask = bwt_offset_mask[offset];
            UNIT_TYPE code = unit & mask;
            uint8_t shift = bwt_offset_shift[offset];
            code >>= shift;
            return decode(code);
        }
};

#endif
