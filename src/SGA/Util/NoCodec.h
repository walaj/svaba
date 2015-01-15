//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// NoCodec - Testing codec for EncodedString
// that does not encode anything - it just returns
// the numerical value of the char 
//
#ifndef NOCODEC_H
#define NOCODEC_H

class NoCodec
{
    public:
        typedef char UNIT_TYPE;
        static const int SYMBOLS_PER_UNIT = 1;

        // Encoded the character b into a value
        inline uint8_t encode(char b) const
        {
            return (uint8_t)b;
        }

        // Decode the value c into a character
        inline char decode(uint8_t c) const
        {
            return (char)c;
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
        inline uint8_t getOffsetShift(size_t /*offset*/) const
        {
            return 0;
        }

        // Store the value v at the i-th encoded position in the data array
        inline void store(UNIT_TYPE* pData, size_t i, char b) const
        {
            UNIT_TYPE& unit = pData[getUnitIndex(i)];
            unit = encode(b);
            
            size_t offset = getUnitOffset(i);
            uint8_t shift = getOffsetShift(offset);
            
            // Clear the currrent position
            uint8_t mask = 0xFF;
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
            (void)offset;
            uint8_t mask = 0xFF;
            UNIT_TYPE code = unit & mask;
            uint8_t shift = 0;
            code >>= shift;
            return decode(code);
        }
};

#endif
