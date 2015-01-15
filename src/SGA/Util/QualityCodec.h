//-----------------------------------------------
// Copyright 2012 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// QualityCodec - Encoder/decoder for a quality
// symbols. Supports downsampling to a small range.
//
#ifndef QUALITYCODEC_H
#define QUALITYCODEC_H

#include <inttypes.h>
#include <iostream>

typedef uint8_t QualityStorageUnit;

template<int BITS>
class QualityCodec
{
    public:
        typedef QualityStorageUnit UNIT_TYPE;
        static const int SYMBOLS_PER_UNIT = 8 / BITS;

        // These values are the lowest and highest quality
        // values that are stored. The other values are interpolated
        // between these range
        static const int PHRED_MIN = 0;
        static const int PHRED_MAX = 40;
        static const int QUALITY_RANGE = PHRED_MAX - PHRED_MIN;
        static const int QUALITY_BINS = 1 << BITS;
        static const int BIN_SIZE = (QUALITY_RANGE + QUALITY_BINS - 1) / QUALITY_BINS; // ceiling

        // Encode the character b into a value
        inline uint8_t encode(char b) const
        {
            // Compute the quality bin to place this value in
            int phred = Quality::char2phred(b);

            // Clamp values
            if(phred < 0)
                phred = 0;
            UNIT_TYPE bin_idx = (phred - PHRED_MIN) / BIN_SIZE;

            if(bin_idx >= QUALITY_BINS)
                bin_idx = QUALITY_BINS - 1;
            //std::cout << "Q: " << b << " phred: " << phred << " c: " << (int)bin_idx << "\n";
            //decode(bin_idx);
            return bin_idx;
        }

        // Decode the value c into a character
        inline char decode(uint8_t c) const
        {
            // Return the midpoint value of the bin
            int phred = (2 * c + 1) * BIN_SIZE / 2;
            char b = Quality::phred2char(phred);
            //std::cout << "c: " << (int)c << " phred: " << phred << " Q: " << b << "\n";
            return b;
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
            return BITS*(SYMBOLS_PER_UNIT - 1 - offset);
        }

        // Return a mask that can be used to select a single value by ANDing a storage unit
        inline uint8_t getOffsetMask(size_t offset) const
        {
            return ((1 << BITS) - 1) << getOffsetShift(offset);
        }

        // Store the value v at the i-th encoded position in the data array
        inline void store(UNIT_TYPE* pData, size_t i, char b) const
        {
            UNIT_TYPE& unit = pData[getUnitIndex(i)];
            size_t offset = getUnitOffset(i);
            uint8_t shift = getOffsetShift(offset);

            // Clear the 
            unsigned char mask = getOffsetMask(offset);


            // Clear position
            unit &= ~mask;

            // Set position
            uint8_t code = encode(b);
            code <<= shift;
            unit |= code; 
            
            /*
            std::cout << "E Unit: " << getUnitIndex(i) << "\n";
            std::cout << "E Offset: " << offset << "\n";
            std::cout << "E Shift: " << (int)shift << "\n";
            std::cout << "E Mask: " << (int)mask << "\n";
            std::cout << "E UnitV: " << (int)unit << "\n";
            */
        }

        // get the character stored at position i
        inline char get(const UNIT_TYPE* pData, size_t i) const
        {
            const UNIT_TYPE& unit = pData[getUnitIndex(i)];
            size_t offset = getUnitOffset(i);
            unsigned char mask = getOffsetMask(offset);
            UNIT_TYPE code = unit & mask;
            uint8_t shift = getOffsetShift(offset);
            code >>= shift;
            
            /*
            std::cout << "D Unit: " << getUnitIndex(i) << "\n";
            std::cout << "D Offset: " << offset << "\n";
            std::cout << "D Shift: " << (int)shift << "\n";
            std::cout << "D Mask: " << (int)mask << "\n";
            std::cout << "D UnitV: " << (int)unit << "\n";
            */
            return decode(code);
        }
};

#endif
