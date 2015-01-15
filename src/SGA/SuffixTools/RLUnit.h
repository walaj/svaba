//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// RLUnit - A run-length encoded unit of the FM-index
//
#ifndef RLUNIT_H
#define RLUNIT_H

//
#define RL_COUNT_MASK 0x1F  //00011111
#define RL_SYMBOL_MASK 0xE0 //11100000
#define RL_FULL_COUNT 31
#define RL_SYMBOL_SHIFT 5
#define RLE_VALIDATE

// A unit of the RLBWT is a pair of a symbol and its count
// The high 3 bits encodes the symbol to store
// The low 5 bits encodes the length of the run
struct RLUnit
{
    RLUnit() : data(0) {}
    RLUnit(char b) : data(1)
    {
        setChar(b);   
    }

    // Returns true if the count cannot be incremented
    inline bool isFull() const
    {
        return (data & RL_COUNT_MASK) == RL_FULL_COUNT;
    }

    inline bool isEmpty() const
    {
        return (data & RL_COUNT_MASK) == 0;
    }

    inline bool isInitialized() const
    {
        return data > 0;
    }

    // Add this run to the AlphaCount
    // Only add up to maxCount symbols. Returns the number
    // of symbols added
    inline size_t addAlphaCount(AlphaCount64& ac, size_t max) const
    {
        size_t count = getCount();
        if(count > max)
            count = max;
        char symbol = getChar();
        ac.add(symbol, count);
        return count;
    }

    // Add this run to the count of base b if it matches
    // Only add up to maxCount symbols. Returns the number
    // of symbols in the current run, up to max
    inline size_t addCount(char b, size_t& base_count, size_t max) const
    {
        size_t run_len = getCount();
        if(run_len > max)
            run_len = max;
        char symbol = getChar();
        if(symbol == b)
            base_count += run_len;
        return run_len;
    }    

    // Subtract this run from the AlphaCount
    // Only subtract up to maxCount symbols. Returns the number
    // of symbols added
    inline size_t subtractAlphaCount(AlphaCount64& ac, size_t max) const
    {
        size_t count = getCount();
        if(count > max)
            count = max;
        char symbol = getChar();
        ac.subtract(symbol, count);
        return count;
    }

    // Subtract this run from the count of base b if it matches
    // Only add up to maxCount symbols. Returns the number
    // of symbols in the current run, up to max
    inline size_t subtractCount(char b, size_t& base_count, size_t max) const
    {
        size_t run_len = getCount();
        if(run_len > max)
            run_len = max;
        char symbol = getChar();
        if(symbol == b)
            base_count -= run_len;
        return run_len;
    }
        
    // 
    inline void incrementCount()
    {
#ifdef RLE_VALIDATE
        assert(!isFull());
#endif
        ++data;
    }

    // 
    inline void decrementCount()
    {
#ifdef RLE_VALIDATE
        assert(!isEmpty());
#endif
        --data;
    }    

    inline uint8_t getCount() const
    {
#ifdef RLE_VALIDATE
        assert((data & RL_COUNT_MASK) != 0);
#endif
        return data & RL_COUNT_MASK;
    }

    // Set the symbol
    inline void setChar(char symbol)
    {
        // Clear the current symbol
        data &= RL_COUNT_MASK;
        
        uint8_t code = BWT_ALPHABET::getRank(symbol);
        code <<= RL_SYMBOL_SHIFT;
        data |= code;
    }

    // Get the symbol
    inline char getChar() const
    {
        uint8_t code = data & RL_SYMBOL_MASK;
        code >>= RL_SYMBOL_SHIFT;
        return BWT_ALPHABET::getChar(code);
    }

    // 
    uint8_t data;

    friend class RLBWTReader;
    friend class RLBWTWriter;
};
typedef std::vector<RLUnit> RLVector;

#endif
