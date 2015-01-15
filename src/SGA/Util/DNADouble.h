//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// DNADouble.h - Simple map from the DNA
// alphabet {A,C,G,T} to a double value
//
#ifndef ALPHADOUBLE_H
#define ALPHADOUBLE_H

#include "Alphabet.h"

class DNADouble
{
    public:
        
        //
        inline DNADouble()
        {
            memset(m_vals, 0, DNA_ALPHABET::size * sizeof(double));
        }

        //
        inline size_t getAlphabetSize() const
        {
            return DNA_ALPHABET::size;
        }

        //
        inline void set(char b, double lp)
        {
            m_vals[DNA_ALPHABET::getBaseRank(b)] = lp;
        }

        // 
        inline double get(char b) const
        {
            return m_vals[DNA_ALPHABET::getBaseRank(b)];
        }

        //
        inline char getMaxBase() const
        {
            char base;
            double val;
            getMax(base, val);
            return base;
        }

        //
        inline double getMaxVal() const
        {
            char base;
            double val;
            getMax(base, val);
            return val;
        }

        //
        inline void getMax(char& base, double& val) const
        {
            double max = -std::numeric_limits<double>::max();
            int maxIdx = 0;
            for(int i = 0; i < DNA_ALPHABET::size; ++i)
            {
                if(m_vals[i] > max)
                {
                    max = m_vals[i];
                    maxIdx = i;
                }
            }
            base = DNA_ALPHABET::getBase(maxIdx);
            val = max;
        }


        //
        inline double getByIdx(const int i) const
        {
            return m_vals[i];
        }

        // Return the base for index i
        static char getBase(size_t i)
        {
            assert(i < DNA_ALPHABET::size);
            return DNA_ALPHABET::getBase(i);
        }

        // Swap the (A,T) and (C,G) entries
        inline void complement()
        {
            double tmp;

            // A,T
            tmp = m_vals[3];
            m_vals[3] = m_vals[0];
            m_vals[0] = tmp;

            // C,G
            tmp = m_vals[2];
            m_vals[2] = m_vals[1];
            m_vals[1] = tmp;
        }

        // Interpret the DNADouble as a log-probabilty vector of some 
        // event occurring and marginalize using a flat prior 
        double marginalize(double prior) const
        {
            double direct = 0.0f;
            for(int i = 0; i < DNA_ALPHABET::size; ++i)
            {
                direct += exp(m_vals[i]);
            }
            double ld = log(direct) + log(prior);
            /*
            double max = -std::numeric_limits<double>::max();
            int max_idx = 0;
            for(int i = 0; i < DNA_ALPHABET::size; ++i)
            {
                if(m_vals[i] > max)
                {
                    max = m_vals[i];
                    max_idx = i;
                }
            }

            double sum = prior;
            for(int i = 0; i < DNA_ALPHABET::size; ++i)
            {
                if(i != max_idx)
                    sum += exp(m_vals[i] - max) * prior;
            }
            double ld = max + log(sum);*/
            return ld;
        }

        // Add VAL to each element
        inline void add(double val)
        {
            m_vals[0] += val;
            m_vals[1] += val;
            m_vals[2] += val;
            m_vals[3] += val;
        }

        // Subtract VAL from each element
        inline void subtract(double val)
        {
            m_vals[0] -= val;
            m_vals[1] -= val;
            m_vals[2] -= val;
            m_vals[3] -= val;
        }

        inline DNADouble& operator+=(const DNADouble& other)
        {
            assert(DNA_ALPHABET::size == 4);
            // Manually unrolled
            m_vals[0] += other.m_vals[0];
            m_vals[1] += other.m_vals[1];
            m_vals[2] += other.m_vals[2];
            m_vals[3] += other.m_vals[3];
            return *this;
        }

        inline DNADouble& operator-=(const DNADouble& other)
        {
            assert(DNA_ALPHABET::size == 4);
            // Manually unrolled
            m_vals[0] -= other.m_vals[0];
            m_vals[1] -= other.m_vals[1];
            m_vals[2] -= other.m_vals[2];
            m_vals[3] -= other.m_vals[3];
            return *this;
        }
        
        
        //
        void printVerbose() const
        {
            std::cout << "AP ";
            for(int i = 0; i < DNA_ALPHABET::size; ++i)
            {
                char b = DNA_ALPHABET::getBase(i);
                std::cout << b << ": " << m_vals[i] << " ";
            }
            std::cout << "\n";
        }

    private:
        double m_vals[DNA_ALPHABET::size];
};

#endif
