//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// KmerDistribution - Histogram of kmer frequencies
//
#ifndef KMERDISTRIBUTION_H
#define KMERDISTRIBUTION_H

#include <vector>
#include <map>
#include <cstddef>
#include <stdio.h>

class KmerDistribution
{
    public:
    
        KmerDistribution();
        ~KmerDistribution() {}
        
        // Returns the proportion of the distribution less than or equal to n
        double getCumulativeProportionLEQ(int n) const;

        // Returns the smallest value n such that the proportion of the data less
        // than n is p. This is the inverse of getCumulativeProportionLEQ.
        size_t getCutoffForProportion(double p) const;

        // Returns the predicted boundary of the left tail which holds erroneous kmers
        int findErrorBoundary() const;
        int findErrorBoundaryByRatio(double ratio) const;

        // Get the mode of the distribution if the first n values are ignored
        int getCensoredMode(size_t n) const;
 
        size_t getTotalKmers() const;
        size_t getNumberWithCount(size_t c) const;

        //
        std::vector<int> toCountVector(int max) const;

        //
        int findFirstLocalMinimum() const;
        void add(int count);
        void print(int max) const; 
        void print(FILE* file, int max) const; 

    private:

        // int -> int map of the number of times a kmer with multiplicty N has been seen
        std::map<int,int> m_data;
};

#endif
