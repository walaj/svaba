//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Pileup - An array of strings containing all the seen
// bases for a given region/read
//
#include "Pileup.h"
#include "Alphabet.h"
#include <math.h>
#include <iostream>

// Calculate the consensus base at this position
// using a simple model where all bases are treated equally
char Pileup::calculateSimpleConsensus() const
{
    assert(!m_data.empty());
    AlphaCount64 ac;
    for(size_t i = 0; i < m_data.size(); ++i)
    {
        ac.increment(m_data[i].base);
    }
    return ac.getMaxBase();
}

// Returns the number of times each base appears in the 
// pileup string
AlphaCount64 Pileup::getAlphaCount() const
{
    AlphaCount64 ac;
    for(size_t i = 0; i < m_data.size(); ++i)
    {
        ac.increment(m_data[i].base);
    }
    return ac;
}

// Calculate the consensus base at this position
// using a simple model where all bases are treated equally
DNADouble Pileup::calculateSimpleAlphaProb() const
{
    DNADouble ap;
    assert(!m_data.empty());
    WARN_ONCE("Fix Pileup::calculateSimpleDNADouble numerical stability");
    assert(false);
    for(size_t i = 0; i < ap.getAlphabetSize(); ++i)
    {
        // Calculate the posterior probability of the data given that
        // b is the true base
        char b = ALPHABET[i];
        double posterior = 0.0f;
        for(size_t i = 0; i < m_data.size(); ++i)
        {
            if(m_data[i].base == b)
                posterior += log(1 - exp(m_data[1].lp));
            else
                posterior += m_data[1].lp;
        }
        ap.set(b, posterior);
    }

    // Calculate the marginal probabilty of the data
    double marginal = 0.0f;
    for(size_t i = 0; i < DNA_ALPHABET_SIZE; ++i)
    {
        char b = ALPHABET[i];
        marginal += exp(ap.get(b));
    }
    

    // Scale the posterior probabilites by the marginal
    marginal = log(marginal);
    for(size_t i = 0; i < ap.getAlphabetSize(); ++i)
    {
        char b = DNA_ALPHABET::getBase(i);
        double lp = ap.get(b);
        //std::cout << "Marginal: " << marginal << " lp: " << lp << " scaled: " << lp - marginal << "\n";
        ap.set(b, lp - marginal);
    }
    return ap;    
}

// Calculate the likelihood of the data given the base call is 
// {A,C,G,T}. 
DNADouble Pileup::calculateLikelihoodNoQuality(double p_error) const
{
    DNADouble ap;
    assert(!m_data.empty());
    double p_correct = 1.0 - p_error;
    
    double log_error = log(p_error);
    double log_correct = log(p_correct);

    for(size_t i = 0; i < ap.getAlphabetSize(); ++i)
    {
        // Calculate the likelihood of the data given b is the true base
        char b = ap.getBase(i);
        double likelihood = 0.0f;
        for(size_t i = 0; i < m_data.size(); ++i)
        {
            if(m_data[i].base == b)
                likelihood += log_correct;
            else
                likelihood += log_error;
        }
        ap.set(b, likelihood);
    }
    return ap;    
}


//
char Pileup::getCount(char base) const
{
    AlphaCount64 ac;
    for(size_t i = 0; i < m_data.size(); ++i)
    {
        ac.increment(m_data[i].base);
    }
    return ac.get(base);
}

//
char Pileup::getBase(size_t idx) const
{
    assert(idx < m_data.size());
    return m_data[idx].base;
}

//
size_t Pileup::getDepth() const
{
    return m_data.size();
}

//
std::string Pileup::toStr() const
{
    std::string out;
    out.reserve(m_data.size());
    for(size_t i = 0; i < m_data.size(); ++i)
        out.append(1, m_data[i].base);
    return out;
}
