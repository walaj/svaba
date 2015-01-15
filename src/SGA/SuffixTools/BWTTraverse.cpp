//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// BWTTraverse - traverse a bwt, outputting all
// strings that are length greater than some threshold
//
#include "BWTTraverse.h"

void BWTTraverse::extract(const BWT* pBWT, unsigned int len)
{
    // Keep a stack of the elements to visit and 
    // a string which is the reverse string corresponding
    // to the current stack elements
    TraverseStack stack;
    std::string rev_str;
    size_t iter_count = 0;
    size_t output_count = 0;

    for(size_t i = 0; i < DNA_ALPHABET_SIZE; ++i)
    {
        char base_char = ALPHABET[i];

        // Initialize the stack and rev_str
        BWTInterval range;
        BWTAlgorithms::initInterval(range, base_char, pBWT);
        AlphaCount ext = BWTAlgorithms::getExtCount(range, pBWT);
        stack.push(TraverseElem(range, ext));
        rev_str.append(1, base_char);

        while(!stack.empty())
        {
            ++iter_count;
            TraverseElem& curr = stack.top();
            curr.goNext();
            bool doPop = !curr.isValid(); 

            // If rev_str is long enough output it and stop the traversal at this depth
            if(rev_str.length() == len)
            {
                if(output_count % 1000000 == 0)
                    std::cerr << "output: " << output_count << "\n";

                doPop = true;
                size_t multiplicity = curr.getRange().size();
                printf(">%zu %u %zu\n%s\n", output_count, len, multiplicity, reverse(rev_str).c_str());
                ++output_count;
            }

            if(doPop)
            {
                // Remove one char from the end of the string 
                // and one element from the stack
                rev_str.erase(rev_str.length() - 1, 1);
                stack.pop();
            }
            else
            {
                // Traverse one level deeper 
                char b = curr.getCurrChar();
                BWTInterval subrange = curr.getRange();
                BWTAlgorithms::updateInterval(subrange, b, pBWT);
                AlphaCount subext = BWTAlgorithms::getExtCount(subrange, pBWT);
                stack.push(TraverseElem(subrange, subext));
                rev_str.append(1, b);
            }
        }
    }

    std::cerr << "Num output: " << output_count << " num iterations: " << iter_count << " count/len: " << (double)iter_count / len << "\n";
}

void BWTTraverse::extractSG(const BWT* pBWT, const BWT* pRevBWT, const unsigned int len)
{
    WARN_ONCE("Skipping bwt[0]");
    // Keep a bool vector marking which entries in pBWT have been visited
    size_t n_elems = pBWT->getBWLen();
    bool_vec visited(n_elems, false);

    size_t count = 0;
    size_t currIdx = 1;
    while(currIdx < n_elems)
    {
        // Invariant: currIdx is the index into pBWT that is the lowest index that has not been visited
        // Left-extend the entry at currIdx into a string of length len
        
        std::string str;
        str.reserve(len);

        // Get the first character of the suffix starting at this position
        char first = pBWT->getF(currIdx);
        str += first;
        BWTInterval range(currIdx, currIdx);

        while(str.length() < len)
        {
            char b = pBWT->getChar(range.lower);
            if(b == '$')
                break;
            str += b;
            BWTAlgorithms::updateInterval(range, b, pBWT);
        }

        // The string was built backwards, reverse it
        str = reverse(str);
        if(str.length() < len)
        {
            visited[currIdx] = true;
        }
        else if(!visited[range.lower])
        {
            markVisited(str, visited, pBWT);
            //std::cout << currIdx << " interval: " << range << " string: " << str << "\n";
            extendLeft(len, str, visited, pBWT, pRevBWT);
            extendRight(len, str, visited, pBWT, pRevBWT, false);
            printf(">%zu %zu 0\n%s\n", count++, str.length(), str.c_str());
        }

        currIdx++;
    }
}


// Extend the sequence by finding a consensus right-ward extension base
void BWTTraverse::extendRight(const unsigned int len, std::string& str, bool_vec& visited, const BWT* pBWT, const BWT* pRevBWT, bool isReverse)
{
    // Extract the last len characters
    unsigned int overlapLen = len - 1;

    // Initialize the range for str
    BWTIntervalPair ranges = BWTAlgorithms::findIntervalPair(pBWT, pRevBWT, str);

    bool done = false;
    while(!done)
    {
        AlphaCount ext_counts = BWTAlgorithms::getExtCount(ranges.interval[1], pRevBWT);

        if(ext_counts.hasUniqueDNAChar())
        {
            char b = ext_counts.getUniqueDNAChar();

            // There is a unique extension from the seed sequence to B
            // Ensure that the reverse is true and that the sequence we are extending to, extends to this one
            std::string joined = str + b;
            std::string back_search = suffix(joined, len);

            // Get the count of bases back to the ending sequence of seed
            // We do this in the opposite direction of extension
            AlphaCount back_count = BWTAlgorithms::calculateExactExtensions(overlapLen, reverse(back_search), pRevBWT, pBWT);

            if(back_count.hasUniqueDNAChar())
            {
                // Assert back is the character we are expecting
                assert(back_count.getUniqueDNAChar() == str[str.length() - len]);
                str = joined;

                // mark the newly visited l-mers
                if(!isReverse)
                {
                    markVisited(back_search, visited, pBWT);
                }
                else
                {
                    markVisited(reverse(back_search), visited, pRevBWT);
                }    
                BWTAlgorithms::updateBothR(ranges, b, pRevBWT);
            }
            else
            {
                done = true;
            }
        }
        else
        {
            done = true;
        }
    }
}

// Extend the string to the left by reversing it and calling extendRight
void BWTTraverse::extendLeft(const unsigned int len, std::string& str, bool_vec& visited, const BWT* pBWT, const BWT* pRevBWT)
{
    str = reverse(str);
    extendRight(len, str, visited, pRevBWT, pBWT, true);
    str = reverse(str);
}

void BWTTraverse::markVisited(const std::string& str, bool_vec& visited, const BWT* pBWT)
{
    BWTInterval range = BWTAlgorithms::findInterval(pBWT, str);
    for(int64_t i = range.lower; i <= range.upper; ++i)
    {
        assert(i >= 0 && i < (int64_t)visited.size());
        visited[i] = true;
    }
}
