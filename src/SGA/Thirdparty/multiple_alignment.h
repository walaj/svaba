//-------------------------------------------------------------------------------
// 
// MultipleAlignment - Class for progressively constructing and managing
// a multiple alignment from a set of pairwise overlaps
//
// Copyright (C) 2011 Jared Simpson (jared.simpson@gmail.com)
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
// ------------------------------------------------------------------------------
#ifndef MULTIPLE_ALIGNMENT_H_
#define MULTIPLE_ALIGNMENT_H_

#include "overlapper.h"
#include <vector>

// Structure to hold a single entry in the multiple alignment
struct MultipleAlignmentElement
{
    // Functions
    MultipleAlignmentElement(const std::string& _name, 
                             const std::string& _sequence, 
                             const std::string& _quality,
                             size_t leading, 
                             size_t trailing);

    // Returns the total number of columns in the multiple alignment
    // All elements in the alignment have the same number of columns
    inline size_t getNumColumns() const
    {
        return leading_columns + padded_sequence.size() + trailing_columns;
    }

    // Returns the symbol for the requested column
    // If the column is in the leading/trailing segment
    // this will return '\0'. Otherwise it will return a base
    // or a gap character
    inline char getColumnSymbol(size_t column_idx) const
    {
        assert(column_idx < getNumColumns());
        if(column_idx < leading_columns || column_idx >= leading_columns + padded_sequence.size()) {
            return '\0';
        } else {
            assert(column_idx - leading_columns < padded_sequence.size());
            return padded_sequence[column_idx - leading_columns];
        }
    }

    // Returns true if there is a base symbol [ACGT] in the given column
    bool hasBaseInColumn(size_t column_idx) const;

    // Returns the quality symbol for the requested column
    // If there is no quality information or the column is in the leading/trailing segment
    // this will return '\0'. Otherwise it will return a quality or a gap character
    inline char getColumnQuality(size_t column_idx) const
    {
        assert(column_idx < getNumColumns());
        if(padded_quality.empty() || 
           column_idx < leading_columns || 
           column_idx >= leading_columns + padded_quality.size()) {
            return '\0';
        } else {
            assert(column_idx - leading_columns < padded_sequence.size());
            return padded_quality[column_idx - leading_columns];
        }
    }

    // Returns the sequence with all padding characters removed
    std::string getUnpaddedSequence() const;

    // Return a substring of this element that can be used to print the full multiple
    // alignment. The leading/trailing columns will be returned as spaces so
    // that all the elements line up.
    std::string getPrintableSubstring(size_t start_column, size_t num_columns) const;

    // Returns the position in the padded string of the base at index idx of
    // the unpadded sequence.
    // Precondition: idx is less than the number of sequence bases
    int getPaddedPositionOfBase(size_t idx) const;

    // Returns a substring of the padded sequence
    // This can have NULL characters if the substring contains leading/trailing columns.
    std::string getPaddedSubstr(size_t start_column, size_t length) const;

    // Returns the column index for the first/last base of the sequence
    size_t getStartColumn() const;
    size_t getEndColumn() const;

    // Insert a new gap before the specified column
    void insertGapBeforeColumn(size_t column_index);
    
    // Extend the length of the trailing columns by n bases
    void extendTrailing(size_t n);

    // Trim the length of the leading column by n bases
    void trimLeading(size_t n);

    // Trim the length of the trailing column by n bases
    void trimTrailing(size_t n);

    // Data
    std::string name;
    std::string padded_sequence;
    std::string padded_quality;

    // The number of columns in the multiple alignment before/after
    // the sequence data
    size_t leading_columns;
    size_t trailing_columns;

};

//
struct SymbolCount
{
    char symbol;
    int count;

    //
    friend std::ostream& operator<<(std::ostream& out, const SymbolCount& a)
    {
        out << a.symbol << ":" << a.count;
        return out;
    }

    // Sorters
    static bool lexicographicOrder(const SymbolCount& a, const SymbolCount& b);
    static bool countOrder(const SymbolCount& a, const SymbolCount& b);
    static bool countOrderDescending(const SymbolCount& a, const SymbolCount& b);
};
typedef std::vector<SymbolCount> SymbolCountVector;

//
class MultipleAlignment
{
    public:
        MultipleAlignment() {}

        // Add the first element to the multiple alignment
        // The quality field is allowed to be empty
        void addBaseSequence(const std::string& name, const std::string& sequence, const std::string& quality);

        // Add a new sequence to the multiple alignment, which overlaps
        // the first sequence added to the multiple alignment (the base sequence).
        // Preconditions are that the base sequence exists and the overlap
        // that is passed in refers to the overlap between the base sequence (as
        // the first set of coordinates) and the incoming sequence. 
        // The quality field is allowed to be empty
        void addOverlap(const std::string& incoming_name,
                        const std::string& incoming_sequence,
                        const std::string& incoming_quality,
                        const SequenceOverlap& reference_incoming_overlap);

        // Add a new sequence to the multiple alignment by extending the last sequence added.
        // This function allows the progressive construction of a multiple alignment for a series
        // of sequences that are laid out into a contig. The SequenceOverlap object must be
        // defined such that the coordinates for the sequence already in the multiple alignment
        // are defined first, and the coordinates for the incoming sequence are defined second.
        // The quality field is allowed to be empty
        void addExtension(const std::string& incoming_name,
                          const std::string& incoming_sequence,
                          const std::string& incoming_quality,
                          const SequenceOverlap& previous_incoming_overlap);

        // Returns an overlap object describing the alignment between the pair of rows
        SequenceOverlap getAlignment(size_t row0, size_t row1) const;

        // Check if the multiple alignment is valid
        // In extension mode after filtering the multiple alignment
        // may not be contiguous. This function checks for this case and
        // return false if so
        bool isValid() const;

        // Returns the index of the first/last column with a symbol of row0 aligned to a symbol of row1.
        // This EXCLUDES leading/trailing gaps. If there is no alignment between the rows, this returns
        // an invalid column index.
        size_t getFirstAlignedBaseColumn(size_t row0, size_t row1) const;
        size_t getLastAlignedBaseColumn(size_t row0, size_t row1) const;

        // Calculate the edit distance between two rows in the multiple alignment
        size_t calculateEditDistanceBetweenRows(size_t row0, size_t row1) const;

        // Calculate a cigar string between two rows
        std::string calculateCigarBetweenRows(size_t row0, size_t row1) const;
    
        // Calculate a new consensus sequence for the base sequence of the multiple alignment
        // A base call is changed only if it has been seen in less than min_call_coverage sequences
        // Leading/trailing bases are trimmed from the consensus sequence if there is less than
        // min_trim_coverage depth at the ends of the base sequence.
        std::string calculateBaseConsensus(int min_call_coverage, int min_trim_coverage);

        // Calculate consensus sequence of the base element that maximizes the likelihood of the multiple alignment
        void calculateBaseConsensusLikelihood(std::string* consensus_sequence, std::string* consensus_quality);
        
        // Calculate the consensus over all columns
        void calculateFullConsensusLikelihood(std::string* consensus_sequence, std::string* consensus_quality);

        // Filter out sequences in the multiple alignment by finding columns with a discrepant
        // consensus base. If two bases in a given column have count more than min_count,
        // then we call a column conflicted and remove the sequences that do not match the base
        // sequence at this column.
        void filterByCount(int min_count);

        // Filter out sequences by the total weight of quality mismatches
        void filterByMismatchQuality(int max_sum_mismatch);

        // Filter out sequences using a vector of bools
        // This is to allow client code to calculate the filtering externally, then
        // apply the vector to the multiple alignment. If vector[i] is TRUE, the sequence is KEPT.
        void filterByVector(const std::vector<bool>& keep_vector);

        // Returns a vector of <symbol,count> pairs for the non-zero symbols of the requested column
        SymbolCountVector getSymbolCountVector(size_t column) const;

        // Returns a formatted string with the number of times each base has been seen for the given column
        std::string getColumnCountString(size_t column) const; 

        // Returns the complete sequence for row r without any padding characters
        std::string getUnpaddedSequence(size_t row) const;

        // Returns a padded substring of a given row 
        std::string getPaddedSubstr(size_t row, size_t start, size_t length) const;

        // Returns the symbol in row r and column c
        // May be the null character if this sequence does not have a base call in this position
        char getSymbol(size_t row, size_t col) const;

        // Calculate the index of the base in the original sequence for a given column
        size_t columnIndexToBaseIndex(size_t row, size_t col) const;

        // Count the length of the homopolymer run in the substring [from, to] inclusive
        // If to is -1, then it is assumed to be until the end of the string
        // if to < from, then the count proceeds backwards, towards the beginning of the string
        size_t countHomopolymer(size_t rowIdx, int from, int to = -1) const;

        // Returns the total number of columns in the multiple alignment.
        // Only valid to call this function if the multiple alignment has been initialized
        size_t getNumColumns() const;

        // Returns the total number of rows (or elements) in the multiple alignment
        size_t getNumRows() const;

        // Return the pileup of bases in the given column
        std::string getPileup(size_t idx) const;

        // Print the alignment to stdout. If the number of columns
        // is greater than max_columns, it will be printed in multiple 
        // segments
        void print(size_t max_columns = 80) const;
    
        // Print a pileup of the base symbol for each column of the alignment
        void printPileup() const;

    private:
     
        // Internal function for performing the addition of a new sequence. Called by
        // addSequenceClipped/addSequenceExtend
        void _addSequence(const std::string& name, 
                          const std::string& sequence, 
                          const std::string& quality, 
                          size_t template_element_index, 
                          const SequenceOverlap& overlap,
                          bool is_extension);
        
        // Calculate the consensus over the range of columns
        void _likelihoodConsensus(size_t start_column, 
                                  size_t end_column, 
                                  std::string* consensus_sequence, 
                                  std::string* consensus_quality);

        // Returns a string with the most frequent base for each column
        // including padding symbols
        std::string getPaddedConsensus() const;

        // Insert a new gap into all sequences in the multiple alignment
        // before the given column
        void insertGapBeforeColumn(size_t column_index);

        // After removing reads from the multiple alignment there
        // may be empty leading or trailing columns. This function
        // removes these
        void trimEmptyColumns();

        // Expand a cigar string by having one symbol per event instead
        // of run length encoding
        std::string expandCigar(const std::string& cigar);

        // Returns the index of a symbol over the alphabet "ACGTN-"
        int symbol2index(char symbol) const;
        
        // Return the probability that a base call is incorrect, given a Phred-33 ascii character
        double quality2prob(char q) const;
        
        // Return a phred+33 scaled quality character for this probability
        char prob2quality(double p) const;

        // Return a vector of counts for each base call for the given column
        std::vector<int> getColumnBaseCounts(size_t idx) const;

        // Data
        std::vector<MultipleAlignmentElement> m_sequences;
        static const size_t m_alphabet_size = 6;
        static const char* m_alphabet;
};

#endif  // MULTIPLE_ALIGNMENT_H_
