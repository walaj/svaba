//-------------------------------------------------------------------------------
// 
// MultipleAlignment - Class for constructing and managing a multiple alignment
// constructed from a set of pairwise overlaps
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
#include "multiple_alignment.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <assert.h>
#include <stdio.h>
#include <limits>
#include <algorithm>

//#define MA_DEBUG 1
//#define MA_DEBUG_CONSENSUS 1

// Initialize static members
const char* MultipleAlignment::m_alphabet = "ACGTN-";

//
// MultipleAlignmentElement
//

//
MultipleAlignmentElement::MultipleAlignmentElement(const std::string& _name, 
                                                   const std::string& _sequence,
                                                   const std::string& _quality,
                                                   size_t leading,
                                                   size_t trailing) : name(_name), 
                                                                      padded_sequence(_sequence),
                                                                      padded_quality(_quality),
                                                                      leading_columns(leading),
                                                                      trailing_columns(trailing)
{
    assert(padded_quality.empty() || padded_quality.size() == padded_sequence.size());
}

//
size_t MultipleAlignmentElement::getStartColumn() const
{
    return leading_columns;
}

//
size_t MultipleAlignmentElement::getEndColumn() const
{
    return getNumColumns() - trailing_columns - 1;
}
    
bool MultipleAlignmentElement::hasBaseInColumn(size_t column_idx) const
{
    char s = getColumnSymbol(column_idx);
    return s != '\0' && s != '-';
}

//
std::string MultipleAlignmentElement::getPaddedSubstr(size_t start_column, size_t length) const
{
    std::string out;
    out.reserve(length);

    for(size_t i = start_column; i < start_column + length; ++i)
        out.push_back(getColumnSymbol(i));

    assert(out.size() == length);
    return out;
}

//
int MultipleAlignmentElement::getPaddedPositionOfBase(size_t idx) const
{
    size_t unpadded_count = 0;
    for(size_t i = 0; i < padded_sequence.size(); ++i) {
        if(padded_sequence[i] != '-') {
            if(unpadded_count == idx)
                return i;
            else
                unpadded_count += 1;
        }
    }
    std::cerr << "Base index out of bounds: " << idx << "\n";
    assert(false);
    return -1;
}

//
void MultipleAlignmentElement::insertGapBeforeColumn(size_t column_index)
{
    // Check if the column to insert the gap falls within the leading columns
    // If so, just increase the number of leading columns to account for the inserted base.
    // If the column index is one greater than the offset, then we want to 
    // insert a gap before the first base. This is equivalent to just
    // extending the offset.
    size_t first_sequence_column_index = leading_columns;
    if(column_index <= first_sequence_column_index) {
        leading_columns += 1;
    }
    else {
        assert(column_index > leading_columns);
        size_t insert_position = column_index - leading_columns;
        if(insert_position < padded_sequence.size()) {
            padded_sequence.insert(insert_position, 1, '-');
            if(!padded_quality.empty())
                padded_quality.insert(insert_position, 1, '-');
        }
        else
            trailing_columns += 1;
    }
}

//
void MultipleAlignmentElement::extendTrailing(size_t n)
{
    trailing_columns += n;
}

//
void MultipleAlignmentElement::trimLeading(size_t n)
{
    assert(n <= leading_columns);
    leading_columns -= n;
}

//
void MultipleAlignmentElement::trimTrailing(size_t n)
{
    assert(n <= trailing_columns);
    trailing_columns -= n;
}

//
std::string MultipleAlignmentElement::getUnpaddedSequence() const
{
    std::string out;
    for(size_t i = 0; i < padded_sequence.size(); ++i) {
        if(padded_sequence[i] != '-')
            out.push_back(padded_sequence[i]);
    }
    return out;
}

//
std::string MultipleAlignmentElement::getPrintableSubstring(size_t start_column, 
                                                            size_t num_columns) const
{
    std::string out;
    size_t i = 0;
    while(i < num_columns) {
        char symbol = getColumnSymbol(start_column + i);
        out.push_back(symbol != '\0' ? symbol : ' ');
        ++i;
    }
    return out;
}

//
// SymbolCount
//

//
bool SymbolCount::lexicographicOrder(const SymbolCount& a, const SymbolCount& b) 
{ 
    return a.symbol < b.symbol; 
}

//
bool SymbolCount::countOrder(const SymbolCount& a, const SymbolCount& b) 
{ 
    return a.count < b.count; 
}

//
bool SymbolCount::countOrderDescending(const SymbolCount& a, const SymbolCount& b) 
{ 
    return a.count > b.count; 
}

//
// MultipleAlignment
//

//
void MultipleAlignment::addBaseSequence(const std::string& name, 
                                        const std::string& sequence, 
                                        const std::string& quality)
{
    m_sequences.push_back(MultipleAlignmentElement(name, sequence, quality, 0, 0));
}

// See header
void MultipleAlignment::addOverlap(const std::string& incoming_name,
                                   const std::string& incoming_sequence,
                                   const std::string& incoming_quality,
                                   const SequenceOverlap& reference_incoming_overlap)
{
    // This function cannot be called before a base element has been added
    assert(!m_sequences.empty());
    _addSequence(incoming_name, incoming_sequence, incoming_quality, 0, reference_incoming_overlap, false);
}

//
void MultipleAlignment::addExtension(const std::string& incoming_name,
                                     const std::string& incoming_sequence,
                                     const std::string& incoming_quality,
                                     const SequenceOverlap& previous_incoming_overlap)
{
    // This function cannot be called before a base element has been added
    assert(!m_sequences.empty());
    _addSequence(incoming_name, incoming_sequence, incoming_quality, m_sequences.size() - 1, previous_incoming_overlap, true);
}

// Adds a new string into the multiple alignment using the overlap
// between the incoming sequence and an existing sequence in the
// multiple alignment to calculate the new padded string
void MultipleAlignment::_addSequence(const std::string& name,
                                     const std::string& sequence,
                                     const std::string& quality, 
                                     size_t template_element_index, 
                                     const SequenceOverlap& overlap,
                                     bool is_extension)
{
    MultipleAlignmentElement* template_element = &m_sequences[template_element_index];

    // Get the padded sequence for the template element
    const std::string& template_padded = template_element->padded_sequence;

    // The output padded sequence for the incoming
    std::string padded_output;
    std::string padded_quality;
    
    // Sanity checks
    assert(quality.empty() || quality.size() == sequence.size());

    // Iterate over the cigar string and the padded sequence of the template element
    // to determine where to insert gap symbols
    size_t cigar_index = 0;
    size_t template_index = template_element->getPaddedPositionOfBase(overlap.match[0].start);
    size_t incoming_index = overlap.match[1].start;

    // If extending the multiple alignment, the first
    // base of the incoming sequence must be aligned.
    if(is_extension) 
        assert(incoming_index == 0);
    
    size_t template_leading = template_element->leading_columns;
    size_t incoming_leading = template_index + template_leading;

    // Expand the cigar for easier parsing
    std::string expanded_cigar = expandCigar(overlap.cigar);
    assert(!expanded_cigar.empty());
    assert(template_index < template_padded.size());
    assert(template_padded[template_index] != '-');

#ifdef MA_DEBUG
    std::cout << "Cigar: " << expanded_cigar << "\n";
    std::cout << "template: " << template_padded.substr(template_index) << "\n";
    std::cout << "incoming: " << sequence.substr(incoming_index) << "\n";
    std::cout << "Template leading: " << template_leading << "\n";
    std::cout << "Incoming leading: " << incoming_leading << "\n";
    std::cout << "Pairwise:\n";
    overlap.printAlignment(template_element->getUnpaddedSequence(), sequence);
#endif

    while(cigar_index < expanded_cigar.size()) {
        // Check if we are in an existing template gap. This must be handled
        // seperately
        bool in_template_gap = template_padded[template_index] == '-';
        if(in_template_gap) {
            // If we are in a incoming sequence insertion
            // (cigar D) then we are adding a base into a known
            // gap. Add the current incoming base to the output
            if(expanded_cigar[cigar_index] == 'I') {
                padded_output.push_back(sequence[incoming_index]);
                if(!quality.empty())
                    padded_quality.push_back(quality[incoming_index]);

                incoming_index += 1;
                cigar_index += 1;
                template_index += 1;
            } else { 
                // This is an insertion that is in some other sequence
                // in the multiple alignment. Add a gap to the padded output
                padded_output.push_back('-');
                if(!quality.empty())
                    padded_quality.push_back('-');
                
                // Increment the template index
                template_index += 1;    
            }
        } else {
            // Not a template gap
            switch(expanded_cigar[cigar_index]) {
                case 'M':
                    padded_output.push_back(sequence[incoming_index]);
                    if(!quality.empty())
                        padded_quality.push_back(quality[incoming_index]);

                    incoming_index += 1;
                    template_index += 1;
                    cigar_index += 1;
                    break;
                case 'I':
                    insertGapBeforeColumn(template_index + template_leading);
                    padded_output.push_back(sequence[incoming_index]);
                    if(!quality.empty())
                        padded_quality.push_back(quality[incoming_index]);

                    incoming_index += 1;
                    cigar_index += 1;
                    template_index += 1; // skip the newly introduced gap
                    break;
                case 'D':
                    padded_output.push_back('-');
                    if(!quality.empty())
                        padded_quality.push_back('-');
                        
                    cigar_index += 1;
                    template_index += 1;
                    break;
                case 'S':
                    cigar_index += 1;
                    break;
                default:
                    std::cerr << "Error: unhandled cigar symbol " << expanded_cigar[cigar_index] << "\n";
                    exit(EXIT_FAILURE);
                    break;
            }
        }
    }

    // Now that the alignment has been built, add the remaining bases of the incoming sequence
    // All other elements of the multiple alignment will be updated to increase the number of
    // trailing columns
    if(is_extension) {
        //printf("Extending the alignment %zu bases\n", sequence.size() - incoming_index);
        padded_output.append(sequence.substr(incoming_index));
        padded_quality.append(quality.substr(incoming_index));

        // Ensure that the incoming sequence ends at least as far as the last base
        // currently in the multiple alignment. Otherwise this sequence is a containment
        // and we can not deal with it.
        size_t incoming_columns = padded_output.size() + incoming_leading;
        assert(incoming_columns >= m_sequences.front().getNumColumns());

        // Extend all other sequences to have the same number of columns as the incoming sequence
        for(size_t i = 0; i < m_sequences.size(); ++i)
            m_sequences[i].extendTrailing(sequence.size() - incoming_index);

    }

    // Calculate the number of unfilled columns of the multiple alignment that come after
    // the padded sequence
    size_t incoming_trailing = template_element->getNumColumns() - padded_output.size() - incoming_leading;

    if(is_extension)
        assert(incoming_trailing == 0);

    MultipleAlignmentElement incoming_element(name, padded_output, padded_quality, 
                                              incoming_leading, incoming_trailing);

    m_sequences.push_back(incoming_element);

    /*
    std::string calculated_cigar = calculateExpandedCigarBetweenRows(template_element_index, m_sequences.size() - 1);
    printf("Incoming cigar: %s\n", expanded_cigar.c_str());
    printf("Calculated cigar: %s\n", calculated_cigar.c_str());
    */
}

SequenceOverlap MultipleAlignment::getAlignment(size_t row0, size_t row1) const
{
    size_t start_column = getFirstAlignedBaseColumn(row0, row1);
    size_t end_column = getFirstAlignedBaseColumn(row0, row1);

    SequenceOverlap overlap;
    overlap.match[0].start = columnIndexToBaseIndex(row0, start_column);
    overlap.match[0].end = columnIndexToBaseIndex(row0, end_column);
    overlap.length[0] = getUnpaddedSequence(row0).length();

    overlap.match[1].start = columnIndexToBaseIndex(row1, start_column);
    overlap.match[1].end = columnIndexToBaseIndex(row1, end_column);
    overlap.length[1] = getUnpaddedSequence(row1).length();

    overlap.cigar = calculateCigarBetweenRows(row0, row1);
    overlap.total_columns = overlap.calculateTotalColumns();
    
    // These fields are unavailable
    overlap.score = -1;
    overlap.edit_distance = -1;

    return overlap;
}


bool MultipleAlignment::isValid() const
{
    size_t total_columns = getNumColumns();
    size_t i = 0;
    while(i < total_columns) {
        if(getPileup(i).empty())
            return false;
        ++i;
    }
    return true;
}

size_t MultipleAlignment::getFirstAlignedBaseColumn(size_t row0, size_t row1) const
{
    size_t num_cols = getNumColumns();
    size_t i = 0;
    while(i < num_cols)
    {
        bool base0 = m_sequences[row0].hasBaseInColumn(i);
        bool base1 = m_sequences[row1].hasBaseInColumn(i);
        if(base0 && base1)
            return i;
        ++i;
    }
    return num_cols; // No alignment between these rows
}

size_t MultipleAlignment::getLastAlignedBaseColumn(size_t row0, size_t row1) const
{
    size_t num_cols = getNumColumns();
    int i = num_cols - 1;
    while(i >= 0)
    {
        bool base0 = m_sequences[row0].hasBaseInColumn(i);
        bool base1 = m_sequences[row1].hasBaseInColumn(i);
        if(base0 && base1)
            return i;
        --i;
    }

    return num_cols; // No alignment between these rows
}

size_t MultipleAlignment::calculateEditDistanceBetweenRows(size_t row0, size_t row1) const
{
    size_t num_columns = getNumColumns();
    assert(row0 < m_sequences.size());
    assert(row1 < m_sequences.size());
    size_t edit_distance = 0;

    for(size_t i = 0; i < num_columns; ++i)
    {
        char row0_symbol = getSymbol(row0, i);
        char row1_symbol = getSymbol(row1, i);

        if(row0_symbol == '\0' || row1_symbol == '\0')
            continue;

        if(row0_symbol != row1_symbol)
            edit_distance += 1;
    }
    return edit_distance;
}

std::string MultipleAlignment::calculateCigarBetweenRows(size_t row0, size_t row1) const
{
    assert(row0 < m_sequences.size());
    assert(row1 < m_sequences.size());
    std::string cigar;
    
    size_t start_column = getFirstAlignedBaseColumn(row0, row1);
    size_t end_column = getLastAlignedBaseColumn(row0, row1);

    for(size_t i = start_column; i <= end_column; ++i)
    {
        char row0_symbol = getSymbol(row0, i);
        char row1_symbol = getSymbol(row1, i);

        // Since we calculated the first/last column above, we can
        // assume every base within the alignment is a valid alignment symbol
        assert(row0_symbol != '\0');
        assert(row1_symbol != '\0');

        // Both in a gap, skip
        if(row0_symbol == '-' && row1_symbol == '-')
            continue;

        if(row0_symbol == '-')
            cigar.push_back('I');
        else if(row1_symbol == '-')
            cigar.push_back('D');
        else
            cigar.push_back('M');
    }
    return Overlapper::compactCigar(cigar);
}

std::string MultipleAlignment::calculateBaseConsensus(int min_call_coverage, int min_trim_coverage)
{
    assert(!m_sequences.empty());
    std::string consensus_sequence;
    MultipleAlignmentElement& base_element = m_sequences.front();
    size_t start_column = base_element.getStartColumn();
    size_t end_column = base_element.getEndColumn();
    
    // This index records the last base in the consensus that had coverage greater than
    // min_trim_coverage. After the consensus calculation the read is trimmed back to this position
    int last_good_base = -1;

    for(size_t c = start_column; c <= end_column; ++c) {
        std::vector<int> counts = getColumnBaseCounts(c);

        char max_symbol = '\0';
        int max_count = -1;
        int total_depth = 0;

#ifdef MA_DEBUG_CONSENSUS
        printf("%zu\t", c);
#endif
        for(size_t a = 0; a < m_alphabet_size; ++a) {
            char symbol = m_alphabet[a];
            total_depth += counts[a];

            if(symbol != 'N' && counts[a] > max_count) {
                max_symbol = symbol;
                max_count = counts[a];
            }
#ifdef MA_DEBUG_CONSENSUS
            printf("%c:%d ", symbol, counts[a]);
#endif
        }

        char base_symbol = base_element.getColumnSymbol(c);
        int base_count = counts[symbol2index(base_symbol)];
        
        // Choose a consensus base for this column. Only change a base
        // if has been seen less than min_call_coverage times and the max
        // base in the column has been seen more times than the base symbol
        char consensus_symbol;
        if(max_count > base_count && base_count < min_call_coverage)
            consensus_symbol = max_symbol;
        else
            consensus_symbol = base_symbol;
        
        // Output a symbol to the consensus. Skip padding symbols and leading
        // bases that are less than the minimum required depth to avoid trimming
        if(consensus_symbol != '-' &&
            (!consensus_sequence.empty() || total_depth >= min_trim_coverage))
                consensus_sequence.push_back(consensus_symbol);

        // Record the position of the last good base
        if(total_depth >= min_trim_coverage) {
            int consensus_index = consensus_sequence.size() - 1;
            if(consensus_index > last_good_base)
                last_good_base = consensus_index;
        }
#ifdef MA_DEBUG_CONSENSUS
        printf("CALL: %c\n", consensus_symbol);
#endif 
    }

    if(last_good_base != -1)
        consensus_sequence.erase(last_good_base + 1);
    else
        consensus_sequence.clear();

    return consensus_sequence;
}

//
void MultipleAlignment::calculateBaseConsensusLikelihood(std::string* consensus_sequence, std::string* consensus_quality)
{
    MultipleAlignmentElement& base_element = m_sequences.front();
    size_t start_column = base_element.getStartColumn();
    size_t end_column = base_element.getEndColumn();
    _likelihoodConsensus(start_column, end_column, consensus_sequence, consensus_quality);
}

//
void MultipleAlignment::calculateFullConsensusLikelihood(std::string* consensus_sequence, std::string* consensus_quality)
{
    _likelihoodConsensus(0, getNumColumns() - 1, consensus_sequence, consensus_quality);
}


void MultipleAlignment::_likelihoodConsensus(size_t start_column, 
                                             size_t end_column, 
                                             std::string* consensus_sequence, 
                                             std::string* consensus_quality)
{
#ifdef MA_DEBUG_CONSENSUS
    print();
#endif
    assert(!m_sequences.empty());

    consensus_sequence->clear();
    consensus_quality->clear();

    // Probability of emitting a gap symbol when there is a true
    // nucleotide in the sequence.
    const double PROBABILITY_GAP = 0.0001;

    for(size_t c = start_column; c <= end_column; ++c) {
        std::vector<int> counts = getColumnBaseCounts(c);

        std::vector<double> likelihoods(m_alphabet_size, 0.0f);

#ifdef MA_DEBUG_CONSENSUS
        std::string pileup;
        std::string quality_pileup;
#endif

        // Update likelihoods using every sequence with a base aligned in this column.
        for(size_t i = 0; i < m_sequences.size(); ++i) {
            char b = m_sequences[i].getColumnSymbol(c);
            char q = m_sequences[i].getColumnQuality(c);
            
            if(b == '\0')
                continue; // no base at this position, skip
            assert(q != '\0');

            // Calculate likelihoods.
            double p_error;
            if(b != '-')
                p_error = quality2prob(q);
            else
                p_error = PROBABILITY_GAP;

            // The probabilty the base call is correct.
            double p_match = 1.0f - p_error;
            double lp_match = log(p_match);

            // The true underlying base is different than base call. 
            // We don't have quality scores for all possible calls so we use p_error / 3
            double p_mismatch = p_error / 3;
            double lp_mismatch = log(p_mismatch);

            // Update likelihoods
            for(size_t a = 0; a < m_alphabet_size; ++a) {
                char a_symbol = m_alphabet[a];
                if(a_symbol == b)
                    likelihoods[a] += lp_match;
                else
                    likelihoods[a] += lp_mismatch;
            }

#ifdef MA_DEBUG_CONSENSUS
            pileup.push_back(b);
            quality_pileup.push_back(q);
#endif
        }

        // Calculate most likely symbol
        char call_symbol = 'N';
        double max_likelihood = -std::numeric_limits<double>::max();
        double sum_probability = 0.0f;

        for(size_t a = 0; a < m_alphabet_size; ++a) {
            if(likelihoods[a] > max_likelihood) {
                call_symbol = m_alphabet[a];
                max_likelihood = likelihoods[a];
            }
            sum_probability += exp(likelihoods[a]);
        }

        double p_consensus_error = 1.0f - (exp(max_likelihood) / sum_probability);
        char quality_symbol = prob2quality(p_consensus_error);

#ifdef MA_DEBUG_CONSENSUS
        // Print stats
        printf("Likelihoods:");
        for(size_t a = 0; a < m_alphabet_size; ++a) {
            char a_symbol = m_alphabet[a];
            printf("%c=%.2lf ", a_symbol, likelihoods[a]);
        }
        printf(" Call: %c Perror: %.8lf Qscore: %c Pileup: %s Quality: %s\n", call_symbol, p_consensus_error, quality_symbol, pileup.c_str(), quality_pileup.c_str());
#endif
        // Output a symbol to the consensus. Skip padding symbols and leading
        // bases that are less than the minimum required depth to avoid trimming
        if(call_symbol != '-') {
            consensus_sequence->push_back(call_symbol);
            consensus_quality->push_back(quality_symbol);
        }
    }    
}

//
void MultipleAlignment::filterByMismatchQuality(int max_sum_mismatch)
{

    assert(!m_sequences.empty());
    MultipleAlignmentElement& base_element = m_sequences.front();
    size_t start_column = base_element.getStartColumn();
    size_t end_column = base_element.getEndColumn();

    int GAP_PENALTY = 30;

    // A vector to record the total phred scores of mismatching bases
    std::vector<int> scores(m_sequences.size(), 0);

    for(size_t c = start_column; c <= end_column; ++c) {
        char base_symbol = base_element.getColumnSymbol(c);
        char base_quality = base_element.getColumnQuality(c);

        for(size_t i = 1; i < m_sequences.size(); ++i) {
            char symbol = m_sequences[i].getColumnSymbol(c);
            char quality = m_sequences[i].getColumnQuality(c);

            if(symbol != '\0' && symbol != base_symbol) {
                // This is a mismatch
                if(symbol == '-' || base_symbol == '-')
                    scores[i] += GAP_PENALTY;
                else
                    scores[i] += std::min(int(quality) - 33, int(base_quality) - 33);
            }
        }
    }

    // A vector to record which sequences pass the filter.
    // keep_vector[i] == 1 means that sequence i should be kept.
    std::vector<bool> keep_vector(m_sequences.size(), 1);

// print(400);
    for(size_t i = 1; i < m_sequences.size(); ++i) {
        if(scores[i] > max_sum_mismatch)
            keep_vector[i] = 0;
    }

    // Erase elements from the vector
    filterByVector(keep_vector);
}

//
void MultipleAlignment::filterByCount(int min_count)
{
    assert(!m_sequences.empty());
    MultipleAlignmentElement& base_element = m_sequences.front();
    size_t start_column = base_element.getStartColumn();
    size_t end_column = base_element.getEndColumn();

    // A vector to record which sequences pass the filter.
    // keep_vector[i] == 1 means that sequence i should be kept.
    std::vector<bool> keep_vector(m_sequences.size(), 1);

    for(size_t c = start_column; c <= end_column; ++c) {
        std::vector<int> counts = getColumnBaseCounts(c);
        char base_symbol = base_element.getColumnSymbol(c);

        // Check that the base sequence has a call in this column
        // If not, we do not filter here
        if(base_symbol == '\0')
            continue;

        int base_count = counts[symbol2index(base_symbol)];
        int num_above_min = 0;

        for(size_t a = 0; a < m_alphabet_size; ++a) {
            if(counts[a] >= min_count)
                num_above_min += 1;
        }

        if(num_above_min > 1 && base_count >= min_count) {
            //printf("Column %zu is conflicted: %s\n", c, getColumnCountString(c).c_str());
            // This column is conflicted. Update keep_vector.
            for(size_t i = 1; i < m_sequences.size(); ++i) {
                char symbol = m_sequences[i].getColumnSymbol(c);

                // Remove this sequence if it overlaps the conflicted column,
                // if it does not match the base sequence, and its base count is greater
                // than the threshold
                if(symbol != '\0' && symbol != base_symbol && counts[symbol2index(symbol)] >= min_count)
                    keep_vector[i] = 0;
            }
        }
    }

    // Erase elements from the vector
    filterByVector(keep_vector);
}

//
void MultipleAlignment::filterByVector(const std::vector<bool>& keep_vector)
{
    std::vector<MultipleAlignmentElement> filtered_sequences;
    for(size_t i = 0; i < m_sequences.size(); ++i) {
        if(keep_vector[i])
            filtered_sequences.push_back(m_sequences[i]);
    }

    m_sequences.swap(filtered_sequences);
    trimEmptyColumns();
}

std::string MultipleAlignment::getUnpaddedSequence(size_t row) const
{
    assert(row < m_sequences.size());
    return m_sequences[row].getUnpaddedSequence();
    
}

//
char MultipleAlignment::getSymbol(size_t row, size_t col) const
{
    assert(row < m_sequences.size());
    return m_sequences[row].getColumnSymbol(col);
}

//
size_t MultipleAlignment::getNumColumns() const
{
    assert(!m_sequences.empty());
    return m_sequences.front().getNumColumns();
}

//
size_t MultipleAlignment::getNumRows() const
{
    return m_sequences.size();
}

// Calculate the index of the base in the original sequence for a given column
size_t MultipleAlignment::columnIndexToBaseIndex(size_t row, size_t col) const 
{ 
    assert(row < m_sequences.size());
    assert(col < getNumColumns());

    // It is only valid to call this function for columns where
    // this row has an actual sequence symbol
    if(!m_sequences[row].hasBaseInColumn(col)) {
        std::cerr << "Error: MultipleAlignment::columnToBaseIndex called on an empty column\n";
        assert(false);
    }

    size_t base_index = 0;
    for(size_t i = 0; i < col; ++i)
    {
        if(m_sequences[row].hasBaseInColumn(i))
            base_index += 1;
    }
    return base_index;
}

//
size_t MultipleAlignment::countHomopolymer(size_t row, int from, int to) const
{
    int num_columns = getNumColumns();
    assert(row < m_sequences.size());
    assert(from < num_columns);
    
    //
    if(to == -1)
        to = num_columns - 1; // inclusive

    // Determine iteration direction
    int step = from <= to ? 1 : -1;
    
    // The first or last column of the multiple alignment was requested
    // This can only be a homopolymer of length 1
    if(from == num_columns - 1 || (from == 0 && step < 0))
        return 1;
    
    // Advance to the next base symbol in the sequence
    char b = '\0';
    int max_position = num_columns - 1;
    while(from >= 0 && from <= max_position) {
        if(!m_sequences[row].hasBaseInColumn(from)) {
            from += step;
        } else {
            b = getSymbol(row, from);
            break;
        }
    }

    // Count the length, skipping padding symbols
    size_t length = 1;
    do
    {
        from += step;
        if(from < 0 || from > max_position)
            return length;

        char s = getSymbol(row, from);

        // skip gaps
        if(s == '-')
            continue;

        // matches the key base
        if(s == b)
            length += 1;
        else // break on mismatch or end of the alignment
            break;
    }
    while(from != to);
    return length;
}

//
std::string MultipleAlignment::getPileup(size_t idx) const
{
    std::string pileup;
    for(size_t i = 0; i < m_sequences.size(); ++i) {
        char symbol = m_sequences[i].getColumnSymbol(idx);
        if(symbol != '\0')
            pileup.push_back(symbol);
    }
    return pileup;
}

// Order the indices of the m_sequences vector by left-coordinate of the alignment on m_sequence[0]
struct SequencePrintOrder
{
    SequencePrintOrder(const std::vector<MultipleAlignmentElement>& sequence) : m_sequences(sequence) {}
    bool operator()(size_t i, size_t j)
    {
        if(i == 0)
            return true;
        if(j == 0)
            return false;
        return m_sequences[i].leading_columns < m_sequences[j].leading_columns;
    }
    private:
        const std::vector<MultipleAlignmentElement>& m_sequences;
};

//
void MultipleAlignment::print(size_t max_columns) const
{
    if(m_sequences.empty())
        return;

    size_t total_columns = m_sequences.front().getNumColumns();

    std::vector<size_t> sequence_index_order;

    // Initialize order
    for(size_t i = 0; i < m_sequences.size(); ++i)
        sequence_index_order.push_back(i);

    // Sort sequences by the first aligned base on sequence[0]
    SequencePrintOrder sorter(m_sequences);
    std::sort(sequence_index_order.begin(), sequence_index_order.end(), sorter);

    // Make a simple consensus to use as a mask
    std::string consensus = getPaddedConsensus();
    assert(consensus.size() == total_columns);

    // Print the multiple alignment in segments
    for(size_t c = 0; c < total_columns; c += max_columns) {
        size_t remaining = total_columns - c;
        size_t slice_size = max_columns < remaining ? max_columns : remaining;

        // Print the consensus
        printf("\t%s\tC\n", consensus.substr(c, slice_size).c_str());

        for(size_t i = 0; i < sequence_index_order.size(); ++i) {
            size_t current_index = sequence_index_order[i];

            // Build the output string
            std::string print_string;
            for(size_t j = c; j < slice_size + c; ++j) {
                char consensus_symbol = consensus[j];
                char row_symbol = m_sequences[current_index].getColumnSymbol(j);

                if(row_symbol == '\0')
                    print_string.push_back(' '); // no base for this column
                else if(row_symbol == '-' || row_symbol != consensus_symbol)
                    print_string.push_back(row_symbol); // always show sequence for row 0, gaps and mismatches
                else
                    print_string.push_back('.'); // mask matches
            }
            // Check if this string is blank, if so don't print it
            if(print_string.find_first_not_of(" ") != std::string::npos)
                printf("\t%s\t%s\n", print_string.c_str(), m_sequences[current_index].name.c_str());
        }
        printf("\n\n");
    }
}

//
void MultipleAlignment::printPileup() const
{
    if(m_sequences.empty())
        return;

    // Get the total number of columns in the alignment
    size_t num_columns = m_sequences.front().getNumColumns();
    size_t num_sequences = m_sequences.size();
    for(size_t i = 0; i < num_columns; ++i) {
        std::string counts_str = getColumnCountString(i);
        std::string pileup;
        std::string quality;
        for(size_t j = 0; j < num_sequences; ++j) {
            // sanity check that the columns are set up correctly
            assert(m_sequences[j].getNumColumns() == num_columns);
            char symbol = m_sequences[j].getColumnSymbol(i);
            if(symbol != '\0')
                pileup.push_back(symbol);

            char quality_symbol =  m_sequences[j].getColumnQuality(i);
            if(quality_symbol != '\0')
                quality.push_back(quality_symbol);

        }
        printf("%zu\t%s\t%s\t%s\n", i, pileup.c_str(), quality.c_str(), counts_str.c_str());
    }
}

std::string MultipleAlignment::getPaddedSubstr(size_t row, size_t start, size_t length) const
{
    assert(row < m_sequences.size());
    return m_sequences[row].getPaddedSubstr(start, length);
}

//
std::string MultipleAlignment::getPaddedConsensus() const
{
    std::string out;
    for(size_t i = 0; i < getNumColumns(); ++i) {
        SymbolCountVector symbol_counts = getSymbolCountVector(i);
        char symbol;
        if(symbol_counts.empty()) {
            symbol = ' ';
        } else {
            std::sort(symbol_counts.begin(), symbol_counts.end(), SymbolCount::countOrderDescending);
            symbol = symbol_counts.front().symbol;
        }
        out.append(1, symbol);
    }
    return out;
}

//
void MultipleAlignment::insertGapBeforeColumn(size_t column_index)
{
    for(size_t i = 0; i < m_sequences.size(); ++i) {
        m_sequences[i].insertGapBeforeColumn(column_index);
    }
}

//
void MultipleAlignment::trimEmptyColumns()
{
    size_t num_columns = getNumColumns();
    size_t leading_empty = 0;
    while(getPileup(leading_empty).empty())
        leading_empty += 1;
   
    size_t trailing_empty = 0;
    while(getPileup(num_columns - 1 - trailing_empty).empty())
        trailing_empty += 1;
   
    size_t num_rows = getNumRows();
    for(size_t i = 0; i < num_rows; ++i) {
        m_sequences[i].trimLeading(leading_empty);
        m_sequences[i].trimTrailing(trailing_empty);
    }
}

//
std::string MultipleAlignment::expandCigar(const std::string& cigar)
{
    std::string out;
    std::stringstream parser(cigar);
    int length;
    char symbol;
    while(parser >> length >> symbol)
        out.append(length, symbol);
    return out;
}

//
int MultipleAlignment::symbol2index(char symbol) const
{
    switch(std::toupper(symbol)) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        case '-':
            return 5;
        default:
            return 4; // all ambiguity codes get index 4
    }

    std::cerr << "Error: Unrecognized symbol in multiple alignment\n";
    exit(EXIT_FAILURE);
    return -1;
}

//
double MultipleAlignment::quality2prob(char q) const
{
    assert(q != '\0');
    return pow(10, -(q - 33)/10.0f);
}

//
char MultipleAlignment::prob2quality(double p) const
{
    int phred = 0;
    double lp = log10(p);

    // Clamp very low log probabilities.
    if(lp < -10.0)
        phred = 60;
    else
        phred = (int)(-10 * lp);

    // Clamp score into range
    phred = std::max(0, phred);
    phred = std::min(phred, 60);
    assert(phred >= 0 && phred <= 60);
    return phred + 33;
}

//
std::vector<int> MultipleAlignment::getColumnBaseCounts(size_t idx) const
{
    std::vector<int> out(m_alphabet_size, 0);
    for(size_t i = 0; i < m_sequences.size(); ++i) {
        char symbol = m_sequences[i].getColumnSymbol(idx);
        if(symbol != '\0') {

            out[symbol2index(symbol)] += 1;
        }
    }
    return out;
}

//
std::string MultipleAlignment::getColumnCountString(size_t column) const
{
    std::vector<int> counts = getColumnBaseCounts(column);
    std::stringstream out;
    for(size_t i = 0; i < m_alphabet_size; ++i) {
        out << m_alphabet[i] << ":" << counts[i] << " ";
    }
    return out.str();
}

//
SymbolCountVector MultipleAlignment::getSymbolCountVector(size_t column) const
{
    SymbolCountVector out;
    std::vector<int> counts = getColumnBaseCounts(column);
    for(size_t i = 0; i < m_alphabet_size; ++i) {
        if(counts[i] > 0) {
            SymbolCount sc = { m_alphabet[i], counts[i] };
            out.push_back(sc);
        }
    }
    return out;
}

