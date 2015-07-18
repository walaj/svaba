//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// RLBWT - Run-length encoded Burrows Wheeler transform
//
#include "RLBWT.h"
#include "Timer.h"
#include "BWTReader.h"
#include "BWTWriter.h"
#include "BWTReader.h"
#include <istream>
#include <queue>
#include <inttypes.h>

// macros
#define OCC(c,i) m_occurrence.get(m_bwStr, (c), (i))
#define PRED(c) m_predCount.get((c))

// Parse a BWT from a file
RLBWT::RLBWT(const std::string& filename, int sampleRate) : m_numStrings(0), 
                                                            m_numSymbols(0), 
                                                            m_largeSampleRate(DEFAULT_SAMPLE_RATE_LARGE),
                                                            m_smallSampleRate(sampleRate)
{
    IBWTReader* pReader = BWTReader::createReader(filename);
    pReader->read(this);
    initializeFMIndex();
    delete pReader;
}

// Construct the BWT from a suffix array
RLBWT::RLBWT(const SuffixArray* pSA, const ReadTable* pRT)
{

      // Set up BWT state
    size_t n = pSA->getSize();
    m_smallSampleRate = DEFAULT_SAMPLE_RATE_SMALL;
    m_largeSampleRate = DEFAULT_SAMPLE_RATE_LARGE;
    m_numStrings = pSA->getNumStrings();
    m_numSymbols = n;

    RLUnit currRun;
    // Set up the bwt string and suffix array from the cycled strings
    for(size_t i = 0; i < n; ++i)
    {

        SAElem saElem = pSA->get(i);
        const SeqItem& si = pRT->getRead(saElem.getID());

        // Get the position of the start of the suffix
        uint64_t f_pos = saElem.getPos();
        uint64_t l_pos = (f_pos == 0) ? si.seq.length() : f_pos - 1;
        char b = (l_pos == si.seq.length()) ? '$' : si.seq.get(l_pos);

        // Add to the current run or append in the new char
        if(currRun.isInitialized())
        {
            if(currRun.getChar() == b && !currRun.isFull())
            {
                currRun.incrementCount();
            }
            else
            {
                // Write out the old run and start a new one
                assert(currRun.isInitialized());
                m_rlString.push_back(currRun);
                currRun = RLUnit(b);
            }        
        }
        else
        {
            // Start a new run
            currRun = RLUnit(b);
        }
    }
    assert(currRun.isInitialized());
    if(currRun.isInitialized())
        m_rlString.push_back(currRun);

    initializeFMIndex();
}

//
void RLBWT::append(char b)
{
    bool increment = false;
    if(!m_rlString.empty())
    {
        RLUnit& lastUnit = m_rlString.back();
        if(lastUnit.getChar() == b && !lastUnit.isFull())
        {
            lastUnit.incrementCount();
            increment = true;
        }
    }

    if(!increment)
    {
        // Must add a new unit to the string
        RLUnit unit(b);
        m_rlString.push_back(unit);
    }
    ++m_numSymbols;
}

// Fill in the FM-index data structures
void RLBWT::initializeFMIndex()
{
    m_smallShiftValue = Occurrence::calculateShiftValue(m_smallSampleRate);
    m_largeShiftValue = Occurrence::calculateShiftValue(m_largeSampleRate);

    // initialize the marker vectors,
    // LargeMarkers are placed every 2048 bases (by default) containing the absolute count
    // of symbols seen up to that point. SmallMarkers are placed every 128 bases with the
    // count over the last 128 symbols. From these relative counts the absolute count
    // every 128 symbols can be interpolated.
    size_t num_large_markers = getNumRequiredMarkers(m_numSymbols, m_largeSampleRate);
    size_t num_small_markers = getNumRequiredMarkers(m_numSymbols, m_smallSampleRate);
    m_largeMarkers.resize(num_large_markers);
    m_smallMarkers.resize(num_small_markers);

    // Fill in the marker values
    // We wish to place markers every sampleRate symbols however since a run may
    // not end exactly on sampleRate boundaries, we place the markers AFTER
    // the run crossing the boundary ends

    // Place a blank markers at the start of the data
    m_largeMarkers[0].unitIndex = 0;
    m_smallMarkers[0].unitCount = 0;

    // State variables for the number of markers placed,
    // the next marker to place, etc
    size_t curr_large_marker_index = 1;
    size_t curr_small_marker_index = 1;

    size_t next_small_marker = m_smallSampleRate;
    size_t next_large_marker = m_largeSampleRate;

    size_t prev_small_marker_unit_index = 0;
    size_t running_total = 0;
    AlphaCount64 running_ac;

    for(size_t i = 0; i < m_rlString.size(); ++i)
    {
        // Update the count and advance the running total
        RLUnit& unit = m_rlString[i];

        char symbol = unit.getChar();
        uint8_t run_len = unit.getCount();
        running_ac.add(symbol, run_len);
        running_total += run_len;

        size_t curr_unit_index = i + 1;
        bool last_symbol = i == m_rlString.size() - 1;

        // Check whether to place a new large marker
        bool place_last_large_marker = last_symbol && curr_large_marker_index < num_large_markers;
        while(running_total >= next_large_marker || place_last_large_marker)
        {
            size_t expected_marker_pos = curr_large_marker_index * m_largeSampleRate;

            // Sanity checks
            // The marker position should always be less than the running total unless 
            // the number of symbols is smaller than the sample rate
            assert(expected_marker_pos <= running_total || place_last_large_marker);
            assert((running_total - expected_marker_pos) <= RL_FULL_COUNT || place_last_large_marker);
            assert(curr_large_marker_index < num_large_markers);
            assert(running_ac.getSum() == running_total);

            LargeMarker& marker = m_largeMarkers[curr_large_marker_index];
            marker.unitIndex = i + 1;
            marker.counts = running_ac;

            next_large_marker += m_largeSampleRate;
            curr_large_marker_index += 1;
            place_last_large_marker = last_symbol && curr_large_marker_index < num_large_markers;
        }    

        // Check whether to place a new small marker
        bool place_last_small_marker = last_symbol && curr_small_marker_index < num_small_markers;
        while(running_total >= next_small_marker || place_last_small_marker)
        {
            // Place markers
            size_t expected_marker_pos = curr_small_marker_index * m_smallSampleRate;

            // Sanity checks
            // The marker position should always be less than the running total unless 
            // the number of symbols is smaller than the sample rate
            assert(expected_marker_pos <= running_total || place_last_small_marker);
            assert((running_total - expected_marker_pos) <= RL_FULL_COUNT || place_last_small_marker);
            assert(curr_small_marker_index < num_small_markers);
            assert(running_ac.getSum() == running_total);
    
            // Calculate the number of rl units that are contained in this block
            if(curr_unit_index - prev_small_marker_unit_index > std::numeric_limits<uint16_t>::max())
            {
                std::cerr << "Error: Number of units in occurrence array block " << curr_small_marker_index 
                          << " exceeds the maximum value.\n";
                exit(EXIT_FAILURE);

            }

            // Calculate the large marker to set the relative count from
            // This is generally the most previously placed large block except it might 
            // be the second-previous in the case that we placed the last large marker.
            size_t large_marker_index = expected_marker_pos >> m_largeShiftValue;
            assert(large_marker_index < curr_large_marker_index); // ensure the last has ben placed
            LargeMarker& prev_large_marker = m_largeMarkers[large_marker_index];

            // Set the 8bit AlphaCounts as the sum since the last large (superblock) marker
            AlphaCount16 smallAC;
            for(size_t j = 0; j < ALPHABET_SIZE; ++j)
            {
                size_t v = running_ac.getByIdx(j) - prev_large_marker.counts.getByIdx(j);
                if(v > smallAC.getMaxValue())
                {
                    std::cerr << "Error: Number of symbols in occurrence array block " << curr_small_marker_index 
                              << " exceeds the maximum value (" << v << " > " << smallAC.getMaxValue() << ")\n";
                    exit(EXIT_FAILURE);
                }
                smallAC.setByIdx(j, v);
            }
            
            // Set the small marker
            SmallMarker& small_marker = m_smallMarkers[curr_small_marker_index];
            small_marker.unitCount = curr_unit_index - prev_large_marker.unitIndex;
            small_marker.counts = smallAC;

            // Update state variables
            next_small_marker += m_smallSampleRate;
            curr_small_marker_index += 1;
            prev_small_marker_unit_index = curr_unit_index;
            place_last_small_marker = last_symbol && curr_small_marker_index < num_small_markers;
        }    
    }

    assert(curr_small_marker_index == num_small_markers);
    assert(curr_large_marker_index == num_large_markers);

    // Initialize C(a)
    m_predCount.set('$', 0);
    m_predCount.set('A', running_ac.get('$')); 
    m_predCount.set('C', m_predCount.get('A') + running_ac.get('A'));
    m_predCount.set('G', m_predCount.get('C') + running_ac.get('C'));
    m_predCount.set('T', m_predCount.get('G') + running_ac.get('G'));
}

// get the number of markers required to cover the n symbols at sample rate of d
size_t RLBWT::getNumRequiredMarkers(size_t n, size_t d) const
{
    // we place a marker at the beginning (with no accumulated counts), every m_sampleRate
    // bases and one at the very end (with the total counts)
    size_t num_markers = (n % d == 0) ? (n / d) + 1 : (n / d) + 2;
    return num_markers;
}

// Print the BWT
void RLBWT::print() const
{
    size_t numRuns = getNumRuns();
    std::string bwt;
    for(size_t i = 0; i < numRuns; ++i)
    {
        const RLUnit& unit = m_rlString[i];
        char symbol = unit.getChar();
        size_t length = unit.getCount();
        for(size_t j = 0; j < length; ++j)
            std::cout << symbol;
        std::cout << " : " << symbol << "," << length << "\n"; 
        bwt.append(length, symbol);
    }
    std::cout << "B: " << bwt << "\n";
}

// Print information about the BWT
void RLBWT::printInfo() const
{
    size_t small_m_size = m_smallMarkers.capacity() * sizeof(SmallMarker);
    size_t large_m_size = m_largeMarkers.capacity() * sizeof(LargeMarker);
    size_t total_marker_size = small_m_size + large_m_size;

    size_t bwStr_size = m_rlString.capacity() * sizeof(RLUnit);
    size_t other_size = sizeof(*this);
    size_t total_size = total_marker_size + bwStr_size + other_size;

    double mb = (double)(1024 * 1024);
    double total_mb = total_size / mb;
    
    printf("\nRLBWT info:\n");
    printf("Large Sample rate: %zu\n", m_largeSampleRate);
    printf("Small Sample rate: %zu\n", m_smallSampleRate);
    printf("Contains %zu symbols in %zu runs (%1.4lf symbols per run)\n", m_numSymbols, m_rlString.size(), (double)m_numSymbols / m_rlString.size());
    printf("Marker Memory -- Small Markers: %zu (%.1lf MB) Large Markers: %zu (%.1lf MB)\n", small_m_size, small_m_size / mb, large_m_size, large_m_size / mb);
    printf("Total Memory -- Markers: %zu (%.1lf MB) Str: %zu (%.1lf MB) Misc: %zu Total: %zu (%lf MB)\n", total_marker_size, total_marker_size / mb, bwStr_size, bwStr_size / mb, other_size, total_size, total_mb);
    printf("N: %zu Bytes per symbol: %lf\n\n", m_numSymbols, (double)total_size / m_numSymbols);
}

// Print the run length distribution of the BWT
void RLBWT::printRunLengths() const
{
    typedef std::map<size_t, size_t> DistMap;
    DistMap rlDist;

    char prevSym = '\0';
    size_t prevRunLen = 0;
    size_t currLen = 0;
    size_t adjacentSingletons = 0;
    size_t numRuns = getNumRuns();
    size_t totalRuns = 0;
    for(size_t i = 0; i < numRuns; ++i)
    {
        const RLUnit& unit = m_rlString[i];
        size_t length = unit.getCount();
        if(unit.getChar() == prevSym)
        {
            currLen += unit.getCount();
        }
        else
        {
            if(prevSym != '\0')
            {
                if(currLen >= 200)
                    rlDist[200]++;
                else
                    rlDist[currLen]++;
                totalRuns++;
            }
            currLen = length;
            prevSym = unit.getChar();
        }

        if(length == 1 && prevRunLen == 1)
        {
            adjacentSingletons += 1;
            prevRunLen = 0;
        }
        prevRunLen = length;
    }
    
    printf("Run length distrubtion\n");
    printf("rl\tcount\tfrac\n");
    double cumulative_mb = 0.0f;
    for(DistMap::iterator iter = rlDist.begin(); iter != rlDist.end(); ++iter)
    {
        cumulative_mb += ((double)iter->second / (1024 * 1024));
        printf("%zu\t%zu\t%lf\t%lf\n", iter->first, iter->second, double(iter->second) / totalRuns, cumulative_mb);
    }
    printf("Total runs: %zu\n", totalRuns);
    printf("Number of adjacent singleton runs: %zu\n", adjacentSingletons);
    printf("Minimal runs: %zu\n", totalRuns - adjacentSingletons / 2);
}

