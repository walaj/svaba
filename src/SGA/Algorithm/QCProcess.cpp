///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// QCProcess - Process to perform quality checks 
// for a sequence work item
//
#include "QCProcess.h"
#include "BWTAlgorithms.h"

//
struct KmerWindow
{
    int64_t getCount(bool bothStrand) const
    {
        int64_t fwdCount = (fwdIntervals.interval[0].isValid()) ? fwdIntervals.interval[0].size() : 0;
        int64_t rcCount = (rcIntervals.interval[0].isValid()) ? rcIntervals.interval[0].size() : 0;
        if (bothStrand)
            return std::min(fwdCount,rcCount);
        else
            return fwdCount + rcCount;
    }

    int start;
    int end;
    BWTIntervalPair fwdIntervals;
    BWTIntervalPair rcIntervals;

    bool isInitialized;
};

//
//
//
QCProcess::QCProcess(QCParameters params) : m_params(params)
{

}

//
QCProcess::~QCProcess()
{

}

//
QCResult QCProcess::process(const SequenceWorkItem& workItem)
{
    QCResult result;

    if(m_params.checkDuplicates) 
    {
        DuplicateCheckResult dupCheckResult = performDuplicateCheck(workItem);
        if(!m_params.substringOnly) 
            result.dupPassed = dupCheckResult == DCR_UNIQUE;
        else
            result.dupPassed = dupCheckResult != DCR_SUBSTRING;
    }
    else
    {
        result.dupPassed = true;
    }
        

    // Only perform the more expensive k-mer test if the dup check succeeded
    if(m_params.checkKmer && result.dupPassed)
        result.kmerPassed = performKmerCheck(workItem);
    else
        result.kmerPassed = true;

    if(result.kmerPassed && result.dupPassed && m_params.checkHPRuns)
        result.hpPassed = performHomopolymerCheck(workItem);
    else
        result.hpPassed = true;

    if(m_params.checkDegenerate && result.dupPassed && result.kmerPassed && result.hpPassed)
        result.degenPassed = performDegenerateCheck(workItem);
    else
        result.degenPassed = true;

    return result;
}

bool QCProcess::performKmerCheck(const SequenceWorkItem& workItem)
{
    // Perform a k-mer filter on the read
    // Each k-mer must be seen at least m times for the read to be kept
    // Naively, this would required k*(l - k + 1) calls to the BWT
    // if each k-mer was computed independently. We improve this by performing
    // a two-stage test. If the kmer k_i is seen at least m times, we extend
    // it by one base. If this k+1-mer is also seen at least m times, it implies
    // that the kmer k_(i+1) is also seen at least m times so it does not need
    // to be recomputed. If the (k+1)-mer is seen less than m times, we recompute
    // kmer k_(i+1) as a final check. Using this optimistic algorithm (checking
    // longer kmers than is required) is significantly faster.
    QCResult result;

    std::string w = workItem.read.seq.toString();

    // Ensure the read is longer than the k-mer length
    if((int)w.size() < m_params.kmerLength)
        return false;

    int k = m_params.kmerLength;
    int n = w.size();
    int nk = n - k + 1;
    int threshold = m_params.kmerThreshold;

    // Are all kmers in the read well-represented?
    bool allSolid = true;

    int i = 0;

    KmerWindow window;
    window.start = 0;
    window.end = 0;
    window.isInitialized = false;
    
    while(i < nk)
    {
        // Determine whether the kmer at position i
        // has been seen more than m times.
        if(window.isInitialized)
        {
            // Extend the current window by 1 base
            // to see if the spanning larger
            // kmer is seen more than m times.
            // If so, we don't need to check this kmer
            // individually
            int next = window.end + k;
            assert(next < (int)w.size());
            char b = w[next];
            char cb = complement(b);

            if(window.fwdIntervals.interval[0].isValid())
                BWTAlgorithms::updateBothR(window.fwdIntervals, b, m_params.pRevBWT);
            if(window.rcIntervals.interval[1].isValid())
                BWTAlgorithms::updateBothR(window.rcIntervals, cb, m_params.pBWT);

            int64_t count = window.getCount(m_params.kmerBothStrand);
            if(count <= threshold)
            {
                // The extended kmer didn't meet the threshold
                // recompute the interval for this kmer
                window.isInitialized = false;
            }
            else
            {
                window.end += 1;
            }
        }
        
        if(!window.isInitialized)
        {
            // initialize the window by computing the
            // BWTIntervals for the kmer starting at
            // i and its reverse complement
            char b = w[i];
            char cb = complement(b);
            BWTAlgorithms::initIntervalPair(window.fwdIntervals, b, m_params.pBWT, m_params.pRevBWT);
            BWTAlgorithms::initIntervalPair(window.rcIntervals, cb, m_params.pRevBWT, m_params.pBWT);

            for(int j = i + 1; j < i + k; ++j)
            {
                // Update intervals rightwards 
                b = w[j];
                cb = complement(b);

                if(window.fwdIntervals.interval[0].isValid())
                    BWTAlgorithms::updateBothR(window.fwdIntervals, b, m_params.pRevBWT);
                if(window.rcIntervals.interval[1].isValid())
                    BWTAlgorithms::updateBothR(window.rcIntervals, cb, m_params.pBWT);
            }

            // record the start/end indices of the kmers spanned by this position
            window.start = i;
            window.end = i;
            window.isInitialized = true;
        }

        // Check the count of this interval
        int64_t count = window.getCount(m_params.kmerBothStrand);
            
        if(count <= threshold)
        {
            allSolid = false;
            break;
        }
        else
        {
            i += 1;
        }
    }

    return allSolid;
}

// Perform duplicate check
// Look up the interval of the read in the BWT. If the index of the read
DuplicateCheckResult QCProcess::performDuplicateCheck(const SequenceWorkItem& workItem)
{
    assert(m_params.pSharedBV != NULL);

    std::string w = workItem.read.seq.toString();
    std::string rc_w = reverseComplement(w);

    // Look up the interval of the sequence and its reverse complement
    BWTIntervalPair fwdIntervals = BWTAlgorithms::findIntervalPair(m_params.pBWT, m_params.pRevBWT, w);
    BWTIntervalPair rcIntervals = BWTAlgorithms::findIntervalPair(m_params.pBWT, m_params.pRevBWT, rc_w);

    // Check if this read is a substring of any other
    // This is indicated by the presence of a non-$ extension in the left or right direction
    AlphaCount64 fwdECL = BWTAlgorithms::getExtCount(fwdIntervals.interval[0], m_params.pBWT);
    AlphaCount64 fwdECR = BWTAlgorithms::getExtCount(fwdIntervals.interval[1], m_params.pRevBWT);

    AlphaCount64 rcECL = BWTAlgorithms::getExtCount(rcIntervals.interval[0], m_params.pBWT);
    AlphaCount64 rcECR = BWTAlgorithms::getExtCount(rcIntervals.interval[1], m_params.pRevBWT);

    if(fwdECL.hasDNAChar() || fwdECR.hasDNAChar() || rcECL.hasDNAChar() || rcECR.hasDNAChar())
    {
        // Substring reads are always removed so no need to update the bit vector
        return DCR_SUBSTRING;
    }

    // Calculate the lexicographic intervals for the fwd and reverse intervals
    BWTAlgorithms::updateBothL(fwdIntervals, '$', m_params.pBWT);
    BWTAlgorithms::updateBothL(rcIntervals, '$', m_params.pBWT);

    // Calculate the canonical index for this string - the lowest
    // value in the two lexicographic index
    assert(fwdIntervals.interval[0].isValid() || rcIntervals.interval[0].isValid());
    int64_t fi = fwdIntervals.interval[0].isValid() ? fwdIntervals.interval[0].lower : std::numeric_limits<int64_t>::max();
    int64_t ri = rcIntervals.interval[0].isValid() ? rcIntervals.interval[0].lower : std::numeric_limits<int64_t>::max();
    int64_t canonicalIdx = std::min(fi, ri);

    // Check if the bit reprsenting the canonical index is set in the shared bit vector
    if(!m_params.pSharedBV->test(canonicalIdx))
    {
        // This read is not a duplicate
        // Attempt to atomically set the bit from false to true
        if(m_params.pSharedBV->updateCAS(canonicalIdx, false, true))
        {
            // Call succeed, return that this read is not a duplicate
            return DCR_UNIQUE;
        }
        else
        {
            // Call failed, some other thread set the bit before
            // this thread. Return that the reead is a duplicate
            return DCR_FULL_LENGTH_DUPLICATE;
        }
    }
    else
    {
        // this read is duplicate
        return DCR_FULL_LENGTH_DUPLICATE;
    }
}

// Perform homopolymer filter
bool QCProcess::performHomopolymerCheck(const SequenceWorkItem& item)
{
    std::string w = item.read.seq.toString();
    size_t k = m_params.hpKmerLength;

    // Skip if the read length is less than the kmer size
    // used for this test
    if(w.size() < k)
        return true;

    size_t maxRunLength = 0;
    size_t maxRunStart = 0;
    
    // Count the longest homopolymer run in the read
    size_t currRunStart = 0;
    size_t currRunLength = 1;
    char prev = w[0];
    char runChar = prev;
    for(size_t i = 1; i < w.size(); ++i)
    {
        if(w[i] == prev)
            currRunLength += 1;

        if(w[i] != prev || i == w.size() - 1)
        {
            if(currRunLength > maxRunLength)
            {
                maxRunLength = currRunLength;
                maxRunStart = currRunStart;
                runChar = prev;
            }
            currRunLength = 1;
            prev = w[i];
            currRunStart = i;
        }
    }

    // If the run length is above the threshold, find
    // the frequency of the covering k-mer and k-mers with
    // different lengths of the homopolymer
    if(maxRunLength >= m_params.hpMinLength && maxRunLength < k / 2)
    {
        // Extract the kmer covering this run
        int hprMiddle = maxRunStart + (maxRunLength / 2);
        int estimatedKmerStart = hprMiddle - (k / 2);
        int trueKmerStart = estimatedKmerStart;
        
        // If the calculated start position of the covering k-mer
        //is before the start of the read, or past the end, shift it
        if(estimatedKmerStart < 0)
            trueKmerStart = 0;

        if(estimatedKmerStart + k > w.size())
            trueKmerStart = w.size() - k;

        std::string prefix = w.substr(trueKmerStart, maxRunStart - trueKmerStart);
        std::string suffix = w.substr(maxRunStart + maxRunLength, trueKmerStart + k - (maxRunStart + maxRunLength));

        // Do not attempt this check if we do not have enough context preceding or following the HPR
        if(prefix.size() < m_params.hpMinContext || suffix.size() < m_params.hpMinContext)
            return true;

        size_t highestCountLength = 0;
        size_t highestCount = 0;
        size_t actualCount = 0;
        for(size_t l = maxRunLength - 2; l <= maxRunLength + 2; ++l)
        {
            std::string composite = prefix + std::string(l, runChar) + suffix;
            size_t count = BWTAlgorithms::countSequenceOccurrences(composite, m_params.pBWT);
            if(l == maxRunLength)
                actualCount = count;

            if(count > highestCount)
            {
                highestCount = count;
                highestCountLength = l;
            }
        }

        double proportion = (double)actualCount / (double)highestCount;
        if(highestCountLength == maxRunLength || actualCount >= m_params.hpHardAcceptCount || proportion >= m_params.hpMinProportion)
        {
            return true;
        }
        else
        {
            if(m_params.verbose > 0)
            {
                printf("Read failed homopolymer filter %s\n", w.c_str());
                printf("Filtered read with poly-%c run. DL: %zu DC: %zu. AL: %zu AC: %zu P: %lf\n", runChar, highestCountLength, highestCount, maxRunLength, actualCount, proportion);
            }
            return false;
        }
    }

    return true;
}

// Check if sequence is composed of predominantely a single base
// Returns true if the sequence is not degenrate
bool QCProcess::performDegenerateCheck(const SequenceWorkItem& item)
{
    std::string w = item.read.seq.toString();
    AlphaCount64 bc;
    for(size_t i = 0; i < w.size(); ++i)
    {
        bc.increment(w[i]);
    }

    size_t maxCount = bc.getMaxCount();
    double prop = (double)maxCount / w.size();
    if(prop > m_params.degenProportion)
    {
        if(m_params.verbose > 0)
            std::cout << "Read " << w << " failed degenerate filter\n";
        return false;
    }
    return true;
}
//
//
//
QCPostProcess::QCPostProcess(std::ostream* pCorrectedWriter,
                             std::ostream* pDiscardWriter) :
                                m_pCorrectedWriter(pCorrectedWriter),
                                m_pDiscardWriter(pDiscardWriter),
                                m_readsKept(0), m_readsDiscarded(0),
                                m_readsFailedKmer(0), m_readsFailedDup(0),
                                m_readsFailedHP(0), m_readsFailedDegen(0)
{

}

//
QCPostProcess::~QCPostProcess()
{
    std::cout << "Reads kept: " << m_readsKept << "\n";
    std::cout << "Reads discarded: " << m_readsDiscarded << "\n";
    std::cout << "Reads failed kmer check: " << m_readsFailedKmer << "\n";
    std::cout << "Reads failed duplicate check: " << m_readsFailedDup << "\n";
    std::cout << "Reads failed homopolymer check: " << m_readsFailedHP << "\n";
    std::cout << "Reads failed degenerate check: " << m_readsFailedDegen << "\n";
}

//
void QCPostProcess::process(const SequenceWorkItem& item, const QCResult& result)
{
    SeqRecord record = item.read;
    if(result.kmerPassed && result.dupPassed && result.hpPassed && result.degenPassed)
    {
        record.write(*m_pCorrectedWriter);
        ++m_readsKept;
    }
    else
    {
        // To be able to rebuild the index after discarding the read, we need to write
        // the rank of the string (its position in the original read file into the read name)
        std::stringstream newID;
        newID << item.read.id << ",seqrank=" << item.idx;
        record.id = newID.str();

        record.write(*m_pDiscardWriter);
        ++m_readsDiscarded;

        if(!result.kmerPassed)
            m_readsFailedKmer += 1;
        else if(!result.dupPassed)
            m_readsFailedDup += 1;
        else if(!result.hpPassed)
            m_readsFailedHP += 1;
        else if(!result.degenPassed)
            m_readsFailedDegen += 1;
    }
}
