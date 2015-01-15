///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// QCProcess - Process to perform quality checks 
// for a sequence work item
//
#ifndef QCPROCESS_H
#define QCPROCESS_H

#include "Util.h"
#include "BWT.h"
#include "SequenceProcessFramework.h"
#include "SequenceWorkItem.h"
#include "BitVector.h"

// Parameters
struct QCParameters
{
    QCParameters() { setDefaults(); }

    // Set reasonable default values for the qc filters
    void setDefaults()
    {
        checkDuplicates = true;
        checkKmer = true;
	kmerBothStrand = false;
        checkHPRuns = true;
        checkDegenerate = true;
        verbose = 0;

        pBWT = NULL;
        pRevBWT = NULL;
        pSharedBV = NULL;

        kmerLength = 27;
        kmerThreshold = 2;

        hpHardAcceptCount = 10;
        hpMinProportion = 0.1f;
        hpKmerLength = 51;
        hpMinLength = 6;
        hpMinContext = 5;

        degenProportion = 0.90;
    }

    const BWT* pBWT;
    const BWT* pRevBWT;
    BitVector* pSharedBV;

    // Control parameters
    bool checkDuplicates;
    bool checkKmer;
    bool kmerBothStrand;
    bool checkHPRuns;
    bool checkDegenerate;
    bool substringOnly;
    int verbose;

    //
    // Kmer frequency filtering parameters
    //
    int kmerLength;
    int kmerThreshold;

    //
    // Homopolymer filter parameters
    //

    // The length a run must be to trigger filtering
    size_t hpMinLength;

    // If the k-mer covering the HP run is seen
    // at least this many times, it is accepted
    size_t hpHardAcceptCount; 

    // The hp run length in the read must be
    // at least this fraction of the dominant
    // run length to be kept
    double hpMinProportion;

    // The context k-mer length for the homopolymer
    // filter.
    size_t hpKmerLength;
    
    // The minimum amount of sequence around the homopolymer
    // to trigger filtering
    size_t hpMinContext;

    //
    // Degenerate filter parameters
    //

    // If the most frequent base in the read
    // makes up more than this proportion, the read
    // is discarded
    double degenProportion;

};

// Results object
class QCResult
{
    public:
        QCResult() : kmerPassed(false), dupPassed(false), hpPassed(false), degenPassed(false) {}

        bool kmerPassed;
        bool dupPassed;
        bool hpPassed;
        bool degenPassed;
};

enum DuplicateCheckResult
{
    DCR_UNIQUE,
    DCR_SUBSTRING,
    DCR_FULL_LENGTH_DUPLICATE
};

// Perform quality checks on the input stream of reads
class QCProcess
{
    public:
        QCProcess(QCParameters params); 
        ~QCProcess();
        QCResult process(const SequenceWorkItem& item);
        
        // Discard reads with low-frequency kmers
        bool performKmerCheck(const SequenceWorkItem& item);

        // Discard reads that are identical to, or a substring of, some other read
        DuplicateCheckResult performDuplicateCheck(const SequenceWorkItem& item);

        // Check whether the sequence has a homopolymer sequencing error. This
        // check finds a kmer covering a homopolymer run then modifies it to check
        // for a shorter/longer run of higher frequency
        bool performHomopolymerCheck(const SequenceWorkItem& item);

        // Check if the sequence read is degenerate
        // A degenerate read is define to be reads consisting almost
        // entirely of a single nucleotide
        bool performDegenerateCheck(const SequenceWorkItem& item);

    private:
        
        const QCParameters m_params;
};

// Write the results from the overlap step to an ASQG file
class QCPostProcess
{
    public:
        QCPostProcess(std::ostream* pCorrectedWriter, std::ostream* pDiscardWriter);
        ~QCPostProcess();

        void process(const SequenceWorkItem& item, const QCResult& result);

    private:

        std::ostream* m_pCorrectedWriter;
        std::ostream* m_pDiscardWriter;

        size_t m_readsKept;
        size_t m_readsDiscarded;
        size_t m_readsFailedKmer;
        size_t m_readsFailedDup;
        size_t m_readsFailedHP;
        size_t m_readsFailedDegen;
};

#endif
