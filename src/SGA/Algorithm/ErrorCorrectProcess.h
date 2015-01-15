///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ErrorCorrectProcess - Wrapper to perform error correction
// for a sequence work item
//
#ifndef CORRECTPROCESS_H
#define CORRECTPROCESS_H

#include "Util.h"
#include "OverlapAlgorithm.h"
#include "SequenceProcessFramework.h"
#include "SequenceWorkItem.h"
#include "MultiOverlap.h"
#include "Metrics.h"
#include "BWTIndexSet.h"
#include "SampledSuffixArray.h"
#include "multiple_alignment.h"

enum ErrorCorrectAlgorithm
{
    ECA_HYBRID, // hybrid kmer/overlap correction
    ECA_KMER, // kmer correction
    ECA_OVERLAP, // overlap correction
    ECA_THREAD // thread the read through a de Bruijn graph
};

enum ECFlag
{
    ECF_NOTCORRECTED,
    ECF_CORRECTED,
    ECF_AMBIGIOUS,
    ECF_DUPLICATE
};

// Parameter object for the error corrector
struct ErrorCorrectParameters
{
    ErrorCorrectAlgorithm algorithm;

    //
    const OverlapAlgorithm* pOverlapper;
    BWTIndexSet indices;

    // Overlap-based corrector params
    int minOverlap;
    int numOverlapRounds;
    double minIdentity;
    int conflictCutoff;
    int depthFilter;

    // k-mer based corrector params
    int numKmerRounds;
    int kmerLength;

    // output options
    bool printOverlaps;
};

class ErrorCorrectResult
{
    public:
        ErrorCorrectResult() : num_prefix_overlaps(0), num_suffix_overlaps(0), kmerQC(false), overlapQC(false) {}

        DNAString correctSequence;
        ECFlag flag;

        // Metrics
        size_t num_prefix_overlaps;
        size_t num_suffix_overlaps;
        bool kmerQC;
        bool overlapQC;
};

//
class ErrorCorrectProcess
{
    public:
        ErrorCorrectProcess(const ErrorCorrectParameters params); 
        ~ErrorCorrectProcess();

        ErrorCorrectResult process(const SequenceWorkItem& item);
        ErrorCorrectResult correct(const SequenceWorkItem& item);

        ErrorCorrectResult kmerCorrection(const SequenceWorkItem& item);
        ErrorCorrectResult overlapCorrection(const SequenceWorkItem& workItem);
        ErrorCorrectResult overlapCorrectionNew(const SequenceWorkItem& workItem);
        ErrorCorrectResult threadingCorrection(const SequenceWorkItem& workItem);
    
    private:

        bool attemptKmerCorrection(size_t i, size_t k_idx, size_t minCount, std::string& readSequence);

        OverlapBlockList m_blockList;
        ErrorCorrectParameters m_params;
};

// Write the results from the overlap step to an ASQG file
class ErrorCorrectPostProcess
{
    public:
        ErrorCorrectPostProcess(std::ostream* pCorrectedWriter, 
                                std::ostream* pDiscardWriter, bool bCollectMetrics);

        ~ErrorCorrectPostProcess();

        void process(const SequenceWorkItem& item, const ErrorCorrectResult& result);
        void writeMetrics(std::ostream* pWriter);

    private:

        void collectMetrics(const std::string& originalSeq, 
                            const std::string& correctedSeq, const std::string& qualityStr);

        std::ostream* m_pCorrectedWriter;
        std::ostream* m_pDiscardWriter;
        bool m_bCollectMetrics;

        ErrorCountMap<char> m_qualityMetrics;
        ErrorCountMap<int64_t> m_positionMetrics;
        ErrorCountMap<char> m_originalBaseMetrics;
        ErrorCountMap<std::string> m_precedingSeqMetrics;

        size_t m_totalBases;
        size_t m_totalErrors;
        size_t m_readsKept;
        size_t m_readsDiscarded;

        size_t m_kmerQCPassed;
        size_t m_overlapQCPassed;
        size_t m_qcFail;
};

#endif
