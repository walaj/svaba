///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// StatsProcess - Compute statistics about the reads
//
#include "StatsProcess.h"
#include "BWTAlgorithms.h"
#include "MultiOverlap.h"

//
//
//
StatsProcess::StatsProcess(const BWT* pBWT, const BWT* pRBWT, int kmerLength, int minOverlap, int branchCutoff, bool bNoOverlap) :
                            m_pBWT(pBWT),
                            m_pRBWT(pRBWT),
                            m_kmerLength(kmerLength),
                            m_minOverlap(minOverlap),
                            m_branchCutoff(branchCutoff),
                            m_bNoOverlap(bNoOverlap),
                            m_pAllOverlapper(NULL)
{
    if (!m_bNoOverlap)
    {
        assert(m_pRBWT);
        m_pAllOverlapper = new OverlapAlgorithm(m_pBWT, m_pRBWT, 0.05, 16, 16, false, m_branchCutoff);
    }
}

//
StatsProcess::~StatsProcess()
{
    if(m_pAllOverlapper)
    {
        delete m_pAllOverlapper;
        m_pAllOverlapper = NULL;
    }
}

//
StatsResult StatsProcess::process(const SequenceWorkItem& workItem)
{
    StatsResult result;

        
    std::string readSequence = workItem.read.seq.toString();

    //
    // Compute the kmer distribution
    //
    int k = m_kmerLength;
    int n = readSequence.size();
    int nk = n - m_kmerLength + 1;

    for(int i = 0; i < nk; ++i)
    {
        std::string kmer = readSequence.substr(i, k);
        int count = BWTAlgorithms::countSequenceOccurrences(kmer, m_pBWT);
        result.kmerCoverage.push_back(count);
    }
    
    //
    // Compute the number of implied errors in the read
    //
    if(!m_bNoOverlap)
    {
        SeqRecord currRead = workItem.read;
        OverlapBlockList blockList;
        OverlapResult overlap_result = m_pAllOverlapper->overlapRead(currRead, m_minOverlap, &blockList);
        
        // Convert the overlap block list into a multi-overlap 
        if(!overlap_result.searchAborted)
        {
            MultiOverlap mo = blockListToMultiOverlap(currRead, blockList);
            int covered = mo.countBasesCovered();
            int wrong = mo.countPotentialIncorrect(m_errorThreshold);
            result.mean_depth = mo.getMeanDepth();
            result.bases_counted += covered;
            result.bases_wrong += wrong;
        }
    }

    return result;
}

//
//
//
StatsPostProcess::StatsPostProcess(bool bPrintKmer) : m_bPrintKmer(bPrintKmer), m_basesCounted(0), m_basesWrong(0), m_depthSum(0.0f), m_numReads(0), m_numPerfect(0)
{
}

//
StatsPostProcess::~StatsPostProcess()
{
    std::cout << "\n*** Stats: \n";
    
    if(m_bPrintKmer)
    {
        int max = 100;
        m_kmerDist.print(max);
    }

    printf("%d out of %d bases are potentially incorrect (%lf)\n", m_basesWrong, m_basesCounted, (double)m_basesWrong/m_basesCounted);
    printf("%d reads out of %d are perfect (%lf)\n", m_numPerfect, m_numReads, (double)m_numPerfect/m_numReads);
    printf("Mean overlap depth: %.2lf\n", m_depthSum / m_numReads);
}

//
void StatsPostProcess::process(const SequenceWorkItem& /*item*/, const StatsResult& result)
{
    for(size_t i = 0; i < result.kmerCoverage.size(); ++i)
        m_kmerDist.add(result.kmerCoverage[i]);

    m_basesCounted += result.bases_counted;
    m_basesWrong += result.bases_wrong;
    m_depthSum += result.mean_depth;
    m_numReads += 1;
    m_numPerfect += (result.bases_wrong == 0 ? 1 : 0);
}
