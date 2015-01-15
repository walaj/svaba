//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// MultiOverlap.h - Data structure containing a set
// of overlaps for a given read
//
#ifndef MULTIOVERLAP_H
#define MULTIOVERLAP_H

#include "Match.h"
#include "Pileup.h"
#include "DNADouble.h"

class MultiOverlap
{
    public:
    struct MOData
    {
        MOData(const std::string& s, const Overlap& o) : seq(s), ovr(o), offset(0), partitionID(0), score(0.0f) {}
        static bool sortOffset(const MOData& a, const MOData& b);
        static bool sortID(const MOData& a, const MOData& b);
        
        // data
        std::string seq;
        Overlap ovr;
        int offset;
        int partitionID;
        double score;

        static bool compareScore(const MOData& a, const MOData& b)
        {
            return a.score > b.score;
        }
    };

    typedef std::vector<MOData> MODVector;

    public:

        MultiOverlap(const std::string& rootID, 
                     const std::string& rootSeq, 
                     const std::string rootQual = "");
        
        ///
        void add(const std::string& seq, const Overlap& ovr);
        void add(const MOData& mod);
        void updateRootSeq(const std::string& newSeq);

        //
        Overlap getOverlap(size_t idx) const;
        size_t getNumBases() const;

        //
        bool isConflicted(size_t cutoff) const;
        std::string simpleConsensus() const;

        // Count the number of bases that are potentially incorrect
        // A base is considered to be incorrect if it is not the majority in the column
        // and it has been seen less than cutoff times
        int countPotentialIncorrect(size_t cutoff) const;

        // Count the number of bases in the read that are covered by an overlap
        int countBasesCovered() const;

        // Partition the multioverlap into groups
        int getPartition(size_t idx) const;
        void setPartition(size_t idx, int p);
        std::string consensusConflict(double p_error, int conflictCutoff);

        // Returns true if each base of the root sequence has enough support from
        // other reads in the multioverlap. The level of support is determined by
        // the quality score.
        bool qcCheck() const;

        // Count the number of prefix and suffix overlaps
        void countOverlaps(size_t& prefix_count, size_t& suffix_count) const;
        double getMeanDepth() const;

        // Calculate the amount of the sequence that is covered
        // by both prefix and suffix overlaps. For instance:
        // Read    ---------------
        // OvrP --------------
        // OvrX             -----------
        //                  xx
        // In this case the function would return 2.
        // If there is no overlap, zero is returned
        // Read    ---------------
        // OvrP ---------
        // OvrX           -------------
        int calculateCoverageOverlap();
        std::string calculateConsensusFromPartition(double p_error);

        char getMODBase(const MOData& mod, int idx) const;

        // IO
        void print(int default_padding = DEFAULT_PADDING, int max_overhang = DEFAULT_MAX_OVERHANG);
        void printMasked();
        void printPileup();
        void printGroups();

    private:

        AlphaCount64 getAlphaCount(int idx) const;
        Pileup getPileup(int idx) const;
        Pileup getPileup(int idx, int numElems) const;
        Pileup getSingletonPileup(int base_idx, int ovr_idx) const;

        void getPartitionedPileup(int idx, Pileup& g0, Pileup& g1) const;
        PileupVector getPartitionedPileup(int idx, int num_parts) const;

        void printRow(int default_padding, int max_overhang, int root_len, 
                      int offset, int overlap_len, int pid, double score,
                      const std::string& seq, const std::string& id);

        // data
        static const int DEFAULT_PADDING = 20;
        static const int DEFAULT_MAX_OVERHANG = 3;

        std::string m_rootID;
        std::string m_rootSeq;
        std::string m_rootQual;

        MODVector m_overlaps;
};

#endif
