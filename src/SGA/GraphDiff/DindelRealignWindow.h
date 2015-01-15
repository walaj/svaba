//-----------------------------------------------
// Copyright 201- Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
//        and Kees Albers (caa@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// DindelRealignWindow - Infer haplotypes using the FM-index
//

#ifndef DINDELREALIGNWINDOW_H
#define DINDELREALIGNWINDOW_H
#include <fstream>
#include "Util.h"
#include "HashMap.h"
#include "multiple_alignment.h"
#include "VCFUtil.h"
#include "BWTIndexSet.h"
#include <iomanip>
#include <list>
#include <set>
#include <map>


// Enums
enum DindelWindowCandHapAlgorithm {LINEAR} ;


// Constants
const size_t DINDEL_HASH_SIZE=8;
const int DINDEL_HMM_BANDWIDTH=8;

// Event types
const int DELETION_NOVEL_SEQUENCE = -1;
const int SNP = -2;
const int MULTINUCLEOTIDE_RUN = -3;
const int INSERTION = -4;
const int LEFTOVERHANG = -5;
const int RIGHTOVERHANG = -6;
const int ADDVARIANT_DEBUG = -9;

// Typedefs
typedef std::string SampleName;

// Forward declares
class DindelRealignWindow;

// Functions
long long int combinations(const int n, const int k);

// Complement a single base
inline char _complement(char base)
{
    switch(base)
    {
        case 'A':
            return 'T';
        case'C':
            return 'G';
        case 'G':
            return 'C';
        case'T':
            return 'A';
        default:
            return 'N';
    }
}

// Add the log-scaled values l1 and l2 using a transform to avoid
// precision errors
inline double addLogs(const double l1, const double l2)
{
    if (l1>l2) {
            double diff=l2-l1;
            return l1+log(1.0+exp(diff));
    } else {
            double diff=l1-l2;
            return l2+log(1.0+exp(diff));
    }
}

inline double log_one_minus_exp_minus_x(double x)
{
    // returns log(1-10^-x)
    if (x == 0.0)
        return -HUGE_VAL;
    else if (x>0.0)
    {
        if (x>10.0)
            return -exp(-x);
        else return log(1.0- exp(-x));
    } else
    {
        assert(1==0);
    }
}

// Convert the type t into a string. Returns true if successful
template <class T>
bool from_string(T& t,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
    std::istringstream iss(s);
    return !(iss >> f >> t).fail();
}

// Convert the string s into a value of type T. Clamp
// the out value at min or max if the input string is out of the range
template <class T>
bool from_string(T& t, T min, T max,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  bool failed = !(iss >> f >> t).fail();

  if (!failed)
  {
      if (t<min) t=min;
      else if (t>max) t=max;
  }
  return failed;
}

// Split a string into tokens based on the given separator
// JS TODO: Util.h has a similar function, can deduplicate
std::vector<std::string> SplitString(const std::string & str, char sep);

// Generic wrapper for aligning two haplotypes to each other
std::string haplotypeAlignment(const std::string& h1, const std::string& h2);

// Align two haplotypes against each other semi-globally and adjust
// cigar strings to give an end-to-end alignment
std::string globalHaplotypeAlignment(const std::string& h1, const std::string& h2);

// Align two haplotypes against each other semi-globally and adjust
// cigar strings to give an end-to-end alignment
SequenceOverlap semiGlobalHaplotypeAlignment(const std::string& h1, const std::string& h2);

// Parse a colon-separated string into a chromosome/start/end triple
void parseRegionString(const std::string & region, std::string & chrom, int & start, int & end);

//
// Minimal information describing a read's alignment on the reference
// These alignments are calculated indirectly by projecting the alignment via a called haplotype.
//
class DindelReadReferenceAlignment
{
    public:
        std::string read_name;
        std::string read_sequence;
        std::string reference_name;
        int reference_start_position;
        std::string cigar;
        bool is_reference_reverse_strand;

        size_t dindel_ref_index;
        size_t dindel_read_index;
};
typedef std::vector<DindelReadReferenceAlignment> DindelReadReferenceAlignmentVector;

//
// VariantPriors - Prior probability values for various variant types
//
class VariantPriors
{
    public:
        VariantPriors();
        double getDefaultProbVariant(const std::string & type) const;
        double m_probSNP, m_probINDEL, m_probMNP;

};

//
// DindelRead - A sequence read used as input to the dindel realignment
//
class DindelRead
{
    public:

        // Constructor
        DindelRead(const SeqRecord & seqRecord, const SampleName & sampleName, double mappingQual, bool isForward);

        // Functions

        const std::string getSequence() const { return m_seqRecord.seq.toString(); }
        const std::string& getQualString() const { return m_seqRecord.qual; }
        const std::vector<double> getBaseQuals() const;
        const std::string getID() const { return m_seqRecord.id; }
        double getMappingQual() const { return m_mappingQual; } //return double(bam->core.qual); }
        double getLogProbNotMapping() const { return m_mappingQual*0.23026; }
        bool isForward() const { return m_isForward; }
        bool mateIsForward() const { assert(false); return false; } //FIXME
        
        bool isUnmapped() const {  return false; } //FIXME
        bool mateIsUnmapped() const { assert(false); return false; } //FIXME
        char getBase(int b) const { return m_seqRecord.seq.get(size_t(b)); }
        int getQual(int b) const { assert((size_t)b < m_seqRecord.qual.length()); return Quality::char2phred(m_seqRecord.qual[b]); }
        int length() const { return int(m_seqRecord.seq.length()); }
        
        //BAM const bam1_t *getBam() const { return bam; }

        void setMappingQuality(double mappingQual) { m_mappingQual=mappingQual; }
        void setRCRead(bool status) { m_rcRead=status; if (status) setupHash(); } // need to recompute hash
        bool getRCRead() const { return m_rcRead; }
        void getLogProbCorrectError(std::vector<double> & lpCorrect, std::vector<double> & lpError) const;
        const std::vector<unsigned int> & getHashKeys();
        const std::string & getSampleName() const { return m_sampleName; }
        const std::string getLibraryName() const { return std::string("NAN"); }//FIXME

    private:
        
        // Functions
        void setupHash();

        // Data
        double m_mappingQual;
        // reverse-complement the read. Intended for unmapped reads. Determined using the status of the mate.
        bool m_rcRead, m_isForward, m_setupHash; 
        SampleName m_sampleName;
        SeqRecord m_seqRecord;
        
        std::vector<unsigned int> m_hashKeys; // stores hashes for the read

};

//
// DindelVariant - A difference from the reference haplotype
//
class DindelVariant
{
    public:

        // Constructos
        DindelVariant();      
        DindelVariant(const std::string & chrom, const std::string & ref, const std::string & alt, int pos);

        // Functions
        const std::string & getChrom() const { return m_chrom; }
        const std::string & getRef() const { return m_ref; }
        const std::string & getAlt() const { return m_alt; }

        int getPos() const { return m_pos; }
        bool operator<(const DindelVariant & v) const { if (m_chrom!=v.m_chrom) return m_chrom<v.m_chrom; else if (m_pos!=v.m_pos) return m_pos<v.m_pos; else return m_id<v.m_id; }
        //bool operator<(const DindelVariant & v) const { return m_id < v.m_id; }

        // checks ref and pos against reference sequence
        //void check(Fasta *pFasta);
        // set min and max positions in haplotype where variant may be positioned without affecting the sequence
        void setHaplotypeUnique(int leftUniquePos, int rightUniquePos) { m_leftUniquePos = leftUniquePos; m_rightUniquePos = rightUniquePos; };
        int getHaplotypeLeftUnique() const;
        int getHaplotypeRightUnique() const;
        int getLengthDifference() const { return m_dlen; };
        //int getIndex() const;
        //void setIndex(int idx);
        std::string getID() const { return m_id; }
        static bool variantFromWindowString(const std::string & varstring, DindelVariant & variant);
        const std::string & getType() const { return m_type; }
        void setPriorProb(double prob) { m_priorProb = prob; }
        double getPriorProb() const { if (m_priorProb==-1.0) { std::cerr << "ERRVAR: "; write(std::cerr); } assert(m_priorProb!=-1.0); return m_priorProb; }
        void write(std::ostream & out) const;

        void setHPLen(int hplen) { m_hplen = hplen;}
        int getHPLen() const { return m_hplen; }

    private:

        // Functions
        void setup();
        void checkSequence(const std::string & seq); // checks whether sequence has no funny characters

        // Data
        std::string m_chrom;
        std::string m_ref;
        std::string m_alt;

        std::string m_id;
        std::string m_type;

        double m_priorProb; // prior probability

        int m_pos;

        // how far can indel be shifted to left and right?
        int m_leftUniquePos, m_rightUniquePos;
        
        // homopolymer length in reference at this position
        int m_hplen; 

        // Dust score in a window around the variant
        double m_dustScore;

        // does the variant change the length of the reference sequence?
        // m_dlen>0 is insertion
        int m_dlen;

        // index of variant in vector of DindelWindow
        int m_idx; 
};

//
// DindelSequenceHash
//
class DindelSequenceHash
{
    public:
        
        // Constructors
        DindelSequenceHash() {};
        DindelSequenceHash(const std::string & sequence);

        // Functions
        const std::set<int> & lookup(const std::string & seq, int pos) const;
        //int lookup(unsigned int key) const;
        const std::list<int> & lookupList(unsigned int key) const;
        
        static inline unsigned int convert(const std::string & seq, int pos)
        {
            unsigned int v=0;
            for (int x=pos, y=0;x<pos+int(DINDEL_HASH_SIZE);x++,y++) {
                v |= (map_char(seq[x]) << (2*y) );
            }
            return v;
        }
    
        static inline unsigned int pushBack(const unsigned int & key, const char & c)
        {
            return (key >> 2) | (map_char(c) << (2*(DINDEL_HASH_SIZE-1)));
        }

        static inline int map_char(const char & c)
        {
            switch (c)
            {
                case 'A': return 0;
                case 'C': return 1;
                case 'G': return 2;
                case 'T' : return 3;
            }
            return 0;
        }
 
        static void getKeys(std::vector<unsigned int> & keys, const std::string & sequence);
        void print() const;

    protected:
        typedef HashMap<unsigned int, std::list<int> > Hash;
        typedef HashMap<unsigned int, int > UniqueHash;

        // Functions
        void makeHash(const std::string & sequence);

        // Data
        Hash m_hash;
        // UniqueHash m_uniqueHash;
        std::set<int> emptySet;
        std::list<int> emptyList;
};

// DindelReferenceMapping
class DindelReferenceMapping
{
public:
    DindelReferenceMapping() { refName = ""; refSeq = ""; refStart = 0; referenceAlignmentScore = 0; isRC = false; probMapCorrect = 0.0; };
    DindelReferenceMapping(const std::string & _refName, 
                           const std::string & _refSeq, 
                           int _refStart, 
                           int _referenceAlignmentScore, 
                           bool _isRC) : refName(_refName),
                                         refSeq(_refSeq),
                                         refStart(_refStart),
                                         referenceAlignmentScore(_referenceAlignmentScore),
                                         isRC(_isRC){};
    std::string refName, refSeq;
    int refStart;
    mutable int referenceAlignmentScore;
    bool isRC;

    // probability that mapping is correct. Can only be computed by estimateMappingProbabilities()
    double probMapCorrect;
    friend bool operator<(const DindelReferenceMapping & a, const DindelReferenceMapping & b)
    {
        if (a.refName != b.refName) return a.refName<b.refName;
        else if (a.refStart != b.refStart) return a.refStart<b.refStart;
        else return a.refSeq.size()<b.refSeq.size();
    }

};

//
// DindelHaplotype
//
class DindelHaplotype
{
    friend class DindelMultiHaplotype;
    public:

        // Constructors
        DindelHaplotype(const std::string & haplotypeSequence, const DindelReferenceMapping & refMapping);
        DindelHaplotype(const DindelHaplotype & haplotype, int copyOptions);
        DindelHaplotype(const DindelHaplotype & haplotype);
        ~DindelHaplotype();

        // Functions
        bool addVariant(const DindelVariant & var);
        const std::vector<DindelVariant> & getVariants() const { return m_variants; }
        const std::string & getSequence() const { return m_seq; }
        bool isReference() const { return m_isReference; } 
        void write(std::ostream & out) const;
        int getHomopolymerLength(int b) const { return m_hplen[b]; }
        
        int getHomopolymerLengthRefPos(int refPos) const 
        { 
            int b=getHapBase(refPos); 
            if (b>=0)
                return m_hplen[b];
            else
                return -1;
        }

        int length() const { return (int) m_seq.length();}
        int getRefBase(int hapBase) const { return m_refPos[hapBase]; };
        int getHapBase(int refPos) const;
        int getRefStart() const { return m_refPos[0]; }
        bool hasVariant(const DindelVariant & variant) const ;
        int getClosestDistance(const DindelVariant & variant, int hapPos1, int hapPos2) const;
        int getClosestDistance(const DindelVariant& variant, int hapPosStartRead, int hapPosEndRead, const DindelRead & read) const;

        const DindelReferenceMapping getReferenceMapping() const { return m_refMapping; }

        // Compute the alignment of the reference sequence onto the haplotype sequence
        // This returns the haplotype as the first element of SequenceOverlap.match
        SequenceOverlap getReferenceToHaplotypeAlignment() const { return m_pMA->getAlignment(1, 0); }
        
        // As above, but with the reference sequence as the first element
        SequenceOverlap getHaplotypeToReferenceAlignment() const { return m_pMA->getAlignment(0, 1); }

    private:

        // 
        void copy(const DindelHaplotype & haplotype, int copyOptions);

        // Do not allow default constructor
        DindelHaplotype();

        // Functions

        // initializes haplotype from reference sequence
        // sets homopolymer length vector
        void alignHaplotype();
        void extractVariants();
        void determineHomopolymerLengths();
        
        // Data
        std::string m_seq;
        std::vector<int> m_hplen; // gives the homopolymer run length for any base in the haplotype sequence going left and right.
        std::vector<int> m_refPos; // position of haplotype base on reference.
        int m_firstAlignedBase;

        // list of VCF4 style variants contained in haplotype
        std::vector<DindelVariant> m_variants;
        HashMap<std::string, std::pair<int, int> > m_variant_to_pos;
  
        bool m_isReference;
        MultipleAlignment *m_pMA;
        DindelReferenceMapping m_refMapping;
};




// Haplotype aligned to multiple positions in the reference sequence

class DindelMultiHaplotype
{
public:
    DindelMultiHaplotype(const std::vector< DindelReferenceMapping > & referenceMappings, const std::string & haplotypeSequence  );
    const std::string & getSequence() const { return m_seq; }
    const std::vector<DindelVariant> & getVariants(int refIdx) const { return m_haplotypes[refIdx].getVariants(); }
    int getNumReferenceMappings() const { return int(m_referenceMappings.size()); }
    const DindelHaplotype & getSingleMappingHaplotype(int refIdx) const { assert(refIdx < (int)m_referenceMappings.size()); return m_haplotypes[refIdx]; }
    int getHomopolymerLength(int b) const { return getSingleMappingHaplotype(0).getHomopolymerLength(b); }
    int getHomopolymerLength(const std::string & chrom, int refPos) const;
    int length() const { return getSingleMappingHaplotype(0).length(); }
    bool isReference() const { return false; } //FIXME allow reference haplotypes?
    double getLogMappingProbability(int refIdx) const { return m_haplotypes[refIdx].m_refMapping.probMapCorrect; }
    
    // stores DindelVariants for each reference location the haplotype has been aligned to.
    const DindelSequenceHash & getHash() const { return m_sequenceHash; }

    const std::vector< DindelReferenceMapping > & getReferenceMappings() const { return m_referenceMappings; }

private:
    void estimateMappingProbabilities();
    std::vector< DindelHaplotype > m_haplotypes;
    std::vector< DindelReferenceMapping > m_referenceMappings;
    std::string m_seq;
    DindelSequenceHash m_sequenceHash;
};


//
// DindelWindow contains a list of variants, the candidate haplotypes generated from those variants, and the relationships
// between the candidate haplotypes and the variants.
// The candidate haplotypes are annotated with respect to the reference sequence.
// The candidate variants are specified as modifications to the reference sequence by any string.replace operation
//
class DindelWindow
{
    public:
        
        // Create window from a set of haplotypes and a reference sequence.
        // Uses MultipleAlignment to annotate the variations in the haplotypes with respect to the reference sequence.
        DindelWindow(const std::vector<std::string> & haplotypeSequences, const std::vector<DindelReferenceMapping>  & referenceMappings);
        DindelWindow(const DindelWindow & window);
        ~DindelWindow();

        const std::vector<DindelMultiHaplotype> & getHaplotypes() const { return m_haplotypes; }
        const std::vector< DindelReferenceMapping > & getReferenceMappings() const { return m_referenceMappings; }
        void addVariant(const DindelVariant& variant, bool addToAll);
        
    private:
        
        //
        void initHaplotypes(const std::vector<std::string> & haplotypeSequences, const std::vector<DindelReferenceMapping> & referenceMappings);
        void doMultipleHaplotypeAlignment();
        void copy(const DindelWindow & window);

        // Data
        
        //parameters
        DindelWindowCandHapAlgorithm m_candHapAlgorithm;
        
        // candidate haplotypes
        std::vector<DindelMultiHaplotype> m_haplotypes;
        std::map<std::string, int> m_hashAltHaps;
        std::vector< DindelReferenceMapping > m_referenceMappings;

        MultipleAlignment *m_pHaplotype_ma;
};

//
// VCFFile - Reader/Writer for VCF files
//
class VCFFile
{
    public:
        class VCFEntry
        {
            friend class VCFFile;
            public:
                
                // Typedefs
                typedef std::map<std::string, std::string> Maps;
                
                // Constructors
                VCFEntry() { m_isEmpty = true; };

                // Functions
                bool isEmpty() { return m_isEmpty; };
                std::string getInfoValue(const std::string & key) const;
                std::vector<std::string> getFilters() const;
                void parseInfoString();
                
                template<class T> 
                bool fromInfoTag(T & val, 
                                 const std::string & tag, 
                                 std::ios_base& (*f)(std::ios_base&));

                void setFilters(const std::string & str);
                void write(std::ostream & out) const;

                // Data
                std::string chrom, ref, alt, id, infoString;
                int pos, qual;
                Maps  info, filters;

            private:

                // Data
                bool m_isEmpty;
        };
        typedef std::vector<VCFEntry> VCFEntryVector;

        // Constructors
        VCFFile();
        VCFFile(const std::string & fileName, const std::string & mode);
        ~VCFFile();

        // Functions
        std::ostream & getOutputStream() { assert(m_mode=="w"); return m_outputFileHandle; }
        VCFEntry getNextEntry();
        void setSamples(const std::vector<std::string> & samples);
        const std::vector<std::string> & getSamples() const { return m_samples; }
        void outputHeader(const std::string & refFile, const std::string & paramString);
        
    private:

        // Data
        std::ofstream m_outputFileHandle;
        std::ifstream m_inputFileHandle;
        std::string m_fileName;
        bool m_isOpen;
        std::string m_mode;
        std::vector<std::string> m_samples;

};


//
// DindelRealignWindowResult
//
class DindelRealignWindowResult
{
    public:

        //
        class HaplotypeProperties
        {
        public:
            HaplotypeProperties(double _logMappingProb, double _haplotypeQual, 
                                double _freq, int _iteration) : logMappingProb(_logMappingProb), qual(_haplotypeQual), 
                                                                freq(_freq), iteration(_iteration) {};
            double logMappingProb;
            double qual;
            double freq;
            int iteration;
        };

        class Inference
        {
            public:

                // Functions
                Inference() :  strandBias(0.0), numReadsForward(0), numReadsReverse(0), 
                               numReadsForwardZeroMismatch(0), numReadsReverseZeroMismatch(0), 
                               numUnmapped(0), numLibraries(0), numReadNames(0), 
                               numRealignedReads(0), numCalledHaplotypes(0) {};
                
                // Output variants in VCF format and read alignments
                void outputAsVCF(const DindelVariant & var, 
                                 const DindelRealignWindowResult & result, 
                                 VCFCollection& out,
                                 const ReadTable* pRefTable) const;
                
                static double computeStrandBias(int numForward, int numReverse);
                
                // these functions will take care of scaling and initialization
                void addDistanceToHistogram(int distance);
                void addAlignLikToHistogram(double logLik); //
                void addMapQToHistogram(double mappingQuality);

                // Data

                // Phred-scaled posterior prob
                double qual;
                //double freq;
                double strandBias;
                int numReadsForward;
                int numReadsReverse;
                int numReadsForwardZeroMismatch;
                int numReadsReverseZeroMismatch;
                int numUnmapped;
                int numLibraries;
                int numReadNames;
                int numRealignedReads; 
                int numCalledHaplotypes;

                // these are properties that only apply to reference-based variants
                std::set<int> readIndex, haplotypeIndex;
                std::vector<HaplotypeProperties> haplotypeProperties;
                std::vector<int> histDistance; // histogram for distances to variant
                std::vector<int> histAlignLik; // histogram for read-haplotype alignment likelihoods
                std::vector<int> histMapQ; // histogram for mapping qualities
                std::string infoStr; // special info string that can be set freely
        };

        class GenotypeCall
        {
            public:

                // Constructor
                GenotypeCall() { qual=0.0; count = 0; called = false;}

                // Data
                double qual;
                int count;
                bool called; // variant was present in called pair of haplotypes
                double gl[3];
        };

        // Typedefs
        typedef std::map<DindelVariant, Inference> VarToInference;
        typedef std::map<DindelVariant, GenotypeCall> VarToGenotypeCall;
        typedef HashMap<std::string, VarToGenotypeCall> SampleToGenotypes;

        // Constructors
        DindelRealignWindowResult() : m_pDindelRealignWindow(NULL){ };
        DindelRealignWindowResult(const DindelRealignWindow & dindelRealignWindow) : m_pDindelRealignWindow(&dindelRealignWindow){ };
    
        // Functions
        void outputVCF(VCFCollection& out, const ReadTable* pRefTable);

        // Data
        std::vector<DindelHaplotype> haplotypes;
        std::vector<double> haplotypeFrequencies;
        // Size is determined by number of reference mappings in m_dindelWindow.
        // Since the haplotypes may map to multiple locations, this keeps track of which proportion of the haplotypes can be attributed to a given reference mapping location
        std::vector<double> weightedRefMapFrequencies;
        std::string outputID;
        
        // integrates same variants across different haplotypes
        VarToInference variantInference; 
        SampleToGenotypes sampleToGenotypes;
        
        // haplotype results
        typedef HashMap<int, Inference> HapIdxToInference;
        HapIdxToInference hapIdxToInference;
        std::vector<int> addOrder;
        int numHapSpecificReads;

        // number of called haplotypes
        int numCalledHaplotypes;
        bool outputGenotypes;
        
    private:
        
        // Data
        const DindelRealignWindow *m_pDindelRealignWindow;
};

//
// CandidateSNP
//
class CandidateSNP
{
    public:
        // Constructor
        CandidateSNP(char refBase)
        {
            m_refBase = refBase;
            baseToReads = BaseToReads(4);
        }

        // Funtions
        char getBase(size_t idx) const { return "ACGT"[idx]; }
        size_t getBaseIdx(char base) const 
        { 
            switch(base) 
            { 
                case 'A' : return 0; 
                case 'C' : return 1; 
                case 'G' : return 2; 
                case 'T' : return 3; 
                case 'N' : throw std::string("N not allowed"); 
            }; 
            return 0; 
        }

        
        void addRead(char refBase, char altBase, size_t readIndex)
        {
            if (m_refBase!=refBase)
            {
                std::cerr << "ERROR: m_refBase: " << m_refBase << " refBase: " << refBase << std::endl;
            }
            assert(m_refBase==refBase);
            baseToReads[getBaseIdx(altBase)].insert(readIndex);
        }

        // Data
        typedef std::vector < std::set<size_t> > BaseToReads;
        BaseToReads baseToReads;
        char m_refBase;
};

//
// MaPosToCandidateSNP
//
class MaPosToCandidateSNP : public HashMap<int, CandidateSNP >
{
    public:
        void addSNP(int readIndex, int refPos, char refBase, char altBase);
};

//
// ReadHaplotypeAlignment
//
class ReadHaplotypeAlignment
{
    public:

        // Constructors
        ReadHaplotypeAlignment() 
        { 
            logLik = 1.0; 
            nmm =-1; 
            postProbLastReadBase=-1.0; 
            isUngapped=false; 
            hapPosLastReadBase=-1;
        }

        ReadHaplotypeAlignment(double _logLik, int _nmm) 
        { 
            logLik = _logLik; 
            nmm = _nmm; 
            postProbLastReadBase=-1.0; 
            isUngapped=false; 
            hapPosLastReadBase=-1;
        }
        
        ReadHaplotypeAlignment(double _logLik, int _nmm, bool _isUngapped) 
        { 
            logLik = _logLik; 
            nmm = _nmm; 
            postProbLastReadBase=-1.0; 
            isUngapped=_isUngapped; 
            hapPosLastReadBase=-1;
        }
       
        // Data
        double logLik;
        
        // posterior probability of last read base
        double postProbLastReadBase; 
        
        // number of mismatches between read and haplotype
        int nmm; 

        // position of last readbase on haplotype. -1 indicates off haplotype (shouldn't happen)
        int hapPosLastReadBase; 

        bool isUngapped;
};

//
// DindelRealignParameters
//
class DindelRealignParameters
{
    friend class DindelRealignAlgorithm;
    friend class DindelRealignWindow; 

    public:

        // Constructors
        DindelRealignParameters();
        DindelRealignParameters(const std::string & paramString);
        DindelRealignParameters(const char * paramString);

        // Functions
        std::string getDefaultParameters() const;
        void initFromString(const std::string & paramString);
        const std::string & getParamString() const { return m_paramString; }
        void checkAndInit();

        // Data
        std::string m_paramString;
        std::string excludeSamplesFile;

        int windowReadBuffer;
        int minVariantSep, haplotypeWidth, minCandidateAlleleCount;
        int maxNumReads; // maximum number for fetch
        int maxNumReadsWindow; // maximum number of reads in window
        int maxNumCandidatesPerWindow;
        double theta_snp;
        double theta_indel;
        double logPriorAddHaplotype;

        VariantPriors variantPriors;
        int maxMappingQuality;

        int addSNPMaxMismatches; // maximum number of mismatches allowed in ungapped alignment for adding candidate SNPs. Note this would include the SNP!
        double addSNPMinMappingQual; // minimum mapping quality of read for adding candidate SNPs
        int addSNPMinBaseQual; // minimum base quality
        size_t addSNPMinCount;
        int addSNPMaxSNPs;
        double minPostProbLastReadBaseForUngapped; // minimum posterior probability of the last base in order to attempt an ungapped alignment for SNP detection

        int genotyping; // genotyping mode? If 1, then each read will be aligned against each haplotype using the HMM

        int showCallReads; // show reads driving the call
        int minNumHaplotypeOverlaps; // in filterHaplotypes, how many reads should overlap the sequence surrounding the variant?

        // derived parameters
        double minLogLikAlignToRef, minLogLikAlignToAlt; // minimum logLik for aligning to a reference and an alternative haplotype
        double addSNPMinLogBaseQual;

        double singleSampleHetThreshold, singleSampleHomThreshold;

        int doEM;
        double EMtol;
        int EMmaxiter;

        int realignMatePairs; // compute a joint likelihood for the alignment of a read pair to the candidate haplotype
        int multiSample;
        int graphDiffStyle; // determines way haplotypes are called in the tumour/child sample
};

// Realigns reads in a given window.
// pOverlapper and dindelReads should correspond 1-1. It is assumed that pOverlapper gives back indices into the dindelReads array.
class DindelRealignWindow
{
    friend class DindelRealignWindowResult;

    public:
        
        // Constructors
        DindelRealignWindow(const DindelWindow *pDindelWindow, 
                            std::vector<DindelRead> & dindelReads, 
                            const DindelRealignParameters & realignParameters);
        
        // Functions
        void run(const std::string & algorithm, std::ostream& out);
        void run(const std::string & algorithm,
                 VCFCollection& out,
                 DindelReadReferenceAlignmentVector* pOutAlignments, 
                 const std::string id,
                 DindelRealignWindowResult * pThisResult,
                 const DindelRealignWindowResult * pPreviousResult,
                 const ReadTable* pRefTable);

        const DindelWindow & getDindelWindow() const { return  m_dindelWindow; }
        const std::vector< std::vector< ReadHaplotypeAlignment > > & getHapReadAlignments() const 
        { 
            return hapReadAlignments; 
        }

    private:

        // Data
        DindelWindow m_dindelWindow;
        std::vector<DindelRead> *m_pDindelReads;
        DindelRealignParameters realignParameters;

        // A set of variant indices that assigned to a called haplotype containing a variant
        std::set<size_t> m_variantReadIndices;

        DindelReadReferenceAlignmentVector m_readReferenceAlignments;

        std::string m_outputID; // to be output in VCF

        // HAPLOTYPE ALIGNMENT BUSINESS
        MaPosToCandidateSNP m_maPosToCandidateSNP;

        // Functions
        void filterHaplotypes();
        void addSNPsToHaplotypes();
        void processHaplotypes(size_t firstHapIdx, size_t lastHapIdx);

        // convert seed positions into a likelihood given the candidate haplotype
        ReadHaplotypeAlignment computeReadHaplotypeAlignment(size_t readIndex, 
                                                             const DindelMultiHaplotype& haplotype, 
                                                             int start, 
                                                             int end, 
                                                             const std::vector<double> & lpError, 
                                                             const std::vector<double> & lpCorrect, 
                                                             bool rcReadSeq);

        void HMMAlignReadAgainstHaplotypes(size_t readIndex, 
                                           size_t firstHap, 
                                           size_t lastHap, 
                                           const std::vector<double> & lpCorrect, 
                                           const std::vector<double> & lpError,
                                           HashMap<std::string, ReadHaplotypeAlignment>& hmm_alignment_cache);

        // HAPLOTYPE FREQUENCY ESTIMATION BUSINESS

        // haplotype frequencies
        class HL
        {
            public:
                HL(double _ll, int _idx) { ll = _ll; idx = _idx;}

                bool operator<(const HL & hl) const 
                { 
                    if (int(ll*10.0)==int(hl.ll*10.0))
                    {
                        if (idx<hl.idx) 
                            return true; 
                        else 
                            return false;
                    } 
                    else
                    {
                        return int(ll*10.0)<int(hl.ll*10.0);
                    }
                }

                double ll;
                int idx;
        };

        // Functions
        void printReadAlignments(int readIdx, std::ostream & out, int offset, bool supportAlt);

        //
        void doEM(const std::vector< std::vector<double> > & hrLik, 
                  const std::vector<int> & calledHaplotypes, 
                  std::vector<double> & haplotypeFrequencies);

        void doEMMultiSample(int numSamples,
                             int maxIter,
                             const std::vector<double> & llHapPairs,
                             const std::vector<int> & allowedHaplotypes,
                             const std::vector<double> & initHaplotypeFrequencies,
                             std::vector<double> & haplotypeFrequencies);

        void addCalledHaplotypeMatePairs(int hapIdx,
                                DindelRealignWindowResult & result,
                                const std::vector<DindelMultiHaplotype> & haplotypes,
                                const std::vector< HashMap<int, double> > & addReads,
                                double minLogLikAlignToAlt,
                                int numReadPairs);
        
        void addCalledHaplotypeSingleRead(int hapIdx,
                                DindelRealignWindowResult & result,
                                const std::vector<DindelMultiHaplotype> & haplotypes,
                                const std::vector< HashMap<int, double> > & addReads,
                                double minLogLikAlignToAlt,
                                int numReads);

        // Project the alignment of a single read from its position on a haplotype to the reference genome.
        void projectReadAlignmentToReference(const std::vector<DindelMultiHaplotype> & haplotypes, 
                                             int readIdx, int hapIdx, int refIdx);

        // Compute the mapping quality of each read onto the reference
        void computeProjectedMappingQuality();

        void setAddReadsMatePairs(int type,
                         int h,
                         const std::vector< std::vector<double> > & hrLik,
                         std::vector< HashMap<int, double> > & addReads,
                         std::vector<double> & addLL,
                         int numReadPairs,
                         int numHaps,
                         const std::vector<int> & added);

        void setAddReadsSingleRead(int type,
                         int h,
                         const std::vector< std::vector<double> > & hrLik,
                         std::vector< HashMap<int, double> > & addReads,
                         std::vector<double> & addLL,
                         int numReads,
                         int numHaps,
                         const std::vector<int> & added);

        void setAddReadsSingleReadMultiSample(int type,
                                      int htest,
                                      const std::vector< std::vector<double> > & hrLik,
                                      std::vector< HashMap<int, double> > & addReads,
                                      std::vector<double> & newLL,
                                      int numReads,
                                      int numHaps,
                                      std::vector<double> & zindNew,
                                      std::vector<double> & hapFreqNew,
                                      const std::vector<int> & added,
                                      const std::vector<double> & hapFreqPrevious,
                                      const std::vector<double> & zindPrevious,
                                      const std::vector<double> & llHapPairs,
                                      const HashMap<std::string, std::list<int> > & sampleToReads);

        void computeAddLLMatePairs(int type,
                          int h,
                          const std::vector< std::vector<double> > & hrLik,
                          std::vector<double> & addLL,
                          int numReadPairs,
                          int numHaps,
                          const std::vector<int> & added);

         void computeAddLLSingleRead(int type,
                          int h,
                          const std::vector< std::vector<double> > & hrLik,
                          std::vector<double> & addLL,
                          int numReads,
                          int numHaps,
                          const std::vector<int> & added);
        void computeHaplotypePairLikelihoods(std::vector<double> & llHapPairs,
                                             HashMap<std::string, std::list<int> > & sampleToReads,
                                             const std::vector< std::vector<double> > & hrLik);
        
        void doReadHaplotypeAlignment(int H, const std::vector<DindelRead> & dReads);
        void showHaplotypeOnlyReadsMatePairs(int h,
                                    const std::vector< std::vector<double> > & hrLik,
                                    int numReadPairs,
                                    int numHaps);

        void showHaplotypeOnlyReadsSingleRead(int h,
                                    const std::vector< std::vector<double> > & hrLik,
                                    int numReads,
                                    int numHaps);
                                
        DindelRealignWindowResult estimateHaplotypeFrequencies(double minLogLikAlignToRef, 
                                                               double minLogLikAlignToAlt, 
                                                               bool capUsingMappingQuality, 
                                                               bool print);
        
        DindelRealignWindowResult estimateHaplotypeFrequenciesModelSelection(double minLogLikAlignToRef, 
                                                                             double minLogLikAlignToAlt, 
                                                                             bool capUsingMappingQuality, 
                                                                             bool print);
        DindelRealignWindowResult estimateHaplotypeFrequenciesModelSelectionMatePairs(double minLogLikAlignToRef, 
                                                                                      double minLogLikAlignToAlt, 
                                                                                      bool capUsingMappingQuality,
                                                                                      const DindelRealignWindowResult * pPreviousResult,
                                                                                      bool print);
    
        DindelRealignWindowResult estimateHaplotypeFrequenciesModelSelectionSingleReads(double minLogLikAlignToRef,
                                                                                        double minLogLikAlignToAlt,
                                                                                        bool capUsingMappingQuality,
                                                                                        const DindelRealignWindowResult * pPreviousResult,
                                                                                        bool print);

        DindelRealignWindowResult estimateHaplotypeFrequenciesModelSelectionSingleReadsMultiSample(double minLogLikAlignToRef,
                                                                                        double minLogLikAlignToAlt,
                                                                                        bool capUsingMappingQuality,
                                                                                        const DindelRealignWindowResult * pPreviousResult,
                                                                                        bool print);

        //
        std::vector< std::vector<double> > fillReadHaplotypeLikelihoods() const;

        double getHaplotypePrior(const DindelHaplotype & h1, const DindelHaplotype & h2) const;

        // this is a function that hides how the read information is actually stored in the various ReadTables
        const DindelRead & getRead(size_t readIndex) const { return m_pDindelReads->at(readIndex); }

        // haplotype-read likelihoods
        
        // first dim is haplotypes, second dim is reads
        std::vector< std::vector< ReadHaplotypeAlignment > > hapReadAlignments; 

        // update hapReadAlignments.
        void computeReadHaplotypeAlignmentsUsingHMM(size_t firstHap, size_t lastHap);
        void addDiploidGenotypes(DindelRealignWindowResult & result, bool useEstimatedHaplotypeFrequencies);
        void addDiploidGenotypes(DindelRealignWindowResult& result, const std::vector<int> & allowedHaplotype, const std::vector< std::vector<double> > & hrLik);
        
        // algorithms
        void algorithm_hmm(VCFCollection& out,
                           DindelReadReferenceAlignmentVector* pOutAlignments,
                           DindelRealignWindowResult * pThisResult,
                           const DindelRealignWindowResult * pPreviousResult,
                           const ReadTable* pRefTable);

        // result
        DindelRealignWindowResult m_result;
    
};

//
// WindowFile
//
class WindowFile
{

    public:
        
        // Constructors
        WindowFile(const std::string & fileName);
        ~WindowFile();

        // Functions
        std::vector<DindelVariant> getNextWindow();
        const std::string & getFilename() const { return m_fileName; }
        int getCurrentLineIndex() const { return m_index; };

    private:

        std::ifstream m_fileHandle;
        bool m_isOpen;
        int m_index;
        std::string m_fileName;
        std::vector<DindelVariant> m_empty;

};

// Add an insertion/deletion into a haplotype
DindelVariant addDeletion(const DindelHaplotype & refHap, const std::string & chrom, 
                          int refPos, int len);
DindelVariant addInsertion(const DindelHaplotype & refHap, const std::string & chrom, 
                           int refPos, const std::string & seq);

//
// CandidateIndelFreq
//
class CandidateIndelFreq
{
    public:

        // Constructors
        CandidateIndelFreq() {};
        CandidateIndelFreq(const std::string & sample, 
                           const std::string & library, 
                           const std::string & fragment) 
        { 
            addObservation(sample, library, fragment); 
        };

        // Functions
        void addObservation(const std::string & sample, const std::string & library, const std::string & fragment);
        int getNumSamples() const { return int(m_samples.size()); };    
        void getStats(int & numLibs, int & numFragments, int & maxFragPerSample) const;

    private:
        class SampleIndels
        {
            public:
                SampleIndels() {};
                void addObservation(const std::string & library, const std::string & fragment);
                typedef HashMap<std::string, int> FragmentToFreq;
                typedef HashMap<std::string, FragmentToFreq > LibToFragment;
                LibToFragment libToFragment;
        };

        HashMap<std::string, SampleIndels> m_samples;
};

typedef std::map<DindelVariant, CandidateIndelFreq> DindelVariantToFreqMap;
typedef std::map<int, DindelVariantToFreqMap > CandidateIndelMap;

#endif  /* DINDELREALIGNWINDOW_H */
