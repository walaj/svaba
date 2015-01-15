//-----------------------------------------------
// Copyright 201- Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
//        and Kees Albers (caa@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// DindelRealignWindow - Infer haplotypes using the FM-index
//

/*
 TODO
 * 1. Fix addSNP. Is currently still based on single reference position. Needs to be changed to position in haplotype multiple alignment.
 * 2. Add read variant coverage statistics for existing variants but supported by reads mapping to new haplotype
      Maybe first store only supporting read idxs and then add statistics after all haplotypes have been called.
 * 3. Store for each variant the vectors [hapqual], [hapfreq], [hap_mapqual] for the haplotypes in which the variant is present
 * 4. Use reference sequence to extract reads? Add reference sequence as haplotype? Maybe that is necessary for accurate genotyping?


 */
#include <string>
#include <sstream>
#include <algorithm>
#include <map>
#include <set>
#include <numeric>
#include <iomanip>
#include "StdAlnTools.h"
#include "MultiAlignment.h"
#include "DindelRealignWindow.h"
#include "BWT.h"
#include "DindelHMM.h"
#include "Profiler.h"
#include "VCFUtil.h"
#include "overlapper.h"
#include "multiple_alignment.h"
#include "HapgenUtil.h"
#include "Quality.h"
#include <cmath>

const int DINDEL_DEBUG=0;
const int QUIET=1;
const int ADDSNPS=0;
const int ALWAYS_REALIGN=1;
const int DEBUG_CALLINDEL=0;
const int REPOSITION_INDEL_WINDOW=1000;
const int SHOWHAPFREQ=0;
const int REPOSITIONVARIANTSSLOW=0; // uses slow code to reposition indels.
const int DINDEL_DEBUG_3 = 0; // useful debugging
const int DINDEL_ADJUST_MAPPINGQUAL = 1;
//#define OVERLAPPER // build overlapper

// used in DindelHaplotype::extractVariants. Variants within this distance from the edge of the reference haplotype will not be annotated
const int MIN_EVENT_DIST_FROM_EDGE = 10;


#include <iostream>
#include <vector>
//#define DEBUG_NEW

long long int combinations(const int n, const int k)
{
    if (n == k || k == 0)
        return 1;
    if (k == 1)
        return n;

    if (n < k)
        return -1;
    else
    {
        int i;
        long long int result, f_n, f_k, f_n_sub_k;

        f_n = 1;
        for (i = 2; i <= n; i++)
            f_n *= i;

        f_k = 1;
        for (i = 2; i <= k; i++)
            f_k *= i;

        f_n_sub_k = 1;
        i = 2;
        while (i <= (n - k))
        {
            f_n_sub_k *= i;
            i++;
        }

        result = f_n / (f_k * f_n_sub_k);
        return result;
    }
}


void parseRegionString(const std::string & region, std::string & chrom, int & start, int & end)
{
        std::string filtered;
        for(size_t x=0;x<region.size();x++) {
            char c=region[x];
            if (c=='-' || c==':') filtered+=' ';
            else filtered+=c;
        }
        std::istringstream is(filtered);
        is >> chrom;
        std::string e; is >> e;
        if (!from_string(start,e,std::dec)) throw std::string("Cannot parse region for start!");
        is >> e;
        if (!from_string(end,e,std::dec)) throw std::string("Cannot parse region end!");

        if (end<start)
        {
            std::string message("Invalid region: start>end : ");
            message+=region;
            throw message;
        }
}

//
std::vector<std::string> SplitString(const std::string & str, char sep)
{
    std::string elem;
    std::vector<std::string> split;

    for (size_t x=0;x<str.size();x++)
    {
        if (str[x]==sep)
        {
            if (!elem.empty()) split.push_back(elem);
            elem = "";
        } else elem += str[x];

    }
    if (!elem.empty()) split.push_back(elem);
    return split;
}

// Convert a local alignment structure into a SequenceOverlap object
SequenceOverlap localAlignmentToSequenceOverlap(const LocalAlignmentResult& local_aln)
{
    SequenceOverlap overlap;
    overlap.match[0].start = local_aln.targetStartIndex;
    overlap.match[0].end = local_aln.targetEndIndex;
    overlap.match[1].start = local_aln.queryStartIndex;
    overlap.match[1].end = local_aln.queryEndIndex;
    overlap.cigar = local_aln.cigar;
    return overlap;
}

SequenceOverlap computeHaplotypeAlignment(const std::string& base_sequence, const std::string& haplotype)
{
    return HapgenUtil::alignHaplotypeToReference(base_sequence, haplotype);
}

void addHaplotypeToMultipleAlignment(MultipleAlignment* pMA,
                                     const std::string& base_sequence, 
                                     const std::string& incoming_name, 
                                     const std::string& incoming_sequence)
{
    SequenceOverlap overlap = computeHaplotypeAlignment(base_sequence, incoming_sequence);
    pMA->addOverlap(incoming_name, incoming_sequence, "", overlap);
}

/*
 *
 * DINDELVARIANT
 *
 *
 */

DindelRead::DindelRead(const SeqRecord & seqRecord, const SampleName & sampleName, double mappingQual, bool isForward) : m_seqRecord(seqRecord)
{
    m_mappingQual = mappingQual;
    m_rcRead = false;
    m_sampleName = sampleName;
    m_setupHash = false;
    m_isForward = isForward;
}

void DindelRead::getLogProbCorrectError(std::vector<double>& lpCorrect, std::vector<double>& lpError) const
{
    lpError=std::vector<double>(this->length());
    lpCorrect=std::vector<double>(this->length());
    for (int b=0;b<this->length();b++)
    {
        double pr=1.0 - exp ( (-2.3026/10.0)*double(this->getQual(b) ));
        double eq=log(.25+.75*pr);
        double uq=log(.75+1e-10-.75*pr);
        lpError[b] = uq;
        lpCorrect[b] = eq;
    }
}


void DindelRead::setupHash()
{
    if (!m_rcRead)
        DindelSequenceHash::getKeys(m_hashKeys, this->getSequence());
    else
    {
        std::string seq = this->getSequence();
        std::string rcSeq=seq;
        size_t rlen=rcSeq.size();
        for (size_t x=0;x<rcSeq.size();x++) rcSeq[rlen-x-1]=_complement(seq[x]);
        DindelSequenceHash::getKeys(m_hashKeys, rcSeq);
    }
    m_setupHash = true;
}

void DindelSequenceHash::getKeys(std::vector<unsigned int>& keys, const std::string& sequence)
{

    if(sequence.size() >= DINDEL_HASH_SIZE)
    {
        keys.clear();
        keys.reserve(sequence.size()-DINDEL_HASH_SIZE);
        unsigned int key=DindelSequenceHash::convert(sequence, 0);
        keys.push_back( key );
        for (size_t x=1;x<sequence.size()-DINDEL_HASH_SIZE;x++)
        {
            key = DindelSequenceHash::pushBack(key, sequence[x+DINDEL_HASH_SIZE]);
            keys.push_back( key );
        }
    }    
}

const std::vector<unsigned int> & DindelRead::getHashKeys()
{
    if (!m_setupHash) setupHash();
    return m_hashKeys;
}

VariantPriors::VariantPriors()
{
    m_probSNP = 0.001;
    m_probINDEL = 0.0001;
    m_probMNP = 0.00001;
}


double VariantPriors::getDefaultProbVariant(const std::string& type) const
{
    if (type == "SNP") return m_probSNP;
    else if (type == "INDEL") return m_probINDEL;
    else if (type == "MNP") return m_probMNP;
    else assert(1==0); // ("Unknown variant type");
}


DindelVariant::DindelVariant()
{
    m_pos = -1;
    m_type = "uninitialized";
}

DindelVariant::DindelVariant(const std::string & chrom, const std::string & ref, const std::string & alt, int pos)
{
    if (ref.empty() || alt.empty() || chrom.empty()) throw std::string("DindelVariant::zero_length_chrom_ref_alt_string");
    //if (pos<0) throw std::string("DindelVariant::negative_variant_position");
    assert(pos>=0);
    m_chrom = chrom;
    m_ref = ref;
    m_alt = alt;
    m_pos = pos;

    m_hplen = -1;
    m_leftUniquePos = -1;
    m_rightUniquePos = -1;
    m_idx = -1;
    m_priorProb=-1;
    
    setup();
}


bool DindelVariant::variantFromWindowString(const std::string & varstring, DindelVariant & variant)
{
    // split string

    std::vector<std::string> vs;
    size_t lastPos = 0;
    for (size_t x = 0; x < varstring.size(); x++)
    {
        if (varstring[x]==',')
        {
            if (x>0)
            {
                vs.push_back(varstring.substr(lastPos, x-lastPos));
                lastPos=x+1;
            }
        }
    }
    if (lastPos<varstring.size())
    {
        vs.push_back(varstring.substr(lastPos, varstring.size()-lastPos));
    }

    if (vs.size()<3) return false;

    int pos;
    std::string chrom, ref, alt;
   
    chrom = vs[0]; if (chrom.empty()) return false;
    if (!from_string<int>(pos, vs[1], std::dec)) return false;
    if (pos<0) return false;

    ref = vs[2]; if (ref.empty()) return false;
    alt = vs[3]; if (alt.empty()) return false;

    for (size_t x=0;x<ref.size();x++) if (ref[x]!='A' && ref[x]!='T' && ref[x] != 'C' && ref[x]!='G' && ref[x] != 'N') return false;
    for (size_t x=0;x<alt.size();x++) if (alt[x]!='A' && alt[x]!='T' && alt[x] != 'C' && alt[x]!='G' && alt[x] != 'N') return false;

    variant = DindelVariant(chrom, ref, alt, pos);

    return true;
}


void DindelVariant::checkSequence(const std::string & seq)
{
    for (size_t x=0;x<seq.size();x++) if (seq[x]!='A' && seq[x] !='C' && seq[x] != 'G' && seq[x] != 'T' && seq[x] != 'N') throw std::string("Variant: sequence "+seq+" not allowed.");

}

void DindelVariant::setup()
{
    // check if there are no commas in alt or ref

    checkSequence(m_alt);
    checkSequence(m_ref);

    m_dlen = m_alt.length()-m_ref.length();

    if (m_dlen==0 && m_alt.size()==1) m_type = "SNP";
    else if (m_dlen==0 && m_alt.size()>1) m_type = "MNP";
    else if (m_dlen!=0) m_type = "INDEL";

    std::stringstream posStr;
    posStr << m_pos;

    m_id = m_chrom+"_"+posStr.str()+"_"+m_ref+"_"+m_alt;
}

int DindelVariant::getHaplotypeLeftUnique() const
{
    if(m_leftUniquePos==-1)
         throw std::string("DindelVariant::call_determineLeftRightUnique_first");
    return m_leftUniquePos;
}

int DindelVariant::getHaplotypeRightUnique() const
{
    if(m_rightUniquePos==-1)
         throw std::string("DindelVariant::call_determineLeftRightUnique_first");
    return m_rightUniquePos;
}

void DindelVariant::write(std::ostream & out) const
{
    out << m_id << " " << m_chrom << " " << m_ref << " " << m_alt << " type: " << m_type << " priorProb: " << m_priorProb << std::endl;
}

/*
 *
 * DINDELHAPLOTYPE
 *
 *
 */


void DindelHaplotype::alignHaplotype()
{
    // semi-globally align haplotypes to the first haplotype (arbitrary)
    const std::string& rootSequence = m_refMapping.refSeq;
    std::string alignSeq;
    if(m_refMapping.isRC)
        alignSeq = reverseComplement(m_seq);
    else
        alignSeq = m_seq;
    
    // The reference sequence is the base sequence of the multiple alignment
    std::stringstream ref_ss;
    ref_ss << "ref=" << m_refMapping.refName << ":" << m_refMapping.refStart;

    m_pMA = new MultipleAlignment();
    m_pMA->addBaseSequence(ref_ss.str(), rootSequence, "");
    SequenceOverlap hap_overlap = computeHaplotypeAlignment(rootSequence, alignSeq);
    m_pMA->addOverlap("haplotype", alignSeq, "", hap_overlap);
    m_firstAlignedBase = hap_overlap.match[1].start;

    if (DINDEL_DEBUG || DINDEL_DEBUG_3)
    {
        std::cout << "DindelHaplotype::alignHaplotype:\n";
        m_pMA->print(500);
        //std::cout << "DindelHaplotype::DindelHaplotype globalAlignmentCigar " << alignSeq << " vs root: " << _ma.expandedCigar << std::endl;
    }
}

void DindelHaplotype::extractVariants()
{
    const int verbose = 0;
    // extract variants from alignment
    int numCols = m_pMA->getNumColumns();

    int refRow = 0;
    int varRow = 1;

    // position of haplotype base on reference.
    m_refPos = std::vector<int>(m_seq.size(), 0);
    m_isReference = false;


    // Compute a mapping from a haplotype base to its aligned position on the reference sequence
    int first_aligned_column = m_pMA->getFirstAlignedBaseColumn(refRow, varRow);
    int last_aligned_column = m_pMA->getLastAlignedBaseColumn(refRow, varRow);
    
    // This calculation assumes the haplotype was completely aligned!
    assert(m_pMA->columnIndexToBaseIndex(varRow, first_aligned_column) == 0);
    int first_aligned_haplotype_base = m_firstAlignedBase + m_pMA->columnIndexToBaseIndex(varRow, first_aligned_column);

    // Fill in the lefthand unaligned haplotype sequence
    for(int i = 0; i < first_aligned_haplotype_base; ++i)
    {
        assert(i < (int)m_refPos.size());
        m_refPos[i] = LEFTOVERHANG;
    }

    for(int i = first_aligned_column; i <= last_aligned_column; ++i)
    {
        char rs = m_pMA->getSymbol(refRow,i);
        char vs = m_pMA->getSymbol(varRow,i);

        // Skip if this is a deleted haplotype base - it has no reference coordinate
        if(vs == '-')
            continue;

        assert(rs != '\0' && vs != '\0');
        size_t hap_idx = m_firstAlignedBase + m_pMA->columnIndexToBaseIndex(varRow, i);
        assert(hap_idx < m_refPos.size());

        if (rs == '-' && vs != '-') 
            m_refPos[hap_idx] = INSERTION;

        if (rs != '-' && vs != '-')
        {
            if(rs != vs)
                m_refPos[hap_idx] = SNP;
            else
                m_refPos[hap_idx] = m_refMapping.refStart + m_pMA->columnIndexToBaseIndex(refRow, i);
        }
    }

    // Fill in the righthand unaligned haplotype sequence
    int last_aligned_haplotype_base = m_firstAlignedBase + m_pMA->columnIndexToBaseIndex(varRow, last_aligned_column);
    for(int i = last_aligned_haplotype_base; i < (int)m_refPos.size(); ++i)
    {
        assert(i < (int)m_refPos.size());
        m_refPos[i] = RIGHTOVERHANG;
    }
    
    if(DINDEL_DEBUG || verbose)
    {
        std::cout << "m_refPos: ";
        for (int i = 0; i < (int)m_refPos.size(); i++) 
            std::cout << "[" << i << " " << m_refPos[i] << "]";
        std::cout << "\n";
    }

    int leftExactMatch = 0;
    int rightExactMatch = 0;

    int eventStart = -1;
    int eventEnd = -1;
    bool isIndel = false;
    bool inLeftExact = true;

    for(int i = first_aligned_column; i <= last_aligned_column; ++i)
    {
        char refSymbol = m_pMA->getSymbol(refRow, i);
        char varSymbol = m_pMA->getSymbol(varRow, i);

        // JS 25/10/11 Keep considering the event to be a variant if the var or ref sequence has a gap here
        bool isVariant = (varSymbol != refSymbol || varSymbol == '-' || refSymbol == '-');

        // Update the counter of the number of exact matches for the leftmost and rightmost bases
        if(isVariant)
        {
            inLeftExact = false; // stop the count of leftmost exact matches
            rightExactMatch  = 0; // reset the rightmost count
        }
        else // this is a good match
        {
            rightExactMatch += 1;
            if(inLeftExact)
                leftExactMatch += 1;
        }

        // Process variants indicated by this portion of the alignment
        if(isVariant)
        {
            if(eventStart == -1)
            {
                // Not currently in a variant event, start a new one
                eventStart = i;
                eventEnd = i;
            }
            else
            {
                // Extend the event
                // NOTE that a SNP next to an indel will be annotated as a single event.
                assert(eventEnd != -1);
                eventEnd += 1;
            }

            // event will be an indel even if it contains single nucleotide mismatches to the reference
            isIndel = isIndel || varSymbol == '-' || refSymbol == '-';
        }
        else
        {
            // Check if this is the end of a variant
            if(eventStart != -1)
            {
                // Extract the substrings of the three sequences describing
                // this event.
                // If the event is an indel, we move the start point
                // of the event backwards until all 3 columns are not padding characters.
                // This is to avoid the following situation:
                // AATATTT--CGT ref
                // AATATTTT-CGT base
                // AATATTTTTCGT var
                // or
                // TTAGGTTTTTTT ref
                // TTAGG-TTTTTT base
                // TTAGG--TTTTT var
                //
                // In the first case, it is a 1 base insertion wrt base but a 2 base insertion wrt to the ref.
                // So we need to move the event to the first column with all Ts to properly capture what is happening.
                // In the second case it is a 1 base deletion wrt base and a 2 base deletion wrt to the ref.
                // Again, we need to move back to the last full column.

                bool okEvent = true;
                if(isIndel && (eventStart == 0 || eventEnd == (int)numCols - 1) )
                    okEvent = false;
                
                int left_distance_to_end = eventStart - first_aligned_column;
                int right_distance_to_end = last_aligned_column - eventEnd;
                if(left_distance_to_end < MIN_EVENT_DIST_FROM_EDGE || right_distance_to_end < MIN_EVENT_DIST_FROM_EDGE) 
                    okEvent = false;
                
                if(okEvent)
                {
                    if(isIndel)
                    {
                        // Move the event leftwards until the reference and haplotype both have a non-gap symbol
                        while(eventStart >= (int)first_aligned_column)
                        {
                            char rc = m_pMA->getSymbol(refRow, eventStart);
                            char vc = m_pMA->getSymbol(varRow, eventStart);
                            if(rc == '-' || vc == '-')
                                eventStart -= 1;
                            else
                                break;
                        }
                        assert(eventStart >= (int)first_aligned_column);
                    }

                    int eventLength = eventEnd - eventStart + 1;
                    std::string refStringPadded = m_pMA->getPaddedSubstr(refRow, eventStart, eventLength);
                    std::string varStringPadded = m_pMA->getPaddedSubstr(varRow, eventStart, eventLength);
                    std::string refString = StdAlnTools::unpad(refStringPadded);
                    std::string varString = StdAlnTools::unpad(varStringPadded);
                    if(verbose > 0)
                    {
                        printf("Ref start: %d col: %d eventStart: %d eventEnd: %d\n", (int) m_refMapping.refStart, (int)i, eventStart, eventEnd);
                        std::cout << "RefString " << refString << "\n";
                        std::cout << "varString " << varString << "\n";
                    }

                    assert(!refString.empty());
                    assert(!varString.empty());

                    // minPos and maxPos are coordinates of window in MultiAlignment where variant can be unambiguously positioned
                    int minPos = eventStart;
                    int maxPos = eventStart + eventLength - 1;

                    // Get the base position in the reference string. This is not necessarily the
                    // same as the eventStart column as the reference may be padded
                    // get positions in haplotype where indel variant may be ambiguously positioned
                    int refBaseOffsetMinPos = m_pMA->columnIndexToBaseIndex(refRow, minPos);
                    
                    const static int HP_COUNT_TOWARDS_END = -1;
                    const static int HP_COUNT_TOWARDS_START = 0;
                    size_t hplen_before = m_pMA->countHomopolymer(refRow, minPos - 1, HP_COUNT_TOWARDS_END) + m_pMA->countHomopolymer(refRow, minPos - 1, HP_COUNT_TOWARDS_START) - 1;
                    size_t hplen_at = m_pMA->countHomopolymer(refRow, minPos, HP_COUNT_TOWARDS_END) + m_pMA->countHomopolymer(refRow, minPos, HP_COUNT_TOWARDS_START) - 1;
                    size_t hplen_after = m_pMA->countHomopolymer(refRow, minPos + 1, HP_COUNT_TOWARDS_END) + m_pMA->countHomopolymer(refRow, minPos + 1, HP_COUNT_TOWARDS_START) - 1;
                    size_t hplen = std::max(std::max(hplen_before, hplen_at), hplen_after);

                    while( minPos > 0 && m_pMA->getSymbol(varRow, minPos)== '-' )
                        minPos--;
                    while( maxPos < numCols - 1 && m_pMA->getSymbol(varRow, maxPos)== '-' )
                        maxPos++;

                    int varBaseOffsetMinPos = m_pMA->columnIndexToBaseIndex(varRow, minPos);
                    int varBaseOffsetMaxPos = m_pMA->columnIndexToBaseIndex(varRow, maxPos);

                    // if the haplotype was aligned to the reverse strand of the reference sequenced, then we need to reverse these coordinates back.
                    if(m_refMapping.isRC)
                    {
                        int tmp = varBaseOffsetMaxPos;
                        varBaseOffsetMaxPos = int(m_seq.length()) - 1 - varBaseOffsetMinPos;
                        varBaseOffsetMinPos = int(m_seq.length()) - 1 - tmp;
                        //assert(varBaseOffsetMinPos<=varBaseOffsetMaxPos);

                        if(!(varBaseOffsetMaxPos<int(m_seq.length())))
                        {
                            std::cout << "varBaseOffsetMinPos: " << varBaseOffsetMinPos 
                                      << " varBaseOffsetMaxPos: " << varBaseOffsetMaxPos 
                                      << " length: " << m_seq.length() 
                                      << " minPos: " << minPos 
                                      << " maxPos: " << maxPos 
                                      << " ma.columnIndexToBaseIndex(varRow, minPos):" << m_pMA->columnIndexToBaseIndex(varRow, minPos) 
                                      << " ma.columnIndexToBaseIndex(varRow, maxPos); " << m_pMA->columnIndexToBaseIndex(varRow, maxPos) << "\n";
                            m_pMA->print(1000);
                        }

                        assert(varBaseOffsetMinPos>=0);
                        assert(varBaseOffsetMaxPos<int(m_seq.length()));
                    }


                    // Use leftmost position for the variant (which can only be change for an insertion or deletion)
                    if (DINDEL_DEBUG || 0) 
                    {
                        std::cout << "VARIANT: " << m_refMapping.refName << " " 
                                  << refString << "/" << varString 
                                  << " pos: " << refBaseOffsetMinPos+m_refMapping.refStart 
                                  << " varBaseOffsetMinPos: " << varBaseOffsetMinPos 
                                  << " varBaseOffsetMaxPos: " << varBaseOffsetMaxPos 
                                  << " m_refMapping.refStart: " << m_refMapping.refStart << std::endl;
                    }

                    DindelVariant var(m_refMapping.refName, refString, varString, refBaseOffsetMinPos+m_refMapping.refStart);
                    var.setHPLen(int(hplen));
                    var.setPriorProb(0.001); //FIXME
                    var.setHaplotypeUnique(varBaseOffsetMinPos, varBaseOffsetMaxPos);
                    
                    m_variants.push_back(var);

                    // this will be used in getClosestDistance
                    std::pair < HashMap<std::string, std::pair<int, int> >::iterator, bool> ins_pair =  
                               m_variant_to_pos.insert( HashMap<std::string, std::pair<int, int> >::value_type ( var.getID(), std::pair<int,int>(varBaseOffsetMinPos, varBaseOffsetMaxPos)));

                    assert (ins_pair.second == true);
                }

                // Reset state
                eventStart = -1;
                eventEnd = -1;
                isIndel = false;
            }
        }
    }
}

DindelHaplotype::DindelHaplotype(const std::string & haplotypeSequence, const DindelReferenceMapping & refMapping)
{
    // initialize
    m_pMA = NULL;
    m_seq = haplotypeSequence;
    m_refMapping = refMapping;
    m_firstAlignedBase = 0;
#ifdef DEBUG_NEW
    std::cout << "refMapping " << m_refMapping.refName << " " 
              << m_refMapping.refStart 
              << " length: " << m_refMapping.refSeq.size() 
              << " score: " << m_refMapping.referenceAlignmentScore << "\n";
#endif

    // determine haplotype sequence properties
    determineHomopolymerLengths();

    // align haplotype to the reference location.
    alignHaplotype();

    // extract the variants from the alignment of haplotype to reference sequence
    extractVariants();
}

DindelHaplotype::~DindelHaplotype()
{
    if (m_pMA != NULL) 
    {
        delete m_pMA;
        m_pMA = NULL;
    }
}

void DindelHaplotype::copy(const DindelHaplotype & haplotype, int copyOptions)
{
   assert(copyOptions == 0);
   m_seq = haplotype.m_seq;
   m_refPos = haplotype.m_refPos;
   m_hplen = haplotype.m_hplen;
   m_variants = haplotype.m_variants;
   m_variant_to_pos = haplotype.m_variant_to_pos;
   m_refMapping = haplotype.m_refMapping;
   m_isReference = haplotype.m_isReference;
   m_firstAlignedBase = haplotype.m_firstAlignedBase;

   if(haplotype.m_pMA != NULL)
       m_pMA = new MultipleAlignment(*haplotype.m_pMA);
   else
       m_pMA = NULL;

   // DO NOT CALL initHaploype()
}

DindelHaplotype::DindelHaplotype(const DindelHaplotype & haplotype, int copyOptions)
{
    copy(haplotype, copyOptions);
}

DindelHaplotype::DindelHaplotype(const DindelHaplotype & haplotype)
{
    copy(haplotype, 0);
}

void DindelHaplotype::determineHomopolymerLengths()
{
    // determine homopolymer length in sequence
    m_hplen = std::vector<int>(m_seq.size(), 1);
    std::vector<int> hp_forward = std::vector<int>(m_seq.size(), 0);
    std::vector<int> hp_reverse = std::vector<int>(m_seq.size(), 0);

    for (size_t r=1;r<m_seq.size();r++)
    {
        int add=0;
        if(m_seq[r]==m_seq[r-1]) 
            add=hp_forward[r-1]+1;
        hp_forward[r] = add;
    }

    for (int r=int(m_seq.size())-2;r>=0;r--)
    {
        int add=0;
        if(m_seq[r]==m_seq[r+1]) 
            add=hp_reverse[r+1]+1;
        hp_reverse[r] = add;
    }

    for (size_t r=0;r<m_seq.size();r++) 
        m_hplen[r] = hp_reverse[r] + hp_forward[r] + 1;
}

void DindelHaplotype::write(std::ostream & out) const
{
    out << "VARIANTS: "; for (size_t x=0;x<m_variants.size();x++) out << " " << m_variants[x].getID(); out << std::endl;
    for (size_t x=0; x<m_seq.size(); x++) 
        out << m_seq[x]; out << std::endl;
    for (size_t x=0; x<m_refPos.size(); x++) 
        out << " " << m_refPos[x]-m_refPos[0]; out << std::endl;
}


int DindelHaplotype::getClosestDistance(const DindelVariant& variant, int hapPos1, int hapPos2) const
{
    HashMap<std::string, std::pair<int, int> >::const_iterator it = this->m_variant_to_pos.find(variant.getID());
    if (it == this->m_variant_to_pos.end())
    {
        std::cout << " query variant: " << variant.getID() << std::endl;
        std::cout << " haplotype variants: "; for (it = m_variant_to_pos.begin();it != m_variant_to_pos.end();it++) std::cout << " " << it->first; std::cout << std::endl;
        std::cout << " h2 variants:        "; for (size_t x=0;x<m_variants.size();x++) std::cout << " " << m_variants[x].getID(); std::cout << std::endl;
        return -1;
    } 
    else 
    {
        int s = it->second.first;
        int e = it->second.second;
    
        // check if there is any overlap at all.
        if (hapPos1>hapPos2)
        {
            int t = hapPos1;
            hapPos1 = hapPos2;
            hapPos2 = t;
        }

        if (hapPos1<=s && hapPos2>=e)
        {
            std::set<int> dists;
            /*
            dists.insert(abs(hapPos1-s));
            dists.insert(abs(hapPos2-s));
            dists.insert(abs(hapPos1-e));
            dists.insert(abs(hapPos2-e));
            */
            // make sure we only count extensions from the left and right
            int d1=hapPos2-e;
            if (d1<0) d1=0;
            int d2=s-hapPos1;
            if (d2<0) d2=0;
            dists.insert( d1 );
            dists.insert( d2 );

            return *dists.begin();
        }
        else
        {
            // no overlap
            return -1;
        }
    }
}


int DindelHaplotype::getClosestDistance(const DindelVariant& variant, int hapPosStartRead, int hapPosEndRead, const DindelRead & read) const
{
    HashMap<std::string, std::pair<int, int> >::const_iterator it = this->m_variant_to_pos.find(variant.getID());

    if (it == this->m_variant_to_pos.end())
    {
        std::cout << " query variant: " << variant.getID() << std::endl;
        std::cout << " haplotype variants: "; for (it = m_variant_to_pos.begin();it != m_variant_to_pos.end();it++) std::cout << " " << it->first; std::cout << std::endl;
        std::cout << " h2 variants:        "; for (size_t x=0;x<m_variants.size();x++) std::cout << " " << m_variants[x].getID(); std::cout << std::endl;
        return -1;
    } else {
        int s = it->second.first;
        int e = it->second.second;

        // check if there is any overlap at all.
        if (hapPosStartRead>hapPosEndRead)
        {
            std::cerr << "Fix me Kees, please." << "\n";
            return -1;
            /*
            int t = hapPos1;
            hapPos1 = hapPos2;
            hapPos2 = t;
            */
        }
        // overlaps with variant?
        if (hapPosStartRead<=e && hapPosEndRead>=s)
        {

            int startRead = s - hapPosStartRead;
            int endRead = e - hapPosStartRead;
            if (endRead>read.length()-1) endRead = read.length()-1;
            if (endRead < 0) endRead = 0;
            if (startRead < 0) startRead = 0;
            if (startRead>read.length()-1) return -1; // ignore this case

            if (DINDEL_DEBUG_3 && 0)
            {
                std::cout << "\n";
                std::cout << "variant: " << variant.getID() << " s: " << s << " e: " << e  << "\n";
                std::cout << "hapPosStartRead: " << hapPosStartRead << " hapPosEndRead: " << hapPosEndRead  << "\n";

                int spacer = 0;
                if (hapPosStartRead<0) spacer = -hapPosStartRead;
                std::cout << std::string(spacer+s, ' '); std::cout << std::string(e-s+1,'*') << "\n";
                std::cout << std::string(spacer, ' ') << m_seq << "\n";

                if (hapPosStartRead>0) spacer = hapPosStartRead;
                std::cout << std::string(spacer, ' ') << read.getSequence() << "\n";
                std::cout << "\n";
            }

            int d = endRead-startRead+1;
            bool allMatch = true;
            for(int b = 0; b < d; b++)
            {
                if (read.getBase(startRead + b) != this->m_seq.at(s+b))
                {
                    allMatch = false;
                    break;
                }
            }

            if (!allMatch) return -1;

            std::set<int> dists;           
            dists.insert(abs( (s+e)/2 - hapPosEndRead) );
            dists.insert(abs( (s+e)/2 - hapPosStartRead) );
            
            // make sure we only count extensions from the left and right
            /*
            int d1=hapPosEndRead-e;
            if (d1<0) d1=0;
            int d2=s-hapPosStartRead;
            if (d2<0) d2=0;
            dists.insert( d1 );
            dists.insert( d2 );
            */
            return *dists.begin();
        }
        else
        {
            // no overlap
            return -1;
        }
    }
}


int DindelHaplotype::getHapBase(int refPos) const
{
    for(size_t x=0; x<this->m_refPos.size(); x++) 
        if (m_refPos[x] == refPos) 
            return x;
    return -1;
}


/*

 * DINDELMULTIHAPLOTYPE

 */

DindelMultiHaplotype::DindelMultiHaplotype(const std::vector< DindelReferenceMapping > & referenceMappings, const std::string & haplotypeSequence)
{
#ifdef DEBUG_NEW
    std::cout << "There are " << referenceMappings.size() << " reference mappings.\n";
#endif
    m_referenceMappings = referenceMappings;
    
    for (size_t i = 0; i < m_referenceMappings.size(); ++i)
        m_haplotypes.push_back(DindelHaplotype(haplotypeSequence, m_referenceMappings[i]));

    m_seq = haplotypeSequence;

    // hash table for realigning reads
    m_sequenceHash = DindelSequenceHash(m_seq);

    estimateMappingProbabilities();
}


void DindelMultiHaplotype::estimateMappingProbabilities()
{
    double sumScore = -HUGE_VAL;
    std::vector<double> scores(m_haplotypes.size(), 0.0);
    
    for(size_t i = 0; i < m_haplotypes.size(); ++i)
    {
        // compute p(H_i|R_i)
        // vars[x] gives the number of differences between the haplotype and the
        const std::vector<DindelVariant> & vars = m_haplotypes[i].getVariants();

        // use the GlobalAlign defaults but may decide to tweak this.
        double vscore = 0.0;
        for (size_t j = 0; j < vars.size(); ++j)
        {
            //std::cout << "\ti: " << i << " " << vars[j].getID() << "\n";
            const std::string & type = vars[j].getType();
            if (type == "SNP") vscore -= 3.0;
            else if (type == "MNP") vscore -= 3.0*double(vars[j].getAlt().length());
            else if (type == "INDEL") vscore -= double( (5 + 2*abs(vars[j].getLengthDifference()) ) );
            else assert(1==0);
        }

        if (DINDEL_ADJUST_MAPPINGQUAL)
	{
            scores[i] = double(m_referenceMappings[i].referenceAlignmentScore) + vscore;
        }
        else
        {
            scores[i] = double(m_referenceMappings[i].referenceAlignmentScore);
            std::cout << "DindelMultiHaplotype::estimateMappingProbabilities NOT adjusting mapping qualities\n";
        }
        sumScore = addLogs(sumScore, scores[i]);
    }

    // normalize and set log mapping probability
    for(size_t i = 0; i < m_haplotypes.size(); ++i)
    {
        m_haplotypes[i].m_refMapping.probMapCorrect = scores[i] - sumScore;
        if (DINDEL_DEBUG_3)
	{
            std::cout << "\tHaplotypes[" << i << "]: " << m_referenceMappings[i].refName << ":" << m_referenceMappings[i].refStart << " refSCore: " << m_referenceMappings[i].referenceAlignmentScore << " scores[i] " << scores[i] << " sumScore: " << sumScore << " probMapCorrect: " << scores[i] - sumScore << std::endl;
            if(m_haplotypes[i].m_refMapping.probMapCorrect>-5.0)
            {
                 const std::vector<DindelVariant> & vars = m_haplotypes[i].getVariants();
		 for (size_t j = 0; j < vars.size(); ++j)
                     std::cout << "\t\t" << vars[j].getID() << "\n";
			     
            }
        }

    }
}

int DindelMultiHaplotype::getHomopolymerLength(const std::string & chrom, int refPos) const
{
    int hp = -1;
    for(size_t i = 0; i < m_haplotypes.size(); ++i)
    {
        int thp = -1;
        if (m_haplotypes[i].m_refMapping.refName == chrom)
            thp = m_haplotypes[i].getHomopolymerLengthRefPos(refPos);

        if (thp >= 0)
        {
            hp = thp;
            break;
        }
    }

    return hp;
}

/*

 DINDELSEQUENCEHASH
 
 */

DindelSequenceHash::DindelSequenceHash(const std::string & sequence)
{
    makeHash(sequence);
}

void DindelSequenceHash::print() const
{
    for (Hash::const_iterator it=m_hash.begin();it!=m_hash.end();it++)
    {
        std::cout << " hash: " << it->first << " => ";
        for (std::list<int>::const_iterator i=it->second.begin();i!=it->second.end();i++) std::cout << " " << *i; std::cout << std::endl;
    }
    
}

/*
int DindelSequenceHash::lookup(unsigned int key) const
{
    UniqueHash::const_iterator iter = m_uniqueHash.find(key);
    if (iter!=m_uniqueHash.end()) return iter->second; else return -1;
}
*/

const std::list<int> & DindelSequenceHash::lookupList(unsigned int key) const
{
    Hash::const_iterator iter = m_hash.find(key);
    if (iter == m_hash.end())
        return emptyList;
    else
        return iter->second;
}

void DindelSequenceHash::makeHash(const std::string & sequence)
{
    for (size_t x=0;x<sequence.size()-DINDEL_HASH_SIZE;x++)
        m_hash[convert(sequence,x)].push_back(x);
}


/*
 *
 * DINDELWINDOW
 *
 *
 */
DindelWindow::DindelWindow(const std::vector<std::string> & haplotypeSequences,
                           const std::vector<DindelReferenceMapping>  & referenceMappings)
{

    m_pHaplotype_ma = NULL;
    m_referenceMappings = referenceMappings;
    
    // filter out duplicate haplotype sequences.
    initHaplotypes(haplotypeSequences, referenceMappings);

    // do multiple alignment of all haplotypes to each other.
    doMultipleHaplotypeAlignment();
}

DindelWindow::DindelWindow(const DindelWindow & window)
{
    copy(window);
}

void DindelWindow::copy(const DindelWindow & window)
{
    m_haplotypes = window.m_haplotypes;
    m_hashAltHaps = window.m_hashAltHaps;
    m_candHapAlgorithm = window.m_candHapAlgorithm;
    m_referenceMappings = window.m_referenceMappings;

    // Deep copy
    if(window.m_pHaplotype_ma != NULL)
        m_pHaplotype_ma = new MultipleAlignment(*window.m_pHaplotype_ma);
    else
        m_pHaplotype_ma = NULL;
}

void DindelWindow::initHaplotypes(const std::vector<std::string> & haplotypeSequences,
                                  const std::vector<DindelReferenceMapping>  & referenceMappings)
{
    for(size_t i = 0; i < haplotypeSequences.size(); ++i)
    {
        if(m_hashAltHaps.find(haplotypeSequences[i]) == m_hashAltHaps.end())
        {
            if (DINDEL_DEBUG_3)
                std::cout << "DindelWindow::MultiHaplotype " << m_haplotypes.size() << "\n";
            m_haplotypes.push_back(DindelMultiHaplotype(referenceMappings, haplotypeSequences[i]));
        }
    }
}

void DindelWindow::doMultipleHaplotypeAlignment()
{
    // globally align haplotypes to the first haplotype (arbitrary)
    std::vector<MAlignData> maVector;

    assert(m_haplotypes.size() >= 1);

    m_pHaplotype_ma = new MultipleAlignment;
    const std::string& rootSequence = m_haplotypes[0].getSequence();
    m_pHaplotype_ma->addBaseSequence("root", rootSequence, "");

    for (size_t h = 0; h < m_haplotypes.size(); ++h)
    {
        std::stringstream ss; 
        ss << "haplotype-" << h;
        addHaplotypeToMultipleAlignment(m_pHaplotype_ma, rootSequence, ss.str(), m_haplotypes[h].getSequence());
    }

    if (DINDEL_DEBUG_3)
    {
        std::cout << "DindelWindow::doMultipleHaplotypeAlignment:\n";
        m_pHaplotype_ma->print(500);

    }
}

DindelWindow::~DindelWindow()
{
    if (m_pHaplotype_ma!=NULL)
        delete m_pHaplotype_ma;
}

//
// DindelRealignWindowResult
//

//
void DindelRealignWindowResult::Inference::outputAsVCF(const DindelVariant & var, 
                                                       const DindelRealignWindowResult & result, 
                                                       VCFCollection& out,
                                                       const ReadTable* pRefTable) const
{
    VCFRecord record;
    std::stringstream hapPropString;
    hapPropString.precision(5);
    hapPropString.setf(std::ios::fixed,std::ios::floatfield);

    //double hapVarQual = 0.0; // this is the quality score of the haplotype this variant was called from.
    double varFreq = 0.0;

    double expNumberOfHaplotypesMappingToThisLocation = 0.0; // NOTE haplotypes that have this variant.
  
    double logProbNoHaplotypeMaps = 0.0;

    bool isEmpty = true;
    bool hasHighMappingQualityVariants = false;

    //double vprob = 0.0;
    double lnprob = 0.0;
    
    // indicates number of haplotypes supporting this variant.
    int numHQHaplotypes = 0;

    // first integrate across different iterations.

    std::map<int, std::vector<double> > logPostNotProbs;
    std::map<int, std::vector<double> > freqs;
    if (DINDEL_DEBUG_3)
        std::cout << "==outputAsVCF: variant: " << var.getID() << "\n";
    
    for(size_t i = 0; i < haplotypeProperties.size(); ++i)
    {
        // probability that at least one haplotype maps with high probability to this position in the reference
        double probHaplotypeDoesNotMap = 1.0-exp(haplotypeProperties[i].logMappingProb);
        if (probHaplotypeDoesNotMap<1e-6) probHaplotypeDoesNotMap=1e-6;
        if (probHaplotypeDoesNotMap>1.0-1e-6) probHaplotypeDoesNotMap=1.0-1e-6;

        int hapIt = haplotypeProperties[i].iteration;
        double hq = haplotypeProperties[i].qual/10.0*log(10.0);
        if (hq<0.0)
            hq = 0.0;

        // note this must be probabilities of NOT being true in log space
        double _lnprob = addLogs(log(probHaplotypeDoesNotMap), -hq+log(1.0-probHaplotypeDoesNotMap));
        logPostNotProbs[hapIt].push_back(_lnprob);

        if (DINDEL_DEBUG_3)
            std::cout << "==outputAsVCF: i: " << i << " hapIt: " << hapIt << " hq: " << hq << " lnprob: " << _lnprob << " haplotypeProperties[i].qual: " << haplotypeProperties[i].qual << " probHaplotypeDoesNotMap: " << probHaplotypeDoesNotMap << "\n";

        freqs[hapIt].push_back(haplotypeProperties[i].freq);

        // additional haplotype stats

        int hmq = int(round(-10.0*log10(probHaplotypeDoesNotMap)));

        if (hmq>0) hasHighMappingQualityVariants = true;

        logProbNoHaplotypeMaps += log(probHaplotypeDoesNotMap);

        expNumberOfHaplotypesMappingToThisLocation += (1.0-probHaplotypeDoesNotMap)*(1.0-exp(-hq));
        // set values in hapPropString
        if(hmq>0)
        {
            if (!isEmpty) hapPropString << ",";
            hapPropString << haplotypeProperties[i].qual <<  ":" << haplotypeProperties[i].freq << ":" << -10.0*log10(probHaplotypeDoesNotMap);
            isEmpty = false;

            if(haplotypeProperties[i].qual>100.0)
                numHQHaplotypes += 1;

        }
    }


    // combine the results from different iterations.
    // sum probabilities within a given iteration, and multiply values of different iterations.
    for (std::map<int, std::vector<double> >::const_iterator it =  logPostNotProbs.begin(); it != logPostNotProbs.end(); it++)
    {
        int iteration = it->first;
        int numHaps = int(it->second.size());

        double s = 0.0;
        for (int h = 0; h < numHaps; h++)
        {
            s += 1.0-exp(it->second[h]);
            varFreq += (1.0-exp(it->second[h]))*freqs[iteration][h];
        }

        double saves = s;
        if (s < 1e-6)
            s = 1e-6;
        
        double notPost = 1.0-s;
        

        lnprob += log(notPost);

        if (DINDEL_DEBUG_3)
            std::cout << "==outputAsVCF: iteration: " << iteration << " notPost: " << notPost << " lnprob: " << lnprob << " saves: " << saves << "\n";

    }
    
    double probNoHaplotypeMaps = exp(logProbNoHaplotypeMaps);
    if (probNoHaplotypeMaps<1e-6) probNoHaplotypeMaps=1e-6;
    int hapMappingQual = int(round(-10.0*log10(probNoHaplotypeMaps)));
    
    double vcfQual = -lnprob/.2303;
    //double vcfQual = -10.0*log10(1.0-exp(lvprob));
    if (hapPropString.str().empty()) 
        hapPropString << "NoMappings";
    
    // FIXME Add this as a parameter?
    if (!hasHighMappingQualityVariants) 
        return;

    int iqual = (vcfQual<0.0) ? 0 : int(round(vcfQual));

    // determine support for filter tag
    int totInHist = 0;
    if (!this->histDistance.empty())
        for (size_t x=0;x<histDistance.size();x++)
            totInHist += histDistance[x];

    // Count homopolymer runs
    std::set<int> hps;
    hps.insert(var.getHPLen());
    int hp = *hps.rbegin();

    double dust_score = HapgenUtil::calculateDustScoreAtPosition(var.getChrom(),
                                                                 var.getPos(),
                                                                 pRefTable);

    // Apply filters
    std::string filter="NoCall";
    if (iqual == 0)
        filter = "NoCall";
    else
    {
        if (totInHist == 0) 
        {
            filter = "NoSupp";
        }
        else
        {
            if (iqual<20)
            {
                if (numHQHaplotypes == 0)
                    filter = "LowQuality";
                else {
                    if (this->qual>100.0 && hapMappingQual<60.0)
                        filter = "AmbiMap";
                    else
                        filter = "LowPosterior";
                }
            }
            else if (iqual>=20)
                filter = "PASS";
        }
    }

    record.refName = var.getChrom();
    record.refPosition = var.getPos();
    record.id = result.outputID.empty() ? "." : result.outputID;
    record.refStr = var.getRef();
    record.varStr = var.getAlt();
    record.quality = iqual;
    record.passStr = filter;
    record.addComment("AF", varFreq);
    record.addComment("NumReads", numRealignedReads);
    record.addComment("NumCalledHaps", numCalledHaplotypes);
    record.addComment("ExpNumHapsWithVarMappingHere", expNumberOfHaplotypesMappingToThisLocation);
    record.addComment("VarQual", this->qual);
    record.addComment("HV", hapPropString.str());
    record.addComment("HMQ", hapMappingQual);
    record.addComment("VarDP", numReadsForward+numReadsReverse);
    record.addComment("NF", numReadsForward);
    record.addComment("NR", numReadsReverse);
    record.addComment("NF0", this->numReadsForwardZeroMismatch);
    record.addComment("NR0", this->numReadsReverseZeroMismatch);
    record.addComment("HSR", result.numHapSpecificReads);
    record.addComment("HPLen", hp);
    record.addComment("Dust", dust_score);
    record.addComment("SB", this->strandBias);
    
    // HistDist
    std::stringstream histdist_ss;
    for (size_t x=0;x<histDistance.size();x++) {
        if (x>0)
            histdist_ss << ",";
        histdist_ss << histDistance[x];
    }
    if (!histdist_ss.str().empty())
        record.addComment("HistDist", histdist_ss.str());

    // HistLik
    std::stringstream histlik_ss;
    for (size_t x=0;x<histAlignLik.size();x++) {
        if (x>0)
            histlik_ss << ",";
        histlik_ss << histAlignLik[x];
    }
    if (!histlik_ss.str().empty())
        record.addComment("HistLik", histlik_ss.str());

    // HistMapQ
    std::stringstream histmap_ss;
    for (size_t x=0;x<histMapQ.size();x++) {
        if (x>0)
            histmap_ss << ",";
        histmap_ss << histMapQ[x];
    }
    if (!histmap_ss.str().empty())
        record.addComment("HistMapQ", histmap_ss.str());


    // output genotyping information
    if (out.samples.size()>0 && result.outputGenotypes)
    {
        record.formatStr = "GT:GQ:GL";
       
        for (size_t x=0;x<out.samples.size();x++)
        {
            std::stringstream sampleString;
            const std::string & sample = out.samples[x];
            const  DindelRealignWindowResult::SampleToGenotypes::const_iterator sample_it=result.sampleToGenotypes.find(sample);
            
            if (sample_it!=result.sampleToGenotypes.end())
            {
                DindelRealignWindowResult::VarToGenotypeCall::const_iterator it_var = sample_it->second.find(var);
                assert(it_var!=sample_it->second.end());
                const DindelRealignWindowResult::GenotypeCall & gc = it_var->second;
                if (gc.qual>5.0) {
                    if (gc.count == 0) sampleString << "0/0";
                    else if (gc.count == 1) sampleString << "0/1";
                    else if (gc.count == 2) sampleString << "1/1";
                } else {                    
                    sampleString << "./.";
                }
                sampleString << ":" << int(gc.qual);
                //out <<  std::iostream::fixed;
                sampleString << std::setprecision(5) << ":" << gc.gl[0] << "," << gc.gl[1] << "," << gc.gl[2];
            } 
            else
            {
                sampleString << "."; // no call
            }
            record.sampleStr.push_back(sampleString.str());
        }
    }

    out.records.push_back(record);

}

double DindelRealignWindowResult::Inference::computeStrandBias(int numForward, int numReverse)
{
    // computes Bayes factor for strand bias
    // Even in the case of no strand bias, the expected frequency is allowed to differ from 0.5 in order to allow for
    // mapper bias.
    int total = numForward+numReverse;
    if (total == 0)
	return 0.0;

    const int N=50;
    double p = (1.0)/double(N);

    double bBias = 0.0;
    double bNoBias = 0.0;

    double nf = numForward;
    double nr = numReverse;
    //double comb=(total<=10) ? double(combinations(total, numForward )) : exp( lgamma(double(total+1))-lgamma(nf+1.0)-lgamma(double(total-numForward+1)));
    double comb;
    if (total<=10) {
       comb = log(double(combinations(total, numForward )));
    } else {
       comb = lgamma(double(total+1))-lgamma(nf+1.0)-lgamma(double(total-numForward+1));
    }

    double tBias = 0.35; // point at which we believe there is a bias
    double lbBias = -HUGE_VAL, lbNoBias = -HUGE_VAL;
    for (int i=0;i<N/2;i++)
    {
        double q = double(i+1)*tBias/(double(N+1)/2.0);
        bBias += pow(q,nf)*pow(1.0-q,nr)*p*comb;
	lbBias = addLogs(lbBias, log(q)*nf+log(1.0-q)*nr+log(p)+comb);

        q = 1.0-double(i+1)*tBias/(double(N+1)*2.0);
        bBias += pow(q,nf)*pow(1.0-q,nr)*p*comb;
	lbBias = addLogs(lbBias, log(q)*nf+log(1.0-q)*nr+log(p)+comb);
    }

    
    for (int i=0;i<N/2;i++)
    {
        double q = 0.5+(0.5-tBias)*double(i+1)/(double(N+1)/2.0);
        bNoBias += pow(q,nf)*pow(1.0-q,nr)*p*comb;
	lbNoBias = addLogs(lbNoBias, log(q)*nf+log(1.0-q)*nr+log(p)+comb);

        q = 0.5-(0.5-tBias)*double(i+1)/(double(N+1)/2.0);
        bNoBias += pow(q,nf)*pow(1.0-q,nr)*p*comb;
	lbNoBias = addLogs(lbNoBias, log(q)*nf+log(1.0-q)*nr+log(p)+comb);
    }
    
     //std::cout << "::compSB nf: " << nf << " nr: " << nr << " comb: " << comb << " bBias: " << bBias << " bNoBias: " << bNoBias << " lbBias: " << lbBias << " lbNoBias: " << lbNoBias << " lgamma(double(total+1)): " << lgamma(double(total+1)) << " lgamma(nf+1.0): " << lgamma(nf+1.0) << " lgamma(double(total-numForward+1)): " << lgamma(double(total-numForward+1)) << std::endl;
   
    //double bf = log(bBias)-log(bNoBias);
    
    return lbBias-lbNoBias;
}

void DindelRealignWindowResult::Inference::addDistanceToHistogram(int distance)
{
    const int D=25;
    if (this->histDistance.empty()) histDistance = std::vector<int>(D,0);
    if (distance>=0) {
        if (distance>D-1) distance = D-1;
        histDistance[distance]++;
    }
}

/// logLik is assumed to be the average per base
void DindelRealignWindowResult::Inference::addAlignLikToHistogram(double logLik)
{
    const int D=25;
    if (this->histAlignLik.empty()) histAlignLik = std::vector<int>(D,0);
    int bin = int(round(logLik* 25.0/(-10.0/50.0)));
    assert(bin>=0);
    if (bin>D-1) bin=D-1;
    histAlignLik[bin]++;
}

void DindelRealignWindowResult::Inference::addMapQToHistogram(double mappingQuality)
{
    const int D=8;
    if (this->histMapQ.empty()) histMapQ = std::vector<int>(D,0);
    int bin = int(round(mappingQuality/10.0)); 
    if (bin>D-1) bin=D-1;
    if (bin>=0) histMapQ[bin]++;
}

void DindelRealignWindowResult::outputVCF(VCFCollection& out, const ReadTable* pRefTable)
{
    VarToInference::const_iterator iter = variantInference.begin();
    for (;iter!=variantInference.end();iter++)
        iter->second.outputAsVCF(iter->first, *this, out, pRefTable);
}

/*
 *
 * DINDELREALIGNWINDOW
 *
 *
 */

DindelRealignWindow::DindelRealignWindow(const DindelWindow* pDindelWindow,
                                         std::vector<DindelRead> & dindelReads,
                                         const DindelRealignParameters & dindelRealignParameters) :
                                         m_dindelWindow(*pDindelWindow),
                                         m_pDindelReads(&dindelReads),
                                         realignParameters(dindelRealignParameters)                                        
{ 
    if (DINDEL_DEBUG) std::cout << "DindelRealignWindow::DindelRealignWindow STARTED" << std::endl;
    if (DINDEL_DEBUG) std::cout << "DindelRealignWindow::DindelRealignWindow DONE" << std::endl;    
}

/*
void DindelRealignWindow::run(const std::string & algorithm,
                              std::ostream& out)
{
    if (algorithm == "hmm")
    {
        algorithm_hmm(out);
    }
    else
    {
        throw std::string("DindelRealignWindow::run::unknown_algorithm");
    }

}
*/
void DindelRealignWindow::run(const std::string & algorithm,
                              VCFCollection& out,
                              DindelReadReferenceAlignmentVector* pOutAlignments, 
                              const std::string id,
                              DindelRealignWindowResult *pThisResult,
                              const DindelRealignWindowResult *pPreviousResult,
                              const ReadTable* pRefTable)
{
    PROFILE_FUNC("DindelRealignWindow::run")

    this->m_outputID = id;
    if (algorithm == "hmm")
    {
        algorithm_hmm(out, pOutAlignments, pThisResult, pPreviousResult, pRefTable);
    }
    else
    {
        throw std::string("DindelRealignWindow::run::unknown_algorithm");
    }

}

ReadHaplotypeAlignment DindelRealignWindow::computeReadHaplotypeAlignment(size_t readIndex, 
                                                                          const DindelMultiHaplotype& haplotype, 
                                                                          int start, 
                                                                          int end, 
                                                                          const std::vector<double> & lpError, 
                                                                          const std::vector<double> & lpCorrect, 
                                                                          bool rcReadSeq)
{
    const DindelRead & read = (*m_pDindelReads)[readIndex];
    // compute the log-likelihood of a single candidate alignment

    const int DRHA=0;

    if (DINDEL_DEBUG &&DRHA) std::cout << std::endl <<  "====> START 1 DindelRealignWindow::computeReadHaplotypeAlignment " << read.getID() << std::endl;
    // scenarios
    /* Capital R's are matches, r's are mismatches
    1. The seeds are too long to detect the possible gapped alignment of read to haplotype
           HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
                     rrrrrrrRRRRRRRRRRRRRRRRRRRRRRRRRRR
        which should be
           HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
                  RRR   rrrrRRRRRRRRRRRRRRRRRRRRRRRRRRR


    2. Deletion in read
           HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
                   RRRRRRR      RRRRRRRRRRRRRRRRRRRRRRRR
    */

    // start and end are inclusive
    // start may be <0, and end may be > haplotype length
    // lpError, lpCorrect

    int dlen = end+1-start-read.length();
    int rlen = read.length();
    int hlen = haplotype.length();
    const double minLpCorrect = realignParameters.addSNPMinLogBaseQual;

    if (DINDEL_DEBUG && DRHA)
    {
	    std::cout << std::endl;
	    std::cout << "HAPLOTYPE: " << haplotype.getSequence() << " isREF " << haplotype.isReference() << std::endl;
	    std::string pad;
	    std::string rseq = read.getSequence();
	    if (start>=0)
	    {
		    pad = std::string(start, ' ');
	    }
	    else {
		   rseq = read.getSequence().substr(-start, read.length());
	    }

	    std::cout << "READ:      " << pad << rseq << std::endl;
	    std::cout << "DLEN: " << dlen << " RLEN: " << rlen << " hlen: " << hlen << " start: " << start << " end: " << end << " read_index: " << readIndex << " " << read.getID() <<  std::endl;
    }

    int nmm=-1; // number of mismatches between read and candidate haplotype
    double logLik;
    bool isUngapped = false;
    if (dlen == 0)
    {
        // possible ungapped alignment
        typedef std::pair<int, char> Mismatch;
        std::list<Mismatch> mismatches;
        bool countMismatch = (realignParameters.addSNPMinMappingQual<read.getMappingQual())?true:false;

        double ll=0.0;
        int b=0;
    	double all_match = 0.0;
	    nmm=0;
        if (!rcReadSeq) {
	    for (int x=start;x<=end;x++,b++)
            {
		bool match = true;
		if (x>=0 && x<=hlen-1)
		{
			match = (read.getBase(b) == haplotype.getSequence()[x])?true:false;
		}
		if (match) ll += lpCorrect[b]; else {
		    if (DINDEL_DEBUG) std::cout << " mm " << b << " lpError: " << lpError[b] << std::endl;
                    ll += lpError[b]; nmm++;
                    if (countMismatch && lpCorrect[b]>minLpCorrect) mismatches.push_back(Mismatch(x, read.getBase(b)));
                }
		all_match += lpCorrect[b];
	    }
	} else {
	    for (int x=start;x<=end;x++,b++)
	    {
	        bool match = true;
	        if (x>=0 && x<=hlen-1)
		{
			match = (_complement(read.getBase(rlen-1-b)) == haplotype.getSequence()[x])?true:false;
		}
	        if (match) ll += lpCorrect[rlen-1-b]; else {
                    ll += lpError[rlen-1-b]; nmm++;
                    if (countMismatch && lpCorrect[rlen-b-1]>minLpCorrect) mismatches.push_back(Mismatch(x, _complement(read.getBase(rlen-1-b))));
                }
	    }
	}

        if (ADDSNPS && nmm<=realignParameters.addSNPMaxMismatches) {
            for (std::list<Mismatch>::const_iterator iter = mismatches.begin(); iter != mismatches.end(); iter++) 
	    {
                // FIXME Need to add wrt haplotype multiple alignment column
		assert(1==0);    
		/*
		int refPos = haplotype.getRefBase(iter->first);
		if (refPos>0)
		{
                    m_posToCandidateSNP.addSNP(readIndex, refPos, haplotype.getSequence()[iter->first], iter->second); // pos, ref, alt resp.
		}
		*/
            }
        }
	logLik = ll;
	isUngapped = true;
    }
    else { assert(1==0); }


    assert(!std::isnan(logLik));
    assert(logLik<0.0);

    logLik -= log(double(DINDEL_HMM_BANDWIDTH)); // added for consistency with the HMM, which has a prior

    // cap likelihood based on mapping quality
    // addLogs(logLik, -read.getMappingQual()*.23026-log(double(DINDEL_HMM_BANDWIDTH)));

    if (DINDEL_DEBUG && DRHA)
    {
	std::cout << "LOGLIK: " << logLik << std::endl;
    }
    if (DINDEL_DEBUG &&DRHA) std::cout << std::endl <<  "====> END DindelRealignWindow::computeReadHaplotypeAlignment" << std::endl;
    return ReadHaplotypeAlignment(logLik, nmm, isUngapped);
}

void DindelRealignWindow::addSNPsToHaplotypes()
{
    typedef std::map<size_t, std::vector<DindelVariant> > FreqCandidates;
    FreqCandidates freqCandidates;
    throw std::string("IMPLEMENT ME");

    for (MaPosToCandidateSNP::const_iterator iter = m_maPosToCandidateSNP.begin(); iter != m_maPosToCandidateSNP.end(); iter++)
    {
        for (size_t baseIdx=0;baseIdx<4;baseIdx++)
        {
            if (iter->second.baseToReads[baseIdx].size()>=realignParameters.addSNPMinCount)
            {
                // FIXME 
                /*
                std::string ref; ref+=iter->second.m_refBase;
                std::string alt; alt+=iter->second.getBase(baseIdx);
                DindelVariant snp(m_dindelWindow.getChrom(), ref, alt, iter->first);
                snp.setPriorProb(realignParameters.variantPriors.getDefaultProbVariant(snp.getType()));
                freqCandidates[iter->second.baseToReads[baseIdx].size()].push_back(snp);
                if(!QUIET) std::cout << "DindelRealignWindow::addSNPsToHaplotypes adding SNP " << m_dindelWindow.getChrom() << " " << ref << " "  <<  alt << " " <<  iter->first << std::endl;

                 */
            }
        }
    }

    // add snps until maximum is reached.
    int numAdded=0;
    for (FreqCandidates::const_reverse_iterator iter=freqCandidates.rbegin();iter!=freqCandidates.rend();iter++)
    {
        bool stop=false;
        for(size_t x=0;x<iter->second.size();x++)
        {
            m_dindelWindow.addVariant(iter->second[x], true);
            numAdded++;
            if (numAdded == realignParameters.addSNPMaxSNPs)
            {
                stop =true;
                break;
            }

        }
        if (stop) break;
    }

}



void DindelRealignWindow::algorithm_hmm(VCFCollection& out, 
                                        DindelReadReferenceAlignmentVector* pOutAlignments, 
                                        DindelRealignWindowResult * pThisResult, 
                                        const DindelRealignWindowResult *pPreviousResult,
                                        const ReadTable* pRefTable)
{
    if (DINDEL_DEBUG) std::cout << "DindelRealignWindow::algorithm_hmm STARTED" << std::endl;

    size_t numHaps = m_dindelWindow.getHaplotypes().size();
    assert(numHaps>=1);
    assert(pThisResult != NULL);

    // process the user-defined haplotypes
    computeReadHaplotypeAlignmentsUsingHMM(0, numHaps-1);

    // now add SNPs identified by the high quality ungapped alignments
    if (ADDSNPS && realignParameters.addSNPMaxSNPs>0)
    {
        addSNPsToHaplotypes();

        // compute the likelihoods for the additional haplotypes
        computeReadHaplotypeAlignmentsUsingHMM(numHaps, m_dindelWindow.getHaplotypes().size()-1);
    }

    DindelRealignWindowResult & result = *pThisResult;

    // estimate haplotype frequencies using user specified max mapping thresholds
    if (realignParameters.doEM)
    {
        if (realignParameters.realignMatePairs)
        {
            result = estimateHaplotypeFrequenciesModelSelectionMatePairs(realignParameters.minLogLikAlignToRef, 
                                                                         realignParameters.minLogLikAlignToAlt, 
                                                                         true, 
                                                                         pPreviousResult, 
                                                                         true);
        }
	    else 
        {
            if (realignParameters.multiSample)
            {
                if (out.samples.size())
                    realignParameters.genotyping = 1;
                else
                    realignParameters.genotyping = 0;
                
                result = estimateHaplotypeFrequenciesModelSelectionSingleReadsMultiSample(realignParameters.minLogLikAlignToRef, 
                                                                                          realignParameters.minLogLikAlignToAlt, 
                                                                                          true, 
                                                                                          pPreviousResult, 
                                                                                          true);
            }
            else
            {
                result = estimateHaplotypeFrequenciesModelSelectionSingleReads(realignParameters.minLogLikAlignToRef,
                                                                               realignParameters.minLogLikAlignToAlt, 
                                                                               true, 
                                                                               pPreviousResult, 
                                                                               true);
            }
        }
    }
    else
    {
        //FIXME There is no alternative yet...
        assert(1==0);
        //result = estimateHaplotypeFrequencies(realignParameters.minLogLikAlignToRef, realignParameters.minLogLikAlignToAlt, true, true);
    }
   
    // output the VCF results. result does marginalization over haplotypes
    result.outputVCF(out, pRefTable);

    // Compute the read-to-reference mapping qualities for the projected variants
    computeProjectedMappingQuality();

    // If the alignment output pointer is not NULL, copy the alignments out.
    if(pOutAlignments != NULL)
        std::copy(m_readReferenceAlignments.begin(), m_readReferenceAlignments.end(), std::back_inserter(*pOutAlignments));
    
    if (DINDEL_DEBUG) 
        std::cout << "DindelRealignWindow::algorithm_hmm DONE" << std::endl;
}

void DindelRealignWindow::computeProjectedMappingQuality()
{
    // First, separate the alignments for each read into vectors
    typedef std::vector<DindelReadReferenceAlignment*> AlignPtrVec;
    typedef std::map<std::string, AlignPtrVec> ReadAlignMap;

    ReadAlignMap read_alignment_map;
    for(size_t i = 0; i < m_readReferenceAlignments.size(); ++i)
    {
        DindelReadReferenceAlignment& a = m_readReferenceAlignments[i];
        AlignPtrVec& vec = read_alignment_map[a.read_name];
        vec.push_back(&a);
    }

    //
    // For each read, calculate the mapping quality of each alignment
    // Calculated as follows:
    //       P(location | R) = P(R, loc)
    //                         ---------
    //                            P(R)
    //
    // Where:
    //       location = the mapping location
    //       R = the read
    //       h = a haplotype
    //
    // And:       
    //   P(R, loc) = sum_h P(R|h) P(h | location) P(location) 
    //   P(R) = sum_location P(R, loc)
    //
    const std::vector<DindelMultiHaplotype>& haplotypes = m_dindelWindow.getHaplotypes();
    const std::vector<DindelReferenceMapping>& references = m_dindelWindow.getReferenceMappings();

    size_t num_locations = references.size();
    size_t num_haplotypes = haplotypes.size();
    double lp_location_prior = log(1.0 / (double)num_locations);

    // Compute a 2 dimensional vector of read likelihoods given haplotypes
    std::vector<std::vector<double> > readLikelihoods = fillReadHaplotypeLikelihoods();

    for(ReadAlignMap::iterator iter = read_alignment_map.begin(); iter != read_alignment_map.end(); ++iter)
    {
        // This is invariant over all the alignments so we can just use the first element of the vector
        assert(!iter->second.empty());
        size_t read_index = iter->second.front()->dindel_read_index;
        assert(read_index < readLikelihoods.size());

        //printf("Read[%zu] %s has %zu alignments\n", read_index, iter->first.c_str(), iter->second.size());
    
        // Compute the joint probability P(R, loc) for each location
        std::vector<double> joint_read_loc_lprobs;

        for(size_t loc = 0; loc < num_locations; ++loc)
        {
            std::vector<double> posteriors;
            for(size_t hap = 0; hap < num_haplotypes; ++hap)
            {
                assert(hap < readLikelihoods[read_index].size());
                double lp_read_given_hap = readLikelihoods[read_index][hap];
                double lp_hap_given_loc = haplotypes[hap].getLogMappingProbability(loc);
                double c = lp_read_given_hap + lp_hap_given_loc + lp_location_prior;
                posteriors.push_back(c);
                //printf("\t P(R|h]) = %lf P(h|loc) = %lf P(loc) = %lf\n", lp_read_given_hap, lp_hap_given_loc, lp_location_prior);
                //printf("\t P(R|h[%zu]) = %lf\n", hap, c);
            }

            // Sum over haplotypes
            double lp_rloc = posteriors[0];
            for(size_t hap = 1; hap < posteriors.size(); ++hap)
                lp_rloc = addLogs(lp_rloc, posteriors[hap]);

            // Store
            joint_read_loc_lprobs.push_back(lp_rloc);
            //printf("\t P(R,loc) = %lf\n", lp_rloc);
        }

        // Compute the normalizing constant P(R)
        double lp_r = joint_read_loc_lprobs[0];
        for(size_t i = 1; i < joint_read_loc_lprobs.size(); ++i)
            lp_r = addLogs(lp_r, joint_read_loc_lprobs[i]);
        //printf("\t P(R) = %lf\n", lp_r);

        // For each alignment, compute the mapping quality
        AlignPtrVec& alignments = iter->second;
        for(size_t i = 0; i < alignments.size(); ++i)
        {
            size_t loc_idx = alignments[i]->dindel_ref_index;
            //const DindelReferenceMapping& mapping = references[loc_idx];
            assert(loc_idx < joint_read_loc_lprobs.size());
            double lp_loc = joint_read_loc_lprobs[loc_idx];
            lp_loc -= lp_r;
            //int mapq = Quality::prob2phred(1.0f - exp(lp_loc));
            
            /*
            printf("Alignment %zu Ref: %s Pos: %d lp_loc: %lf p: %lf MQ: %d\n", loc_idx, mapping.refName.c_str(), 
                                                                         mapping.refStart, lp_loc, exp(lp_loc), mapq);
            */

        }
    }
}

void MaPosToCandidateSNP::addSNP(int readIndex, int refPos, char refBase, char altBase)
{

   //std::cout << " MaPosToCandidateSNP::size(): " << this->size() << " readIndex: " << readIndex << " refPos: " << refPos << " refBase: " << refBase << " altBase: " << altBase << std::endl;

    MaPosToCandidateSNP::iterator iter = this->find(refPos);
    if (iter == this->end())
    {
        std::pair<MaPosToCandidateSNP::iterator, bool> pib  = this->insert(MaPosToCandidateSNP::value_type(refPos, CandidateSNP(refBase)));
        iter = pib.first;
    }

    iter->second.addRead(refBase, altBase, readIndex);

}





void DindelRealignWindow::HMMAlignReadAgainstHaplotypes(size_t readIndex, size_t firstHap, size_t lastHap, 
                                                        const std::vector<double> & lpCorrect, const std::vector<double> & lpError,
                                                        HashMap<std::string, ReadHaplotypeAlignment>& hmm_alignment_cache)
{
    const std::vector<DindelMultiHaplotype>& haplotypes = m_dindelWindow.getHaplotypes();
    DindelRead& read = m_pDindelReads->at(readIndex);

    for (size_t h=firstHap;h<=lastHap;++h)
    {
        const DindelMultiHaplotype & haplotype = haplotypes[h];

        // Look up the read in the alignment cache using the haplotype index, the read sequence and read quality string as a key
        std::stringstream cache_key;
        cache_key << h << ":" << read.getSequence() << ":" << read.getQualString();
        ReadHaplotypeAlignment rha_hmm;

        HashMap<std::string, ReadHaplotypeAlignment>::iterator iter = hmm_alignment_cache.find(cache_key.str());
        if(iter != hmm_alignment_cache.end())
        {
            // Use the cached alignment
            rha_hmm = iter->second;
        }
        else
        {
            DindelHMM hmm(read, haplotype);
            rha_hmm = hmm.getAlignment();
            hmm_alignment_cache[cache_key.str()] = rha_hmm;
        }

        if (DINDEL_DEBUG)
        {
            std::cout << " HMM read " << readIndex 
                      << " haplotype " <<  h  
                      << " loglik: " << rha_hmm.logLik 
                      << " postProb: " << rha_hmm.postProbLastReadBase 
                      << " hapPosLastReadBase: " << rha_hmm.hapPosLastReadBase << std::endl;
        }

        // attempt ungapped alignment for detecting SNPs
        if (rha_hmm.postProbLastReadBase>this->realignParameters.minPostProbLastReadBaseForUngapped)
        {
            int end = rha_hmm.hapPosLastReadBase;
            int start = end-read.length()+1;
            bool rcReadSeq=read.getRCRead();
            ReadHaplotypeAlignment rha_ung = computeReadHaplotypeAlignment(readIndex, haplotype, start, end, lpError, lpCorrect, rcReadSeq);

            if (rha_ung.nmm<=2)
            {
                rha_hmm.isUngapped=true;
                rha_hmm.nmm = rha_ung.nmm;
            }

            if(DINDEL_DEBUG) 
                std::cout << " HMM read " << readIndex << " haplotype " << h  << " ungapped loglik: " << rha_ung.logLik << " nmm: " << rha_ung.nmm << std::endl;
        }
    
        hapReadAlignments[h][readIndex]=rha_hmm;
    }
}

void DindelRealignWindow::computeReadHaplotypeAlignmentsUsingHMM(size_t firstHap, size_t lastHap)
{
    PROFILE_FUNC("DindelRealignWindow::computeReadHaplotypeAlignmentsUsingHMM")

    if (firstHap>lastHap) return;
    if (DINDEL_DEBUG)
        std::cout << "DindelRealignWindow::computeReadHaplotypeAlignmentsUsingHMM firstHap " << firstHap << " " << lastHap << std::endl;
    // store information whether reads were realigned against all haplotypes or just the reference haplotype?
    const std::vector<DindelMultiHaplotype> haplotypes = m_dindelWindow.getHaplotypes();

    assert(haplotypes.size()>lastHap);
    for (size_t h=firstHap;h<=lastHap;h++)
    {
        if (haplotypes[h].isReference())
            hapReadAlignments.push_back( std::vector<ReadHaplotypeAlignment>(m_pDindelReads->size(), ReadHaplotypeAlignment(realignParameters.minLogLikAlignToRef,-1)));
        else
            hapReadAlignments.push_back( std::vector<ReadHaplotypeAlignment>(m_pDindelReads->size(), ReadHaplotypeAlignment(realignParameters.minLogLikAlignToAlt,-1)));
    }
    
    // We cache the results of the HMM alignment to save time when the read set contains many duplicate reads (like 1000 genomes).
    HashMap<std::string, ReadHaplotypeAlignment> hmm_alignment_cache;

    for (size_t r=0;r<m_pDindelReads->size();r++)
    {
        const DindelRead & read = (*m_pDindelReads)[r];
        if  (DINDEL_DEBUG) std::cout << "\n*****\nDindelRealignWindow::computeReadHaplotypeAlignmentsUsingHMM reads[" << r << "]: " << read.getID() << std::endl;
    
        // only realign if it has not yet been found to map to the reference
        std::vector<double> lpCorrect, lpError;
        read.getLogProbCorrectError(lpCorrect, lpError);

        HMMAlignReadAgainstHaplotypes(r, firstHap, lastHap, lpCorrect, lpError, hmm_alignment_cache);
    }
}

double DindelRealignWindow::getHaplotypePrior(const DindelHaplotype & h1, const DindelHaplotype & h2) const
{
   const std::vector<DindelVariant> & v1 = h1.getVariants();
   const std::vector<DindelVariant> & v2 = h2.getVariants();
   std::set<DindelVariant> uniqueVariants;
   for (size_t x=0;x<v1.size();x++) uniqueVariants.insert(v1[x]);
   for (size_t x=0;x<v2.size();x++) uniqueVariants.insert(v2[x]);

   double logPrior = 0.0;
  
   for (std::set<DindelVariant>::const_iterator it1=uniqueVariants.begin();it1!=uniqueVariants.end();it1++)
       logPrior += log(it1->getPriorProb());

   return logPrior;  
}


void DindelRealignWindow::printReadAlignments(int readIdx, std::ostream & out, int offset = 10, bool supportAlt = false)
{

    // print alignments against each candidate haplotype

   const std::vector<DindelMultiHaplotype> & haplotypes = m_dindelWindow.getHaplotypes();
   if (hapReadAlignments.size() != haplotypes.size()) throw std::string("DindelRealignWindow::printReadAlignments incomplete haplotype alignments.");
   DindelRead & read = (*m_pDindelReads)[readIdx];
   

   out << "ALIGNMENTS for read " << read.getID() << " sample: " << read.getSampleName() << " seq: " << read.getSequence() << " ref lik: " << hapReadAlignments[0][readIdx].logLik <<  std::endl;

   int hstart = 0;
   if (supportAlt) hstart = 0;
   bool printRef = false;
   for (int h=int(haplotypes.size())-1;h>=hstart;h--)
   {
       if (supportAlt && h!=0 && (hapReadAlignments[h][readIdx].logLik<-10.0 || (h!=0 && hapReadAlignments[h][readIdx].logLik <= hapReadAlignments[0][readIdx].logLik))) continue;
       if (h!=0) printRef = true;
       if (h==0 && !printRef) break;
       std::string pad(offset, ' ');


       // print variant annotations in the haplotype
       for(int refIdx = 0; refIdx < haplotypes[h].getNumReferenceMappings(); ++refIdx)
       {
            const DindelHaplotype & haplotype = haplotypes[h].getSingleMappingHaplotype(refIdx);

	       out << pad;
	       for (int b=0;b<haplotype.length();b++)
	       {
		   int rb = haplotype.getRefBase(b);

		   if (rb<0) 
		   {
		       if (rb == SNP || rb == MULTINUCLEOTIDE_RUN) out << "S";
		       else if (rb == INSERTION) out << "I";
		       else out << -rb;
		   }
		   else
		   {
		       if (b==0) out << " ";
		       else
		       {
			   // check for deletion
			   int prb = haplotype.getRefBase(b-1);
			   if (prb>=0 && rb>=0 && rb-prb>1) 
			   {
			       out << "D";
			   }
			   else
			       out << " ";
		       }
		   }
	       }
	       out << std::endl;
       }

       // print haplotype
       std::ostringstream os;
       os << h;
       out << os.str() << std::string(offset-os.str().length(),' ') << haplotypes[h].getSequence() << std::endl;


       // print read alignment
       const ReadHaplotypeAlignment & rha = hapReadAlignments[h][readIdx];

       int end = rha.hapPosLastReadBase;
       int start = end-read.length()+1;
       std::string readseq = read.getSequence();

       std::string rpad = "";
       if (start>0)
       {
          rpad = std::string(start, ' ');
       }

       out << pad << rpad;
       for (int b=start,r=0;b<=end && r<read.length();b++,r++)
       {
           if (b>=0)
           {
               if (b<haplotypes[h].length())
               {
                   // on haplotype
                   if (haplotypes[h].getSequence()[b]==readseq[r]) out << "."; else out << readseq[r];
               } 
           }
       }

       // output read haplotype alignment statistics. Still on same line

       out << " logLik: " << rha.logLik << " postProbLastReadBase: " << rha.postProbLastReadBase << " newpos: " << haplotypes[h].getSingleMappingHaplotype(0).getRefStart()+start << " start: " << start << " end: " << end;

       out << std::endl;
       out << std::endl;
   }
}

void DindelRealignWindow::doEM(const std::vector< std::vector<double> > & hrLik, const std::vector<int> & calledHaplotypes, std::vector<double> & haplotypeFrequencies)
{
    PROFILE_FUNC("DindelRealignWindow::doEM")

    // first remove uncalled haplotypes
    size_t nh=0;
    for(size_t h=0;h < calledHaplotypes.size(); ++h) if (calledHaplotypes[h]) nh++;
    size_t nr = hrLik.size();


    std::vector<double> rl(nh*nr,0.0); // read given haplotype likelihoods1
    
    size_t hidx = 0;
    for (size_t h = 0; h < calledHaplotypes.size(); ++h) if (calledHaplotypes[h])
    {
        
        for (size_t r=0;r<nr;r++)
        {
            rl[hidx*nr+r]=hrLik[r][h];
        }
        ++hidx;
    }

    if (DINDEL_DEBUG) std::cout << "doEM: nh: " << nh << " nr: " << nr << std::endl;
    

    std::vector<double> z(nh*nr,0.0); // expectations of read-haplotype indicator variables
    std::vector<double> pi(nh); // log of haplotype frequencies
    std::vector<double> nk(nh,0.0); // counts for each haplotype

    std::vector<double> hapFreqs=nk;

    // initialize frequencies
    for (size_t h=0;h<nh;h++) pi[h]=log(1.0/double(nh));

    // initialize expectations of indicator variables
    for (size_t h=0;h<nh;h++) for (size_t r=0;r<nr;r++)
    {
        z[h*nr+r]=0.5;
    }


    bool converged=false;
    double tol=this->realignParameters.EMtol;

    double eOld=-HUGE_VAL, eNew;

    int iter=0;
    while (!converged)
    {

        // compute expectation of indicator variables
        for (size_t h=0;h<nh;h++) nk[h]=0.0;

        for (size_t r=0;r<nr;r++)
        {
            double lognorm=-HUGE_VAL;
            // compute responsibilities
            for (size_t h=0;h<nh;h++)
            {
                z[h*nr+r]=pi[h]+(rl[h*nr+r]);
                lognorm=addLogs(lognorm, z[h*nr+r]);
            }
            // normalize and compute counts
            for (size_t h=0;h<nh;h++)
            {
                z[nr*h+r]-=lognorm;
                z[nr*h+r]=exp(z[nr*h+r]);

                nk[h]+=z[nr*h+r];
            }
        }


       // compute frequencies
       double zh = 0.0;
       for (size_t h=0;h<nh;h++) zh += nk[h];
       for (size_t h=0;h<nh;h++) pi[h] = log(nk[h]/zh); 

       eNew = 0.0;
       for (size_t r=0;r<nr;r++)
       {
           double t = -HUGE_VAL;
           for (size_t h=0;h<nh;h++)
           {
               // compute responsibilities
               t = addLogs(t, pi[h]+rl[h*nr+r]);
           }
           eNew += t;
       }
 

       if (DINDEL_DEBUG) std::cout << " EM iter: " << iter << " " << eNew << " eOld-eNew: " <<  eOld-eNew << std::endl;
       //
       //assert (eOld-eNew<1e-10);
       if (eOld-eNew > 1e-10) 
       {
           std::cout << "!!EM ERROR: " << std::endl;
           std::cout << "!!nr: " << nr << " nh: " << nh << std::endl;
           for(size_t h = 0; h < nh; h++) for (size_t r = 0; r < nr; r++)
           {
               std::cout << "!rl[" << h << "," << r << ",]: " << rl[nr*h+r] << std::endl;
           }
       }

       converged=(fabs(eOld-eNew))<tol || iter>realignParameters.EMmaxiter;

       eOld=eNew;


       iter++;
    }

    haplotypeFrequencies = std::vector<double>(calledHaplotypes.size(),0.0);

    hidx = 0;
    for(size_t h=0;h < calledHaplotypes.size();h++)
    {
        if (calledHaplotypes[h]) haplotypeFrequencies[h] = exp(pi[hidx++]); else haplotypeFrequencies[h] = 0.0;
    }

    if(DINDEL_DEBUG_3 || SHOWHAPFREQ)
    {
        std::cout  << "DindelRealignWindow::doEM haplotype Frequencies: ";
        hidx = 0;
        for (size_t h=0;h<calledHaplotypes.size();h++) std::cout << " [ " << h << " " << haplotypeFrequencies[h] << "];";
        std::cout << std::endl;
    }
       
}

void DindelRealignWindow::doEMMultiSample(int numSamples,
                                          int maxIter,
                                          const std::vector<double> & llHapPairs,
                                          const std::vector<int> & allowedHaplotypes,
                                          const std::vector<double> & initHaplotypeFrequencies,
                                          std::vector<double> & haplotypeFrequencies)
{
    PROFILE_FUNC("DindelRealignWindow::doEMMultiSample")

    // llHapPairs is ordered by sample then by haplotype pair h1,h2 with h2>=h1
    size_t ns = (size_t) numSamples;
    size_t nh = allowedHaplotypes.size();
    size_t numPairs = nh*(nh+1)/2;
    // check if the dimensions are correct
    assert (llHapPairs.size() == ns*numPairs);

  
    std::vector<double> z(numPairs*ns,0.0); // expectations of read-haplotype indicator variables
    std::vector<double> pi(nh); // log of haplotype frequencies
    std::vector<double> nk(nh,0.0); // counts for each haplotype

    std::vector<double> hapFreqs=nk;

    // initialize frequencies
    for (size_t h = 0;h<allowedHaplotypes.size();h++)
    {
        if (allowedHaplotypes[h])
            pi[h] = (initHaplotypeFrequencies[h]<1e-16)?log(1e-16):log(initHaplotypeFrequencies[h]);
    }
    
    bool converged=false;
    double tol=this->realignParameters.EMtol;

    double eOld=-HUGE_VAL, eNew;

    int iter=0;
    while (!converged)
    {

        // compute expectation of indicator variables
        for (size_t h=0;h<nh;h++) nk[h]=0.0;

        for (size_t s=0;s<ns;s++)
        {
            size_t sampleOffset = s*numPairs;
            double lognorm=-HUGE_VAL;

            // compute responsibilities
            
            for (size_t h1=0, gt = 0;h1<nh;h1++)
                for (size_t h2=h1;h2<nh;h2++,gt++)
                {
                    if (allowedHaplotypes[h1] && allowedHaplotypes[h2])
                    {
                        double lr = pi[h1]+pi[h2]+llHapPairs[sampleOffset+gt];
                        z[sampleOffset+gt]= lr;
                        lognorm=addLogs(lognorm, lr);
                        if (h1!=h2)
                            lognorm=addLogs(lognorm, lr); // do twice for off-diagonal terms
                    }
                }

            // normalize and compute counts
            for (size_t h1=0, gt = 0;h1<nh;h1++)
                for (size_t h2=h1;h2<nh;h2++,gt++)
                {
                    if (allowedHaplotypes[h1] && allowedHaplotypes[h2])
                    {
                        double _z = exp(z[sampleOffset + gt] - lognorm);
                        nk[h1] += _z;
                        if (h1!=h2)
                            nk[h2] += _z;
                    }
                }
       }


       // compute frequencies
       double zh = 0.0;
       for (size_t h=0;h<nh;h++)
           zh += nk[h];
       for (size_t h=0;h<nh;h++)
       {
           if (allowedHaplotypes[h])
               pi[h] = log(nk[h]/zh);
       }
       eNew = 0.0;
       for (size_t s=0;s<ns;s++)
       {
           double t = -HUGE_VAL;
           size_t sampleOffset = s*numPairs;
           for (size_t h1=0, gt = 0;h1<nh;h1++)
                for (size_t h2=h1;h2<nh;h2++,gt++)
                {
                    if (allowedHaplotypes[h1] && allowedHaplotypes[h2])
                    {
                        double lr = pi[h1]+pi[h2]+llHapPairs[sampleOffset+gt];
                        t = addLogs(t, lr);
                        if (h1!=h2)
                            t = addLogs(t, lr);
                    }
                }
           eNew += t;
       }


       if (DINDEL_DEBUG_3) std::cout << " EM MultiSample iter: " << iter << " " << eNew << " eOld-eNew: " <<  eOld-eNew << std::endl;
       //
       //assert (eOld-eNew<1e-10);
       if (eOld-eNew > 1e-10)
       {
           std::cout << "!!EM MultiSample ERROR: " << std::endl;
       }

       converged=(fabs(eOld-eNew))<tol || iter>maxIter;

       eOld=eNew;

       iter++;
    }

    haplotypeFrequencies = std::vector<double>(allowedHaplotypes.size(),0.0);

    for(size_t h=0;h < allowedHaplotypes.size();h++)
    {
        if (allowedHaplotypes[h]) haplotypeFrequencies[h] = exp(pi[h]); else haplotypeFrequencies[h] = 0.0;
    }

    if(DINDEL_DEBUG_3 || SHOWHAPFREQ)
    {
        std::cout  << "DindelRealignWindow::doEMMultiSample haplotype Frequencies: ";        
        for (size_t h=0;h<allowedHaplotypes.size();h++) std::cout << " [ " << h << " " << haplotypeFrequencies[h] << "];";
        std::cout << std::endl;
    }

}


void DindelRealignWindow::addDiploidGenotypes(DindelRealignWindowResult& result, const std::vector<int> & allowedHaplotype, const std::vector< std::vector<double> > & hrLik)
{
  
   const std::vector<DindelMultiHaplotype> & haplotypes = m_dindelWindow.getHaplotypes();
   const std::vector<DindelRead> & reads = *m_pDindelReads;
   if (DINDEL_DEBUG) std::cerr << "GENOTYPING START" << std::endl;
   if (hapReadAlignments.size() != haplotypes.size()) throw std::string("DindelRealignWindow::addDiploidGenotypes incomplete haplotype alignments.");
   int numHaps = int(haplotypes.size());

   typedef HashMap<std::string, std::list<int> > SampleToReads;
   SampleToReads sampleToReads;

   for (size_t r=0;r<reads.size();r++)
   {
       const std::string & sampleName = reads[r].getSampleName();
       sampleToReads[sampleName].push_back(r);
   }

   typedef std::map<DindelVariant, HashMap<int, int> > VariantToHaplotype;
   VariantToHaplotype variantToHaplotype;

   for (int i=0;i<numHaps;i++)
   {
       for (int refIdx = 0; refIdx<haplotypes[i].getNumReferenceMappings(); refIdx++)
       {
           const std::vector<DindelVariant> & vars = haplotypes[i].getVariants(refIdx);
           for (size_t x=0;x<vars.size();x++) variantToHaplotype[vars[x]][i] = refIdx;
       }
   }


   for (SampleToReads::const_iterator sample_it = sampleToReads.begin(); sample_it != sampleToReads.end(); sample_it++)
   {
       const std::string & sample = sample_it->first;

       // infer likelihoods for each pair of haplotypes
       std::vector<double> llHapPairs(numHaps*numHaps);

       double max_ll=-HUGE_VAL;
       int h1max=-1, h2max=-1;

       for (int h1=0;h1<numHaps;h1++) 
           for (int h2=h1;h2<numHaps;h2++) 
           {
               double ll=0.0;
               for (std::list<int>::const_iterator _r = sampleToReads[sample].begin(); _r != sampleToReads[sample].end(); ++_r)
               {
                   int r = *_r;
                   double lh1 = hrLik[r][h1];
                   double lh2 = hrLik[r][h2];
                   //if (h1==0 && h2==1) {
                    //       std::cerr << "r: " << r << " lh1: " << lh1 << " lh2: " << lh2 << std::endl;
                   //}
                   ll += log(.5) + addLogs(lh1, lh2);
               }

               llHapPairs[h1*numHaps+h2]=ll;
               llHapPairs[h2*numHaps+h1]=ll;

               if (allowedHaplotype[h1] && allowedHaplotype[h2])
               {
                   // only use allowed haplotypes to find the most likely pair of haplotypes.
                   if (ll>max_ll)
                   {
                       max_ll = ll;
                       h1max = h1;
                       h2max = h2;
                   }
               }
           }

       assert(h1max!=-1 && h2max !=-1);


       // genotype call for most likely pair of haplotypes
       result.sampleToGenotypes[sample] = DindelRealignWindowResult::VarToGenotypeCall();
       DindelRealignWindowResult::VarToGenotypeCall & gc = result.sampleToGenotypes[sample];
       
       // initialize all variant genotypes to NOT called (ref/ref)
       for (VariantToHaplotype::const_iterator it = variantToHaplotype.begin();it!=variantToHaplotype.end();it++)
           gc[it->first] = DindelRealignWindowResult::GenotypeCall();

       // use referenceMappings to determine how many haplotypes were mapped to a given reference location.
       // This is a crude way to determine homozygous from heterozygous positions (as in two or one haplotype mapped to a reference location for this sample)
       DindelRealignWindowResult::VarToGenotypeCall::iterator it1;
       for (it1=gc.begin();it1!=gc.end();it1++)
       {
           // initialize quality
           it1->second.qual = 123456.0;
           it1->second.called = true;
           for (int x=0;x<3;x++) it1->second.gl[x] = -HUGE_VAL;
       }

       for (int refIdx = 0; refIdx < int( m_dindelWindow.getReferenceMappings().size() ); refIdx++)
       {

           // set genotypes based on best haplotype pair
           const std::vector<DindelVariant> & v1 = haplotypes[h1max].getVariants(refIdx);
           const std::vector<DindelVariant> & v2 = haplotypes[h2max].getVariants(refIdx);



           // gc contains genotypes for maximum likelihood pair of haplotypes
           for (size_t x=0;x<v1.size();x++) gc[v1[x]].count += 1;
           for (size_t x=0;x<v2.size();x++) gc[v2[x]].count += 1;

          

           // always allow the reference haplotype for the mapping location so that we don't get -inf likelihoods for genotypes involving the reference
           for (int h1=0;h1<numHaps;h1++) // if (allowedHaplotype[h1] || haplotypes[h1].getVariants(refIdx).size()==0)
               for (int h2=h1;h2<numHaps;h2++) // if (allowedHaplotype[h2] || haplotypes[h2].getVariants(refIdx).size()==0)
               {
                   double ll_this_pair = llHapPairs[h1*numHaps + h2];

                   if (true)// (h1!=h1max || h2!=h2max)) NOT NECESSARY?
                   {
                       // different pair of haplotypes
                       // only consider likelihood difference if it is higher, otherwise genotype quality will be zero.
                        const std::vector<DindelVariant> & _v1 = haplotypes[h1].getVariants(refIdx);
                        const std::vector<DindelVariant> & _v2 = haplotypes[h2].getVariants(refIdx);
                        DindelRealignWindowResult::VarToGenotypeCall _gc;

                        for (size_t x=0;x<_v1.size();x++) _gc[_v1[x]].count += 1;
                        for (size_t x=0;x<_v2.size();x++) _gc[_v2[x]].count += 1;

                        // update genotype qualities for called variants
                        for (it1 = gc.begin();it1 != gc.end();it1++)
                        {
                            DindelRealignWindowResult::VarToGenotypeCall::const_iterator it2;
                            const DindelVariant & variant = it1->first;
                            it2 = _gc.find(variant);
                            //bool countIsDifferent = true;
                            if (it2 == _gc.end())
                            {
                                // variant not called in this pair of haplotypes
                                int count = 0;
                                it1->second.gl[count] = addLogs(it1->second.gl[count], ll_this_pair);
                            } else
                            {
                                it1->second.gl[it2->second.count] = addLogs(it1->second.gl[it2->second.count], ll_this_pair);
                            }

                        }
                   }
               }


           for (it1 = gc.begin();it1 != gc.end();it1++)
           {
               int count = it1->second.count;
               double qual = 1000.0;
               for (int c=2;c>0;c--)
               {
                   double d = it1->second.gl[count]-it1->second.gl[c];
                   if (c!=count && d<qual)
                   {
                       qual = d;
                   }
               }
               if (qual<0.0) qual = 0.0;
               it1->second.qual = 10.0*qual/log(10.0);
               // prevent non-ref genotype calls when there is little evidence
               if (count != 0 && fabs(it1->second.gl[count]-it1->second.gl[0])<0.05)
                   it1->second.count = 0; // assign ref 
           }
       }
   }
   if (DINDEL_DEBUG) std::cerr << "GENOTYPING END" << std::endl;

}

void DindelRealignWindow::computeHaplotypePairLikelihoods(std::vector<double> & llHapPairs,
                                                          HashMap<std::string, std::list<int> > & sampleToReads,
                                                          const std::vector< std::vector<double> > & hrLik)
{
   const std::vector<DindelRead> & reads = *m_pDindelReads;
   if (reads.size()==0)
       return;
   
   typedef HashMap<std::string, std::list<int> > SampleToReads;
   sampleToReads.clear();
   for (size_t r=0;r<reads.size();r++)
   {
       const std::string & sampleName = reads[r].getSampleName();
       sampleToReads[sampleName].push_back(r);
   }

   int numSamples = int (sampleToReads.size());

   int numAH = int(hrLik[0].size());
   int numGT = numAH*(numAH+1)/2;
   llHapPairs = std::vector<double>(numSamples*numGT,-HUGE_VAL);


   // first compute all haplotype pair likelihoods
   std::vector<double> tmpLL(numAH*numAH,0.0);
   int sample_idx = 0;
   for (SampleToReads::const_iterator strIt = sampleToReads.begin(); strIt != sampleToReads.end(); strIt++, sample_idx++)
   {
       // infer likelihoods for non-candidate haplotype/cand haplotype genotypes

       for (int x=0;x<numAH*numAH;x++) tmpLL[x]=0.0;

       for (std::list<int>::const_iterator rit = strIt->second.begin(); rit != strIt->second.end(); rit++)
       {
           int r = *rit;
   
           for (int h1=0;h1<numAH;h1++)
               for (int h2=h1;h2<numAH;h2++)
               {
                   double lh1 = hrLik[r][h1];//addLogs(hapReadAlignments[allowedHaplotypes[h1]][r].logLik, mq);
                   double lh2 = hrLik[r][h2];//addLogs(hapReadAlignments[allowedHaplotypes[h2]][r].logLik, mq);
                   //if (h1==0 && h2==1) {
                    //       std::cerr << "r: " << r << " lh1: " << lh1 << " lh2: " << lh2 << std::endl;
                   //}
                   double ll = log(.5) + addLogs(lh1, lh2);
                   tmpLL[h1*numAH+h2]+=ll;
               }
               if (DEBUG_CALLINDEL) {
                       std::cerr << "AFTER READ " << r << std::endl;
                       for (int x=0;x<numAH*numAH;x++) std::cerr << "    TMPLIK[" << x << "]:" << tmpLL[x] << std::endl;
               }
       }

       int gt = 0;
       for (int h1=0;h1<numAH;h1++)
           for (int h2=h1;h2<numAH;h2++)
           {
               llHapPairs[sample_idx*numGT+gt] = addLogs(llHapPairs[sample_idx*numGT+gt], tmpLL[h1*numAH+h2]);
               gt++;
           }
   }
}

void DindelRealignWindow::addCalledHaplotypeMatePairs(int hapIdx,
                                             DindelRealignWindowResult & result,
                                             const std::vector<DindelMultiHaplotype> & haplotypes,
                                             const std::vector< HashMap<int, double> > & addReads,
                                             double minLogLikAlignToAlt,
                                             int numReadPairs)
{

        if (DINDEL_DEBUG_3)
            std::cout << "Adding haplotype " << hapIdx << "\n";

        const std::vector<DindelRead> & reads = (*m_pDindelReads);
   
        
        bool indelAdded = false;

        // store haplotype calling results.
        std::pair<DindelRealignWindowResult::HapIdxToInference::iterator, bool> hapit_pair = result.hapIdxToInference.insert(DindelRealignWindowResult::HapIdxToInference::value_type(hapIdx,DindelRealignWindowResult::Inference()));
        DindelRealignWindowResult::Inference & hapInf = hapit_pair.first->second;

        hapInf.numReadsForward=0; hapInf.numReadsReverse=0; hapInf.numReadsForwardZeroMismatch=0; hapInf.numReadsReverseZeroMismatch=0; hapInf.numUnmapped = 0; hapInf.numLibraries = 0; hapInf.numReadNames = 0;
        // store order of adding haplotypes.
        result.addOrder.push_back(hapIdx);

        // determine which variants have not been called yet
        for (int refIdx = 0; refIdx < haplotypes[hapIdx].getNumReferenceMappings(); ++refIdx)
        {
            const std::vector<DindelVariant> & vars = haplotypes[hapIdx].getVariants(refIdx);

            for (size_t x = 0; x < vars.size(); x++)
            {
                if (DINDEL_DEBUG_3 && haplotypes[hapIdx].getLogMappingProbability(refIdx)>-5.0)
                    std::cout << "  adding var[" << x << "]: hapIdx: " << hapIdx << " refIdx: " << refIdx << " " << vars[x].getID() << std::endl;

                DindelRealignWindowResult::VarToInference::iterator vit = result.variantInference.find(vars[x]);
                if (vit == result.variantInference.end())
                {
                    std::pair<DindelRealignWindowResult::VarToInference::iterator, bool> it_pair = result.variantInference.insert(DindelRealignWindowResult::VarToInference::value_type(vars[x],hapInf));
                    vit = it_pair.first;
                    vit->second.numRealignedReads = int(reads.size());
                    vit->second.qual = 0.0;
                }

                DindelRealignWindowResult::Inference & varInf = vit->second;

                HashMap<std::string, int> libraries, readnames;
                varInf.haplotypeIndex.insert(hapIdx);

                if (vars[x].getType()=="INDEL") indelAdded = true;

                //for (std::set<int>::const_iterator r = addReads[hapIdx].begin(); r != addReads[hapIdx].end(); r++) // read pair indices
                for(HashMap<int, double>::const_iterator hit = addReads[hapIdx].begin(); hit != addReads[hapIdx].end(); hit++)
                {
                    int r = hit->first;
                    double dH = hit->second;
                    for (int mate = 0; mate <= 1; ++mate)
                    {
                        int read_idx = r + mate*numReadPairs;

                        if (varInf.readIndex.find(read_idx) == varInf.readIndex.end())
                        {

                            // add variant convering statistics
                            if (hapReadAlignments[hapIdx][read_idx].logLik>minLogLikAlignToAlt &&
                                hapReadAlignments[hapIdx][read_idx].nmm<=2 &&
                                hapReadAlignments[hapIdx][read_idx].nmm>=0 &&
                                hapReadAlignments[hapIdx][read_idx].isUngapped)
                            {
                                int dist = haplotypes[hapIdx].getSingleMappingHaplotype(refIdx).getClosestDistance(vars[x],
                                                                                                                   hapReadAlignments[hapIdx][read_idx].hapPosLastReadBase-reads[read_idx].length()+1,
                                                                                                                   hapReadAlignments[hapIdx][read_idx].hapPosLastReadBase,
                                                                                                                   reads[read_idx]);
                                if(dist != -1)
                                {

                                    // read overlaps the variant.


                                    if (DINDEL_DEBUG_3 && 0)
                                    {
                                        std::cout << "\t likelihood hapIdx: " << hapIdx << " read_idx: " << read_idx << " lik: " << hapReadAlignments[hapIdx][read_idx].logLik << "\n";
                                    }

                                    // use it only once if it overlaps. Maybe in some read-haplotype alignments this read doesn't support the variant...
                                    varInf.readIndex.insert(read_idx);

                                    // add quality in phred scale
                                    // note that this quality uses only reads that overlap this position of the variant in the haplotype
                                    varInf.qual += dH/.23026;

                                    varInf.addDistanceToHistogram(dist);
                                    // update statistics using supporting read
                                    if (reads[read_idx].isForward()) varInf.numReadsForward++; else varInf.numReadsReverse++;
                                    if (hapReadAlignments[hapIdx][read_idx].nmm==0)
                                    {
                                        if (reads[read_idx].isForward()) varInf.numReadsForwardZeroMismatch++; else varInf.numReadsReverseZeroMismatch++;
                                    }
                                    if ((*m_pDindelReads)[read_idx].isUnmapped()) varInf.numUnmapped++;
                                    varInf.addAlignLikToHistogram(hapReadAlignments[hapIdx][read_idx].logLik/double(reads[read_idx].length()));
                                    varInf.addMapQToHistogram(reads[read_idx].getMappingQual());
                                    libraries[reads[read_idx].getLibraryName()]=1;
                                    readnames[reads[read_idx].getID()]=1;
                                }

                            }
                        }
                    }
                } // loop over supporting reads addReads[hapIdx]

                varInf.strandBias = DindelRealignWindowResult::Inference::computeStrandBias(varInf.numReadsForward, varInf.numReadsReverse);
                varInf.numLibraries = int(libraries.size());
                varInf.numReadNames = int(readnames.size());
            }
        } // for (int refIdx = 0; refIdx < haplotypes[hapIdx].getNumReferenceMappings(); ++refIdx)

        if (indelAdded && realignParameters.showCallReads!=0)
        {
            if(DINDEL_DEBUG_3)
            {
                std::cout << "DindelRealignWindow::estimateHaplotypeFrequenciesModelSelectionReadPairs Adding haplotype " << hapIdx << " with INDEL " << std::endl;
                std::cout << "  variants: ";
            }
            for (int refIdx = 0; refIdx < haplotypes[hapIdx].getNumReferenceMappings(); ++refIdx)
            {
                if(DINDEL_DEBUG_3)
                    std::cout << "   refmapping[" << refIdx << "]: ";
                const std::vector<DindelVariant> & vars = haplotypes[hapIdx].getVariants(refIdx);
                for (size_t x=0;x<vars.size();x++) std::cout << " " << vars[x].getID(); std::cout << std::endl;
            }

            for (HashMap<int, double>::const_iterator hit = addReads[hapIdx].begin(); hit != addReads[hapIdx].end();hit++)
            {
                int r = hit->first;
                for (int mate = 0; mate <= 1; ++mate)
                    printReadAlignments(r+mate*numReadPairs, std::cout);
            }
        }

    
}

void DindelRealignWindow::addCalledHaplotypeSingleRead(int hapIdx,
                                                       DindelRealignWindowResult & result,
                                                       const std::vector<DindelMultiHaplotype> & haplotypes,
                                                       const std::vector< HashMap<int, double> > & addReads,
                                                       double minLogLikAlignToAlt,
                                                       int numReads)
{

    (void) numReads;
    if (DINDEL_DEBUG_3)
        std::cout << "Adding haplotype " << hapIdx << "\n";

    const std::vector<DindelRead> & reads = (*m_pDindelReads);
    bool indelAdded = false;

    // store haplotype calling results.
    std::pair<DindelRealignWindowResult::HapIdxToInference::iterator, bool> hapit_pair = 
            result.hapIdxToInference.insert(DindelRealignWindowResult::HapIdxToInference::value_type(hapIdx,DindelRealignWindowResult::Inference()));

    DindelRealignWindowResult::Inference & hapInf = hapit_pair.first->second;

    hapInf.numReadsForward=0; 
    hapInf.numReadsReverse=0; 
    hapInf.numReadsForwardZeroMismatch=0; 
    hapInf.numReadsReverseZeroMismatch=0; 
    hapInf.numUnmapped = 0; 
    hapInf.numLibraries = 0; 
    hapInf.numReadNames = 0;
    
    // store order of adding haplotypes.
    result.addOrder.push_back(hapIdx);

    // determine which variants have not been called yet
    for (int refIdx = 0; refIdx < haplotypes[hapIdx].getNumReferenceMappings(); ++refIdx)
    {
        const std::vector<DindelVariant> & vars = haplotypes[hapIdx].getVariants(refIdx);

        for (size_t x = 0; x < vars.size(); x++)
        {
            if (DINDEL_DEBUG_3 && haplotypes[hapIdx].getLogMappingProbability(refIdx)>-5.0)
                std::cout << "  adding var[" << x << "]: hapIdx: " << hapIdx << " refIdx: " << refIdx << " " << vars[x].getID() << std::endl;

            DindelRealignWindowResult::VarToInference::iterator vit = result.variantInference.find(vars[x]);
            if (vit == result.variantInference.end())
            {
                std::pair<DindelRealignWindowResult::VarToInference::iterator, bool> it_pair = 
                        result.variantInference.insert(DindelRealignWindowResult::VarToInference::value_type(vars[x],hapInf));

                vit = it_pair.first;
                vit->second.numRealignedReads = int(reads.size());
                vit->second.qual = 0.0;
            }

            DindelRealignWindowResult::Inference & varInf = vit->second;

            HashMap<std::string, int> libraries, readnames;
            varInf.haplotypeIndex.insert(hapIdx);

            if (vars[x].getType()=="INDEL") 
                indelAdded = true;

            //for (std::set<int>::const_iterator r = addReads[hapIdx].begin(); r != addReads[hapIdx].end(); r++) // read pair indices
            for(HashMap<int, double>::const_iterator hit = addReads[hapIdx].begin(); hit != addReads[hapIdx].end(); hit++)
            {
                int r = hit->first;
                double dH = hit->second;
                {
                    int read_idx = r;

                    if (varInf.readIndex.find(read_idx) == varInf.readIndex.end())
                    {

                        // add variant convering statistics
                        if (hapReadAlignments[hapIdx][read_idx].logLik>minLogLikAlignToAlt &&
                            hapReadAlignments[hapIdx][read_idx].nmm<=2 &&
                            hapReadAlignments[hapIdx][read_idx].nmm>=0 &&
                            hapReadAlignments[hapIdx][read_idx].isUngapped)
                        {
                            int dist = haplotypes[hapIdx].getSingleMappingHaplotype(refIdx).getClosestDistance(vars[x],
                                                                                                               hapReadAlignments[hapIdx][read_idx].hapPosLastReadBase-reads[read_idx].length()+1,
                                                                                                               hapReadAlignments[hapIdx][read_idx].hapPosLastReadBase,
                                                                                                               reads[read_idx]);
                            if(dist != -1)
                            {
                                // read overlaps the variant.
                                if (DINDEL_DEBUG_3)
                                {
                                    std::cout << "\t likelihood hapIdx: " << hapIdx << 
                                                 " read_idx: " << read_idx << " lik: " << 
                                                 hapReadAlignments[hapIdx][read_idx].logLik << "\n";
                                }

                                // use it only once if it overlaps. Maybe in some read-haplotype alignments this read doesn't support the variant...
                                varInf.readIndex.insert(read_idx);

                                // add quality in phred scale
                                // note that this quality uses only reads that overlap this position of the variant in the haplotype
                                varInf.qual += dH/.23026;

                                varInf.addDistanceToHistogram(dist);

                                // update statistics using supporting read
                                if (reads[read_idx].isForward()) 
                                    varInf.numReadsForward++; 
                                else 
                                    varInf.numReadsReverse++;

                                if (hapReadAlignments[hapIdx][read_idx].nmm==0)
                                {
                                    if (reads[read_idx].isForward()) 
                                        varInf.numReadsForwardZeroMismatch++; 
                                    else 
                                        varInf.numReadsReverseZeroMismatch++;
                                }

                                if ((*m_pDindelReads)[read_idx].isUnmapped()) 
                                    varInf.numUnmapped++;
                                varInf.addAlignLikToHistogram(hapReadAlignments[hapIdx][read_idx].logLik/double(reads[read_idx].length()));
                                varInf.addMapQToHistogram(reads[read_idx].getMappingQual());
                                libraries[reads[read_idx].getLibraryName()]=1;
                                readnames[reads[read_idx].getID()]=1;
                                
                                // Project the alignment of read to the haplotype onto the reference
                                projectReadAlignmentToReference(haplotypes, read_idx, hapIdx, refIdx);
                            }

                        }
                    }
                }
            } // loop over supporting reads addReads[hapIdx]

            varInf.strandBias = DindelRealignWindowResult::Inference::computeStrandBias(varInf.numReadsForward, varInf.numReadsReverse);
            varInf.numLibraries = int(libraries.size());
            varInf.numReadNames = int(readnames.size());
        }
    } // for (int refIdx = 0; refIdx < haplotypes[hapIdx].getNumReferenceMappings(); ++refIdx)

    if (indelAdded && realignParameters.showCallReads!=0)
    {
        if(DINDEL_DEBUG_3)
        {
            std::cout << "DindelRealignWindow::estimateHaplotypeFrequenciesModelSelectionReadPairs Adding haplotype " << hapIdx << " with INDEL " << std::endl;
            std::cout << "  variants: ";
        }
        for (int refIdx = 0; refIdx < haplotypes[hapIdx].getNumReferenceMappings(); ++refIdx)
        {
            if(DINDEL_DEBUG_3)
                std::cout << "   refmapping[" << refIdx << "]: ";
            const std::vector<DindelVariant> & vars = haplotypes[hapIdx].getVariants(refIdx);
            for (size_t x=0;x<vars.size();x++) std::cout << " " << vars[x].getID(); std::cout << std::endl;
        }

        for (HashMap<int, double>::const_iterator hit = addReads[hapIdx].begin(); hit != addReads[hapIdx].end();hit++)
        {
            int r = hit->first;                
            printReadAlignments(r, std::cout);
        }
    }
}

void DindelRealignWindow::projectReadAlignmentToReference(const std::vector<DindelMultiHaplotype> & haplotypes,
                                                          int readIdx, int hapIdx, int refIdx)
{
    PROFILE_FUNC("DindelRealignWindow::projectReadAlignmentToReference")
    const DindelRead& read = getRead(readIdx);
    std::string read_sequence = read.getSequence();

    const DindelHaplotype& aligned_haplotype = haplotypes[hapIdx].getSingleMappingHaplotype(refIdx);
    const DindelReferenceMapping& reference_mapping = aligned_haplotype.getReferenceMapping();
    const std::string& reference = reference_mapping.refSeq;

    std::string haplotype = aligned_haplotype.getSequence();
    if(reference_mapping.isRC)
    {
        haplotype = reverseComplement(haplotype);
        read_sequence = reverseComplement(read_sequence);
    }
    
    // Query the haplotype to get the reference-to-haplotype alignment
    SequenceOverlap ref2hap_overlap = aligned_haplotype.getReferenceToHaplotypeAlignment();

    // Align the read to the haplotype using overlapper
    SequenceOverlap read2hap_overlap = Overlapper::computeOverlap(haplotype, read_sequence);

    /*
    std::cout << "ref2hap:  " << ref2hap_overlap.cigar << "\n";
    std::cout << "read2hap: " << read2hap_overlap.cigar << "\n";
    std::cout << "Compacted: " << StdAlnTools::compactCigar(StdAlnTools::expandCigar(read2hap_overlap.cigar)) << "\n";
    hap2ref_overlap.printAlignment(haplotype, reference);
    read2hap_overlap.printAlignment(haplotype, read_sequence);
    */

    MultipleAlignment projector_ma;
    projector_ma.addBaseSequence("haplotype", haplotype, "");
    projector_ma.addOverlap("reference", reference, "", ref2hap_overlap);
    projector_ma.addOverlap("read", read_sequence, "", read2hap_overlap);
    //projector_ma.print(500);

    // Calculate the alignment of the read onto the reference sequence
    size_t REF_ROW = 1;
    size_t READ_ROW = 2;
    size_t num_columns = projector_ma.getNumColumns();
    size_t read_offset = 0;

    // Find the bounds of the alignment of the read on the reference
    size_t read_align_start = 0;
    size_t read_align_end = num_columns - 1;
    size_t ref_bases_pre_skipped = ref2hap_overlap.match[1].start;
    size_t read_bases_pre_skipped = 0;
    size_t read_bases_post_skipped = 0;

    // Find the start
    while(read_align_start < num_columns)
    {
        char ref_symbol = projector_ma.getSymbol(REF_ROW, read_align_start);
        char read_symbol = projector_ma.getSymbol(READ_ROW, read_align_start);
        bool is_read_base = read_symbol != '\0' && read_symbol != '-';
        bool is_ref_base = ref_symbol != '\0' && ref_symbol != '-';

        if(is_read_base && is_ref_base)
            break; // start found

        if( is_read_base && ref_symbol == '\0')
            read_bases_pre_skipped += 1;
        if( is_ref_base && read_symbol == '\0')
            ref_bases_pre_skipped += 1;

        read_align_start += 1;
    }

    // Find the end
    bool end_found = false;
    size_t i = read_align_start;
    while(i < num_columns)
    {
        char ref_symbol = projector_ma.getSymbol(REF_ROW, i);
        char read_symbol = projector_ma.getSymbol(READ_ROW, i);
        bool is_read_base = read_symbol != '\0' && read_symbol != '-';

        if( (ref_symbol == '\0' || read_symbol == '\0') && !end_found)
        {
            end_found = true;
            read_align_end = i - 1;
        }
        
        if( is_read_base && ref_symbol == '\0')
            read_bases_post_skipped += 1;

        i += 1;
    }

    // Fill in the cigar.
    std::string expanded_cigar;
    for(size_t i = read_align_start; i <= read_align_end; ++i)
    {
        char ref_symbol = projector_ma.getSymbol(REF_ROW, i);
        char read_symbol = projector_ma.getSymbol(READ_ROW, i);

        // Within the range of columns there should be a proper alignment between the read and reference
        // with no skipped bases.
        assert(ref_symbol != '\0');
        assert(read_symbol != '\0');

        if(read_symbol == '-' && ref_symbol != '-')
            expanded_cigar.push_back('D');
        else if(read_symbol != '-' && ref_symbol == '-')
            expanded_cigar.push_back('I');
        else if(read_symbol != '-' && ref_symbol != '-')
            expanded_cigar.push_back('M');
    }

    // Add softclipping to the cigar if the read was not completely aligned to the haplotype
    size_t start_clip = read2hap_overlap.match[1].start + read_bases_pre_skipped;
    if(start_clip > 0)
        expanded_cigar.insert(0, start_clip, 'S');
    
    size_t end_clip = read_sequence.size() - (read2hap_overlap.match[1].end + 1) + read_bases_post_skipped;
    if(end_clip > 0)
        expanded_cigar.append(end_clip, 'S');
    
    DindelReadReferenceAlignment drra;
    drra.dindel_ref_index = refIdx;
    drra.dindel_read_index = readIdx;
    drra.cigar = StdAlnTools::compactCigar(expanded_cigar);
    drra.read_name = read.getID();
    drra.reference_name = reference_mapping.refName;
    drra.reference_start_position = reference_mapping.refStart + read_offset + ref_bases_pre_skipped;

    // DindelReads are on the same strand as the haplotype. 
    // Write the read sequence field as the original sequencing strand
    drra.read_sequence = read.isForward() ? read.getSequence() : reverseComplement(read.getSequence());

    // Calculate the the read's strand with respect to the refernce
    bool rc_to_ref = read.isForward() ? reference_mapping.isRC : !reference_mapping.isRC;
    drra.is_reference_reverse_strand = rc_to_ref;
    m_readReferenceAlignments.push_back(drra);
}   

void DindelRealignWindow::setAddReadsMatePairs(int type,
                                      int h,
                                      const std::vector< std::vector<double> > & hrLik,
                                      std::vector< HashMap<int, double> > & addReads,
                                      std::vector<double> & addLL,
                                      int numReadPairs,
                                      int numHaps,
                                      const std::vector<int> & added)
{

    if (type == 0)
    {
       for (int r=0;r<numReadPairs;r++) 
       {
           //double max_ll = -HUGE_VAL;
           //for (int h1=0;h1<numHaps;h1++) if (h1!=h && hrLik[r][h1]>max_ll) max_ll = hrLik[r][h1];
           double dH = (hrLik[r][h] - ( - (*m_pDindelReads)[r].getMappingQual()*.23026-2.0*log(double(DINDEL_HMM_BANDWIDTH))));
           if (dH>2.0) addReads[h][r] = dH;
           if (dH>0.0) addLL[h] += dH;

       }
    } else
    {

       for (int r=0;r<numReadPairs;r++)
       {

           //bool print = false;
           if (!added[h])
           {
               // find the called haplotype that gives the best likelihood.
               double min_ll_diff = HUGE_VAL;
               for (int h1=0;h1<numHaps;h1++) if (added[h1])
               {
                   double tdiff = hrLik[r][h]-hrLik[r][h1];
                   if (tdiff<min_ll_diff)
                   {
                       min_ll_diff = tdiff;
                   }
               }

               double diff = min_ll_diff;
               if (diff > 0.0)
               {
                   addLL[h] += diff;
                   if (diff > 2.0)
                   {
                       addReads[h][r] = diff;
                   }
               }
           }
       }
    }
}

void DindelRealignWindow::setAddReadsSingleRead(int type,
                                      int h,
                                      const std::vector< std::vector<double> > & hrLik,
                                      std::vector< HashMap<int, double> > & addReads,
                                      std::vector<double> & addLL,
                                      int numReads,
                                      int numHaps,
                                      const std::vector<int> & added)
{

    if (type == 0)
    {
       for (int r=0;r<numReads;r++)
       {
           //double max_ll = -HUGE_VAL;
           //for (int h1=0;h1<numHaps;h1++) if (h1!=h && hrLik[r][h1]>max_ll) max_ll = hrLik[r][h1];
           double dH = (hrLik[r][h] - ( - (*m_pDindelReads)[r].getMappingQual()*.23026-log(double(DINDEL_HMM_BANDWIDTH))));
           if (dH>2.0) addReads[h][r] = dH;
           if (dH>0.0) addLL[h] += dH;

       }
    } else
    {

       for (int r=0;r<numReads;r++)
       {

           //bool print = false;
           if (!added[h])
           {
               // find the called haplotype that gives the best likelihood.
               double min_ll_diff = HUGE_VAL;
               for (int h1=0;h1<numHaps;h1++) if (added[h1])
               {
                   double tdiff = hrLik[r][h]-hrLik[r][h1];
                   if (tdiff<min_ll_diff)
                   {
                       min_ll_diff = tdiff;
                   }
               }

               double diff = min_ll_diff;
               if (diff > 0.0)
               {
                   addLL[h] += diff;
                   if (diff > 2.0)
                   {
                       addReads[h][r] = diff;
                   }
               }
           }
       }
    }
}


void DindelRealignWindow::setAddReadsSingleReadMultiSample(int type,
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
                                      const HashMap<std::string, std::list<int> > & sampleToReads)
{
    (void) hrLik;
    (void) numReads;
    (void) type;

    std::vector<int> allowed(numHaps, 0);
    std::vector<double> hapFreqs(numHaps,0.0);

    int numAdded = 0;
    for (size_t h = 0; h < added.size(); h++)
    {
        allowed[h] = added[h];
        if (added[h]!=0)
            numAdded++;
    }

    double scale = double(numAdded)/double(numAdded+1);

    for (size_t h = 0; h < added.size(); h++)
        hapFreqs[h] = hapFreqPrevious[h]*scale; //FIXME. Initialize with previous frequency estimates

    assert(htest<numHaps);
    assert(added[htest]==0);

    hapFreqs[htest] = 1.0-scale;
    allowed[htest] = 1;

    int numSamples = int(sampleToReads.size());
    assert (int(zindNew.size()) == numSamples );
    
    std::vector<double> logPairFreqs;
    std::vector<int> allowedGenotypes;

    int GT = int(llHapPairs.size())/numSamples;
    assert ( GT * numSamples == int(llHapPairs.size()));

    // estimate haplotype frequencies
    if (DINDEL_DEBUG_3)
    {
        std::cout << "setAddReadsSingleReadMultiSample [allowed, initFreq]: ";
        for (int h = 0; h < numHaps; h++)
            std::cout << " [" << allowed[h] << "," << hapFreqs[h] << "]";
        std::cout << "\n";
    }
    doEMMultiSample(numSamples, 15, llHapPairs, allowed, hapFreqs, hapFreqs);


    // compute improvement to the likelihood

    int gt = 0;
    for (int h1 = 0; h1 < numHaps; h1++)
        for (int h2 = h1; h2 < numHaps; h2++)
        {
            if (allowed[h1] && allowed[h2])
            {
                allowedGenotypes.push_back(gt);
                logPairFreqs.push_back(log(hapFreqs[h1])+log(hapFreqs[h2]));
            }
            gt++;
        }

    assert(gt == GT);
    
    double z = 0.0;
    int s = 0;
    bool addToReads = true;
    for (HashMap<std::string, std::list<int> >::const_iterator sit = sampleToReads.begin(); sit != sampleToReads.end(); sit++, s++)
    {
        double zind = -HUGE_VAL;
        int idx = s*GT;
        for (size_t x = 0; x < allowedGenotypes.size(); x++)
        {
            gt = allowedGenotypes[x];
            double pf = logPairFreqs[x];
            zind = addLogs(zind, llHapPairs[idx+gt]+pf);
        }

        double diff = zind - zindPrevious[s];
        if (DINDEL_DEBUG_3)
            std::cout << "\tzind[" << s << "]: " << zind << " zindPrevious: " << zindPrevious[s] << "\n";
        if (addToReads && diff>0.0)
        {
            if (numAdded == 0) {
                for (std::list<int>::const_iterator rit = sit->second.begin(); rit != sit->second.end(); rit++)
                {
                    int r = *rit;
                    double min = HUGE_VAL;
                    for (int h = 0; h < numHaps; h++)
                        if (h!=htest && hrLik[r][h]<min)
                            min = hrLik[r][h];
                    double rdiff = hrLik[r][htest]-min;
                    if (rdiff>2.0)
                        addReads[htest][r] = rdiff;
                }
            } else {
                for (std::list<int>::const_iterator rit = sit->second.begin(); rit != sit->second.end(); rit++)
                {
                    int r = *rit;
                    double max = -HUGE_VAL;
                    for (int h = 0; h < numHaps; h++)
                    {
                        if (added[h] && hrLik[r][h]>max)
                            max = hrLik[r][h];
                    }

                    double rdiff = hrLik[r][htest]-max;
                    if (rdiff>2.0)
                        addReads[htest][r] = rdiff;
                }
            }
        }

        zindNew[s] = zind;
        
        z += zind;
    }

    newLL[htest] = z;
    hapFreqNew = hapFreqs;
    
}

void DindelRealignWindow::computeAddLLMatePairs(int type,
                                       int h,
                                       const std::vector< std::vector<double> > & hrLik,
                                       std::vector<double> & addLL,
                                       int numReadPairs,
                                       int numHaps,
                                       const std::vector<int> & added)
{

    if (type == 0)
    {
       for (int r=0;r<numReadPairs;r++)
       {
           double dH = (hrLik[r][h] - ( - (*m_pDindelReads)[r].getMappingQual()*.23026-2.0*log(double(DINDEL_HMM_BANDWIDTH))));
           if (dH>0.0) addLL[h] += dH;
       }
    } else
    {

       for (int r=0;r<numReadPairs;r++)
       {

           //bool print = false;
           if (!added[h])
           {
               // find the called haplotype that gives the best likelihood.
               double min_ll_diff = HUGE_VAL;
               for (int h1=0;h1<numHaps;h1++) if (added[h1])
               {
                   double tdiff = hrLik[r][h]-hrLik[r][h1];
                   if (tdiff<min_ll_diff)
                   {
                       min_ll_diff = tdiff;
                   }
               }

               double diff = min_ll_diff;
               if (diff > 0.0)
               {
                   addLL[h] += diff;
               }
           }
       }
    }
}

void DindelRealignWindow::computeAddLLSingleRead(int type,
                                       int h,
                                       const std::vector< std::vector<double> > & hrLik,
                                       std::vector<double> & addLL,
                                       int numReads,
                                       int numHaps,
                                       const std::vector<int> & added)
{

    if (type == 0)
    {
       for (int r=0;r<numReads;r++)
       {
           double dH = (hrLik[r][h] - ( - (*m_pDindelReads)[r].getMappingQual()*.23026-log(double(DINDEL_HMM_BANDWIDTH))));
           if (dH>0.0) addLL[h] += dH;
       }
    } else
    {

       for (int r=0;r<numReads;r++)
       {

           //bool print = false;
           if (!added[h])
           {
               // find the called haplotype that gives the best likelihood.
               double min_ll_diff = HUGE_VAL;
               for (int h1=0;h1<numHaps;h1++) if (added[h1])
               {
                   double tdiff = hrLik[r][h]-hrLik[r][h1];
                   if (tdiff<min_ll_diff)
                   {
                       min_ll_diff = tdiff;
                   }
               }

               double diff = min_ll_diff;
               if (diff > 0.0)
               {
                   addLL[h] += diff;
               }
           }
       }
    }
}


void DindelRealignWindow::doReadHaplotypeAlignment(int H, const std::vector<DindelRead> & dReads)
{
    const std::vector<DindelMultiHaplotype> & haplotypes = m_dindelWindow.getHaplotypes();

    // globally align haplotypes to the first haplotype (arbitrary)
    
    for (int h = H; h == H; ++h)
    {
        std::cout << "ALIGNING EVERYTHING AGAINST HAPLOTYPE " << h << "\n";
        MultipleAlignment ma;
        const std::string rootSequence = haplotypes[h].getSequence();
        std::stringstream h_ss;
        h_ss << "haplotype-" << h;
        ma.addBaseSequence("root", rootSequence, "");

        for(size_t r = 0; r < dReads.size(); ++r)
        {
            std::stringstream r_ss;
            r_ss << "read-" << r << "("  << dReads[r].getID() << ")";
            SequenceOverlap overlap = Overlapper::computeOverlap(rootSequence, dReads[r].getSequence());
            ma.addOverlap(r_ss.str(), dReads[r].getSequence(), "", overlap);
        }
        ma.print(100000);
    }
}

void DindelRealignWindow::showHaplotypeOnlyReadsMatePairs(int h,
                                       const std::vector< std::vector<double> > & hrLik,
                                       int numReadPairs,
                                       int numHaps)
{
   for (int r=0;r<numReadPairs;r++)
   {
       // find the called haplotype that gives the best likelihood.
       double min_ll_diff = HUGE_VAL;
       for (int h1=0;h1<numHaps;h1++) if (h1!=h)
       {
           double tdiff = hrLik[r][h]-hrLik[r][h1];
           if (tdiff<min_ll_diff)
               min_ll_diff = tdiff;
       }

       double diff = min_ll_diff;
       if (diff > 2.0)
       {
           {
               const ReadHaplotypeAlignment & rha = hapReadAlignments[h][r];

               std::cout << " read[" << r << "]: hap[" << h << "] min_diff: " << diff 
                         << " logLik: " << rha.logLik << " nmm: " << rha.nmm 
                         << " isUngapped: " << rha.isUngapped << " " 
                         << (*m_pDindelReads)[r].getSequence();

               for (int i = 0; i < numHaps; i++)
                   std::cout << " lik[" << i << "]: " << hapReadAlignments[i][r].logLik;

               std::cout << "\n";
           }

           {
               const ReadHaplotypeAlignment & rha = hapReadAlignments[h][r+numReadPairs];

               std::cout << " read[" << r << "]: hap[" << h << "] min_diff: " << diff 
                         << " logLik: " << rha.logLik << " nmm: " << rha.nmm 
                         << " isUngapped: " << rha.isUngapped << " " 
                         << (*m_pDindelReads)[r+numReadPairs].getSequence();

               for (int i = 0; i < numHaps; i++)
                   std::cout << " lik[" << i << "]: " << hapReadAlignments[i][r+numReadPairs].logLik;
               std::cout << "\n";
           }

           std::vector<DindelRead> dReads;
           dReads.push_back(this->m_pDindelReads->at(r));
           dReads.push_back(this->m_pDindelReads->at(r+numReadPairs));
           for (int k = 0; k < numHaps; k++) {
               std::cout << "Aligning reads to haplotype " << k << " for HAPLOTYPE-ONLY reads for haplotype " << h << "\n";
               doReadHaplotypeAlignment(k, dReads);
           }
       }
   }


}

void DindelRealignWindow::showHaplotypeOnlyReadsSingleRead(int h,
                                       const std::vector< std::vector<double> > & hrLik,
                                       int numReads,
                                       int numHaps)
{
   for (int r=0;r<numReads;r++)
   {

       //bool print = false;
       // find the called haplotype that gives the best likelihood.
       double min_ll_diff = HUGE_VAL;
       for (int h1=0;h1<numHaps;h1++) if (h1!=h)
       {
           double tdiff = hrLik[r][h]-hrLik[r][h1];
           if (tdiff<min_ll_diff)
           {
               min_ll_diff = tdiff;
           }
       }

       double diff = min_ll_diff;
       if (diff > 2.0)
       {
           {
               const ReadHaplotypeAlignment & rha = hapReadAlignments[h][r];

               std::cout << " read[" << r << "]: hap[" << h << "] min_diff: " << diff << " logLik: " << rha.logLik << " nmm: " << rha.nmm << " isUngapped: " << rha.isUngapped << " " << (*m_pDindelReads)[r].getSequence();
               for (int i = 0; i < numHaps; i++)
                   std::cout << " lik[" << i << "]: " << hapReadAlignments[i][r].logLik;
               std::cout << "\n";
           }

           std::vector<DindelRead> dReads;
           dReads.push_back(this->m_pDindelReads->at(r));
           for (int k = 0; k < numHaps; k++) {
               std::cout << "Aligning reads to haplotype " << k << " for HAPLOTYPE-ONLY reads for haplotype " << h << "\n";
               doReadHaplotypeAlignment(k, dReads);
           }
       }
   }


}


DindelRealignWindowResult DindelRealignWindow::estimateHaplotypeFrequenciesModelSelectionMatePairs(double minLogLikAlignToRef, 
                                                                                                   double minLogLikAlignToAlt,
                                                                                                   bool capUsingMappingQuality,
                                                                                                   const DindelRealignWindowResult *pPreviousResult,
                                                                                                   bool print = false)
{
   assert(realignParameters.realignMatePairs);
   (void) pPreviousResult;
   if (DINDEL_DEBUG)
       std::cout << "DindelRealignWindow::estimateHaplotypeFrequencies START" << std::endl;
   
   const std::vector<DindelMultiHaplotype> & haplotypes = m_dindelWindow.getHaplotypes();
   if (hapReadAlignments.size() != haplotypes.size()) throw std::string("DindelRealignWindow::estimateHaplotypeFrequencies incomplete haplotype alignments.");

   // init result
   DindelRealignWindowResult result(*this);

// construct read X haplotype vector
   int numReads = int(m_pDindelReads->size());

   // it is assumed that the first half are the first read of each fragment, and the second half their mates
   assert(numReads%2==0);
   int numReadPairs = numReads/2;

   int numHaps = (int) haplotypes.size();

   std::vector< std::vector<double> > hrLik( numReadPairs, std::vector<double>(numHaps, realignParameters.minLogLikAlignToAlt));
   
   for (int r=0;r<numReadPairs;r++)
   {
       // Currently we do not take into account the insert size distribution.
       // For tandem repeats and long insertions and deletions this would be desirable.
       for (int h=0;h<numHaps;h++)
       {
           double matePairLogLik = hapReadAlignments[h][r].logLik+hapReadAlignments[h][r+numReadPairs].logLik;

           if (haplotypes[h].isReference())
               hrLik[r][h]=addLogs(matePairLogLik, minLogLikAlignToRef);
           else
               hrLik[r][h]=addLogs(matePairLogLik, minLogLikAlignToAlt);

           if (capUsingMappingQuality) hrLik[r][h] = addLogs(hrLik[r][h], -(*m_pDindelReads)[r].getMappingQual()*.23026-2.0*log(double(DINDEL_HMM_BANDWIDTH)));

	   if (DINDEL_DEBUG) std::cout << "mpl: " << hapReadAlignments[h][r].logLik << " " << hapReadAlignments[h][r+numReadPairs].logLik << " " << " hrLik[r][h]: " << hrLik[r][h] << std::endl;
           
       }
       // std::cout << " hrLik[" << r << "]: " << hrLik[r][0] << " " << hrLik[r][1] << " 1>0 " << int(hrLik[r][0]<hrLik[r][1]) << "\n";
   }
   

   if (DEBUG_CALLINDEL)
   {
       for (int r=0;r<numReads;r++)
       {
            printReadAlignments(r, std::cout,10,true);
       }
   }

   // now add haplotypes until we have obtained a maximal set of hapltoypes

   std::vector<int> added(haplotypes.size(),0);
   std::vector<int> addedInIteration(haplotypes.size(),-1);


   size_t numAdded = 0;
   bool done = false;

   // If we are doing normal/tumour style calling the set of called haplotypes may be initialized to have haplotypes already 'added' from the normal sample
   /*
   if(pPreviousResult != NULL)
   {
       if(this->realignParameters.graphDiffStyle == 1)
       {
           numAdded = 1;
           assert(pPreviousResult->addOrder.size()>=1);
           assert(haplotypes.size() == pPreviousResult->haplotypes.size());
           added[ pPreviousResult->addOrder[0] ] = 1;
           result.addOrder.push_back( pPreviousResult->addOrder[0] );
       } else
       {
           assert(false);
       }
    }
   */


   std::vector<double> addLL(numHaps, 0.0); // loglikelihood gained by adding haplotypes
   std::vector<double> haplotypeQual(numHaps, 0.0); // loglikelihood gained by adding haplotypes
   std::vector<double> haplotypeFrequencies(numHaps, 0.0);

   std::vector< HashMap<int, double> > addReads(numHaps); // read pairs supporting a haplotype

   int iteration = 0;
   bool stopLoop = false;
   while (!stopLoop)
   {
       if(DINDEL_DEBUG_3)
           std::cout << "Starting iteration " << iteration << "\n";
       if (numAdded == haplotypes.size()) break; // always keep at least one haplotype

       // compute for each haplotype the penalty for removing it.

       // only itialize kept haplotypes.
       for (int h=0;h<numHaps;h++) if (!added[h])
       {
          addLL[h]=0.0;
          addReads[h].clear();
       }

       int addType = (numAdded == 0)?0:1; // 1 means quality will depend on difference to existing haplotypes
       for (int h = 0; h < numHaps; h++)
           setAddReadsMatePairs(addType, h, hrLik, addReads, addLL, numReadPairs, numHaps, added);
       

       std::multiset<HL> orderedAddCandidates; // ordered by number of new variants added


       for (int h=0;h<numHaps;h++) if (!added[h])
       {
           // if we would add haplotype h, what would be the increase in the likelihood?
           // FIXME very simple prior..
           double logPrior = this->realignParameters.logPriorAddHaplotype;
           double qual = logPrior + addLL[h];
           if (DINDEL_DEBUG_3)
               std::cout << "estimateHaplotypeFrequenciesModelSelectionMatePairs: estimate haplotype " << h << " logPrior  " << logPrior << " haplotype qual: " << qual << std::endl;
      
           haplotypeQual[h] = qual;
           orderedAddCandidates.insert(HL(qual, h));
       }
       
       done = true;
       std::multiset<HL>::const_reverse_iterator iter = orderedAddCandidates.rbegin();
       if (iter != orderedAddCandidates.rend())
       {
           if (iter->ll>=5.0)
           {
                // this haplotype will be added
                int hapIdx = iter->idx;
                double ll = iter->ll ; // use addLL because haplotypeQual was only set for the called haplotypes
                done = false;
                added[hapIdx]=1;
             
                // see which other haplotypes yield a similar improvement to the likelihood
                std::vector<int> idxSimilar;
                iter++;
                while (iter != orderedAddCandidates.rend())
                {
                    if (DINDEL_DEBUG_3)
                        std::cout << "iter->idx: " << iter->idx << " ll: " << iter->ll << "\n";
                    if (ll - iter->ll <= 20.0 && iter->ll>5.0)
                    {
                        idxSimilar.push_back(iter->idx);
                        if (DINDEL_DEBUG_3)
                            std::cout << "Haplotype " << iter->idx << " has a similar likelihood.\n";
                    }
                    iter++;
                }

                // also call these haplotypes and renormalize scores
                std::vector<int> idxCallHaplotypes;
                idxCallHaplotypes.push_back(hapIdx);

                if (idxSimilar.size()!=0)
                {
                    // chech which ones of these would be explained away by adding hapIdx
                    std::vector<double> testAddLL(addLL.size(),0.0);
                    for (int h = 0; h < numHaps; h++)
                        computeAddLLMatePairs(1, h, hrLik, testAddLL, numReadPairs, numHaps, added);

                    double logZ = -HUGE_VAL;
                    for(size_t x = 0; x < idxSimilar.size();++x)
                    {
                        int h = idxSimilar[x];
                        // only call haplotype h in this iteration if testAddLL[h] is small,
                        // which means that h and hapIdx explain more or less the same set of reads
                        if (DINDEL_DEBUG_3)
                            std::cout << " haplotype " << h << " : " << testAddLL[h] << "\n";
                        if (testAddLL[h] <= 10.0)
                        {
                            idxCallHaplotypes.push_back(h);
                            if (DINDEL_DEBUG_3)
                                std::cout << "Also adding haplotype " << h << " with likelihood " << testAddLL[h] << "\n";
                        }
                    }


                    
                    for(size_t x = 0; x < idxCallHaplotypes.size();++x)
                        logZ = addLogs(logZ, this->realignParameters.logPriorAddHaplotype + addLL[idxCallHaplotypes[x]]);

                    for(size_t x = 0; x < idxCallHaplotypes.size();++x)
                    {
                        int _idx = idxCallHaplotypes[x];
                        double t = 1.0-exp(this->realignParameters.logPriorAddHaplotype + addLL[_idx] - logZ);
                        if (t<1e-100) t = 1e-100;
                        double renormQual = -log(t);
                        if (DINDEL_DEBUG_3)
                            std::cout << " original HapQual: " << haplotypeQual[_idx] << " renormQual: " << renormQual;
                        if (renormQual<haplotypeQual[_idx])
                            haplotypeQual[_idx] = renormQual;
                        if (DINDEL_DEBUG_3)
                            std::cout << " new qual: " << haplotypeQual[_idx] << "\n";


                    }
                }

                for(size_t x = 0; x < idxCallHaplotypes.size();++x)
                {
                    int _idx = idxCallHaplotypes[x];
                    added[_idx] = 1;
                    addedInIteration[_idx] = iteration; // store in which iteration the haplotype was added.
                    addCalledHaplotypeMatePairs(_idx, result, haplotypes, addReads, minLogLikAlignToAlt, numReadPairs);
                    numAdded++;

                    if (DINDEL_DEBUG_3)
                    {
                        std::cout << "Haplotype[" << _idx << "] ONLY reads\n";
                        showHaplotypeOnlyReadsMatePairs(_idx, hrLik, numReadPairs, numHaps);
                    }

                }

           }
       }


       if(done)
       {
           /*
           if (pPreviousResult != NULL)
           {
              if(this->realignParameters.graphDiffStyle == 1)
              {
                  // set values for first haplotype
                  int h = result.addOrder[0];
                  setAddReads(1, h, hrLik, addReads, addLL, numReadPairs, numHaps, added);
                  double qual = this->realignParameters.logPriorAddHaplotype + addLL[h];
                  haplotypeQual[h] = qual;
                  addCalledHaplotype(h, result, haplotypes, addReads, minLogLikAlignToAlt, numReadPairs);
              }
           }
           */
           stopLoop = true;
       }


       iteration++;
   }

   if (DINDEL_DEBUG_3)
       std::cout << "Number of haplotypes added: " << numAdded << "\n";

   // estimate haplotype frequencies for called haplotypes using EM.

   doEM(hrLik, added, haplotypeFrequencies);

   if (!QUIET) std::cout << "DindelRealignWindow::estimateHaplotypeFrequencies Added " << numAdded << " out of " << haplotypes.size() << " haplotypes." << std::endl;

   if (realignParameters.showCallReads && print)
   {
       for (int h=0;h<numHaps;h++)
       {
               std::cout << "DindelRealignWindow::estimateHaplotypeFrequencies final penalty haplotype[" <<h << "]: " << haplotypeQual[h] << std::endl;
       }

   }


   const std::vector< DindelReferenceMapping > & refMappings = this->m_dindelWindow.getReferenceMappings();
   result.weightedRefMapFrequencies = std::vector<double>(refMappings.size(), 0.0);

   for (int hapIdx = 0; hapIdx < numHaps; ++hapIdx)
   {
       double hq = haplotypeQual[hapIdx];
       if (hq<0.)
           hq = 0.0;
       double probHap = 1.0-exp(-hq);
       double hapFreq = haplotypeFrequencies[hapIdx];

       for (int refIdx = 0; refIdx < haplotypes[hapIdx].getNumReferenceMappings(); ++refIdx)
       {
           double mapProb = exp(haplotypes[hapIdx].getLogMappingProbability(refIdx));
           result.weightedRefMapFrequencies[refIdx] += hapFreq * probHap * mapProb;
       }
   }

   if (DINDEL_DEBUG_3)
   {
       for (size_t x = 0; x < result.weightedRefMapFrequencies.size();x++) {
           std::cout << " refMapping #" << x << " : weighted hapFreq " << result.weightedRefMapFrequencies[x] << "\n";
       }

   }

   // Add haplotype properties for all variants
   for (int hapIdx = 0; hapIdx < numHaps; ++hapIdx)
   {
       double hapQual = haplotypeQual[hapIdx]/.23026;
       double hapFreq = haplotypeFrequencies[hapIdx];

        for (int refIdx = 0; refIdx < haplotypes[hapIdx].getNumReferenceMappings(); ++refIdx)
        {
            const std::vector<DindelVariant> & vars = haplotypes[hapIdx].getVariants(refIdx);

            for (size_t x = 0; x < vars.size(); x++)
            {
                DindelRealignWindowResult::VarToInference::iterator vit = result.variantInference.find(vars[x]);
                if (vit == result.variantInference.end())
                {
                    // should set quality and freq to zero.
                    DindelRealignWindowResult::Inference varInf;
                    varInf.haplotypeProperties.push_back(DindelRealignWindowResult::HaplotypeProperties(haplotypes[hapIdx].getLogMappingProbability(refIdx), hapQual, (result.weightedRefMapFrequencies[refIdx]<1e-6)?0:hapFreq/result.weightedRefMapFrequencies[refIdx], addedInIteration[hapIdx]));
                    varInf.numCalledHaplotypes = numAdded;
                    result.variantInference[vars[x]] = varInf;
                }
                else
                {
                    DindelRealignWindowResult::Inference & varInf = vit->second;
                    varInf.haplotypeProperties.push_back(DindelRealignWindowResult::HaplotypeProperties(haplotypes[hapIdx].getLogMappingProbability(refIdx), hapQual, (result.weightedRefMapFrequencies[refIdx]<1e-6)?0:hapFreq/result.weightedRefMapFrequencies[refIdx], addedInIteration[hapIdx]));
                    varInf.numCalledHaplotypes = numAdded;
                }

            }
        }

   }



  // assign result


   result.haplotypeFrequencies = haplotypeFrequencies;
   size_t totReads = 0; for (size_t x=0;x<addReads.size();x++) totReads += addReads[x].size();
   result.numHapSpecificReads = (int) totReads;
   result.numCalledHaplotypes = (int) numAdded; // -1 makes sure we don't count the reference haplotype. CHANGE when we are doing proper genotyping, because then we align each read to each haplotype

   result.outputID = this->m_outputID;

   return result;

}

DindelRealignWindowResult DindelRealignWindow::estimateHaplotypeFrequenciesModelSelectionSingleReads(double minLogLikAlignToRef,
                                                                                                   double minLogLikAlignToAlt,
                                                                                                   bool capUsingMappingQuality,
                                                                                                   const DindelRealignWindowResult *pPreviousResult,
                                                                                                   bool print = false)
{
   
   (void) pPreviousResult;
   if (DINDEL_DEBUG)
       std::cout << "DindelRealignWindow::estimateHaplotypeFrequencies START" << std::endl;

   const std::vector<DindelMultiHaplotype> & haplotypes = m_dindelWindow.getHaplotypes();
   if (hapReadAlignments.size() != haplotypes.size()) throw std::string("DindelRealignWindow::estimateHaplotypeFrequencies incomplete haplotype alignments.");

   // init result
   DindelRealignWindowResult result(*this);

   // construct read X haplotype vector
   int numReads = int(m_pDindelReads->size());
   int numHaps = (int) haplotypes.size();

   std::vector< std::vector<double> > hrLik(numReads, std::vector<double>(numHaps, realignParameters.minLogLikAlignToAlt));

   for (int r=0;r<numReads;r++)
   {
       // Currently we do not take into account the insert size distribution.
       // For tandem repeats and long insertions and deletions this would be desirable.
       for (int h=0;h<numHaps;h++)
       {
           double readLogLik = hapReadAlignments[h][r].logLik;
           if (haplotypes[h].isReference())
               hrLik[r][h]=addLogs(readLogLik, minLogLikAlignToRef);
           else
               hrLik[r][h]=addLogs(readLogLik, minLogLikAlignToAlt);

           hrLik[r][h] = readLogLik;

           if (capUsingMappingQuality) 
               hrLik[r][h] = addLogs(hrLik[r][h], -(*m_pDindelReads)[r].getMappingQual()*.23026-log(double(DINDEL_HMM_BANDWIDTH)));

       }
       // std::cout << " hrLik[" << r << "]: " << hrLik[r][0] << " " << hrLik[r][1] << " 1>0 " << int(hrLik[r][0]<hrLik[r][1]) << "\n";
   }


   if (DEBUG_CALLINDEL)
   {
       for (int r=0;r<numReads;r++)
       {
            printReadAlignments(r, std::cout,10,true);
       }
   }

   // now add haplotypes until we have obtained a maximal set of hapltoypes

   std::vector<int> added(haplotypes.size(),0);
   std::vector<int> addedInIteration(haplotypes.size(),-1);


   size_t numAdded = 0;
   bool done = false;

   // If we are doing normal/tumour style calling the set of called haplotypes may be initialized to have haplotypes already 'added' from the normal sample
   /*
   if(pPreviousResult != NULL)
   {
       if(this->realignParameters.graphDiffStyle == 1)
       {
           numAdded = 1;
           assert(pPreviousResult->addOrder.size()>=1);
           assert(haplotypes.size() == pPreviousResult->haplotypes.size());
           added[ pPreviousResult->addOrder[0] ] = 1;
           result.addOrder.push_back( pPreviousResult->addOrder[0] );
       } else
       {
           assert(false);
       }
    }
   */


   std::vector<double> addLL(numHaps, 0.0); // loglikelihood gained by adding haplotypes
   std::vector<double> haplotypeQual(numHaps, 0.0); // loglikelihood gained by adding haplotypes
   std::vector<double> haplotypeFrequencies(numHaps, 0.0);

   std::vector< HashMap<int, double> > addReads(numHaps); // read pairs supporting a haplotype

   int iteration = 0;
   bool stopLoop = false;
   while (!stopLoop)
   {
       if(DINDEL_DEBUG_3)
           std::cout << "Starting iteration " << iteration << "\n";
       if (numAdded == haplotypes.size()) break; // always keep at least one haplotype

       // compute for each haplotype the penalty for removing it.

       // only itialize kept haplotypes.
       for (int h=0;h<numHaps;h++) if (!added[h])
       {
          addLL[h]=0.0;
          addReads[h].clear();
       }

       int addType = (numAdded == 0)?0:1; // 1 means quality will depend on difference to existing haplotypes
       for (int h = 0; h < numHaps; h++)
           setAddReadsSingleRead(addType, h, hrLik, addReads, addLL, numReads, numHaps, added);


       std::multiset<HL> orderedAddCandidates; // ordered by number of new variants added


       for (int h=0;h<numHaps;h++) if (!added[h])
       {
           // if we would add haplotype h, what would be the increase in the likelihood?
         // FIXME very simple prior..
           double logPrior = this->realignParameters.logPriorAddHaplotype;

           double qual = logPrior + addLL[h];
           if (DINDEL_DEBUG_3)
               std::cout << "estimateHaplotypeFrequenciesModelSelectionMatePairs: estimate haplotype " << h << " logPrior  " << logPrior << " haplotype qual: " << qual << std::endl;

           haplotypeQual[h] = qual;

           orderedAddCandidates.insert(HL(qual, h));
       }


       /*
       if (realignParameters.showCallReads)
       {
           for (size_t h=0;h<haplotypes.size();h++)
           {

               std::cout << "Haplotype[" << h << "]: ";
               const std::vector<DindelVariant> & vars = haplotypes[h].getVariants();
               for (size_t x=0;x<vars.size();x++) std::cout << " " << vars[x].getID();
               std::cout << std::endl;
           }
       }
       */


       done = true;
       std::multiset<HL>::const_reverse_iterator iter = orderedAddCandidates.rbegin();
       if (iter != orderedAddCandidates.rend())
       {

           if (iter->ll>=5.0)
           {
                // this haplotype will be added
                int hapIdx = iter->idx;
                double ll = iter->ll ; // use addLL because haplotypeQual was only set for the called haplotypes
                done = false;
                added[hapIdx]=1;

                // see which other haplotypes yield a similar improvement to the likelihood
                std::vector<int> idxSimilar;
                iter++;
                while (iter != orderedAddCandidates.rend())
                {
                    if (DINDEL_DEBUG_3)
                        std::cout << "iter->idx: " << iter->idx << " ll: " << iter->ll << "\n";
                    if (ll - iter->ll <= 20.0 && iter->ll>5.0)
                    {
                        idxSimilar.push_back(iter->idx);
                        if (DINDEL_DEBUG_3)
                            std::cout << "Haplotype " << iter->idx << " has a similar likelihood.\n";
                    }
                    iter++;
                }


                // also call these haplotypes and renormalize scores
                std::vector<int> idxCallHaplotypes;
                idxCallHaplotypes.push_back(hapIdx);


                if (idxSimilar.size()!=0)
                {
                    // chech which ones of these would be explained away by adding hapIdx
                    std::vector<double> testAddLL(addLL.size(),0.0);
                    for (int h = 0; h < numHaps; h++)
                        computeAddLLSingleRead(1, h, hrLik, testAddLL, numReads, numHaps, added);

                    double logZ = -HUGE_VAL;
                    for(size_t x = 0; x < idxSimilar.size();++x)
                    {
                        int h = idxSimilar[x];
                        // only call haplotype h in this iteration if testAddLL[h] is small,
                        // which means that h and hapIdx explain more or less the same set of reads
                        if (DINDEL_DEBUG_3)
                            std::cout << " haplotype " << h << " : " << testAddLL[h] << "\n";
                        if (testAddLL[h] <= 10.0)
                        {
                            idxCallHaplotypes.push_back(h);
                            if (DINDEL_DEBUG_3)
                                std::cout << "Also adding haplotype " << h << " with likelihood " << testAddLL[h] << "\n";
                        }
                    }



                    for(size_t x = 0; x < idxCallHaplotypes.size();++x)
                        logZ = addLogs(logZ, this->realignParameters.logPriorAddHaplotype + addLL[idxCallHaplotypes[x]]);

                    for(size_t x = 0; x < idxCallHaplotypes.size();++x)
                    {
                        int _idx = idxCallHaplotypes[x];
                        double t = 1.0-exp(this->realignParameters.logPriorAddHaplotype + addLL[_idx] - logZ);
                        if (t<1e-100) t = 1e-100;
                        double renormQual = -log(t);
                        if (DINDEL_DEBUG_3)
                            std::cout << " original HapQual: " << haplotypeQual[_idx] << " renormQual: " << renormQual;
                        if (renormQual<haplotypeQual[_idx])
                            haplotypeQual[_idx] = renormQual;
                        if (DINDEL_DEBUG_3)
                            std::cout << " new qual: " << haplotypeQual[_idx] << "\n";


                    }
                }

                for(size_t x = 0; x < idxCallHaplotypes.size();++x)
                {
                    int _idx = idxCallHaplotypes[x];
                    added[_idx] = 1;
                    addedInIteration[_idx] = iteration; // store in which iteration the haplotype was added.
                    addCalledHaplotypeSingleRead(_idx, result, haplotypes, addReads, minLogLikAlignToAlt, numReads);
                    numAdded++;

                    if (DINDEL_DEBUG_3)
                    {
                        std::cout << "Haplotype[" << _idx << "] ONLY reads\n";
                        showHaplotypeOnlyReadsSingleRead(_idx, hrLik, numReads, numHaps);
                    }

                }

           }
       }


       if(done)
       {
           /*
           if (pPreviousResult != NULL)
           {
              if(this->realignParameters.graphDiffStyle == 1)
              {
                  // set values for first haplotype
                  int h = result.addOrder[0];
                  setAddReads(1, h, hrLik, addReads, addLL, numReadPairs, numHaps, added);
                  double qual = this->realignParameters.logPriorAddHaplotype + addLL[h];
                  haplotypeQual[h] = qual;
                  addCalledHaplotype(h, result, haplotypes, addReads, minLogLikAlignToAlt, numReadPairs);
              }
           }
           */
           stopLoop = true;
       }


       iteration++;
   }

   if (DINDEL_DEBUG_3)
       std::cout << "Number of haplotypes added: " << numAdded << "\n";

   // estimate haplotype frequencies for called haplotypes using EM.

   doEM(hrLik, added, haplotypeFrequencies);

   if (!QUIET) std::cout << "DindelRealignWindow::estimateHaplotypeFrequencies Added " << numAdded << " out of " << haplotypes.size() << " haplotypes." << std::endl;

   if (realignParameters.showCallReads && print)
   {
       for (int h=0;h<numHaps;h++)
       {
               std::cout << "DindelRealignWindow::estimateHaplotypeFrequencies final penalty haplotype[" <<h << "]: " << haplotypeQual[h] << std::endl;
       }

   }

   // compute the expected sum of haplotype frequencies for each mapping location
   const std::vector< DindelReferenceMapping > & refMappings = this->m_dindelWindow.getReferenceMappings();
   result.weightedRefMapFrequencies = std::vector<double>(refMappings.size(), 0.0);

   for (int hapIdx = 0; hapIdx < numHaps; ++hapIdx)
   {
       double hq = haplotypeQual[hapIdx];
       if (hq<0.)
           hq = 0.0;
       double probHap = 1.0-exp(-hq);
       double hapFreq = haplotypeFrequencies[hapIdx];

       for (int refIdx = 0; refIdx < haplotypes[hapIdx].getNumReferenceMappings(); ++refIdx)
       {
           double mapProb = exp(haplotypes[hapIdx].getLogMappingProbability(refIdx));
           result.weightedRefMapFrequencies[refIdx] += hapFreq * probHap * mapProb;
       }
   }

   if (DINDEL_DEBUG_3)
   {
       for (size_t x = 0; x < result.weightedRefMapFrequencies.size();x++) {
           std::cout << " refMapping #" << x << " : weighted hapFreq " << result.weightedRefMapFrequencies[x] << "\n";
       }

   }

   // Add haplotype properties for all variants
   for (int hapIdx = 0; hapIdx < numHaps; ++hapIdx)
   {
       double hapQual = haplotypeQual[hapIdx]/.23026;
       double hapFreq = haplotypeFrequencies[hapIdx];

        for (int refIdx = 0; refIdx < haplotypes[hapIdx].getNumReferenceMappings(); ++refIdx)
        {
            const std::vector<DindelVariant> & vars = haplotypes[hapIdx].getVariants(refIdx);

            for (size_t x = 0; x < vars.size(); x++)
            {
                DindelRealignWindowResult::VarToInference::iterator vit = result.variantInference.find(vars[x]);
                if (vit == result.variantInference.end())
                {
                    // should set quality and freq to zero.
                    DindelRealignWindowResult::Inference varInf;
                    varInf.haplotypeProperties.push_back(DindelRealignWindowResult::HaplotypeProperties(haplotypes[hapIdx].getLogMappingProbability(refIdx), hapQual, (result.weightedRefMapFrequencies[refIdx]<1e-6)?0:hapFreq/result.weightedRefMapFrequencies[refIdx], addedInIteration[hapIdx]));
                    varInf.numCalledHaplotypes = numAdded;
                    result.variantInference[vars[x]] = varInf;
                }
                else
                {
                    DindelRealignWindowResult::Inference & varInf = vit->second;
                    varInf.haplotypeProperties.push_back(DindelRealignWindowResult::HaplotypeProperties(haplotypes[hapIdx].getLogMappingProbability(refIdx), hapQual, (result.weightedRefMapFrequencies[refIdx]<1e-6)?0:hapFreq/result.weightedRefMapFrequencies[refIdx], addedInIteration[hapIdx]));
                    varInf.numCalledHaplotypes = numAdded;
                }

            }
        }

   }

  

  // assign result


   result.haplotypeFrequencies = haplotypeFrequencies;
   size_t totReads = 0; for (size_t x=0;x<addReads.size();x++) totReads += addReads[x].size();
   result.numHapSpecificReads = (int) totReads;
   result.numCalledHaplotypes = (int) numAdded; // -1 makes sure we don't count the reference haplotype. CHANGE when we are doing proper genotyping, because then we align each read to each haplotype

   result.outputID = this->m_outputID;

   return result;

}

DindelRealignWindowResult DindelRealignWindow::estimateHaplotypeFrequenciesModelSelectionSingleReadsMultiSample(double minLogLikAlignToRef,
                                                                                                   double minLogLikAlignToAlt,
                                                                                                   bool capUsingMappingQuality,
                                                                                                   const DindelRealignWindowResult *pPreviousResult,
                                                                                                   bool print = false)
{
   (void) pPreviousResult;
   if (DINDEL_DEBUG_3)
       std::cout << "DindelRealignWindow::estimateHaplotypeFrequenciesModelSelectionSingleReadsMultiSample START" << std::endl;

   const std::vector<DindelMultiHaplotype> & haplotypes = m_dindelWindow.getHaplotypes();
   if (hapReadAlignments.size() != haplotypes.size()) throw std::string("DindelRealignWindow::estimateHaplotypeFrequencies incomplete haplotype alignments.");

   // init result
   DindelRealignWindowResult result(*this);

// construct read X haplotype vector
   int numReads = int(m_pDindelReads->size());
   int numHaps = (int) haplotypes.size();

   std::vector< std::vector<double> > hrLik(numReads, std::vector<double>(numHaps, realignParameters.minLogLikAlignToAlt));

   for (int r=0;r<numReads;r++)
   {
       // Currently we do not take into account the insert size distribution.
       // For tandem repeats and long insertions and deletions this would be desirable.
       for (int h=0;h<numHaps;h++)
       {
           double readLogLik = hapReadAlignments[h][r].logLik;

           if (haplotypes[h].isReference())
               hrLik[r][h]=addLogs(readLogLik, minLogLikAlignToRef);
           else
               hrLik[r][h]=addLogs(readLogLik, minLogLikAlignToAlt);

           hrLik[r][h] = readLogLik;

           if (capUsingMappingQuality) 
               hrLik[r][h] = addLogs(hrLik[r][h], -(*m_pDindelReads)[r].getMappingQual()*.23026-log(double(DINDEL_HMM_BANDWIDTH)));

       }
       // std::cout << " hrLik[" << r << "]: " << hrLik[r][0] << " " << hrLik[r][1] << " 1>0 " << int(hrLik[r][0]<hrLik[r][1]) << "\n";
   }

   std::vector<double> llHapPairs;
   HashMap<std::string, std::list<int> > sampleToReads;

   computeHaplotypePairLikelihoods(llHapPairs, sampleToReads, hrLik);

   int numSamples = int(sampleToReads.size());

   // now add haplotypes until we have obtained a maximal set of hapltoypes

   std::vector<int> added(haplotypes.size(),0);
   std::vector<int> addedInIteration(haplotypes.size(),-1);


   size_t numAdded = 0;
   bool done = false;

   std::vector<double> newLL(numHaps, 0.0); // new loglikelihood after adding haplotype h
   std::vector<double> addLL(numHaps, 0.0); // loglikelihood difference after adding haplotype h
   std::vector<double> haplotypeQual(numHaps, 0.0); // loglikelihood gained by adding haplotypes

   std::vector<double> haplotypeFrequencies(numHaps, 1.0);
   
   
   std::vector< HashMap<int, double> > addReads(numHaps); // read pairs supporting a haplotype


   std::vector< std::vector<double> > zindNew(numHaps, std::vector<double>(numSamples, 0.0));

  
   std::vector<double> zindPrevious(numSamples, -100.0);


   // initial haplotype frequencies are uniform
   std::vector<double> hapFreqPrevious(numHaps, 1.0);


   int iteration = 0;
   bool stopLoop = false;

   double logLikPrevious = 0.0; //


   while (!stopLoop)
   {
       if(DINDEL_DEBUG_3)
           std::cout << "Starting iteration " << iteration << "\n";
       if (numAdded == haplotypes.size()) break; // always keep at least one haplotype

       
       // only itialize kept haplotypes.
       for (int h=0;h<numHaps;h++) if (!added[h])
       {
          addLL[h]=0.0;
          addReads[h].clear();
       }

       int addType = -1;
       std::vector<double> hapFreqNew(numHaps, 0.0);
       
       for (int h = 0; h < numHaps; h++) if (!added[h])
           setAddReadsSingleReadMultiSample(addType, h, hrLik, addReads, newLL, numReads, numHaps, zindNew[h], hapFreqNew, added, hapFreqPrevious, zindPrevious, llHapPairs, sampleToReads);

       if (iteration == 0)
       {
           // make sure the addLL is positive, as haplotypes are only called when addLL > 5
           double min_lik = *(std::min_element(newLL.begin(), newLL.end()));
           for (int h = 0; h < numHaps; h++)
               addLL[h] = newLL[h]-min_lik;
           min_lik = *(std::min_element(addLL.begin(), addLL.end()));
           if (min_lik < 10.0)
               for (int h = 0; h < numHaps; h++)
                   addLL[h] += 100.0;
       }
       else
       {
           for (int h = 0; h < numHaps; h++)
               addLL[h] = newLL[h] - logLikPrevious;
       }


       std::multiset<HL> orderedAddCandidates; // ordered by number of new variants added

       for (int h=0;h<numHaps;h++) if (!added[h])
       {
           // if we would add haplotype h, what would be the increase in the likelihood?
         // FIXME very simple prior..
           double logPrior = this->realignParameters.logPriorAddHaplotype;

           double qual = logPrior + addLL[h];
           if (DINDEL_DEBUG_3)
               std::cout << "estimateHaplotypeFrequenciesModelSelectionMatePairs: estimate haplotype " << h << " logPrior  " << logPrior << " haplotype qual: " << qual << std::endl;

           haplotypeQual[h] = qual;

           orderedAddCandidates.insert(HL(qual, h));
       }


       done = true;
       std::multiset<HL>::const_reverse_iterator iter = orderedAddCandidates.rbegin();
       if (iter != orderedAddCandidates.rend())
       {

           if (iter->ll>=5.0)
           {
                // this haplotype will be added
                int hapIdx = iter->idx;
                double ll = iter->ll ; // use addLL because haplotypeQual was only set for the called haplotypes
                done = false;
                added[hapIdx]=1;

                // see which other haplotypes yield a similar improvement to the likelihood
                std::vector<int> idxSimilar;
                iter++;
                while (iter != orderedAddCandidates.rend())
                {
                    if (DINDEL_DEBUG_3)
                        std::cout << "iter->idx: " << iter->idx << " ll: " << iter->ll << "\n";
                    if (ll - iter->ll <= 20.0 && iter->ll>5.0)
                    {
                        idxSimilar.push_back(iter->idx);
                        if (DINDEL_DEBUG_3)
                            std::cout << "Haplotype " << iter->idx << " has a similar likelihood.\n";
                    }
                    iter++;
                }


                // also call these haplotypes and renormalize scores
                std::vector<int> idxCallHaplotypes;
                idxCallHaplotypes.push_back(hapIdx);


                if (idxSimilar.size()!=0)
                {
                    // chech which ones of these would be explained away by adding hapIdx
                    std::vector<double> testNewLL(addLL.size(),0.0);
                    std::vector< HashMap<int, double> > testAddReads(numHaps); // read pairs supporting a haplotype
                    std::vector<double> testZindNew(numSamples, 0.0);
                    std::vector<double> testHapFreqNew(numHaps, 0.0);
                    
                    for (int h = 0; h < numHaps; h++)
                        if (!added[h])
                            setAddReadsSingleReadMultiSample(-1, h, hrLik, testAddReads, testNewLL, numReads, numHaps, testZindNew, testHapFreqNew, added, hapFreqPrevious, zindPrevious, llHapPairs, sampleToReads);

                    // idea is that if difference with original newLL is zero for those haplotypes with similar score, the two haplotypes would have similar posteriors
                    for (int h = 0; h < numHaps; h++)
                        testNewLL[h] -= newLL[h];
                    
                    double logZ = -HUGE_VAL;
                    for(size_t x = 0; x < idxSimilar.size();++x)
                    {
                        int h = idxSimilar[x];
                        // only call haplotype h in this iteration if testNewLL[h] is small,
                        // which means that h and hapIdx explain more or less the same set of reads
                        if (DINDEL_DEBUG_3)
                            std::cout << " haplotype " << h << " : " << testNewLL[h] << "\n";
                        if (testNewLL[h] <= 10.0)
                        {
                            idxCallHaplotypes.push_back(h);
                            if (DINDEL_DEBUG_3)
                                std::cout << "Also adding haplotype " << h << " with likelihood " << testNewLL[h] << "\n";
                        }
                    }



                    for(size_t x = 0; x < idxCallHaplotypes.size();++x)
                        logZ = addLogs(logZ, this->realignParameters.logPriorAddHaplotype + addLL[idxCallHaplotypes[x]]);

                    for(size_t x = 0; x < idxCallHaplotypes.size();++x)
                    {
                        int _idx = idxCallHaplotypes[x];
                        double t = 1.0-exp(this->realignParameters.logPriorAddHaplotype + addLL[_idx] - logZ);
                        if (t<1e-100) t = 1e-100;
                        double renormQual = -log(t);
                        if (DINDEL_DEBUG_3)
                            std::cout << " original HapQual: " << haplotypeQual[_idx] << " renormQual: " << renormQual;
                        if (renormQual<haplotypeQual[_idx])
                            haplotypeQual[_idx] = renormQual;
                        if (DINDEL_DEBUG_3)
                            std::cout << " new qual: " << haplotypeQual[_idx] << "\n";


                    }
                }

                logLikPrevious = -HUGE_VAL;
                for(size_t x = 0; x < idxCallHaplotypes.size();++x)
                {
                    int _idx = idxCallHaplotypes[x];
                    added[_idx] = 1;
                    addedInIteration[_idx] = iteration; // store in which iteration the haplotype was added.
                    addCalledHaplotypeSingleRead(_idx, result, haplotypes, addReads, minLogLikAlignToAlt, numReads);
                    numAdded++;
                    // update logLikPrevious
                    if (newLL[_idx]>logLikPrevious) {
                        logLikPrevious = newLL[_idx];
                        zindPrevious = zindNew[_idx];
                    }
                    
                    if (DINDEL_DEBUG_3 && 0)
                    {
                        std::cout << "Haplotype[" << _idx << "] ONLY reads\n";
                        showHaplotypeOnlyReadsSingleRead(_idx, hrLik, numReads, numHaps);
                    }

                }

           }
       }


       if(done)
           stopLoop = true;

       iteration++;
   }

    if (DINDEL_DEBUG_3)
        std::cout << "Number of haplotypes added: " << numAdded << "\n";

    // estimate haplotype frequencies for called haplotypes using EM. 
    doEMMultiSample(numSamples, 25, llHapPairs, added, haplotypeFrequencies, haplotypeFrequencies);

    if (realignParameters.genotyping) {
        this->addDiploidGenotypes(result, added, hrLik);
        result.outputGenotypes = true;
    } else {
        result.outputGenotypes = false;
    }

    if (!QUIET) std::cout << "DindelRealignWindow::estimateHaplotypeFrequencies Added " << numAdded << " out of " << haplotypes.size() << " haplotypes." << std::endl;

    if (realignParameters.showCallReads && print)
    {
        for (int h=0;h<numHaps;h++)
               std::cout << "DindelRealignWindow::estimateHaplotypeFrequencies final penalty haplotype[" <<h << "]: " << haplotypeQual[h] << std::endl;
    }

    // compute the expected sum of haplotype frequencies for each mapping location
    const std::vector< DindelReferenceMapping > & refMappings = this->m_dindelWindow.getReferenceMappings();
    result.weightedRefMapFrequencies = std::vector<double>(refMappings.size(), 0.0);

    for (int hapIdx = 0; hapIdx < numHaps; ++hapIdx)
    {
        double hq = haplotypeQual[hapIdx];
        if (hq<0.)
            hq = 0.0;
        double probHap = 1.0-exp(-hq);
        double hapFreq = haplotypeFrequencies[hapIdx];

        for (int refIdx = 0; refIdx < haplotypes[hapIdx].getNumReferenceMappings(); ++refIdx)
        {
            double mapProb = exp(haplotypes[hapIdx].getLogMappingProbability(refIdx));
            result.weightedRefMapFrequencies[refIdx] += hapFreq * probHap * mapProb;
        }
    }

    if (DINDEL_DEBUG_3)
    {
        for (size_t x = 0; x < result.weightedRefMapFrequencies.size();x++) {
            std::cout << " refMapping #" << x << " : weighted hapFreq " << result.weightedRefMapFrequencies[x] << "\n";
        }
    }

    // Add haplotype properties for all variants
    for (int hapIdx = 0; hapIdx < numHaps; ++hapIdx)
    {
        double hapQual = haplotypeQual[hapIdx]/.23026;
        double hapFreq = haplotypeFrequencies[hapIdx];

        for (int refIdx = 0; refIdx < haplotypes[hapIdx].getNumReferenceMappings(); ++refIdx)
        {
            const std::vector<DindelVariant> & vars = haplotypes[hapIdx].getVariants(refIdx);

            for (size_t x = 0; x < vars.size(); x++)
            {
                DindelRealignWindowResult::VarToInference::iterator vit = result.variantInference.find(vars[x]);
                if (vit == result.variantInference.end())
                {
                    // should set quality and freq to zero.
                    DindelRealignWindowResult::Inference varInf;
                    varInf.haplotypeProperties.push_back(DindelRealignWindowResult::HaplotypeProperties(haplotypes[hapIdx].getLogMappingProbability(refIdx), hapQual, (result.weightedRefMapFrequencies[refIdx]<1e-6)?0:hapFreq/result.weightedRefMapFrequencies[refIdx], addedInIteration[hapIdx]));
                    varInf.numCalledHaplotypes = numAdded;
                    result.variantInference[vars[x]] = varInf;
                }
                else
                {
                    DindelRealignWindowResult::Inference & varInf = vit->second;
                    varInf.haplotypeProperties.push_back(DindelRealignWindowResult::HaplotypeProperties(haplotypes[hapIdx].getLogMappingProbability(refIdx), hapQual, (result.weightedRefMapFrequencies[refIdx]<1e-6)?0:hapFreq/result.weightedRefMapFrequencies[refIdx], addedInIteration[hapIdx]));
                    varInf.numCalledHaplotypes = numAdded;
                }

            }
        }
    }



    // assign result
    result.haplotypeFrequencies = haplotypeFrequencies;
    size_t totReads = 0; for (size_t x=0;x<addReads.size();x++) totReads += addReads[x].size();
    result.numHapSpecificReads = (int) totReads;
    result.numCalledHaplotypes = (int) numAdded; // -1 makes sure we don't count the reference haplotype. CHANGE when we are doing proper genotyping, because then we align each read to each haplotype
    result.outputID = this->m_outputID;

    return result;

}

std::vector< std::vector<double> > DindelRealignWindow::fillReadHaplotypeLikelihoods() const
{
   const std::vector<DindelMultiHaplotype> & haplotypes = m_dindelWindow.getHaplotypes();
   int numReads = int(m_pDindelReads->size());
   int numHaps = (int) haplotypes.size();

   std::vector< std::vector<double> > hrLik(numReads, std::vector<double>(numHaps, realignParameters.minLogLikAlignToAlt));

   for (int r=0;r<numReads;r++)
   {
       for (int h=0;h<numHaps;h++)
       {
           double readLogLik = hapReadAlignments[h][r].logLik;

           if (haplotypes[h].isReference())
               hrLik[r][h]=addLogs(readLogLik, realignParameters.minLogLikAlignToRef);
           else
               hrLik[r][h]=addLogs(readLogLik, realignParameters.minLogLikAlignToAlt);

           hrLik[r][h] = readLogLik;
           hrLik[r][h] = addLogs(hrLik[r][h], -(*m_pDindelReads)[r].getMappingQual()*.23026-log(double(DINDEL_HMM_BANDWIDTH)));

       }
       // std::cout << " hrLik[" << r << "]: " << hrLik[r][0] << " " << hrLik[r][1] << " 1>0 " << int(hrLik[r][0]<hrLik[r][1]) << "\n";
   }
   return hrLik;
}

WindowFile::WindowFile(const std::string& fileName)
{
    m_isOpen = false;
    m_fileName = fileName;
    m_fileHandle.open(fileName.c_str());
    if (!m_fileHandle.is_open())
    {
        std::string msg = "WindowFile::WindowFile::cannot_open_file_";
        msg.append(fileName);
        throw msg;
    }
    m_isOpen = true;
    m_index = 0;

    m_empty = std::vector<DindelVariant>();
}

std::vector<DindelVariant> WindowFile::getNextWindow()
{
    std::string line;
    std::getline(m_fileHandle, line);
    if (line.empty()) return m_empty;
    
    // update line counter
    m_index++;

    std::istringstream is(line);
  
    std::vector<DindelVariant> variants;

    while (!is.eof()) {
        std::string varstring;
        is >> varstring;
    
        DindelVariant var;
        bool success=DindelVariant::variantFromWindowString(varstring, var);
        if (!success) 
        {
            std::cout << "WindowFile::getNextWindow Could not parse variant " + varstring + " in line " << m_index+1 << " of window file " << m_fileName << std::endl;
        }
        else
        {
            variants.push_back(var);
        }
    }

    return variants;
}

WindowFile::~WindowFile()
{
    if (m_isOpen)
    {
        m_fileHandle.close();
    }

}
/*
 *
 * REALIGNALGORITHM
 *
 */

DindelRealignParameters::DindelRealignParameters()
{
     initFromString(getDefaultParameters());
}

DindelRealignParameters::DindelRealignParameters(const std::string & paramString)
{
    initFromString(getDefaultParameters());
    // override with user parameters
    initFromString(paramString);
}

DindelRealignParameters::DindelRealignParameters(const char * paramString)
{
    initFromString(getDefaultParameters());
    // override with user parameters
    std::string params(paramString);
    initFromString(params);
}

std::string DindelRealignParameters::getDefaultParameters() const
{
     std::string paramString = "genotyping:0,maxNumReads:100000,maxNumReadsWindow:100000,showCallReads:0,minNumHaplotypeOverlaps:0,maxNumCandidatesPerWindow:32,windowReadBuffer:500,minVariantSep:10,haplotypeWidth:60,minCandidateAlleleCount:0,probSNP:0.001,probINDEL:0.0001,probMNP:0.00001";
     paramString += ",maxMappingQuality:80,addSNPMaxSNPs:0,addSNPMaxMismatches:3,addSNPMinMappingQual:30,addSNPMinBaseQual:20,addSNPMinCount:2,minPostProbLastReadBaseForUngapped:0.95";
     paramString += ",singleSampleHetThreshold:20,singleSampleHomThreshold:20,EMtol:0.0001,EMmaxiter:200,doEM:1,realignMatePairs:0,multiSample:0,priorAddHaplotype:0.0001,graphDiffStyle:1";
     return paramString; 
}      


void DindelRealignParameters::checkAndInit()
{
    bool print = false;
    this->minLogLikAlignToAlt=-double(this->maxMappingQuality)*.2302585-log(double(DINDEL_HMM_BANDWIDTH)); // note last term accounts for base prior in HMM.
    this->minLogLikAlignToRef=this->minLogLikAlignToAlt;
    this->addSNPMinLogBaseQual=log(pow(10.0, -this->addSNPMinBaseQual/10.0));

    std::ostringstream os;

    os <<  "probSNP:" << variantPriors.m_probSNP;
    os << ",probINDEL:" << variantPriors.m_probINDEL;
    os << ",probMNP:" << variantPriors.m_probMNP;
    os << ",minVariantSep:" <<  minVariantSep;
    os << ",windowReadBuffer:" <<  windowReadBuffer;
    os << ",haplotypeWidth:" <<haplotypeWidth;
    os << ",minCandidateAlleleCount:" <<minCandidateAlleleCount;
    os << ",genotyping:" << genotyping;
    os << ",showCallReads:" << showCallReads;
    os << ",minNumHaplotypeOverlaps:" << minNumHaplotypeOverlaps;
    os << ",maxMappingQuality:" << maxMappingQuality;
    os << ",maxMappingQuality:" << maxMappingQuality;
    os << ",maxNumReads:" << maxNumReads;
    os << ",addSNPMinMappingQual:" << addSNPMinMappingQual;
    os << ",addSNPMinBaseQual:" << addSNPMinBaseQual;
    os << ",addSNPMinCount:" << addSNPMinCount;
    os << ",addSNPMaxSNPs:" << addSNPMaxSNPs;
    os << ",minPostProbLastReadBaseForUngapped:" << minPostProbLastReadBaseForUngapped;
    os << ",maxNumCandidatesPerWindow:" << maxNumCandidatesPerWindow;
    os << ",maxNumReadsWindow:" << maxNumReadsWindow;
    os << ",singleSampleHetThreshold:" << int(singleSampleHetThreshold);
    os << ",singleSampleHomThreshold:" << int(singleSampleHomThreshold);
    os << ",EMtol:" << EMtol;
    os << ",EMmaxiter:" << int(EMmaxiter);
    os << ",doEM:" << int(doEM);
    os << ",realignMatePairs:" << int(realignMatePairs);
    os << ",multiSample:" << int(multiSample);
    os << ",priorAddHaplotype:" << exp(logPriorAddHaplotype);
    os << ",graphDiffStyle:" << int(graphDiffStyle);
    if (print) 
    {
        std::cout << "Realignment parameters:\n";
        this->m_paramString = os.str();
        for (size_t x=0;x<m_paramString.size();x++)
        {
        char c = m_paramString[x];
        if (c==':') std::cout << "\t";
        else if (c== ',') std::cout << "\n";
        else std::cout << c;
        }
        std::cout << std::endl;
    }

}

void DindelRealignParameters::initFromString(const std::string& paramString)
{
    assert(!paramString.empty());

    // check that it contains no spaces
    for (size_t x=0;x<paramString.size();x++) if (paramString[x]==' ' || paramString[x]=='\t')
    {
        std::cout << "Parameter string cannot contain spaces or tabs." << std::endl;
        exit(1);
    }

    std::vector<std::string> pairs = SplitString(paramString,',');

    for (size_t p=0;p<pairs.size();p++)
    {
        std::vector<std::string> pair = SplitString(pairs[p],':');

        if (pair.size()==1)
        {
            // don't have any yet
        std::cout << "Unrecognized flag" << std::endl;
        exit(1);
            
    }
        else if (pair.size()==2)
        {
            const std::string & k = pair[0];
            const std::string & v = pair[1];
            bool fail=false;
            double minProbVariant = 1e-10;
            double maxProbVariant = 1.0;

            if (k == "probSNP") { if (!from_string<double>(variantPriors.m_probSNP,minProbVariant, maxProbVariant, v, std::dec)) fail = true; }
            else if (k == "probINDEL") { if (!from_string<double>(variantPriors.m_probINDEL,minProbVariant, maxProbVariant, v, std::dec)) fail = true; }
            else if (k == "probMNP") { if (!from_string<double>(variantPriors.m_probMNP,minProbVariant, maxProbVariant, v, std::dec)) fail = true; }
            else if (k == "minVariantSep") { if (!from_string<int>(minVariantSep,1, 100, v, std::dec)) fail = true; }
            else if (k == "windowReadBuffer") { if (!from_string<int>(windowReadBuffer, 10, 10000, v, std::dec)) fail = true; }
            else if (k == "haplotypeWidth") { if (!from_string<int>(haplotypeWidth, 20, 100, v, std::dec)) fail = true; }
            else if (k == "minCandidateAlleleCount") { if (!from_string<int>(minCandidateAlleleCount,0, 1000, v, std::dec)) fail = true; }
            else if (k == "maxMappingQuality") { if (!from_string<int>(maxMappingQuality,40,100,v, std::dec)) fail = true; }
            else if (k == "addSNPMaxMismatches") { if (!from_string<int>(addSNPMaxMismatches,1,10, v, std::dec)) fail = true; }
            else if (k == "addSNPMinMappingQual") { if (!from_string<double>(addSNPMinMappingQual,5.,50., v, std::dec)) fail = true; }
            else if (k == "addSNPMinBaseQual") { if (!from_string<int>(addSNPMinBaseQual,0, 20, v, std::dec)) fail = true; }
            else if (k == "addSNPMinCount") { if (!from_string<size_t>(addSNPMinCount,1, 20, v, std::dec)) fail = true; }
            else if (k == "addSNPMaxSNPs") { if (!from_string<int>(addSNPMaxSNPs,1, 6, v, std::dec)) fail = true; }
            else if (k == "minPostProbLastReadBaseForUngapped") { if (!from_string<double>(minPostProbLastReadBaseForUngapped,0.0, 1.0, v, std::dec)) fail = true; }
            else if (k == "genotyping") { if (!from_string<int>(genotyping,0, 1, v, std::dec)) fail = true; }
            else if (k == "maxNumReads") { if (!from_string<int>(maxNumReads,0, 200000, v, std::dec)) fail = true; }
            else if (k == "showCallReads") { if (!from_string<int>(showCallReads,0, 1, v, std::dec)) fail = true; }
            else if (k == "minNumHaplotypeOverlaps") { if (!from_string<int>(minNumHaplotypeOverlaps,0, 100, v, std::dec)) fail = true; }
            else if (k == "excludeSamplesFile") { if (!v.empty()) { excludeSamplesFile = v; } else { fail = true; } }
            else if (k == "maxNumCandidatesPerWindow") { if (!from_string<int>(maxNumCandidatesPerWindow,0, 64, v, std::dec)) fail = true; }
            else if (k == "maxNumReadsWindow") { if (!from_string<int>(maxNumReadsWindow,0, 100000, v, std::dec)) fail = true; }
            else if (k == "singleSampleHetThreshold") { if (!from_string<double>(singleSampleHetThreshold,0, 1000, v, std::dec)) fail = true; }
            else if (k == "singleSampleHomThreshold") { if (!from_string<double>(singleSampleHomThreshold,0, 1000, v, std::dec)) fail = true; }
            else if (k == "EMtol") { if (!from_string<double>(EMtol,0, 1, v, std::dec)) fail = true; }
            else if (k == "EMmaxiter") { if (!from_string<int>(EMmaxiter,0, 10000, v, std::dec)) fail = true; }
            else if (k == "doEM") { if (!from_string<int>(doEM,0, 1, v, std::dec)) fail = true; }
            else if (k == "realignMatePairs") { if (!from_string<int>(realignMatePairs,0, 1, v, std::dec)) fail = true; }
            else if (k == "multiSample") { if (!from_string<int>(multiSample,0, 1, v, std::dec)) fail = true; }
            else if (k == "priorAddHaplotype")  { double tmp; if (!from_string<double>(tmp,1e-10, 1.0, v, std::dec)) fail = true; else logPriorAddHaplotype=log(tmp); }
            else if (k == "graphDiffStyle")  { if (!from_string<int>(graphDiffStyle,0, 1, v, std::dec)) fail = true; }
            else throw std::string("Unrecognized parameter: "+k);

            if (fail) throw std::string("Cannot determine value for parameter " + k + " from "+v);
      }
        else throw std::string("Incorrectly formatted parameter string: two consecutive colons");
    }

    checkAndInit();
}

/*


 VCF File


 */


std::vector<std::string> VCFFile::VCFEntry::getFilters() const
{
    std::vector<std::string> vFilters;
    Maps::const_iterator it = filters.begin();
    for (;it!=filters.end();it++) vFilters.push_back(it->first);
    return vFilters;
}

std::string VCFFile::VCFEntry::getInfoValue(const std::string& key) const
{
    Maps::const_iterator it = info.find(key);
    if (it == info.end()) return "-"; else return it->second;
}

void VCFFile::VCFEntry::parseInfoString()
{
    std::vector<std::string> infos = SplitString(this->infoString,';');

    for (size_t x=0;x<infos.size();x++)
    {
        std::vector<std::string> kv = SplitString(infos[x],'=');
        if (kv[0].empty()) throw std::string("VCFEntry: cannot parse info tag "+kv[0]);
        if (kv.size()==1) info[kv[0]]="TAG";
        else if (kv.size()==2) info[kv[0]]=kv[1];
        else throw std::string("VCFEntry: cannot parse info string");
    }

}

void VCFFile::VCFEntry::setFilters(const std::string & str)
{
     std::vector<std::string> f = SplitString(str,';');
     for (size_t x=0;x<f.size();x++) filters[f[x]] = "PRESENT";
}


template<class T> bool VCFFile::VCFEntry::fromInfoTag(T & val, const std::string & tag, std::ios_base& (*f)(std::ios_base&))
{
    std::string tagv = getInfoValue(tag);
    if (tagv=="-") return false;
    if (!from_string(val, tagv, f)) std::cout << "Candidate VCF File: could not parse INFO tag " << tag << std::endl;
    return true;
}

void VCFFile::VCFEntry::write(std::ostream & out) const
{
    out << chrom << " " << ref << " " << alt << " " << infoString << std::endl;

}

VCFFile::VCFFile()
{
    m_isOpen = false;
    m_mode = "null";
}

VCFFile::VCFFile(const std::string& fileName, const std::string & mode)
{
    m_isOpen = false;
    m_fileName = fileName;
    assert(mode == "w" || mode == "r");

    // open the file

    if (mode == "w")
    {
        m_outputFileHandle.open(fileName.c_str());
    if (m_outputFileHandle.is_open()) m_isOpen = true;

    }
    else if (mode == "r")
    {
        m_inputFileHandle.open(fileName.c_str());
    if (m_inputFileHandle.is_open()) m_isOpen = true;

    }

    if (!m_isOpen)
    {
        std::string msg = "VCFFile::VCFFile: Cannot open file: '";
        msg.append(fileName);
        msg.append("'");
        if(fileName.empty()) msg.append("EMPTY FILENAME");
        throw msg;
    }
    m_mode = mode;
}

void VCFFile::setSamples(const std::vector<std::string> & samples)
{
    m_samples = samples;
}

void VCFFile::outputHeader(const std::string & refFile, const std::string& /*paramString*/)
{
    assert(m_isOpen && m_mode == "w");
    m_outputFileHandle << "##fileformat=VCFv4.1" << std::endl;
    m_outputFileHandle << "##source=SGA/Dindel" << std::endl;
    m_outputFileHandle << "##reference=" << refFile << std::endl;

    m_outputFileHandle << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=NumReads,Number=1,Type=Integer,Description=\"Number of reads passed to Dindel\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=NumCalledHaps,Number=1,Type=Integer,Description=\"Number of haplotypes called by Dindel\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=ExpNumHapsWithVarMappingHere,Number=1,Type=Float,Description=\"Expected number of haplotypes mapping to this location\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=VarQual,Number=1,Type=Float,Description=\"Variant quality\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=HV,Number=.,Type=String,Description=\"Haplotype property string\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=HMQ,Number=1,Type=Integer,Description=\"Haplotype mapping quality\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=VarDP,Number=1,Type=Integer,Description=\"Number of reads containing the variant\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=NF,Number=1,Type=Integer,Description=\"Number of reads preferentially aligning to variant haplotype on reverse strand\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=NR,Number=1,Type=Integer,Description=\"Number of reads preferentially aligning to variant haplotype on forward strand\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=NF0,Number=1,Type=Integer,Description=\"Number of perfect reads preferentially aligning to variant haplotype on reverse strand\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=NR0,Number=1,Type=Integer,Description=\"Number of perfect reads preferentially aligning to variant haplotype on forward strand\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=HSR,Number=1,Type=Integer,Description=\"Number of haplotype-specific (informative) reads\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=HPLen,Number=1,Type=Integer,Description=\"Homopolymer length\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand bias\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=HPLen,Number=1,Type=Integer,Description=\"Homopolymer length\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=Dust,Number=1,Type=Float,Description=\"Dust score for a 64bp window centered at the variant start\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=HistDist,Number=25,Type=Integer,Description=\"Histogram of variant distance from the end of the reads\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=HistLik,Number=25,Type=Integer,Description=\"Histogram of read likelihoods\">" << std::endl;
    m_outputFileHandle << "##INFO=<ID=HistMapQ,Number=8,Type=Integer,Description=\"Histogram of mapping quality\">" << std::endl;

    m_outputFileHandle << "##FILTER=<ID=Candidate,Description=\"Variant is candidate for testing\">" << std::endl;
    m_outputFileHandle << "##FILTER=<ID=NoCall,Description=\"No call made\">" << std::endl;
    m_outputFileHandle << "##FILTER=<ID=LowQuality,Description=\"Low quality variant\">" << std::endl;

    m_outputFileHandle << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    m_outputFileHandle << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" << std::endl;
    m_outputFileHandle << "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">" << std::endl;
    m_outputFileHandle << "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Log-likelihoods\">" << std::endl;
    m_outputFileHandle << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    if (m_samples.size()>0)
    {
        m_outputFileHandle << "\tFORMAT";
        for (size_t x=0;x<m_samples.size();x++) m_outputFileHandle << "\t" << m_samples[x];
    }

    m_outputFileHandle << std::endl;
}

VCFFile::VCFEntry VCFFile::getNextEntry()
{
    assert(m_mode=="r");


    if (m_inputFileHandle.eof()) return VCFFile::VCFEntry();

    std::string line;

    bool foundEntry = false;
    while (!foundEntry)
    {
        if (m_inputFileHandle.eof()) 
            return VCFFile::VCFEntry();
        std::getline(m_inputFileHandle,line);

        if (line.size()<1) 
            return VCFFile::VCFEntry();
        if (line[0]!='#') 
            foundEntry=true; 
        else 
            line.clear();
    }
    
    std::istringstream is(line);

    VCFEntry e;

    int idx = 0;
    while (!is.eof())
    {
        std::string field;
        is >> field;
        switch(idx)
        {
            case 0: e.chrom=field; break;
            case 1:
                if (!from_string<int>(e.pos, field, std::dec)) throw std::string("VCF: cannot convert position in entry "+line);
                break;
            case 2: e.id = field; break;
            case 3: e.ref = field; break;
            case 4: e.alt = field; break;
            case 5:
                float qual;
                if (!from_string<float>(qual, field, std::dec)) throw std::string("VCF: cannot convert quality in entry "+line);
                e.qual = int(qual);
                break;
            case 6: e.setFilters(field); break;
            case 7: e.infoString = field; e.parseInfoString(); break;
        }
        idx++;
    }
    e.m_isEmpty = false;
    
    if (idx<7) 
        throw std::string("VCF: line did not have all required fields " + line);

    return e;
}

VCFFile::~VCFFile()
{
   if (m_mode == "w")
    m_outputFileHandle.close();
   else if (m_mode == "r")
    m_inputFileHandle.close();
}





DindelVariant addInsertion(const DindelHaplotype & refHap, const std::string & chrom, int refPos, const std::string & seq)
{
    // insertion is added BEFORE refPos (see extractIndelsFromBamLeftJustified why)

    // determine left and right boundaries

    const std::string & refSeq = refHap.getSequence();
    int refseqlen = refSeq.length();

    int relPos = refPos-refHap.getRefStart();
    if (relPos<=0) throw std::string("relpos error");

    int len = (int) seq.length();

    int left=relPos-REPOSITION_INDEL_WINDOW;
    //while (left>=0 && refSeq.substr(left,len) == seq) left-=len;
    if (left<0) left=0;

    int right=relPos+REPOSITION_INDEL_WINDOW;
    if (right>refseqlen) right=refseqlen;
    //while (right+len<=refseqlen && refSeq.substr(right,len)==seq) right+=len;

    std::string eir = refSeq.substr(left, right-left+1);
    std::string alt = eir;
    alt.insert(relPos-left, seq);

    if (DINDEL_DEBUG)
    {
        std::cout << "\n ***** INS \n ";
        std::cout << "refPos: " << refPos << " relPos: " << relPos << " left: " << left << " right: " << right << " len: " << len << " seq: " << seq << "\n";
        std::cout << refSeq.substr(relPos-25,25) << "|" << refSeq.substr(relPos,25) << "\n";
        std::cout << "eir: " << eir << "\n";
        std::cout << "alt: " << alt << "\n";
    }


    // shift to the left as far as possible

    int pos=relPos-left;
    std::string newSeq=seq;
    for (int p=relPos-left-1;p>=0;p--)
    {
        std::string _alt = eir;
        std::string nseq =alt.substr(p,len);
       _alt.insert(p, nseq);
       if (_alt == alt) {
           pos=p;
           newSeq=nseq;
       }

    }
    int newPos = pos+left-1; // -1 is for conversion to VCF4 convention

    std::string ref; ref += refSeq[newPos];
    alt=ref+newSeq;

    if (DINDEL_DEBUG)
        std::cout << " ref: " << ref << " alt: " << alt << " pos: " << newPos << std::endl;
    return DindelVariant(chrom, ref, alt, newPos+refHap.getRefStart());
    
}

DindelVariant addDeletion(const DindelHaplotype & refHap, const std::string & chrom, int refPos, int len)
{
    const std::string & refSeq = refHap.getSequence();
    int refseqlen = refSeq.length();
    int relPos = refPos-refHap.getRefStart();
    assert(relPos>0);

    assert(relPos+len<=int(refSeq.length()));
    std::string seq = refSeq.substr(relPos, len);

    int left=relPos-REPOSITION_INDEL_WINDOW;
    //while (left>=0 && refSeq.substr(left,len) == seq) left-=len;
    if (left<0) left=0;

    int right=relPos+REPOSITION_INDEL_WINDOW;
    if (right>refseqlen) right=refseqlen;
 
    std::string eir = refSeq.substr(left, right-left+1);
    std::string alt = eir;
    alt.erase(relPos-left, len);

    if (DINDEL_DEBUG)
    {
        std::cout << "\n ***** DEL \n ";
        std::cout << "refPos: " << refPos << " relPos: " << relPos << " left: " << left << " right: " << right << " len: " << len << "\n";
        std::cout << refSeq.substr(relPos-25,25) << "|" << refSeq.substr(relPos,25) << "\n";
        std::cout << "eir: " << eir << "\n";
        std::cout << "alt: " << alt << "\n";
    }

    // shift to the left as far as possible

    int pos=relPos-left;
    std::string newSeq=eir.substr(pos, len);
    for (int p=relPos-left-1;p>=0;p--)
    {
        std::string _alt = eir;
        std::string nseq =eir.substr(p,len);
       _alt.erase(p, nseq.length());
       if (_alt == alt) {
           pos=p;
           newSeq=nseq;
       }

    }
    int newPos = pos+left-1;

    std::string ref; ref += refSeq[newPos]; 
    alt=ref;
    ref+=newSeq;
    assert(newSeq[0]==refSeq[newPos+1]);
    if (DINDEL_DEBUG)
    {
        std::cout << "ref: " << ref << " alt: " << alt << " newPos: " << newPos << std::endl;
        std::cout <<"\n";
    }
    return DindelVariant(chrom, ref, alt, newPos+refHap.getRefStart());
}
