#ifndef SNOWMAN_ALIGNED_CONTIG_H__
#define SNOWMAN_ALIGNED_CONTIG_H__

#include "SeqLib/BamRecord.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/GenomicRegion.h"
#include <unordered_map>
#include <string>
#include <set>
#include "BreakPoint.h"

#define MAX_CONTIG_SIZE 5000000

// forward declare
class AlignedContig;

  /*! This class contains a single alignment fragment from a contig to
   * the reference. For a multi-part mapping of a contig to the reference,
   * an object of this class represents just a single fragment from that alignment.
   */
  class AlignmentFragment {

  public:
    
    friend class AlignedContig;
    
    /*! Construct an AlignmentFragment from a BWA alignment
     * @param const reference to an aligned sequencing read
     * @param flip If the contig sequence was flipped (rev of BAM record), need to track this. This flipping occurs in AlignedContig::AlignedContig
     */
    AlignmentFragment(const SeqLib::BamRecord &talign, bool flip);
    
    // sort AlignmentFragment objects by start position
    bool operator < (const AlignmentFragment& str) const { return (start < str.start); }

    void indelCigarMatches(const std::unordered_map<std::string, SeqLib::CigarMap>& cmap);
    
    // print the AlignmentFragment
    friend std::ostream& operator<<(std::ostream &out, const AlignmentFragment& c); 

    BreakEnd makeBreakEnd(bool left);
    
    /*! @function
     * @abstract Parse an alignment frag for a breakpoint
     * @param reference to a BreakPoint to be created
     * @return boolean informing whether there was a remaining indel break
     */
    bool parseIndelBreak(BreakPoint &bp);
    
    // check whether the alignment fragement overlaps with the given windows
    bool checkLocal(const SeqLib::GenomicRegion& window);
    
    const std::vector<BreakPoint>& getIndelBreaks() const { return m_indel_breaks; }
    
    // write the alignment record to a BAM file
    void writeToBAM(SeqLib::BamWriter& bw) const;
    
    std::vector<AlignmentFragment> secondaries;

    void SetIndels(const AlignedContig * c);

    private:
    
    SeqLib::BamRecord m_align; // BWA alignment to reference

    int sub_n = 0; // number of sub optimal alignments
    
    std::vector<BreakPoint> m_indel_breaks; // indel variants on this alignment
    
    SeqLib::Cigar m_cigar; //cigar oriented to assembled orientation 
    
    size_t idx = 0; // index of the cigar where the last indel was taken from 
    
    int break1 = -1; // 0-based breakpoint 1 on contig 
    int break2 = -1; // 0-based breakpoint 2 on contig 
    int gbreak1 = -1; // 0-based breakpoint 1 on reference chr
    int gbreak2 = -1; // 0-based breakpoint 1 on reference chr 
    
    int start; // the start position of this alignment on the reference. 
    
    bool local = false; // boolean to note whether this fragment aligns to same location is was assembled from 
    
    int di_count = 0; // number of indels

    int num_align = 0;
  };

typedef std::vector<AlignmentFragment> AlignmentFragmentVector;


#endif
