#pragma once

#include <unordered_map>
#include <string>
#include <set>

#include "SeqLib/BamRecord.h"

#include "BreakPoint.h"

using SeqLib::BamRecordPtr;
using SeqLib::Cigar;
using SeqLib::CigarMap;
using SeqLib::GenomicRegion;
using SeqLib::BamHeader;

// forward declare
class AlignedContig;

  /*! This class contains a single alignment fragment from a contig to
   * the reference. For a multi-part mapping of a contig to the reference,
   * an object of this class represents just a single fragment from that alignment.
   */
class AlignmentFragment {
  
  friend class BreakPoint;
  friend class AlignedContig;
  
public:
  
  AlignmentFragment() = delete; // prevent default construction
  
  /*! Construct an AlignmentFragment from a BWA alignment
   * @param flip If the contig sequence was flipped (rev of BAM record), need to track this. This flipping occurs in AlignedContig::AlignedContig
   */
  // talign is an alignment of the contig to reference
  AlignmentFragment(BamRecordPtr &talign,
		    bool flip,
		    const GenomicRegion& local_region,
		    const SvabaSharedConfig* sc_);
  
  // sort AlignmentFragment objects by start position
  bool operator < (const AlignmentFragment& str) const; 
  
  void indelCigarMatches(const std::unordered_map<std::string, CigarMap>& cmap);
  
  // print the AlignmentFragment
  std::string printToAlignmentsFile() const;
  
  void fillRearrangementBreakEnd(bool left, BreakEnd& b);
    
  // check whether the alignment fragement overlaps with the given windows
  bool checkLocal() const;
  
  //const std::vector<BreakPoint>& getIndelBreaks() const { return m_indel_breaks; }
  
  // write the alignment record to a BAM file
  //void writeToBAM(SeqLib::BamWriter& bw) const;
  
  void SetIndels();
  
  BamRecordPtr m_align; // BWA alignment of contig to reference
  
  std::vector<AlignmentFragment> secondaries;    
  
  BreakPointPtrVector m_indel_breaks; // indel variants on this alignment
  
  Cigar m_cigar; //cigar oriented to assembled orientation
  std::string m_seq; // seequence, orientated to assembled orientation
  
  int break1 = -1; // 0-based breakpoint 1 on contig 
  int break2 = -1; // 0-based breakpoint 2 on contig 
  int gbreak1 = -1; // 0-based breakpoint 1 on reference chr
  int gbreak2 = -1; // 0-based breakpoint 1 on reference chr 
  
  const SvabaSharedConfig* sc; //pointer to allow sort later
  
  GenomicRegion region_;
  
  int start = -1; // where on the contig does this alignment fragment start
  
  int num_align = -1;

  bool flipped = false;
};

typedef std::vector<AlignmentFragment> AlignmentFragmentVector;
