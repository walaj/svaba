#ifndef ALIGNED_CONTIG_H
#define ALIGNED_CONTIG_H

#include <algorithm>

#include "SnowTools/GenomicRegion.h"
#include "SnowTools/GenomicRegionCollection.h"
#include "SnowTools/SnowUtils.h"

#include "BamToolsUtils.h"
#include "BreakPoint.h"

using SnowTools::GRC;

typedef std::vector<BamTools::CigarOp> CigarOpVec;
typedef std::unordered_map<std::string, size_t> CigarMap;

class AlignedContig;
typedef std::unordered_map<std::string, AlignedContig> ContigMap;
typedef std::vector<AlignedContig> AlignedContigVec;

/*! This class contains a single alignment fragment from a contig to
 * the reference. For a multi-part mapping of a contig to the reference,
 * an object of this class represents just a single fragment from that alignment.
 */
struct AlignmentFragment {

  friend AlignedContig;

  /*! Construct an AlignmentFragment from a BWA alignment
   * @param const reference to a BamAlignment
   * @param const reference to a GenomicRegion window where this contig was assembled from
   */
  AlignmentFragment(const BamTools::BamAlignment &talign, const SnowTools::GenomicRegion &window, const CigarMap &nmap, const CigarMap &tmap);
  
  //! sort AlignmentFragment objects by start position
  bool operator < (const AlignmentFragment& str) const { return (start < str.start); }

  //! print the AlignmentFragment
  friend std::ostream& operator<<(std::ostream &out, const AlignmentFragment& c); 

  /*! @function
   * @abstract Parse an alignment frag for a breakpoint
   * @param reference to a BreakPoint to be created
   * @return boolean informing whether there was a remaining indel break
   */
  bool parseIndelBreak(BreakPoint &bp);

  /*! @function check if breakpoints match any cigar strings direct from the reads
   * This is important in case the assembly missed a read (esp normal) that indicates
   * that there is an indel. Basically this function says that if assembly calls a breakpoint
   * and it agrees with a read alignment breakpoint, combine the info.
   * @param const CigarMap reference containing hash with key=chr_breakpos_indeltype, val=tumor count
   * @param const CigarMap reference containing hash with key=chr_breakpos_indeltype, val=normal count
   */
  void indelCigarMatches(const CigarMap &nmap, const CigarMap &tmap);  

  BPVec indel_breaks; /**< indel variants on this alignment */

  CigarOpVec cigar; /**< cigar oriented to assembled orientation */

  private:

  BamTools::BamAlignment align; /**< BWA alignment to reference */

  size_t idx = 0; // index of the cigar where the last indel was taken from 

  int break1 = -1; // 0-based breakpoint 1 on contig 
  int break2 = -1; /**< 0-based breakpoint 2 on contig */
  int gbreak1 = -1; /**< 0-based breakpoint 1 on reference chr */
  int gbreak2 = -1; /**< 0-based breakpoint 1 on reference chr */
  
  size_t start; /**< the start position of this alignment on the reference. */

  bool local = false; /**< boolean to note whether this fragment aligns to same location is was assembled from */

  std::string m_name; // name of the entire contig
  
  std::string m_seq; // sequence of the entire contig

};

// define a way to order the contigs by start
/*struct AlignmentOrdering {
  inline bool operator() (const AlignmentFragment& struct1, const AlignmentFragment& struct2) {
    return (struct1.start < struct2.start);
  }
  };*/

//! vector of AlignmentFragment objects
typedef std::vector<AlignmentFragment> AlignmentFragmentVector;

/*! Contains the mapping of an aligned contig to the reference genome,
 * along with pointer to all of the reads aligned to this contig, and a 
 * store of all of the breakpoints associated with this contig
 */
class AlignedContig {

  friend AlignmentFragment;

 public:  

  /*! Constructor which parses an alignment record from BWA (a potentially multi-line SAM record)
   * @param const reference to a string representing a SAM alignment (contains newlines if multi-part alignment)
   * @param const pointer to a BamReader, which is used to convert chr ids to strings (e.g. X, Y)
   * @param const reference to a SnowTools::GenomicRegion window specifying where in the reference this contig was assembled from.
   */
  AlignedContig(const std::string &sam, const BamTools::BamReader * reader, const SnowTools::GenomicRegion &twindow, 
		const CigarMap &nmap, const CigarMap &tmap);

  std::string samrecord; /**< the original SAM record */

  SnowTools::GenomicRegion window; /**< reference window from where this contig was assembled */

  /*! @function Determine if this contig has identical breaks and is better than another.
   * @param const reference to another AlignedContig
   * @return bool bool returning true iff this contig has identical info has better MAPQ, or equal MAPQ but longer */
  bool isWorse(const AlignedContig &ac) const;
  
  /*! @function
    @abstract  Get whether the query is on the reverse strand
    @param  b  pointer to an alignment
    @return    boolean true if query is on the reverse strand
  */
  void addAlignment(const BamTools::BamAlignment &align, const SnowTools::GenomicRegion &window, 
		    const CigarMap &nmap, const CigarMap &tmap);

  //! add a discordant cluster that maps to same regions as this contig
  void addDiscordantCluster(DiscordantCluster dc) { m_dc.push_back(dc); } 

  /*! @function loop through the vector of DiscordantCluster objects
   * associated with this contig and print
   */
  std::string printDiscordantClusters() const;

  //! return the name of the contig
  std::string getContigName() const { assert(m_align.size()); return m_align[0].align.Name; }

  /*! @function get the maximum mapping quality from all alignments
   * @return int max mapq
   */
  int getMaxMapq() const { 
    int m = -1;
    for (auto& i : m_align)
      if (i.align.MapQuality > m)
	m = i.align.MapQuality;
    return m;
  }

  /*! @function get the minimum mapping quality from all alignments
   * @return int min mapq
   */
  int getMinMapq() const { 
    int m = 1000;
    for (auto& i : m_align)
      if (i.align.MapQuality < m)
	m = i.align.MapQuality;
    return m;
  }

  /*! @function loop through all of the breakpoints and
   * calculate the split read support for each. Requires 
   * alignedReads to have been run first (will error if not run).
   */
  void splitCoverage();
  
  /*! Checks if any of the indel breaks are in a blacklist. If so, mark the
   * breakpoints of the indels for skipping. That is, hasMinmal() will return false;
   * @param grv The blaclist regions
   */
  void blacklist(GRC& grv);

  /*! @function if this is a Discovar contig, extract
   * the tumor and normal read support
   * @param reference to int to fill for normal support
   * @param reference to int to fill for tumor support
   * @return boolean reporting if this is a discovar name
   */
  bool parseDiscovarName(size_t &tumor, size_t &normal);

  /*! @function dump the contigs to a fasta
   * @param ostream to write to
   */
  void printContigFasta(std::ofstream &os) const;

  /*! @function set the breakpoints on the reference by combining multi-mapped contigs
   */
  void setMultiMapBreakPairs();

  /*! @function align reads to contig and modify their tags to show contig mapping
   * Currently this function will attempt a SmithWaterman alignment for all reads
   * that don't have an exact mapping to the contig.
   * @param bav Vector of read smart pointers to align. Modifies their SW tag
   */
  void alignReadsToContigs(ReadVec &bav);

  //! return the contig sequence as it came off the assembler
  std::string getSequence() const { assert(m_seq.length()); return m_seq; }

  //! detemine if the contig contains a subsequence
  bool hasSubSequence(const std::string& subseq) const { 
    return (m_seq.find(subseq) != std::string::npos);
  }
  
  //! print this contig
  friend std::ostream& operator<<(std::ostream &out, const AlignedContig &ac);

  /*! @function query if this contig contains a potential variant (indel or multi-map)
   * @return true if there is multimapping or an indel
   */
  bool hasVariant() const;

  /*! @function retrieves all of the breakpoints by combining indels with global mutli-map break
   * @return vector of ind
   */
  std::vector<BreakPoint> getAllBreakPoints() const;
  std::vector<const BreakPoint*> getAllBreakPointPointers() const ;

  std::vector<BreakPoint> m_local_breaks; // store all of the multi-map BreakPoints for this contigs 

  BreakPoint m_global_bp;  // store the single spanning BreakPoing for this contig e

  ReadVec m_bamreads; // store smart pointers to all of the reads that align to this contig 

  bool m_skip = false; // flag to specify that we should minimally process and simply dump to contigs_all.sam 

  AlignmentFragmentVector m_align; // store all of the individual alignment fragments 

 private:

  bool m_hasvariant = false; // flag to specify whether this alignment has some potential varaint (eg indel)

  bool m_tried_align_reads = false; // flag to specify whether we tried to align reads.

  std::string m_seq = ""; // sequence of contig as it came off of assembler

  std::vector<DiscordantCluster> m_dc; // collection of all discordant clusters that map to same location as this contig

};

struct PlottedRead {

  int pos;
  std::string seq;
  std::string info;

  bool operator<(const PlottedRead& pr) const {
    return (pos < pr.pos);
  }

};

typedef std::vector<PlottedRead> PlottedReadVector;

struct PlottedReadLine {

  std::vector<PlottedRead*> read_vec;
  int available = 0;
  int contig_len = 0;

  void addRead(PlottedRead *r) {
    read_vec.push_back(r);
    available = r->pos + r->seq.length() + 5;
  }

  bool readFits(PlottedRead &r) {
    return (r.pos >= available);
  }

  friend std::ostream& operator<<(std::ostream& out, const PlottedReadLine &r) {
    int last_loc = 0;
    for (auto& i : r.read_vec) {
      assert(i->pos - last_loc >= 0);
      out << std::string(i->pos - last_loc, ' ') << i->seq;
      last_loc = i->pos + i->seq.length();
    }
    int name_buff = r.contig_len - last_loc;
    assert(name_buff < 10000);
    out << std::string(max(name_buff, 5), ' ');
    for (auto& i : r.read_vec) { // add the data
      out << i->info << ",";
    }
    return out;
  }

};

typedef std::vector<PlottedReadLine> PlottedReadLineVector;



#endif


