#pragma once

#include <algorithm>

#include "SeqLib/BWAAligner.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BamWriter.h"

#include "BreakPoint.h"
#include "DiscordantCluster.h"
#include "AlignmentFragment.h"
#include "SvabaRead.h"

class R2CDatabase;  // forward decl; full def in R2CDatabase.h

/*! Contains the mapping of an aligned contig to the reference genome,
 * along with pointer to all of the reads aligned to this contig, and a 
 * store of all of the breakpoints associated with this contig
 */
class AlignedContig {
  
  friend class AlignmentFragment;
  friend class BreakPoint;
  
 public:  
  
  AlignedContig() = delete; // prevent default construction
  
  // make an AlignedContig from a set of contig alignments
  AlignedContig(BamRecordPtrVector& bav,
		const GenomicRegion& region,
		const SvabaSharedConfig* sc_);
  
  // Return as a genomic region vector
  SeqLib::GenomicRegionVector getAsGenomicRegionVector() const;

  bool checkLocal(const SeqLib::GenomicRegion& window);
  
  // apply repeat filter to each indel break
  //void assessRepeats();
  
  // Loop through the vector of DiscordantCluster objects
  // associated with this contig and print
  std::string printDiscordantClusters(const BamHeader& h) const;
  
  // return the name of the contig
  std::string getContigName() const;

  // Debug accessors for compile-time trace (SvabaDebug.h)
  size_t getFragCount() const { return m_frag_v.size(); }
  bool hasGlobalBP() const { return m_global_bp != nullptr; }
  size_t getLocalBreakCount() const { return m_local_breaks.size(); }
  size_t getIndelBreakCount() const {
    size_t n = 0;
    for (const auto& f : m_frag_v) n += f.m_indel_breaks.size();
    return n;
  }

  // Return the max mapping quality from all alignments
  int getMaxMapq() const;
    
  // Return the min mapping quality from all alignments
  int getMinMapq() const;
  
  // Loop through all of the breakpoints and
  // calculate the split read support for each. Requires 
  // alignedReads to have been run first (will error if not run).
  void splitCoverage();
  
  // Checks if any of the indel breaks are in a blacklist. If so, mark the
  // breakpoints of the indels for skipping. That is, hasMinmal() will return false;
  void blacklist(SeqLib::GRC& grv);
  
  // Dump the contigs to a fasta
  //void printContigFasta(std::ofstream &os) const;
  
  // Set the breakpoints on the reference by combining multi-mapped contigs
  void setMultiMapBreakPairs();
  
  //! return the contig sequence as it came off the assembler
  //std::string getSequence() const; 

  // SvABA2.0 v4: emit this contig's r2c info directly into a SQLite
  // database (one row in `contigs`, plus one row per r2c-aligned read in
  // `reads`). Replaces the older printToR2CTsv() emitter — instead of
  // building a tab-separated string we'd then write to a per-thread
  // gzip stream, we bind values directly into prepared statements on
  // the per-thread R2CDatabase. Saves the build-string-then-parse-string
  // round trip (which was ~half the per-record cost of the TSV path)
  // and lets queries run on the file directly via sqlite3 / sql.js.
  //
  // The R2CDatabase passed in is the worker-local instance owned by
  // svabaThreadUnit::r2c_db_; per-thread isolation means no cross-thread
  // SQLite locking. Postprocess merges the per-thread .db files into
  // ${ID}.r2c.db via R2CDatabase::merge_from() (ATTACH + INSERT).
  void writeToR2cDb(R2CDatabase& db, const SeqLib::BamHeader& h) const;
  
  // Return if this contig contains a potential variant (indel or multi-map)
  bool hasVariant() const;
  
  // Write all of the contig alignment records to a BAM file
  //void writeToBAM(SeqLib::BamWriter& bw) const;
  
  // Write all of the sequencing reads as aligned to contig to a BAM file
  //void writeAlignedReadsToBAM(SeqLib::BamWriter& bw);

  // Remove indels that map extremely close to rearrangement break points 
  //void filterIndelsAtMultiMapSites(size_t buff);

  // returns whether any of the alignments are local to m_window
  //bool hasLocal() const;

  // Retrieves all of the breakpoints by combining indels with global mutli-map break
  BreakPointPtrVector getAllBreakPoints() const;

  //std::vector<BreakPoint> getAllBreakPointsSecondary() const;

  //  void refilterComplex();

  std::vector<const BreakPoint*> getAllBreakPointPointers() const ;

  void addDiscordantCluster(DiscordantClusterMap& dmap);
  
  std::pair<int, int> getCoverageAtPosition(int pos) const;

  // add a new read aligned to this contig
  // this is a svabaRead (read to genome), but with an r2c in it
  void AddAlignedRead(svabaReadPtr& br);

  // return number of bam reads
  size_t NumBamReads() const { return m_bamreads.size(); }

 private:

  // int insertion_against_contig_read_count = 0;

  // int deletion_against_contig_read_count = 0;

  // store all of the reads aligned to contig
  // these are the same alignments as the BAM, but have an r2c as well
  svabaReadPtrVector m_bamreads; 

  // coverage of each base in contig, whether it has alignment 
  std::vector<int> aligned_coverage; 

  AlignmentFragmentVector m_frag_v; // store all of the individual alignment fragments 

  std::vector<BreakPointPtr> m_local_breaks; // store all of the multi-map BreakPoints for this contigs 

  //std::vector<BreakPoint> m_local_breaks_secondaries; // store all of the multi-map BreakPoints for this contigs 

  BreakPointPtr m_global_bp;  // store the single spanning BreakPoing for this contig e

  //size_t m_index_of_stored_seq = 0; // which alignment did mseq come from?
  
  //std::vector<BreakPoint> m_global_bp_secondaries;  // store the single spanning BreakPoing for this contig e

  std::string m_seq; // sequence of contig as it came off of assembler
  
  std::vector<DiscordantCluster> m_dc; // collection of all discordant clusters that map to same location as this contig

  const SvabaSharedConfig* sc;
  
};

typedef std::unordered_map<std::string, AlignedContig> ContigMap;
typedef std::vector<AlignedContig> AlignedContigVec; 
