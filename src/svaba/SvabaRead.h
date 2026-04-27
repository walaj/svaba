#pragma once

#include "SeqLib/BamRecord.h" 

#include <vector>
#include <string_view>
#include <unordered_map>

using SeqLib::BamRecordPtr;
using SeqLib::BamRecordPtrVector;

/** Store information about a read to contig alignment */
struct r2c {

  int32_t start_on_contig = 0;  // start pos on contig (from r.Position())
  int32_t end_on_contig   = 0;  // end pos on contig (from r.PositionEnd())
  int32_t start_on_read    = 0; // start pos on read (from r.AlignmentPosition())
  int32_t end_on_read     = 0;  // end pos on read (from r.AlignmentPosition())
  bool rc = false;    // reverse complement wrt contig?
  SeqLib::Cigar cig; // cigar of read to contig
  bool supports_var = false; // does this support a variant?
  bool is_split = false; // is this a split read?
  int left_or_right = 0; //-1 read aligns on left of contig, 1 on right
  bool supports_discordant = false; // true if this is part of a discordant pair that supports the break

  // SvABA2.0: NM tag on the r2c alignment (edit distance to contig).
  // -1 sentinel means "not set" (legacy paths that filled r2c without
  // going through AddAlignment). Surfaced in alignments.txt.gz so the
  // human-readable dump shows per-read edit distance to the contig.
  int32_t nm = -1;

  void AddAlignment (const BamRecordPtr& b) {
    start_on_contig = b->Position();
    end_on_contig = b->PositionEnd();
    start_on_read = b->AlignmentPosition();
    cig = b->GetCigar();
    int _nm = 0;
    if (b->GetIntTag("NM", _nm)) nm = _nm;
  }
  
  friend std::ostream& operator<<(std::ostream& out, const r2c& a);
};

typedef std::unordered_map<std::string, r2c> R2CMap;

class svabaRead;
typedef std::shared_ptr<svabaRead> svabaReadPtr;
typedef std::vector<svabaReadPtr> svabaReadPtrVector;

class svabaRead : public SeqLib::BamRecord {

 public:
  
  svabaRead();
  
  svabaRead(const SeqLib::BamRecord& r,
	    std::string_view prefix);

  // Delete the copy constructor                                                                                                                                                      
  svabaRead(const svabaRead&) = delete;                                                                                                                                               
                                                                                                                                                                                      
  // Optionally also delete copy assignment                                                                                                                                           
  svabaRead& operator=(const svabaRead&) = delete;                                                                                                                                    
                                                                                                                                                                                      
  // Still allow move operations:                                                                                                                                                     
  svabaRead(svabaRead&&) = default;                                                                                                                                                   
  svabaRead& operator=(svabaRead&&) = default;  
  
  std::string CorrectedSeq() const;

  /// True if BFC correction or quality trimming changed the read sequence
  /// relative to the original BAM record. When false, the original BAM's
  /// CIGAR/NM is a valid native alignment (assuming the same aligner).
  /// Compares the stored seq_corrected against the BAM 4-bit encoding
  /// without constructing a second string for the original.
  bool CorrectedSeqChanged() const;

  std::string Prefix() const;

  void SetPrefix(const std::string_view pref) { p = pref; }
  
  void SetCorrectedSeq(const std::string_view nseq);  
  
  std::string UniqueName() const;

  int GetDD() const { return dd; }

  void SetDD(int d) { dd = d; }

  bool Tumor() const { return p[0] == 't'; }

  int SeqLength() const; 

  void AddR2C(const std::string& contig_name, const r2c& r);

  r2c GetR2C(const std::string& contig_name) const;

  /// Does this read have any r2c alignments?
  bool HasR2C() const { return !m_r2c.empty(); }

  /// Trim the read based on quality score and store seq in char
  void QualityTrimRead();

  int CorrectedSeqLength() const;

  // discordant read status
  // < 0 is bad discordant read (see DiscordantRealigner.h)
  // == 0 not discordant
  // 1 = good
  int dd = 0;

  // SvABA2.0: post-BFC re-alignment of the *corrected* read sequence to the
  // reference, populated by SvabaRegionProcessor before assembly. The point
  // is to give BreakPoint::splitCoverage's "r2c better than native" gate
  // an apples-to-apples comparison: both sides use the corrected read
  // sequence and the same BWA-MEM parameters that svaba uses internally,
  // rather than mixing svaba-corrected r2c against the input BAM's
  // pre-correction CIGAR/NM (which was the source of an asymmetric gate
  // letting reads with mirror r2c indels through as variant supporters).
  //
  // Sentinel: corrected_native_nm == -1 means "not populated" (e.g. read
  // with to_assemble == false, or empty corrected sequence). In that case
  // splitCoverage falls back to GetCigar()/GetIntTag("NM",...) on the
  // original BAM record.
  SeqLib::Cigar corrected_native_cig;
  int32_t       corrected_native_nm = -1;

  friend class svabaBamWalker;

  bool to_assemble = true;
  
 private:

  std::string p; // prefix for file ID (e.g. t001)

  bool train = false;
  
  std::string seq_corrected; // quality trimmed and/or error corrected

  // store the r2c alignment information. key is contig name
  R2CMap m_r2c; 

};

typedef std::vector<svabaRead> svabaReadVector;

