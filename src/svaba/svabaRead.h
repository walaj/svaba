#ifndef SVABA_READ_H
#define SVABA_READ_H

#include "SeqLib/BamRecord.h" 

#include <vector>
#include <unordered_map>

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

  void AddAlignment (const SeqLib::BamRecord& b) {
    start_on_contig = b.Position();
    end_on_contig = b.PositionEnd();
    start_on_read = b.AlignmentPosition();
    cig = b.GetCigar();
  }
 
  friend std::ostream& operator<<(std::ostream& out, const r2c& a);
};

typedef std::unordered_map<std::string, r2c> R2CMap;

class svabaRead : public SeqLib::BamRecord {

 public:
  
  svabaRead();
  
  svabaRead(const SeqLib::BamRecord r, const std::string& prefix);

  std::string Seq() const;

  std::string Prefix() const;

  void SetSeq(const std::string& nseq);
  
  std::string SR() const;

  int GetDD() const { return dd; }

  void SetDD(int d) { dd = d; }

  bool Tumor() const { return p[0] == 't'; }

  int SeqLength() const; 

  void AddR2C(const std::string& contig_name, const r2c& r);

  r2c& GetR2C(const std::string& contig_name) const;

 private:

  SeqLib::BamRecord r;

  SeqPointer<char> seq;

  char p[4]; // prefix for file ID (e.g. t001)
  
  int dd = 0; // discordant read status 0 

  SeqPointer<R2CMap> m_r2c; // store the r2c alignment information. key is contig name

};

typedef std::vector<svabaRead> svabaReadVector;

#endif
