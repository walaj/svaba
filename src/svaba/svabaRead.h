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
  //bool supports_var = false; // does this support a variant?
  bool is_split = false; // is this a split read?

  std::unordered_map<std::string, bool> supports_var; // does this r2c support a variant (for a particular variant, 
  // since one contig can have multiple variants

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

  int SeqLength() const { return strlen(seq.get()); }

  void AddR2C(const std::string& breakpoint_name, const r2c& r) {
    if (!m_r2c) 
      m_r2c = SeqPointer<R2CMap>(new R2CMap());
    m_r2c->insert({breakpoint_name, r});
  }

  r2c& GetR2C(const std::string& contig_name) const {
    assert(m_r2c);
    R2CMap::const_iterator ff = m_r2c->find(contig_name);
    assert(ff != m_r2c->end());
    return m_r2c->find(contig_name)->second;
  }

 private:

  SeqLib::BamRecord r;

  SeqPointer<char> seq;

  char p[4]; // prefix for file ID (e.g. t001)
  
  int dd = 0; // discordant read status 0 

  SeqPointer<R2CMap> m_r2c; // store the r2c alignment information. key is contig name

};

typedef std::vector<svabaRead> svabaReadVector;

#endif
