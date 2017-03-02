#ifndef SVABA_READ_H
#define SVABA_READ_H

#include "SeqLib/BamRecord.h" 

#include <vector>

class svabaRead : public SeqLib::BamRecord {

 public:
  
  svabaRead() : seq(nullptr) {}
  
  svabaRead(const SeqLib::BamRecord r, const std::string& prefix);
  
  std::string Seq() const;

  void SetSeq(const std::string& nseq);
  
  std::string SR() const;

  int GetDD() const { return dd; }

  void SetDD(int d) { dd = d; }

  bool Tumor() const { return p[0] == 't'; }

  int SeqLength() const { return strlen(seq.get()); }

  //void Reassign(const svabaRead& s);

 private:

  SeqLib::BamRecord r;

  SeqPointer<char> seq;

  char p[4];
  
  int dd = 0;

};

typedef std::vector<svabaRead> svabaReadVector;

#endif
