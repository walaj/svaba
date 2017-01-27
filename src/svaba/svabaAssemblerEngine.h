#ifndef SVABA_ASSEMBLER_ENGINE_H__
#define SVABA_ASSEMBLER_ENGINE_H__

#include "Util.h"
//#include "contigs.h"
#include "SGUtil.h"
#include "ReadTable.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/UnalignedSequence.h"

class svabaAssemblerEngine
{
 public:

  svabaAssemblerEngine() {}
  
  svabaAssemblerEngine(const std::string& id, double er, size_t mo, size_t rl) : m_id(id), m_error_rate(er), m_min_overlap(mo), m_readlen(rl) {}
  
  bool hasRepeat(const std::string& seq);
  
  void fillReadTable(SeqLib::BamRecordVector& r);
  
  void fillReadTable(const std::vector<std::string>& r);
  
  bool performAssembly(int num_assembly_rounds);
  
  //void doAssembly(ReadTable *pRT, ContigVector &contigs, int pass);
  void doAssembly(ReadTable *pRT, SeqLib::UnalignedSequenceVector &contigs, int pass);
  
  void setToWriteASQG() { m_write_asqg = true; }
  
  SeqLib::UnalignedSequenceVector getContigs() const { return m_contigs; }
  //ContigVector getContigs() const { return m_contigs; }
  
  void clearContigs() { m_contigs.clear(); }

  ReadTable* removeDuplicates(ReadTable* pRT);

  void calculateSeedParameters(int read_len, const int minOverlap, int& seed_length, int& seed_stride) const;

 private:

  void print_results(const SeqLib::UnalignedSequenceVector& cc) const;

  // void remove_exact_dups(ContigVector& cc) const;
  void remove_exact_dups(SeqLib::UnalignedSequenceVector& cc) const;

  void write_asqg(const StringGraph * oGraph, std::stringstream& asqg_stream, std::stringstream& hits_stream, int pass) const;
  
  std::string m_id;
  double m_error_rate;
  size_t m_min_overlap;
  size_t m_readlen;
  
  size_t numBubbleRounds = 3;
  float divergence = 0.00; //0.05
  float gap_divergence = 0.00;
  int maxEdges = 128;
  int numTrimRounds = 0; //
  int trimLengthThreshold = -1; // doesn't matter
  bool bPerformTR = false; // transitivie edge reducetion
  bool bValidate = false;
  int resolveSmallRepeatLen = -1; 
  int maxIndelLength = 20;
  bool bExact = true;
  std::string outVariantsFile = ""; // dummy
  
  bool m_write_asqg = false;
  
  ReadTable m_pRT;
  
  SeqLib::BamRecordVector m_reads;
  
  //ContigVector m_contigs;
  SeqLib::UnalignedSequenceVector m_contigs;

};

#endif
