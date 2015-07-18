#ifndef SNOWMAN_ASSEMBLER_ENGINE_H__
#define SNOWMAN_ASSEMBLER_ENGINE_H__

#include "Util.h"
#include "ReadTable.h"
#include "contigs.h"

#include "SnowTools/BamRead.h"

class SnowmanAssemblerEngine
{
 public:
  
  SnowmanAssemblerEngine(const std::string& id, double er, size_t mo, size_t rl) : m_id(id), m_error_rate(er), m_min_overlap(mo), m_readlen(rl) {}
    
    void fillReadTable(SnowTools::BamReadVector& r);

    bool performAssembly();

    void doAssembly(ReadTable *pRT, ContigVector &contigs, int pass);

    void setToWriteASQG() { m_write_asqg = true; }

    ContigVector getContigs() const { return m_contigs; }

 private:

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
    
    SnowTools::BamReadVector m_reads;
    
    ContigVector m_contigs;

};

#endif
