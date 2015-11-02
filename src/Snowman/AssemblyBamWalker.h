#ifndef SNOWMAN_ASSEMBLY_BAM_WALKER_H__
#define SNOWMAN_ASSEMBLY_BAM_WALKER_H__

#include "SnowTools/BamWalker.h"
#include "SnowTools/AlignedContig2.h"
#include "SnowmanBamWalker.h"

/** Walk along a BAM file generated from a de novo assembly,
 * realigned to the genome (preferably by BWA-MEM).
 */
class AssemblyBamWalker: public SnowTools::BamWalker 
{
 public:
  
  //SnowTools::Bam5
  
  /** Construct a new read-only walker to move along the assembly bam
   * @param in File path of the assembly BAM
   */
  AssemblyBamWalker(const std::string& in) : SnowTools::BamWalker(in) {}

    /** Move along a BAM file generated from Discovar and make the AlignedContigs
     */
    void walkDiscovar();
      
    void sendThread();

    int numThreads = 1;


    faidx_t * findex = nullptr;

    SnowmanBamWalker twalk, nwalk;

    std::string tbam, nbam;

    std::shared_ptr<hts_idx_t> tindex, nindex;

};

bool runAC(SnowTools::BamReadVector& brv, faidx_t * f, std::shared_ptr<hts_idx_t> pt, std::shared_ptr<hts_idx_t> pn,
	   const std::string& t, const std::string& n, const SnowTools::GenomicRegionVector& regions);

class AssemblyWalkerWorkItem {

 private:
  
  SnowTools::BamReadVector m_brv;
  int m_number;  
  faidx_t * m_findex;
  std::shared_ptr<hts_idx_t> m_tindex, m_nindex;  
  SnowTools::GenomicRegionVector m_regions;
  std::string m_t, m_n;


 public:
 AssemblyWalkerWorkItem(SnowTools::BamReadVector& brv, int number, faidx_t *f, std::shared_ptr<hts_idx_t> pt, std::shared_ptr<hts_idx_t> pn, SnowTools::GenomicRegionVector r,
			const std::string& t, const std::string& n)  
   : m_brv(brv), m_number(number), m_findex(f), m_tindex(pt), m_nindex(pn), m_regions(r), m_t(t), m_n(n) {}
    ~AssemblyWalkerWorkItem() {}
    
    int getNumber() { return m_number; }
    
    bool run() { return runAC(m_brv, m_findex, m_tindex, m_nindex, m_t, m_n, m_regions); }
    
};



#endif
