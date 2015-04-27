#ifndef SNOWMAN_ASSEMBLY_BAM_WALKER_H__
#define SNOWMAN_ASSEMBLY_BAM_WALKER_H__

#include "SnowTools/BamWalker.h"

/** Walk along a BAM file generated from a de novo assembly,
 * realigned to the genome (preferably by BWA-MEM).
 */
class AssemblyBamWalker: public SnowTools::BamWalker 
{
 public:
  
  /** Construct a new read-only walker to move along the assembly bam
   * @param in File path of the assembly BAM
   */
  AssemblyBamWalker(const std::string& in) : SnowTools::BamWalker(in) {}
      
};


#endif
