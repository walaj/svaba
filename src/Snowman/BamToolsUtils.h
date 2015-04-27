#ifndef SNOWMAN_BAMTOOLS_UTILS_H__
#define SNOWMAN_BAMTOOLS_UTILS_H__

#include <string>
#include <vector>
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/algorithms/Sort.h"

using BamTools::BamAlignment;
typedef std::vector<BamTools::CigarOp> CigarOpVec;


namespace BamToolsUtils {

  std::vector<std::string> GetStringTag(const BamAlignment& a, const std::string tag);

  void SmartAddTag(BamAlignment &a, const std::string tag, const std::string val);

  CigarOpVec stringToCigar(const std::string& val);
  
  std::vector<int> GetIntTag(const BamAlignment& a, const std::string tag);

  void parseTags(const std::string& val, BamAlignment &a);

  std::string cigarToString(const CigarOpVec &cig);

  void flipCigar(CigarOpVec &cig);
}

#endif
