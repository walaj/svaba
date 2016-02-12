#ifndef SNOWTOOLS_DBSNP_FILTER_H__
#define SNOWTOOLS_DBSNP_FILTER_H__

#include <cstring>
#include <vector>
#include <iostream>
#include <sstream>
#include <unordered_set>

#include "SnowTools/GenomicRegionCollection.h"
#include "BreakPoint2.h"

namespace SnowTools {

  class DBSnpSite: public GenomicRegion {

  public:

    DBSnpSite(const std::string& tchr, const std::string& pos, const std::string& rs, const std::string& ref, const std::string& alt);

    std::string m_rs;
    std::string m_ref;
    std::string m_alt;

    friend std::ostream& operator<<(std::ostream& out, const DBSnpSite& d);
    
  };

  typedef GenomicRegionCollection<DBSnpSite> DBC;

  class DBSnpFilter {
    
  public:

    DBSnpFilter() {}

    DBSnpFilter(const std::string& db); 

    /** Test whether the variant overlaps a DBSnp site 
     * If it does, fill the BreakPoint rs field
     */
    bool queryBreakpoint(BreakPoint& bp);

    bool queryRead(const BamRead& r) const;

    bool queryHash(const std::string& r) const;
    
    friend std::ostream& operator<<(std::ostream& out, const DBSnpFilter& d);

  private:
    
    // initialize here once
    std::stringstream cig;

    DBC m_sites;
    std::unordered_set<std::string> m_hash;
    
  };

}

#endif
