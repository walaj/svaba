#pragma once

#include <vector>
#include <cstring>
#include <iostream>
#include <unordered_set>
#include <functional> //for std::hash

#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/BamHeader.h"

#include "BreakPoint.h"

class DBSnpSite: public SeqLib::GenomicRegion {

  public:

  DBSnpSite(const std::string& tchr, const std::string& pos, const std::string& rs, const std::string& ref, const std::string& alt, const SeqLib::BamHeader& h);

    //std::string m_rs;
    //std::string m_ref;
    //std::string m_alt;

    friend std::ostream& operator<<(std::ostream& out, const DBSnpSite& d);
    
  };

typedef SeqLib::GenomicRegionCollection<DBSnpSite> DBC;

  class DBSnpFilter {
    
  public:

    DBSnpFilter() {}

    DBSnpFilter(const std::string& db, const SeqLib::BamHeader& h); 

    /** Test whether the variant overlaps a DBSnp site 
     * If it does, fill the BreakPoint rs field
     */
    bool queryBreakpoint(BreakPoint& bp);

    bool queryRead(const SeqLib::BamRecord& r) const;

    bool queryHash(const std::string& r) const;
    
    friend std::ostream& operator<<(std::ostream& out, const DBSnpFilter& d);

  private:

    std::hash<std::string> hasher;
    
    // initialize here once
    std::stringstream cig;

    DBC m_sites;
    std::unordered_set<std::string> m_hash;

    std::unordered_set<size_t> m_int_hash;
    
  };

