#ifndef SNOWTOOLS_PONFILTER_H__
#define SNOWTOOLS_PONFILTER_H__

#include <string>
#include <unordered_map>
#include <cstdlib>

namespace SnowTools {

  class PONFilter {
    
  public: 

    PONFilter() {}
    
    PONFilter(const std::string& file);

    friend std::ostream& operator<<(std::ostream& out, const PONFilter& p);

    bool count(const std::string& s) const { return m_map.count(s); }

    int NSamps(const std::string& s) const { 
      std::unordered_map<std::string, size_t>::const_iterator it = m_map.find(s);
      if (it == m_map.end())
	return 0;
      return it->second; 
    }

  private:

    std::unordered_map<std::string, size_t> m_map;

  };

}

#endif
