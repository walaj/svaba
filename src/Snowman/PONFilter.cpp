#include "PONFilter.h"

#include <regex>
#include <sstream>

#include "gzstream.h"
#include "SeqLib/SeqLibUtils.h"

using namespace SeqLib;

  std::ostream& operator<<(std::ostream& out, const PONFilter& p) {
    size_t max_samples = 0;
    for (auto& i : p.m_map)
      if (i.second > max_samples)
	max_samples = i.second;
    out << "Indel PON Num sites: " << AddCommas(p.m_map.size()) << " Max Samples Found " << AddCommas(max_samples);
    return out;
  }

  PONFilter::PONFilter(const std::string& file) {

    // import the pon
    igzstream izp(file.c_str());
    if (!izp) {
      std::cerr << "Can't read file " << file << std::endl;
      exit(EXIT_FAILURE);
    }

    std::string pval;
    while (std::getline(izp, pval, '\n')) {

      std::istringstream gg(pval);
      std::string tval;
      std::string key;

      int sample_count_total = 0;

      size_t c = 0;
      while (std::getline(gg, tval, '\t')) {
	++c;
	if (c == 1 && tval.length()) {
	  key = tval.substr(1,tval.length() - 1);
	  if (key.at(0) == 'T') // only accumulate normal
	    break;
	}
	else if (tval.length())
	  try { 
	    sample_count_total += (stoi(tval) > 0 ? 1 : 0); 
	  } catch(...) { 
	    std::cerr << "stoi error in PON read with val " << tval << " on line " << pval << std::endl; 
	  }
	//else if (tval.length() && c==2)
	//	try { read_count_total += stoi(tval); } catch(...) { std::cerr << "stoi error in PON read with val " << tval << " on line " << pval << std::endl; }
	
      }
      
      // trim it down
      std::regex regc("(.*?_.*?)_.*");
      std::smatch smatchr;
      if (!std::regex_search(key, smatchr, regc)) {
	std::cerr << "regex failed on " << key << std::endl; 
      } else {
	key = smatchr[1].str();
      }

      if (sample_count_total > 1) {
	m_map[key] = sample_count_total;
      }

    }

    
  }
