#include "SnowUtils.h"
#include <boost/algorithm/string.hpp>
#include <regex>
  
// get an integer tag that might be separted by "x"
std::vector<int> SnowUtils::GetIntTag(const Read& a, const std::string tag) {
  
  std::vector<int> out;
  std::string tmp;
  
  r_get_Z_tag(a, tag.c_str(), tmp);
  assert(tmp.length());

  if (tmp.find("x") != std::string::npos) {
    std::istringstream iss(tmp);
    std::string line;
    while (std::getline(iss, line, 'x')) {
      try { out.push_back(stoi(line)); } catch (...) { std::cerr << "Failed to read parsed int tag " << tag << " for value " << tmp << " with line " << line << std::endl; std::exit(EXIT_FAILURE); }
      }
  } else {
    try { out.push_back(stoi(tmp)); } catch (...) { std::cerr << "Failed to read int tag " << tag << " for value " << tmp << std::endl; std::exit(EXIT_FAILURE); }
  }

  assert(out.size());
  return out;
  
}

// get a string tag that might be separted by "x"
std::vector<std::string> SnowUtils::GetStringTag(const Read& a, const std::string tag) {
  
  std::vector<std::string> out;
  std::string tmp;
  
  r_get_Z_tag(a, tag.c_str(), tmp);
  assert(tmp.length());

  if (tmp.find("x") != std::string::npos) {
    std::istringstream iss(tmp);
    std::string line;
    while (std::getline(iss, line, 'x')) {
      out.push_back(line);
    }
  } else {
    out.push_back(tmp);
  }
  
  assert(out.size());
  return out;
  
}

// add a tag that might already be there, separete by 'x'
void SnowUtils::SmartAddTag(Read &a, const std::string tag, const std::string val) {
  
  std::string tmp;
  r_get_Z_tag(a, tag.c_str(), tmp);

  if (tmp.length()) {
    tmp += "x"  + val;
    r_remove_tag(a, tag.c_str());
    r_add_Z_tag(a, tag.c_str(), tmp);
    //a->EditTag(tag,"Z",tmp + "x" + val);
  } else { // normal with no x
    r_add_Z_tag(a, tag.c_str(), val);
    //a->AddTag(tag,"Z", val);
  }

}

CigarOpVec SnowUtils::stringToCigar(const std::string& val) {
   
   std::string v = val;
   std::vector<std::string> str_vec; // #2: Search for tokens
   boost::split(str_vec, v, boost::is_any_of("0123456789"), boost::token_compress_on ); // SplitVec == { "hello abc","ABC","aBc goodbye" }
   
   std::vector<std::string> len_vec; // #2: Search for tokens
   boost::split(len_vec, v, boost::is_any_of("MIDSHPN"), boost::token_compress_on ); // SplitVec == { "hello abc","ABC","aBc goodbye" }
   
   assert(len_vec.size() == str_vec.size());
   assert(len_vec.size());
   
   CigarOpVec cigop;
   for (size_t kk = 0; kk < (str_vec.size() - 1); kk++) // first strvec is empty and last len_vec is empty (due to token orderingin cigar)
     cigop.push_back(CigarOp(str_vec[kk+1].at(0), std::stoi(len_vec[kk])));
   
   assert(cigop.size());
   
   return cigop;
}

void SnowUtils::flipCigar(CigarOpVec &cig) {

   CigarOpVec new_cig;
   for (CigarOpVec::const_iterator it = cig.end() - 1; it != cig.begin() - 1; it--) 
     new_cig.push_back(*it);
   cig = new_cig;
 }


std::string SnowUtils::cigarToString(const CigarOpVec &cig) {
  std::stringstream cigstring;
  for (auto& i : cig)
    cigstring << i.Length << i.Type;
  return cigstring.str();
}

void SnowUtils::parseTags(const std::string& val, BamTools::BamAlignment &a) {

  std::regex reg_xp("^XA:Z:(.*)");
  std::regex reg_nm("^NM:[A-Za-z]:(.*)");

  std::smatch match;
  if (val.find("XP") != std::string::npos)
    if (std::regex_search(val, match, reg_xp)) {
      std::string xp = match[1].str();
      a.AddTag("XP","Z",xp);
      return;
    } 
  
  if (val.find("NM") != std::string::npos)
    if (std::regex_search(val, match, reg_nm)) {
      int nm = std::stoi(match[1].str());
      a.AddTag("NM","i",nm);
    }

}
