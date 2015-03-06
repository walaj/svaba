#include "SnowUtils.h"
  
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
