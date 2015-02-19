#include "SnowUtils.h"

  
// get an integer tag that might be separted by "x"
std::vector<int> SnowUtils::GetIntTag(const BamAlignmentUP& a, const std::string tag) {
  
  std::vector<int> out;
  std::string tmp;
  
  assert(a->GetTag(tag, tmp));
  if (tmp.find("x") != std::string::npos) {
    std::istringstream iss(tmp);
    std::string line;
    while (std::getline(iss, line, 'x')) {
      try { out.push_back(stoi(line)); } catch (...) { std::cerr << "Failed to read tag " << tag << " for value " << tmp << std::endl; std::exit(EXIT_FAILURE); }
      }
  } else {
    try { out.push_back(stoi(tmp)); } catch (...) { std::cerr << "Failed to read tag " << tag << " for value " << tmp << std::endl; std::exit(EXIT_FAILURE); }
  }

  assert(out.size());
  return out;
  
}

// get a string tag that might be separted by "x"
std::vector<std::string> SnowUtils::GetStringTag(const BamAlignmentUP& a, const std::string tag) {
  
  std::vector<std::string> out;
  std::string tmp;
  
  assert(a->GetTag(tag, tmp));
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
void SnowUtils::SmartAddTag(BamAlignmentUP &a, const std::string tag, const std::string val) {

  if (a->HasTag(tag)) {
    std::string tmp;
    a->GetTag(tag, tmp);
    a->EditTag(tag,"Z",tmp + "x" + val);
  } else {
    a->AddTag(tag,"Z", val);
  }

}
