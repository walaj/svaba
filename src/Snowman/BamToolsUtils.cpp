#include "BamToolsUtils.h"
#include <regex>
#include <sstream>

#define BT_r_get_Z_tag(b, t, v) if (!(b).GetTag(t, v)) v = ""
#define BT_r_remove_tag(b, t) (b).RemoveTag(t);
#define BT_r_add_Z_tag(b, t, v) ((b).AddTag(t, "Z", v))
#define BT_r_add_int32_tag(b, t, v) (b).AddTag(t, "i", v)


// get a string tag that might be separted by "x"
std::vector<std::string> BamToolsUtils::GetStringTag(const BamTools::BamAlignment& a, const std::string tag) {
  
  std::vector<std::string> out;
  std::string tmp;
  
  BT_r_get_Z_tag(a, tag.c_str(), tmp);
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
void BamToolsUtils::SmartAddTag(BamTools::BamAlignment &a, const std::string tag, const std::string val) {
  
  std::string tmp;
  BT_r_get_Z_tag(a, tag.c_str(), tmp);

  if (tmp.length()) {
    tmp += "x"  + val;
    //BT_r_remove_tag(a, tag.c_str());
    //BT_r_add_Z_tag(a, tag.c_str(), tmp);
    a.EditTag(tag,"Z",tmp + "x" + val);
  } else { // normal with no x
    //BT_r_add_Z_tag(a, tag.c_str(), val);
    a.AddTag(tag,"Z", val);
  }

}

CigarOpVec BamToolsUtils::stringToCigar(const std::string& val) {
   
   std::string v = val;
   std::vector<std::string> str_vec; // #2: Search for tokens
   //boost::split(str_vec, v, boost::is_any_of("0123456789"), boost::token_compress_on ); // SplitVec == { "hello abc","ABC","aBc goodbye" }
   
   std::vector<std::string> len_vec; // #2: Search for tokens
   //boost::split(len_vec, v, boost::is_any_of("MIDSHPN"), boost::token_compress_on ); // SplitVec == { "hello abc","ABC","aBc goodbye" }

   // move through and get the values
   std::string len= "";
   std::string type = "";

   for (size_t i = 0; i < val.length(); ++i)
     {
       if (val.at(i) == 'M' ||  
	   val.at(i) == 'I' ||
	   val.at(i) == 'D' ||
	   val.at(i) == 'S' ||
	   val.at(i) == 'H' ||
	   val.at(i) == 'P' ||
	   val.at(i) == 'N') {
	 len_vec.push_back(len);
	 len = "";
	 str_vec.push_back(std::string(&val.at(i)));
       } else {
	 len.append(1, val.at(i));

       }
       
     }

   assert(len_vec.size() == str_vec.size());
   assert(len_vec.size());
   
   CigarOpVec cigop;
   for (size_t kk = 0; kk < (str_vec.size()); kk++) // first strvec is empty and last len_vec is empty (due to token orderingin cigar)
     try {
       cigop.push_back(BamTools::CigarOp(str_vec[kk].at(0), std::stoi(len_vec[kk])));
     } catch (...) {
       std::cerr << "stoi failure on: " << len_vec[kk] << " from cigar " << val << std::endl;

     }
   
   //debug
   //std::cout << "  Intiial cigar is " << val << "  final cigar is ";
   //for (auto& i : cigop)
   //  std::cout << i.Length << i.Type << "-";
   //std::cout << std::endl;

   assert(cigop.size());
   
   return cigop;
}

void BamToolsUtils::flipCigar(CigarOpVec &cig) {

   CigarOpVec new_cig;
   for (CigarOpVec::const_iterator it = cig.end() - 1; it != cig.begin() - 1; it--) 
     new_cig.push_back(*it);
   cig = new_cig;
 }


std::string BamToolsUtils::cigarToString(const CigarOpVec &cig) {
  std::stringstream cigstring;
  for (auto& i : cig)
    cigstring << i.Length << i.Type;
  return cigstring.str();
}

void BamToolsUtils::parseTags(const std::string& val, BamTools::BamAlignment &a) {

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

// get an integer tag that might be separted by "x"
std::vector<int> BamToolsUtils::GetIntTag(const BamTools::BamAlignment& a, const std::string tag) {
  
  std::vector<int> out;
  std::string tmp;
  
  BT_r_get_Z_tag(a, tag.c_str(), tmp);
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

bool BamToolsUtils::parseXP2BamAlignment(const std::string& xp, std::vector<BamTools::BamAlignment>& a)
{

  if (!xp.length())
    return true;
  
  // loop through alll of the pieces
  std::istringstream iss(xp);
  std::string val;

  while (std::getline(iss, val, ';'))
    {
      std::regex reg("([A-Z0-9]+),(\\+|\\-)([0-9]+),([A-Z0-9]+),([0-9]+),([0-9]+);?");
      std::smatch match;

      BamTools::BamAlignment b;
      
      bool regpass = false;
      try {
	regpass = std::regex_search(val, match, reg);
      } catch (...) {
	std::cerr << "Regexp error on value: " << val << " xp tag " << xp << std::endl;
      }
      
      if (regpass) {
	
	std::string chr = match[1].str();
	if (chr == "X")
	  chr = "23";
	else if (chr == "Y")
	  chr = "24";
	else if (chr == "M")
	  chr = "25";
	else if (chr.at(0) == 'G')
	  chr = "26";
	try { 
	  b.RefID = std::stoi(chr) - 1; 
	} catch(...) { 
	  std::cerr << "Error converting chr " << chr << " to int " << std::endl; 
	  return false;
	}
	
	// set the stirng
	if (match[2].str() == "+")
	  b.AlignmentFlag = 0;
	else
	  b.AlignmentFlag = 16;
	
	b.CigarData = BamToolsUtils::stringToCigar(match[4].str());
	
	try { 
	  b.Position = std::stoi(match[3].str()); 
	} catch (...) {
	  std::cerr << "Error converting XP position " << match[3].str() << " to int " << std::endl; 
	  return false;
	}
	
	try { b.MapQuality = std::stoi(match[5].str());} catch(...) { std::cerr << "Error converting XP mapq " << match[5].str() << " to int " << std::endl; }
      } else {
	std::cerr << "Failed to parse xp tag " << xp << " on alignment " << val << std::endl;
	return false;
      }
      a.push_back(b);
    }   
  return true;
}
  
