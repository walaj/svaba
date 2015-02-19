#include "MiniRules.h"
#include "VariantBamReader.h"
#include <regex>
#include "gzstream.h"

using namespace std;
using namespace BamTools;

  // define what is a valid condition
static const unordered_map<string,bool> valid = 
  { 
  {"duplicate",     true},
  {"supplementary", true},
  {"qcfail",        true},
  {"hardclip",      true},
  {"fwd_strand",    true},
  {"rev_strand",    true},
  {"mate_fwd_strand",  true},
  {"mate_rev_strand",  true},
  {"mapped",           true},
  {"mate_mapped",      true},
  {"isize", true},
  {"clip",  true},
  {"phred", true},
  {"length",   true},
  {"nm",    true},
  {"mapq",  true},
  {"all",   true},
  {"ff", true},
  {"fr", true},
  {"rr", true},
  {"rf", true},
  {"ic", true},
  {"discordant", true},
  {"seq", true},
  {"nbases", true},
  {"ins", true},
  {"del", true}
};


bool MiniRules::isValid(BamAlignment &a) {

  for (auto& it : m_abstract_rules)
    if (it.isValid(a)) 
       return true; // it is includable in at least one. 
      
  return false;

}

// check whether a BamAlignment (or optionally it's mate) is overlapping the regions
// contained in these rules
bool MiniRules::isOverlapping(BamAlignment &a) {

  // if this is a whole genome rule, it overlaps
  if (m_whole_genome)
    return true;

  // TODO fix a.MatePosition + a.Length is using wrong length

  // check whether a read (or maybe its mate) hits a rule
  GenomicIntervalVector grv;
  if (m_tree.count(a.RefID) == 1) // check that we have a tree for this chr
    m_tree[a.RefID].findOverlapping(a.Position, a.Position + a.Length, grv);
  if (m_tree.count(a.MateRefID) == 1 && m_applies_to_mate) // check that we have a tree for this chr
    m_tree[a.MateRefID].findOverlapping (a.MatePosition, a.MatePosition + a.Length, grv);
  return grv.size() > 0;
  
}

// checks which rule a read applies to (using the hiearchy stored in m_regions).
// if a read does not satisfy a rule it is excluded.
string MiniRulesCollection::isValid(BamAlignment &a) {

  if (m_regions.size() == 0) {
    cerr << "Empty MiniRules" << endl;
    exit(EXIT_FAILURE);
  }

  size_t which_region = 0;
  size_t which_rule = 0;
  
  // find out which rule it is a part of
  // lower number rules dominate

  for (auto it : m_regions) {
    which_rule = 0;
    bool rule_hit = false;
    if (it->isOverlapping(a)) // read overlaps a region
      for (auto& jt : it->m_abstract_rules) { // loop rules in that region
	if (jt.isValid(a)) {
	  rule_hit = true;
	  break;
	}
	which_rule++;
      }

    // found a hit in a rule
    if (rule_hit)
      break;
    // didnt find hit, move it up one
    which_region++;
  }
  
  // isn't in a rule or it never satisfied one. Remove
  if (which_region >= m_regions.size())
    return ""; 

  string out = "rg" + to_string(++which_region) + "rl" + to_string(++which_rule);
  return out; 
  
}

// convert a region BED file into an interval tree map
void MiniRules::setIntervalTreeMap(string file) {
  
  m_region_file = file;
  GenomicRegionVector grv = GenomicRegion::regionFileToGRV(file, pad);
  m_grv = GenomicRegion::mergeOverlappingIntervals(grv); 
  sort(m_grv.begin(), m_grv.end());

  // set the width
  for (auto& it : m_grv)
    m_width += it.width();
 
  size_t grv_size = m_grv.size();
  if (grv_size == 0) {
    cerr << "Warning: No regions dected in file: " << file << endl;
    return;
  }

  m_tree = GenomicRegion::createTreeMap(m_grv);
  return;
}

// constructor to make a MiniRulesCollection from a rules file.
// This will reduce each individual BED file and make the 
// GenomicIntervalTreeMap
MiniRulesCollection::MiniRulesCollection(string file) {

  // parse the rules file
  vector<string> region_files;
  ifstream iss_file(file.c_str());
  char delim = '%';
  if (iss_file) {
    string temp;
    file = "";
    while(getline(iss_file, temp)) {
      file += temp + "\n";
    }
    iss_file.close();
    delim = '\n';
  }
  istringstream iss_rules(file.c_str());
  
  // loop through the rules file and grab the rules
  string line;
  int level = 1;

  // define a default rule set
  vector<AbstractRule> all_rules;

  // default a default rule
  AbstractRule rule_all;
  
  while(getline(iss_rules, line, delim)) {

    //exclude comments and empty lines
    bool line_empty = line.find_first_not_of("\t\n ") == string::npos;
    bool line_comment = false;
    if (!line_empty)
      line_comment = line.at(0) == '#';
    
    if (!line_comment && !line_empty) {

      // check that it doesn't have too many rules on it
      if (count(line.begin(), line.end(), '@') > 1) {
	cerr << "ERROR: Every line must start with region@, global@ or specify a rule, and only one region/rule per line" << endl;
	cerr << "  If separating lines in -r flag, separate with %. If in file, use \\n" << endl;
	cerr << "  Offending line: " << line << endl;
	exit(EXIT_FAILURE);
      }

      //////////////////////////////////
      // its a rule line, get the region
      //////////////////////////////////
      if (line.find("region@") != string::npos) {
	
	// check that the last one isnt empty. 
	// if it is, add the global to it
	if (m_regions.size() > 0)
	  if (m_regions.back()->m_abstract_rules.size() == 0)
	    m_regions.back()->m_abstract_rules.push_back(rule_all);

	// start a new MiniRule set
	MiniRules * mr = new MiniRules();
	
	// add the defaults
	//mr->m_abstract_rules = all_rules;

	// check if the mate aplies
	if (line.find(";mate") != string::npos) {
	  mr->m_applies_to_mate = true;
	}
	// check if we should pad 
	regex reg_pad(".*?;pad\\[([0-9]+)\\].*");
	smatch pmatch;
	if (regex_search(line,pmatch,reg_pad))
	  try { mr->pad = stoi(pmatch[1].str()); } catch (...) { cerr << "Cant read pad value for line " << line << ", setting to 0" << endl; }
	  

	if (line.find("@WG") != string::npos) {
	  mr->m_whole_genome = true;
        } else {
	  regex file_reg("region@(.*?)(;|$)");
	  smatch match;
	  if (regex_search(line,match,file_reg))
	    mr->setIntervalTreeMap(match[1].str());
	  else {
	    cerr << "Could not parse line: " << line << " to grab region " << endl;
	    exit(EXIT_FAILURE);
	  }
	}
	mr->m_level = level++;
	m_regions.push_back(mr);
      }
      ////////////////////////////////////
      // its a global rule
      ///////////////////////////////////
      else if (line.find("global@") != string::npos) {
	rule_all.parseRuleLine(line);
      }
      ////////////////////////////////////
      // its an rule
      ////////////////////////////////////
      else {
	AbstractRule ar = rule_all;

	// parse the line
	ar.parseRuleLine(line);
	m_regions.back()->m_abstract_rules.push_back(ar);

	// check for "discordant" shortcut
	regex  regex_disc( ".*?discordant\\[([0-9]+),([0-9]+)\\]($|;)");
	smatch omatch;
	if (regex_search(line, omatch, regex_disc)) {
	  bool isneg = line.find("!discordant[") != string::npos;
	  try {
	    // fill in the isize condition 
	    m_regions.back()->m_abstract_rules.back().isize.min = stoi(omatch[1].str());
	    m_regions.back()->m_abstract_rules.back().isize.max = stoi(omatch[2].str());
	    m_regions.back()->m_abstract_rules.back().isize.inverted = !isneg;
	    m_regions.back()->m_abstract_rules.back().isize.every = false;
	    // use the template to set a sequene of orientation rules
	    AbstractRule aro = ar;
	    if (isneg)
	      aro.fr.ff.setOff();
	    else
	      aro.fr.ff.setOn();
	    aro.fr.na = false;
	    m_regions.back()->m_abstract_rules.push_back(aro);
	    // set another rr rule
	    aro = ar;
	    if (isneg)
	      aro.fr.rr.setOff();
	    else
	      aro.fr.rr.setOn();
	    aro.fr.na = false;
	    m_regions.back()->m_abstract_rules.push_back(aro);
	    // set another rf rule
	    aro = ar;
	    if (isneg)
	      aro.fr.rf.setOff();
	    else
	      aro.fr.rf.setOn();
	    aro.fr.na = false;
	    m_regions.back()->m_abstract_rules.push_back(aro);
	    // set another ic rule
	    aro = ar;
	    if (isneg)
	      aro.fr.ic.setOff();
	    else
	      aro.fr.ic.setOn();
	    aro.fr.na = false;
	    m_regions.back()->m_abstract_rules.push_back(aro);
	  } catch (...) {
	    cerr << "Caught error trying to parse for discordant " << " on line " << line << " match[1] " << omatch[1].str() << " match[2] " << omatch[2].str() << endl;     
	    exit(EXIT_FAILURE);
	  }
	} // end discorant regex
	  
      }


    } //end comment check
  } // end \n parse

  // check that the last one isnt empty. 
  // if it is, add the global to it
  if (m_regions.size() > 0)
    if (m_regions.back()->m_abstract_rules.size() == 0)
      m_regions.back()->m_abstract_rules.push_back(rule_all);
  
  
  
}

// print the MiniRulesCollectoin
ostream& operator<<(ostream &out, const MiniRulesCollection &mr) {

  cout << "----------MiniRulesCollection-------------" << endl;
  for (auto it : mr.m_regions)
    out << (*it);
  cout << "------------------------------------------" << endl;
  return out;

}

// print a MiniRules information
ostream& operator<<(ostream &out, const MiniRules &mr) {
  
  string file_print = mr.m_whole_genome ? "WHOLE GENOME" : SnowUtils::getFileName(mr.m_region_file);
  out << "--Region: " << file_print;;
  if (!mr.m_whole_genome) {
    out << " --Size: " << SnowUtils::AddCommas<int>(mr.m_width); 
    out << " --Pad: " << mr.pad;
    out << " --Include Mate: " << (mr.m_applies_to_mate ? "ON" : "OFF") << endl;
  } else {
    out << endl;
  }
  for (auto it : mr.m_abstract_rules) 
    out << it << endl;
  
  return out;
}

// merge all of the intervals into one and send to a bed file
void MiniRulesCollection::sendToBed(string file) {

  ofstream out(file);
  if (!out) {
    cerr << "Cannot write BED file: " << file << endl;
    return;
  }
  out.close();

  GenomicRegionVector merged = sendToGrv();
  // send to BED file
  GenomicRegion::sendToBed(merged, file);
  return;
}

// parse a rule line looking for flag values
void FlagRule::parseRuleLine(string line) {

  istringstream iss(line);
  string val;
  while (getline(iss, val, ';')) {
    regex reg_dup("^!?dup.*");
    regex reg_sup("^!?supp.*");
    regex reg_qc("^!?qcfail$");
    regex reg_fs("^!?fwd_strand$");
    regex reg_hc("^!?hardclip.*");
    regex reg_rs("^!?rev_strand$");
    regex reg_mf("^!?mate_fwd_strand$");
    regex reg_mr("^!?mate_rev_strand$");
    regex reg_mp("^!?mapped$");
    regex reg_mm("^!?mate_mapped$");
    regex reg_ff("^!?ff$");
    regex reg_fr("^!?fr$");
    regex reg_rf("^!?rf$");
    regex reg_rr("^!?rr$");
    regex reg_ic("^!?ic$");

    if (dup.parseRuleLine(val, reg_dup))   na = false;
    if (supp.parseRuleLine(val, reg_sup))  na = false;
    if (qcfail.parseRuleLine(val, reg_qc)) na = false;
    if (hardclip.parseRuleLine(val, reg_hc))        na = false;
    if (fwd_strand.parseRuleLine(val, reg_fs))      na = false;
    if (mate_rev_strand.parseRuleLine(val, reg_mr)) na = false;
    if (mate_fwd_strand.parseRuleLine(val, reg_mf)) na = false;
    if (mate_mapped.parseRuleLine(val, reg_mm))     na = false;
    if (mapped.parseRuleLine(val, reg_mp)) na = false;
    if (ff.parseRuleLine(val, reg_ff))     na = false;
    if (fr.parseRuleLine(val, reg_fr))     na = false;
    if (rf.parseRuleLine(val, reg_rf))     na = false;
    if (rr.parseRuleLine(val, reg_rr))     na = false;
    if (ic.parseRuleLine(val, reg_ic))     na = false;
  }

}

// modify the rules based on the informaiton provided in the line
void AbstractRule::parseRuleLine(string line) {

  // get everything but the global keyword, if there is one
  regex reg_noname("global@(.*)");
  smatch nnmatch;
  string noname;
  if (regex_search(line, nnmatch, reg_noname)) {
    noname = nnmatch[1].str();
  } else {
    noname = line;
  }

  // check that the conditoins are valid
  istringstream iss_c(noname);
  string tmp;
  while (getline(iss_c, tmp, ';')) {
    regex reg(".*?!?([a-z_]+).*");
    smatch cmatch;
    if (regex_search(tmp, cmatch, reg)) {
      if (valid.count(cmatch[1].str()) == 0) {
	cerr << "Invalid condition of: " << tmp << endl;
	exit(EXIT_FAILURE);
      }
    } else {
      cerr << "Invalid condition of: " << tmp << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  // check for every/none flags
  if (noname.find("all") != string::npos) 
    setEvery();
  if (noname.find("!all") != string::npos)
    setNone();

  // modify the ranges if need to
  isize.parseRuleLine(noname);
  mapq.parseRuleLine(noname);
  len.parseRuleLine(noname);
  clip.parseRuleLine(noname);
  phred.parseRuleLine(noname);
  nbases.parseRuleLine(noname);
  ins.parseRuleLine(noname);
  del.parseRuleLine(noname);
  nm.parseRuleLine(noname);

  // parse the line for flag rules (also checks syntax)
  fr.parseRuleLine(noname);

  parseSeqLine(noname);
  
}

// parse for range
void Range::parseRuleLine(string line) {
  
  istringstream iss(line);
  string val;
  while (getline(iss, val, ';')) {
    
    string i_reg_str = "!" + pattern + ":?\\[(.*?),(.*?)\\]";
    string   reg_str = pattern + ":?\\[(.*?),(.*?)\\]";
    
    string n_reg_str = pattern + ":?!all";
    string a_reg_str = pattern + ":?all";
    
    regex ireg(i_reg_str);
    regex  reg(reg_str);
    regex nreg(n_reg_str);
    regex areg(a_reg_str);
    
    smatch match;
    if (regex_search(val, match, areg)) {
      setEvery();
    } else if (regex_search(val, match, ireg)) {
      try {
	min = stoi(match[1].str());
	max = stoi(match[2].str());
	inverted = true;
	every = false; none = false;
	return;
      } catch (...) {
	cerr << "Caught error trying to parse inverted for " << pattern << " on line " << line << " match[1] " << match[1].str() << " match[2] " << match[2].str() << endl;     
	exit(EXIT_FAILURE);
      }
    } else if (regex_search(val, match, reg)) {
      try {
	min = stoi(match[1].str());
	max = stoi(match[2].str());
	inverted = false;
	every = false; none = false;
	return;
      } catch (...) {
	cerr << "Caught error trying to parse for " << pattern << " on line " << line << " match[1] " << match[1].str() << " match[2] " << match[2].str() << endl;     
	exit(EXIT_FAILURE);
      }
    }
    
  } // end getline
}

// main function for determining if a read is valid
bool AbstractRule::isValid(BamAlignment &a) {
  
  // check if its keep all or none
  if (isEvery())
    return true;

  // check if is discordant
  bool isize_pass = isize.isValid(abs(a.InsertSize));

  if (!isize_pass) {
    return false;
  }

  // check for valid mapping quality
  if (!mapq.isEvery())
    if (!mapq.isValid(a.MapQuality)) 
      return false;

  // check for valid flags
  if (!fr.isValid(a))
    return false;

  // check the CIGAR
  if (!ins.isEvery() || !del.isEvery()) {
    int imin = 1000; int imax = -1;
    int dmin = 1000; int dmax = -1;
    for (auto& i : a.CigarData) {
      if (i.Type == 'I') {
	imin = (i.Length < imin) ? i.Length : imin;
	imax = (i.Length > imax) ? i.Length : imax;
      } else if (i.Type == 'D') {
	dmin = (i.Length < dmin) ? i.Length : dmin;
	dmax = (i.Length > dmax) ? i.Length : dmax;
      }
    }
    if (ins.isValid(imax))
      return false;
    if (del.isValid(dmax))
      return false;
  }
  
  // if we dont need to because everything is pass, just just pass it
  bool need_to_continue = !nm.isEvery() || !clip.isEvery() || !len.isEvery() || !nbases.isEvery();
  if (!need_to_continue)
    return true;

  // now check if we need to build char if all we want is clip
  unsigned clipnum = 0;
  if (!clip.isEvery()) {
    clipnum = VariantBamReader::getClipCount(a);
    if (nm.isEvery() && len.isEvery() && !clip.isValid(clipnum)) // if clip fails, its not going to get better by trimming. kill it now before building teh char data
      return false;
  }

  if (a.Name == "") {// only build once 
    a.BuildCharData();
  }

  // check for valid NM
  if (!nm.isEvery()) {
    uint32_t nm_val;
    if (!a.GetTag("NM",nm_val))
      nm_val = 0;
    int nm2 = nm_val;
    if (!nm.isValid(nm2))
      return false;
  }

  // trim the read, then check length
  int new_len = a.QueryBases.length();
  int new_clipnum = clipnum;

  if (!phred.isEvery()) {
    string trimmed_bases, trimmed_quals;
    if (!a.GetTag("TS", trimmed_bases)) {
      trimmed_bases = a.QueryBases;
      trimmed_quals = a.Qualities;
      VariantBamReader::qualityTrimRead(phred.min, trimmed_bases, trimmed_quals); 
      a.AddTag("TS", "Z", trimmed_bases);
    } 
    new_len = trimmed_bases.length();
    new_clipnum = max(0, static_cast<int>(clipnum - (a.Length - new_len)));
    
    if (new_len == 0)
      return false;

    // check the N
    if (!nbases.isEvery()) {
      size_t n = count(trimmed_bases.begin(), trimmed_bases.end(), 'N');
      if (!nbases.isValid(n))
	return false;
    }

  }

  // check the N if we didn't do phred trimming
  if (!nbases.isEvery() && phred.isEvery()) {
    size_t n = count(a.QueryBases.begin(), a.QueryBases.end(), 'N');
    if (!nbases.isValid(n))
      return false;
  }

  // check for valid length
  if (!len.isValid(new_len))
    return false;

  // check for valid clip
  if (!clip.isValid(new_clipnum))
    return false;

  if (atm_file.length()) {
    bool m = ahomatch(a.QueryBases);
    if ( (!m && !atm_inv) || (m && atm_inv) )
      return false;
  }

  return true;
}

bool FlagRule::isValid(BamAlignment &a) {
  
  if (isEvery())
    return true;
  
  if (!dup.isNA()) 
    if ((dup.isOff() && a.IsDuplicate()) || (dup.isOn() && !a.IsDuplicate()))
      return false;
  if (!supp.isNA()) 
    if ((supp.isOff() && !a.IsPrimaryAlignment()) || (supp.isOn() && a.IsPrimaryAlignment()))
      return false;
  if (!qcfail.isNA())
    if ((qcfail.isOff() && a.IsFailedQC()) || (qcfail.isOn() && !a.IsFailedQC()))
      return false;
  if (!mapped.isNA())
    if ( (mapped.isOff() && a.IsMapped()) || (mapped.isOn() && !a.IsMapped()) )
      return false;
  if (!mate_mapped.isNA())
    if ( (mate_mapped.isOff() && a.IsMateMapped()) || (mate_mapped.isOn() && !a.IsMateMapped()) )
      return false;

  // check for hard clips
  if (!hardclip.isNA())  {// check that we want to chuck hard clip
    if (a.CigarData.size() > 1) { // check that its not simple
      bool ishclipped = false;
      for (auto& cig : a.CigarData)
	if (cig.Type == 'H') {
	  ishclipped = true;
	  break;
	}
      if ( (ishclipped && hardclip.isOff()) || (!ishclipped && hardclip.isOn()) )
	return false;
    }
  }

  // check for orientation
  // check first if we need to even look for orientation
  bool ocheck = !ff.isNA() || !fr.isNA() || !rf.isNA() || !rr.isNA() || !ic.isNA();
  if ( ocheck ) {

    bool first = a.Position < a.MatePosition;
    bool bfr = (first && (!a.IsReverseStrand() && a.IsMateReverseStrand())) || (!first &&  a.IsReverseStrand() && !a.IsMateReverseStrand());
    bool brr = a.IsReverseStrand() && a.IsMateReverseStrand();
    bool brf = (first &&  (a.IsReverseStrand() && !a.IsMateReverseStrand())) || (!first && !a.IsReverseStrand() &&  a.IsMateReverseStrand());
    bool bff = !a.IsReverseStrand() && !a.IsMateReverseStrand();
      
    bool bic = a.MateRefID != a.RefID;

    // its FR and it CANT be FR (off) or its !FR and it MUST be FR (ON)
    // orienation not defined for inter-chrom, so exclude these with !ic
    if (!bic) { // PROCEED IF INTRA-CHROMOSOMAL
      if ( (bfr && fr.isOff()) || (!bfr && fr.isOn())) 
	return false;
      // etc....
      if ( (brr && rr.isOff()) || (!brr && rr.isOn())) 
	return false;
      if ( (brf && rf.isOff()) || (!brf && rf.isOn())) 
	return false;
      if ( (bff && ff.isOff()) || (!bff && ff.isOn())) 
	return false;
    }
    if ( (bic && ic.isOff()) || (!bic && ic.isOn()))
      return false;
      
  }
  return true;
  
}

// define how to print
ostream& operator<<(ostream &out, const AbstractRule &ar) {

  out << "  Rule: " << ar.name << " -- ";;
  if (ar.isEvery()) {
    out << "  KEEPING ALL" << endl;
  } else if (ar.isNone()) {
    out << "  KEEPING NONE" << endl;  } else {
    if (!ar.isize.isEvery())
      out << "isize:" << ar.isize << " -- " ;
    if (!ar.mapq.isEvery())
      out << "mapq:" << ar.mapq << " -- " ;
    if (!ar.len.isEvery())
      out << "length:" << ar.len << " -- ";
    if (!ar.clip.isEvery())
      out << "clip:" << ar.clip << " -- ";
    if (!ar.phred.isEvery())
      out << "phred:" << ar.phred << " -- ";
    if (!ar.nm.isEvery())
      out << "nm:" << ar.nm << " -- ";
    if (!ar.nbases.isEvery())
      out << "nbases:" << ar.nbases << " -- ";
    if (!ar.ins.isEvery())
      out << "ins:" << ar.ins << " -- ";
    if (!ar.del.isEvery())
      out << "del:" << ar.del << " -- ";

    if (ar.atm_file != "")
      out << "matching on " << ar.atm_count << " subsequences from file " << ar.atm_file << " -- ";
    out << ar.fr;
  }
  return out;
}

// define how to print
ostream& operator<<(ostream &out, const FlagRule &fr) {

  if (fr.isEvery()) {
    out << "Flag: ALL";
    return out;
  } 

  string keep = "Flag ON: ";
  string remo = "Flag OFF: ";

  if (fr.dup.isOff())
    remo += "duplicate,";
  if (fr.dup.isOn())
    keep += "duplicate,";

  if (fr.supp.isOff())
    remo += "supplementary,";
  if (fr.supp.isOn())
    keep += "supplementary,";

  if (fr.qcfail.isOff())
    remo += "qcfail,";
  if (fr.qcfail.isOn())
    keep += "qcfail,";

  if (fr.hardclip.isOff())
    remo += "hardclip,";
  if (fr.hardclip.isOn())
    keep += "hardclip,";

  if (fr.ic.isOff())
    remo += "ic,";
  if (fr.ic.isOn())
    keep += "ic,";

  if (fr.ff.isOff())
    remo += "ff,";
  if (fr.ff.isOn())
    keep += "ff,";

  if (fr.fr.isOff())
    remo += "fr,";
  if (fr.fr.isOn())
    keep += "fr,";

  if (fr.rr.isOff())
    remo += "rr,";
  if (fr.rr.isOn())
    keep += "rr,";

  if (fr.rf.isOff())
    remo += "rf,";
  if (fr.rf.isOn())
    keep += "rf,";

  if (fr.mapped.isOff())
    remo += "mapped,";
  if (fr.mapped.isOn())
    keep += "mapped,";

  if (fr.mate_mapped.isOff())
    remo += "mate_mapped,";
  if (fr.mate_mapped.isOn())
    keep += "mate_mapped,";


  /*

  // get the strings
  for (auto it : fr.flags) {
    if (it.second.isNA())
      na += it.first + ",";
    else if (it.second.isOn())
      keep += it.first + ",";
    else if (it.second.isOff())
      remo += it.first + ",";
    else // shouldn't get here
      exit(1); 
  }
  */
  if (!fr.isEvery())
    out << keep << " -- " << remo;

  return out;
}

// define how to print
ostream& operator<<(ostream &out, const Range &r) {
  if (r.isEvery())
    out << "all";
  else
    out << (r.inverted ? "NOT " : "") << "[" << r.min << "," << r.max << "]";
  return out;
}

// convert a MiniRulesCollection into a GRV
GenomicRegionVector MiniRulesCollection::sendToGrv() const {

  // make a composite
  GenomicRegionVector comp;
  for (auto it : m_regions)
    comp.insert(comp.begin(), it->m_grv.begin(), it->m_grv.end()); 
  
  // merge it down
  GenomicRegionVector merged = GenomicRegion::mergeOverlappingIntervals(comp);

  return merged;
}

bool Flag::parseRuleLine(string &val, regex &reg) {

  smatch match;
  if (regex_search(val, match, reg)) {
    //auto ff = flags.find(match[1].str()); 
    if (val.at(0) == '!') { // it is a val in flags and is off
      setOff();
      return true;
    } else  { // is in a val in flags and is on
      setOn();
      return true;
    } //else if (ff == flags.end() && valid.count(match[1].str()) == 0) { // its not anything and its bad
      //cerr << "Not a valid condition: " << match[1].str() << " on val " << val << endl;
      //exit(EXIT_FAILURE);
    //}
  }

  return false;
  
}

// check if a string contains a substring using Aho Corasick algorithm
bool AbstractRule::ahomatch(const string& seq) {

  // make into Ac strcut
  AC_TEXT_t tmp_text = {seq.c_str(), static_cast<unsigned>(seq.length())};
  ac_automata_settext (atm, &tmp_text, 0);

  // do the check
  AC_MATCH_t * matchp;  
  matchp = ac_automata_findnext(atm);

  if (matchp) 
    return true;
  else 
    return false;
  
  
}

// add the aho subsequences by reading in the sequence file
void AbstractRule::parseSeqLine(string line) {

  // get the sequence file out
  regex reg("^!?seq\\[(.*)\\].*");
  smatch match;
  if (regex_search(line, match, reg)) {
    line = match[1].str();
    atm_file = line;
  } else {
    return;
  }
  
  // open the sequence file
  igzstream iss(line.c_str());
  if (!iss) {
    cerr << "ERROR: Cannot read the sequence file: " << line << endl;
    exit(EXIT_FAILURE);
  }

  // should it be inverted?
  string inv;
  if (line.at(0) == '!') {
    atm_inv = true;
    inv = " -- Inverted -- ";
  }

  // initialize it
  atm = ac_automata_init();
  
  // make the Aho-Corasick key
  cout << "...generating Aho-Corasick key"  << inv << " from file " << atm_file << endl;
  string pat;
  size_t count = 0;
  while (getline(iss, pat, '\n')) {
    count++;
    AC_PATTERN_t tmp_pattern;
    tmp_pattern.astring = pat.c_str();
    tmp_pattern.length = static_cast<unsigned>(pat.length());
    ac_automata_add(atm, &tmp_pattern);
  }
  ac_automata_finalize(atm);
  cout << "Done generating Aho-Corasick key of size " << count << endl;  

  atm_count = count;
  return;

}
