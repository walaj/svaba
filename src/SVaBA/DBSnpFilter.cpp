#include "DBSnpFilter.h"
#include "gzstream.h"

using namespace SeqLib;

DBSnpSite::DBSnpSite(const std::string& tchr, const std::string& pos, const std::string& rs, const std::string& ref, const std::string& alt, const BamHeader& h) {
    
    // make the genomic region
    try { 
      GenomicRegion gr(tchr, pos, pos, h);
      chr = gr.chr; 
      pos1 = std::stoi(pos);
    } catch (...) {
      //std::cerr << "DBSnpSite: Error trying to convert " << tchr << ":" << pos << " to number" << std::endl;
    }

    //m_rs = rs;
    //m_ref = ref;
    //m_alt = alt;

    if (ref.length() == 0 || alt.length() == 0) // || m_rs.length() == 0)
      std::cerr << "DBSnpSite: Is the VCF formated correctly for this entry? Ref " << ref << " ALT " << alt << " rs " << rs << std::endl;

    // insertion
    if (ref.length() == 1)
      pos2 = pos1 + 1;
    // deletion
    else
      pos2 = pos1 + ref.length() + 1;
    
  }

DBSnpFilter::DBSnpFilter(const std::string& db, const BamHeader& h) {
    
    // read in the file
    if (!read_access_test(db)) {
      std::cerr << std::endl << "**** Cannot read DBSnp database " << db << "   Expecting a VCF file" << std::endl;
      return;
    }

    // read it in
    igzstream in(db.c_str());
    if (!in) {
      std::cerr << std::endl << "**** Cannot read DBSnp database " << db << "   Expecting a VCF file" << std::endl;
      return;
    }

    std::string line;
    while (std::getline(in, line)) {

      if (line.find("#") != std::string::npos || line.length() == 0)
	continue;

      std::istringstream thisline(line);
      std::string val;
      int this_count = -1;
      std::string chr, pos, rs, ref, alt;
      while (std::getline(thisline, val, '\t')) {
	++this_count;
	switch (this_count) { 
	case 0: chr = val; break;
	case 1: pos = val; break;
	case 2: rs = val; break;
	case 3: ref = val; break;
	case 4: alt = val; break;
	}
      }

      DBSnpSite db(chr, pos, rs, ref, alt, h);

      // for now reject SNP sites
      if (ref.length() + alt.length() > 2) {
	m_sites.add(db);

	// make the hash
	cig.str(std::string());
	//cig << db.chr << "_" << db.pos1 << "_" << (db.pos2-db.pos1) << (db.m_ref.length() == 1 ? "I" : "D");
	cig << db.chr << "_" << db.pos1;
	//m_hash.insert(cig.str());
	
	m_int_hash.insert(hasher(cig.str()));
	//std::cerr << line << " hash " << cig.str() << std::endl;
	
      }
    }
    
    // build the tree
    m_sites.CreateTreeMap();

  }

  std::ostream& operator<<(std::ostream& out, const DBSnpFilter& d) {
    out << "DBSnpFilter with a total of " << AddCommas<size_t>(d.m_sites.size());
    return out;
  }

  std::ostream& operator<<(std::ostream& out, const DBSnpSite& d) {
    //out << d.chr << ":" << d.pos1 << "-" << d.pos2 << "\t" << d.m_rs << " REF " << d.m_ref << " ALT " << d.m_alt;
    out << d.chr << ":" << d.pos1 << "-" << d.pos2;;
    return out;
  }

  bool DBSnpFilter::queryHash(const std::string& h) const {
    //return m_hash.count(h);
    return m_int_hash.count(hasher(h));
  }

  bool DBSnpFilter::queryBreakpoint(BreakPoint& bp) {
    
    std::vector<int32_t> sub, que;
    GenomicRegion gr = bp.b1.gr;
    gr.Pad(2);
    GRC subject(gr);
    GRC out = subject.FindOverlaps(m_sites, sub, que, true); // true = ignore_strand
    
    if (que.size()) {
      //bp.rs = m_sites[sub[0]].m_rs;
      bp.rs = "D"; 
      //for (auto& j : que) {
      //bp.rs += m_sites[j].m_rs + "_";
      //}
      //bp.rs.pop_back(); // just drop the last comma
      return true;
    }
    return false;
  }
