#include "DBSnpFilter.h"
#include "gzstream.h"
#include "SvabaLogger.h"

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

DBSnpFilter::DBSnpFilter(const std::string& db,
                         const SeqLib::BamHeader& header,
                         SvabaLogger& logger)
{

  logger.log(true, true "...loading the DBsnp database ", db); 
  // First, try opening the file (gzipped or not)
  igzstream in(db.c_str());
  if (!in) {
    logger.log(true,true, //toerr, tolog
      "ERROR: Cannot open DBSNP database file '", db, "' for reading."
    );
    throw std::runtime_error("Cannot open DBSnp database: " + db);
  }

  std::string line;
  std::ostringstream cig;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    std::istringstream  iss(line);
    std::string        chr, pos, rs, ref, alt;
    if (!(iss >> chr >> pos >> rs >> ref >> alt)) {
      logger.log(
        true, true,
        "WARNING: malformed VCF line in ", db, ": '", line, "'. Skipping."
      );
      continue;
    }

    // we only care about indels (ref+alt length > 2)
    if (ref.size() + alt.size() <= 2)
      continue;

    DBSnpSite site(chr, pos, rs, ref, alt, header);
    m_sites.add(site);

    // build a simple chr_pos hash
    cig.str("");
    cig << chr << "_" << site.pos1;
    m_int_hash.insert(m_hasher(cig.str()));
  }

  // finalize our index
  m_sites.CreateTreeMap();
  logger.log(true,true, "Loaded ", m_sites.size(), " indel sites from DBSNP file '", db, "'.")'
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
