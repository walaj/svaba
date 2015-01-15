#ifndef CONTIG_SNOW_H
#define CONTIG_SNOW_H

#include "contigs.h"

/*void Contig::print() const { 
    std::cout << "Name: " << m_name << std::endl;
    for (unsigned i = 0; i < m_read_names.size(); i++)
      std::cout << "Read: " << m_read_names[i] << std::endl;
      }*/

/*bool Contig::writeContig( const std::string read ) { 
    unsigned i = 0; 
    bool found = false;
    while (i < m_read_names.size() && !found) {
       found = m_read_names[i].compare(read) == 0;
       i++;
    }
    return found;
    }*/

/*void Contig::addRead(const std::string read, const unsigned align, const std::string seq) {
    m_reads.push_back(read);
    m_alignments.push_back(align);
    m_read_seqs.push_back(seq);
    }
*/

/*bool Contig::readInContig( const std::string read ) { 

  bool out = m_name_map.find(read) != m_name_map.end();
  return out;
//   unsigned i = 0; 
//   bool found = false;
//   while (i < m_read_names.size() && !found) {
//     found = m_read_names[i].compare(read) == 0;
//     i++;
//   }
//   return found;
}*/

void Contig::addRead(BamTools::BamAlignment read, const int align, bool isInContig) {

    if (isInContig) 
      m_contig_read_count++;

    // add the tag
    //read.AddTag("CN", "Z", m_name);
    //read.AddTag("AL", "Z", std::to_string(align));
    //std::string tmp; 
    //read.GetTag("JW", tmp);

    m_bamalignments.push_back(read);
    //m_name_map.insert(std::pair<std::string, unsigned>(tmp, 0));

}

void Contig::printBamAlignments() const {
  std::string tmp, ts, al; 
    for (std::vector<BamTools::BamAlignment>::const_iterator it = m_bamalignments.begin(); it != m_bamalignments.end(); it++) {
      it->GetTag("JW", tmp);
      it->GetTag("TS", ts);
      it->GetTag("AL", al);
      std::cout << tmp << " Pos: " << it->Position << " TS: " << ts << " AlignmentFlag: " << it->AlignmentFlag << " AL: " << al << std::endl;
    }
}

/*
std::string Contig::printAlignments() const {
    std::stringstream sstream;
    for (unsigned i=0; i < m_reads.size(); i++) {
      sstream << m_reads[i] << "\t" << m_name << "\t" << m_alignments[i] << 
	"\t" << m_read_seqs[i] << std::endl;
    }
    return sstream.str();
}
*/

#endif
