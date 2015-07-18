#ifndef CONTIG_H
#define CONTIG_H

#include <vector>
#include <string>
#include <unordered_map>
#include "EncodedString.h"

using namespace std;

typedef unordered_map<string, unsigned> StringMap;

class Contig {

 public:

  Contig() {}

  Contig(const string name, const string seq) :
    m_name(name), m_seq(seq) {}

  ~Contig() {}

  size_t getLength() const { return m_seq.length(); }
  
  //void clearReads() {
  //  m_bamalignments.clear();
  //}
  
  // remove unneccesary tags to save space
  /*  void clearFinalTags() {

    vector<BamTools::BamAlignment>::iterator it = m_bamalignments.begin();
    for (; it != m_bamalignments.end(); it++) {
      it->RemoveTag("J2");
      it->RemoveTag("RP");
      it->RemoveTag("HP");
      it->RemoveTag("TS");

      string jw; 
      it->GetTag("JW", jw);
      jw.erase(1,jw.length()-1);
      it->EditTag("JW", "Z", jw);
      
      //SVBamReader::clearFinalTags(it);
      }
      }*/

  string getID() const { return m_name; }

  string getSeq() const { return m_seq.toString(); }

  //void addRead(BamTools::BamAlignment read, const int align, bool isInContig);

  //void addRead(BamTools::BamAlignment read) { m_bamalignments.push_back(read); }

  //vector<BamTools::BamAlignment> getBamAlignments() const { return m_bamalignments; }

  //void printBamAlignments() const;

  //int getReadCount() const {
  //  return m_bamalignments.size();
  //}

  //int getContigReadCount() const {
  //  return m_contig_read_count;
  //}

  /*int getContigTumorReadCount() const {
    int out = 0;
    for (vector<BamTools::BamAlignment>::const_iterator it = m_bamalignments.begin(); it != m_bamalignments.end(); it++) {
      string tmp;
      string al;
      it->GetTag("JW", tmp);
      it->GetTag("AL", al);
      if (tmp.at(0) == 't' && al.compare("-1") != 0)
	out++;
    }
    return out;
    }*/


  bool operator < (const Contig &c) const { 
    return m_seq.length() < c.getLength();
  }

   int tcount = 0;
   int ncount = 0;
 
 private: 
   string m_name;
   //string m_seq;
   DNAEncodedString m_seq;
   //vector<BamTools::BamAlignment> m_bamalignments;
   int m_contig_read_count = 0; // count of reads that actually built contig

};

typedef vector<Contig> ContigVector;

#endif
