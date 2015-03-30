#include "DiscordantCluster.h"

void DiscordantCluster::addRead(string name) {
  unordered_map<string, bool>::iterator ff = qnames.find(name);
  if (ff == qnames.end())
    qnames.insert(pair<string, bool>(name, true));
  return;
}

DiscordantCluster::DiscordantCluster(ReadVec &this_reads, ReadVec &all_reads) {

  if (this_reads.size() == 0)
    return;
  if (all_reads.size() == 0)
    return;

  // check the orientations, fill the reads
  bool rev = r_is_rstrand(this_reads[0]); //this_reads[0]->IsReverseStrand();
  bool mrev = r_is_mrstrand(this_reads[0]); //this_reads[0]->IsMateReverseStrand();

  id = r_qname(this_reads[0]); //this_reads[0]->Name;

  //cout << "rev " << rev << " mrev " << mrev << endl;
  for (auto& i : this_reads) {

    //cout << "flag " << r_flag(i) << " frev " << r_is_rstrand(i) << " mrev " << r_is_mrstrand(i) << endl;
    assert(rev == r_is_rstrand(i) && mrev == r_is_mrstrand(i)); //i->IsReverseStrand() && mrev == i->IsMateReverseStrand());
    string tmp;
    r_get_Z_tag(i, "SR", tmp);

    assert(id.length());
    r_add_Z_tag(i, "DC", id);

    assert(tmp.length());
    reads[tmp] = i;
    qnames[r_qname(i)] = true;
    if (r_qname(i) < id)
      id = r_qname(i);
    if (tmp.at(0) == 't')
      tcount++;
    else
      ncount++;
  }
  
  addMateReads(all_reads);

  assert(reads.size() == mates.size());
  assert(reads.size() > 0);

  // set the mapq
  reads_mapq = getMeanMapq(false);
  mates_mapq = getMeanMapq(true);
  
  // set the regions
  reg1 = GenomicRegion(-1,500000000,-1); // read region
  for (auto& i : reads) {
    reg1.strand = r_strand(i.second); //(!i.second->IsReverseStrand()) ? '+' : '-';
    reg1.chr = r_id(i.second); //i.second->RefID;
    if (r_pos(i.second) < reg1.pos1)
      reg1.pos1 = r_pos(i.second); //i.second->Position;
    int endpos = r_endpos(i.second);
    if (endpos > reg1.pos2)
      reg1.pos2 = endpos;
  }

  reg2 = GenomicRegion(-1,500000000,-1); // mate region
  for (auto& i : mates) {
    reg2.strand = r_strand(i.second); //(!i.second->IsReverseStrand()) ? '+' : '-';
    reg2.chr = r_id(i.second); //i.second->RefID;
    if (r_pos(i.second) < reg2.pos1)
      reg2.pos1 = r_pos(i.second);
    int endpos = r_endpos(i.second);
    if (endpos > reg2.pos2)
      reg2.pos2 = endpos;
  }

  cluster = toRegionString();

  // add the tags
  //debug
  //for (auto& i : reads)
  //  SnowUtils::SmartAddTag(i.second, "DC", cluster);
  //for (auto& i : mates)
  //  SnowUtils::SmartAddTag(i.second, "DC", cluster);
}

/*
DiscordantCluster::DiscordantCluster(ReadVec &this_reads, ReadVec &all_reads) {

  if (this_reads.size() == 0)
    return;
  if (all_reads.size() == 0)
    return;

  // check the orientations, fill the reads
  bool rev = r_is_rstrand(this_reads[0]); //this_reads[0]->IsReverseStrand();
  bool mrev = r_is_mrstrand(this_reads[0]); //this_reads[0]->IsMateReverseStrand();

  id = r_qname(this_reads[0]); //this_reads[0]->Name;

  for (auto& i : this_reads) {
    assert(rev == i->IsReverseStrand() && mrev == i->IsMateReverseStrand());
    string tmp;
    i->GetTag("SR", tmp);
    reads[tmp] = i;
    qnames[i->Name] = true;
    if (i->Name < id)
      id = i->Name;
    if (tmp.at(0) == 't')
      tcount++;
    else
      ncount++;
  }
  
  addMateReads(all_reads);

  assert(reads.size() == mates.size());
  assert(reads.size() > 0);

  // set the mapq
  reads_mapq = getMeanMapq(false);
  mates_mapq = getMeanMapq(true);
  
  // set the regions
  reg1 = GenomicRegion(-1,500000000,-1); // read region
  for (auto& i : reads) {
    reg1.strand = (!i.second->IsReverseStrand()) ? '+' : '-';
    reg1.chr = i.second->RefID;
    if (i.second->Position < reg1.pos1)
      reg1.pos1 = i.second->Position;
    if (i.second->GetEndPosition() > reg1.pos2)
      reg1.pos2 = i.second->GetEndPosition();
  }

  reg2 = GenomicRegion(-1,500000000,-1); // mate region
  for (auto& i : mates) {
    reg2.strand = (!i.second->IsReverseStrand()) ? '+' : '-';
    reg2.chr = i.second->RefID;
    if (i.second->Position < reg2.pos1)
      reg2.pos1 = i.second->Position;
    if (i.second->GetEndPosition() > reg2.pos2)
      reg2.pos2 = i.second->GetEndPosition();
  }

  cluster = toRegionString();

  // add the tags
  for (auto& i : reads)
    SnowUtils::SmartAddTag(i.second, "DC", cluster);
  for (auto& i : mates)
    SnowUtils::SmartAddTag(i.second, "DC", cluster);
}
*/

void DiscordantCluster::addMateReads(ReadVec &bav) { // BamAlignmentUPVector &bav) {

  for (auto& i : bav) {
    string sr;
    if (qnames.count(r_qname(i))) {
      string tmp;
      r_get_Z_tag(i, "SR", tmp);
      r_add_Z_tag(i, "DC", id);
      if (reads.count(tmp) == 0) // only add if this is a mate read
	mates[tmp] = i;
    }
  }

}

double DiscordantCluster::getMeanMapq(bool mate) const {
  double mean = 0;
  vector<int> tmapq;
  if (mate) {
    for (auto& i : mates)
      tmapq.push_back(r_mapq(i.second));
  } else {
    for (auto& i : reads)
      tmapq.push_back(r_mapq(i.second));
  }

  if (tmapq.size() > 0)
    mean = accumulate(tmapq.begin(), tmapq.end(), 0.0) / tmapq.size();
  return mean;
}

double DiscordantCluster::getMeanMapq() const {
  double mean = 0;
  vector<int> tmapq;
  for (auto& i : mates)
    tmapq.push_back(r_mapq(i.second));
  for (auto& i : reads)
    tmapq.push_back(r_mapq(i.second));

  if (tmapq.size() > 0)
    mean = accumulate(tmapq.begin(), tmapq.end(), 0.0) / tmapq.size();
  return mean;
}

/*
DiscordantCluster::DiscordantCluster(string tcluster) {
  cluster = tcluster;
  
  // parse out the GenomicRegions
  stringstream iss(tcluster);
  string region_string1, region_string2;
  getline(iss, region_string1, '_');
  getline(iss, region_string2, '_');
  reg1 = GenomicRegion(region_string1);
  reg2 = GenomicRegion(region_string2);
}
*/

string DiscordantCluster::toRegionString() const {
  int pos1 = (reg1.strand == '+') ? reg1.pos2 : reg1.pos1;
  int pos2 = (reg2.strand == '+') ? reg2.pos2 : reg2.pos1;
  stringstream ss;
  ss << reg1.chr+1 << ":" << pos1 << "(" << reg1.strand << ")" << "-" << 
    reg2.chr+1 << ":" << pos2 << "(" << reg2.strand << ")";
  return ss.str();
    
}

// define how to print this to stdout
std::ostream& operator<<(std::ostream& out, const DiscordantCluster& dc) {
  
  out << dc.toRegionString() << " Tcount: " << dc.tcount << 
    " Ncount: "  << dc.ncount << " Mean MAPQ: " 
      << dc.reads_mapq << " Mean Mate MAPQ: " << dc.mates_mapq;
  /*for (auto& i : dc.reads) {
    string tmp;
    i.second->GetTag("SR",tmp);
    out << "   " << i.second->RefID << ":" << i.second->Position << (i.second->IsReverseStrand() ? "(-)" : "(+)") << " - " << i.second->MateRefID << ":" << i.second->MatePosition <<  (i.second->IsMateReverseStrand() ? "(-)" : "(+)") << " " << tmp << endl;
    }*/
  return out;
}


// define how to print to file
string DiscordantCluster::toFileString(bool with_read_names /* false */) const { 

  string sep = "\t";

  // add the reads names (currently off)
  string reads_string = "x";
  if (with_read_names) {
    for (auto& i : reads) {
      string tmp;
      r_get_Z_tag(i.second, "SR", tmp);
      //i.second->GetTag("SR", tmp);
      reads_string += tmp + ",";
    }
    
    //debug
    if (reads_string.length() == 0)
      reads_string = "x";
    else
      reads_string = reads_string.substr(0,reads_string.length() - 1);
  }

  stringstream out;
  out << reg1.chr+1 << sep << reg1.pos1 << sep << reg1.strand << sep 
      << reg2.chr+1 << sep << reg2.pos1 << sep << reg2.strand << sep 
      << tcount << sep << ncount << sep << reads_mapq << sep 
      << mates_mapq << sep << reads_string;

  return (out.str());
  /*
  string sep = ",";
  stringstream ss;

  // take the position closest to the break
  int pos1 = (reg1.strand=='+') ? reg1.pos2 : reg1.pos1;
  int pos2 = (reg2.strand=='+') ? reg2.pos2 : reg2.pos1;
  double meanmapq = getMeanMapq();

  //chr1, pos1, str1, chr2, pos2, str2, tum, norm, contig
  ss << reg1.chr+1 << sep << pos1 << sep << reg1.strand << sep <<
    reg2.chr+1 << sep << pos2 << sep << reg2.strand << sep <<
    tcount << sep << ncount << sep << contig << sep << meanmapq;
  return ss.str(); 
  */
}

// define how to sort theses
bool DiscordantCluster::operator<(const DiscordantCluster &b) const {

  if (reg1.chr < b.reg1.chr)
    return true;
  if (reg1.pos1 < b.reg1.pos1)
    return true;
  return false;
  
}
