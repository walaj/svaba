#include "DiscordantCluster.h"

void DiscordantCluster::addRead(string name) {
  unordered_map<string, bool>::iterator ff = qnames.find(name);
  if (ff == qnames.end())
    qnames.insert(pair<string, bool>(name, true));
  return;
}

DiscordantCluster::DiscordantCluster(BamAlignmentUPVector &this_reads, BamAlignmentUPVector &all_reads) {

  if (this_reads.size() == 0)
    return;
  if (all_reads.size() == 0)
    return;

  // check the orientations, fill the reads
  bool rev = this_reads[0]->IsReverseStrand();
  bool mrev = this_reads[0]->IsMateReverseStrand();

  id = this_reads[0]->Name;

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

void DiscordantCluster::addMateReads(BamAlignmentUPVector &bav) {

  for (auto& i : bav) {
    if (qnames.count(i->Name) == 1) {
      string tmp;
      i->GetTag("SR", tmp);
      if (reads.count(tmp) == 0)
	mates[tmp] = i;;
    }
  }

}

double DiscordantCluster::getMeanMapq(bool mate) const {
  double mean = 0;
  vector<int> tmapq;
  if (mate) {
    for (auto& i : mates)
      tmapq.push_back(i.second->MapQuality);
  } else {
    for (auto& i : reads)
      tmapq.push_back(i.second->MapQuality);
  }

  if (tmapq.size() > 0)
    mean = accumulate(tmapq.begin(), tmapq.end(), 0.0) / tmapq.size();
  return mean;
}

double DiscordantCluster::getMeanMapq() const {
  double mean = 0;
  vector<int> tmapq;
  for (auto& i : mates)
    tmapq.push_back(i.second->MapQuality);
  for (auto& i : reads)
    tmapq.push_back(i.second->MapQuality);

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
string DiscordantCluster::toFileString() const { 
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
