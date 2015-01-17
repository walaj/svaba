#include "BreakPoint.h"

// send breakpoint to a string
string BreakPoint::toString() const {
    stringstream out;
    
    out << cname << " " << refID1 + 1 << ":" << pos1 << "(" << strand1 << ") to " <<
       refID2 + 1 << ":" << pos2 << "(" << strand2 << ")" << " Span: " << span << " MAPQ: " << 
       mapq1 << "/" << mapq2 << " homology: " << 
       homology << " insertion: " << insertion << " NumDups: " << num_dups << " Nsplit: " << 
       nsplit << " Tsplit: " << tsplit; // << " isBest: " << isBest;
    return out.str();
}

// test whether are same break 
bool BreakPoint::sameBreak(BreakPoint &bp) const {
    return bp.refID1 == refID1 && bp.pos1 == pos1 && bp.refID2 == refID2 && bp.pos2 == pos2;
}

// order them
void BreakPoint::order() {

    if (refID1 < refID2)
      return;
    if (refID1 == refID2 && pos1 < pos2)
      return;

    unsigned tmppos1 = pos1;  pos1 = pos2;  pos2 = tmppos1;
    unsigned tmprefID1 = refID1;  refID1 = refID2;  refID2 = tmprefID1;
    unsigned tmpcpos1 = cpos1;  cpos1 = cpos2;  cpos2 = tmpcpos1;
    unsigned tmpmapq1 = mapq1;  mapq1 = mapq2;  mapq2 = tmpmapq1;
    unsigned tmpstrand1 = strand1;  strand1 = strand2;  strand2 = tmpstrand1;
    unsigned tmptsplit1 = tsplit1;  tsplit1 = tsplit2;  tsplit2 = tmptsplit1;
    unsigned tmpnsplit1 = nsplit1;  nsplit1 = nsplit2;  nsplit2 = tmpnsplit1;

}

// return whether a bp is good to move on
bool BreakPoint::isGoodSomatic(int mapq, size_t tsplit_cutoff, size_t nsplit_cutoff) const {

  mapq = 0;

  // set whether there is enoguh somatic read support
  bool issom = min(tsplit1, tsplit2) >= tsplit_cutoff; // split
  issom = issom || (dc.tcount >= tsplit_cutoff && dc.getMeanMapq() >= 0); //discordant

  // make sure that there are no germline reads
  bool nogerm = min(nsplit1, nsplit2) <= nsplit_cutoff; //split
  bool good_enough = static_cast<double>(min(tsplit1, tsplit2)) / static_cast<double>(max(nsplit1, nsplit2)) >= 5 && max(nsplit1, nsplit2) < 2;
  nogerm = nogerm && (dc.ncount == 0 || (dc.ncount == 1 && dc.tcount >= 10) );
  nogerm = nogerm || good_enough;

  // if it was aligned to the right spot, forget mapq, it's OK
  //int tmp_mapq1 = local1 ? 100 : mapq1; 
  //int tmp_mapq2 = local2 ? 100 : mapq2;
  
  int tmp_mapq1 = mapq1;
  int tmp_mapq2 = mapq2;

  // if the matching length is good enough, keep it
  tmp_mapq1 = matchlen1 > 250 ? 100 : tmp_mapq1;
  tmp_mapq2 = matchlen2 > 250 ? 100 : tmp_mapq2;

  // make sure the mapq quality is high enough
  bool qual = (min(tmp_mapq1, tmp_mapq2) >= mapq) || dc.tcount >= 1;

  // if it has split and discordant support, forget about mapq if it's high enoguh
  qual = qual || (min(tsplit1, tsplit2) >= 1 && dc.tcount >= 3 && max(mapq1, mapq2) >= 30);

  //bool enough_tum_reads = static_cast<double>(tall)/static_cast<double>(nall) >= 2 && max(tsplit1, tsplit2) > 0 && max(nsplit1, nsplit2) == 0 && tall >= 5;
  
  bool passall = issom && nogerm && qual; 

  return passall;
}

// return whether a bp is good to move on
bool BreakPoint::isGoodGermline(int mapq, size_t allsplit) const {

  mapq = 0;

  // set whether there is enoguh somatic read support
  bool isger = min(nsplit1, nsplit2) >= allsplit; // split
  int total_disc_count = dc.ncount + dc.tcount;
  isger = isger || (dc.ncount >= 1 && total_disc_count >= 12 && dc.getMeanMapq() > 20); //discordant

  // is this close enoguht to be somatic?
  bool good_enough = static_cast<double>(min(tsplit1, tsplit2)) / static_cast<double>(max(nsplit1, nsplit2)) >= 5 && max(nsplit1, nsplit2) < 2;

  // be extra extra careful for germline inter-chrom
  bool inter_ok = (refID1 == refID2) || 
    ( (mapq1 == 60 && mapq2 == 60) && (nm1 < 5 && nm2 < 5)  && (matchlen1 > 70 && matchlen2 > 70)  && total_disc_count >= 6 && dc.getMeanMapq() > 20);

  int tmp_mapq1 = mapq1;
  int tmp_mapq2 = mapq2;

  // if the matching length is good enough, keep it
  tmp_mapq1 = matchlen1 > 250 ? 100 : tmp_mapq1;
  tmp_mapq2 = matchlen2 > 250 ? 100 : tmp_mapq2;

  // make sure the mapq quality is high enough
  bool qual = min(tmp_mapq1, tmp_mapq2) >= mapq || dc.ncount >= 1;

  // if it has split and discordant support, forget about mapq, it's OK
  qual = qual || (min(nsplit1+tsplit, nsplit2+tsplit) >= 2 && total_disc_count >= 3);

  //bool enough_tum_reads = static_cast<double>(tall)/static_cast<double>(nall) >= 2 && max(tsplit1, tsplit2) > 0 && max(nsplit1, nsplit2) == 0 && tall >= 5;
  
  bool passall = isger && qual && inter_ok && !good_enough; 

  return passall;

}

// print to file
void BreakPoint::printToFile(ofstream &of, const BamAlignmentVector &bamreads) {
  string sep = ",";
  
  int discordant_tum = dc.tcount;
  int discordant_norm = dc.ncount;
  string contig_name = cname;
  if (discovar) {
    discordant_tum = disco_tum;
    discordant_norm = disco_norm;
    contig_name = "discovar_" + cname;
  } 
  
  // set the evidence string
  bool isdisc = (dc.tcount + dc.ncount) != 0;
  //bool issplit = (tsplit + nsplit) != 0;
  bool issplit = dc.contig != ""; 
  if (discovar)
    evidence = "DSCVR";
  else if ( isdisc && issplit )
    evidence = "ASDIS";
  else if ( isdisc )
    evidence = "DSCRD";
  else if (num_align == 2)
    evidence = "ASSMB";
  else if (num_align > 2)
    evidence = "COMPL";
  
  confidence = "";
  // low split, high span, no discordant. Suspicous
  if ((tsplit + nsplit) < 6 && span > 1500 && !hasDiscordant()) 
    confidence = "NODISC";
  // assembly break, low mapq, no discordant
  else if ((tsplit + nsplit) != 0 && mapq1 != 60 && mapq2 != 60 && !hasDiscordant())
    confidence = "LOWMAPQ";
  // discordant only
  else if ((tsplit + nsplit) == 0 && hasDiscordant() && (dc.tcount + dc.ncount) < 6)
    confidence = "WEAKDISC";
  // assembly only, low evidence
  else if ( (tsplit + nsplit) <= 3 && !hasDiscordant() && span > 1500)
    confidence = "WEAKASSEMBLY";
  // disc is 6+, has assembly + disc, assembly is 3+ with 60/60 mapq
  else 
    confidence = "PASS";

  assert(span > -2);

  string supporting_reads = "";
  unordered_map<string, bool> reads;
  //add the discordant reads
  for (unordered_map<string, bool>::const_iterator it = dc.qnames.begin(); it != dc.qnames.end(); it++) {
    if (reads.find(it->first) == reads.end())
      if (it->first.size() > 0)
	reads.insert(pair<string, bool>(it->first.substr(1, it->first.size()-1), true));
  }
  //add the contig reads
  //BamAlignmentVector bamreads = (*contigs)[cname].m_bamreads;
  for (BamAlignmentVector::const_iterator it = bamreads.begin(); it != bamreads.end(); it++) {
    if (reads.find(it->Name) == reads.end())
      reads.insert(pair<string, bool>(it->Name, true));
  }
   
  // print reads to a string
  for (unordered_map<string, bool>::const_iterator it = reads.begin(); it != reads.end(); it++) 
    supporting_reads = supporting_reads + "_" + it->first;
  if (supporting_reads.size() > 0)
    supporting_reads = supporting_reads.substr(1, supporting_reads.size() - 1); // remove first _

  of << evidence << sep << refID1+1 << sep << pos1 << sep << strand1 << sep 
     << refID2+1 << sep << pos2 << sep << strand2 << sep
     << mapq1 <<  sep << mapq2 << sep 
     << nsplit << sep << tsplit << sep
     << discordant_norm << sep << discordant_tum << sep
    //<< nall << sep << tall << sep 
     << homology << sep << insertion << sep << contig_name << sep
     << span << sep << num_dups << sep << num_align << sep 
     << confidence << sep << supporting_reads << endl;
}

// make a breakpoint from a discordant cluster 
BreakPoint::BreakPoint(DiscordantCluster tdc) {

  dc = tdc;
  pos1 = (tdc.reg1.strand == '+') ? tdc.reg1.pos2 : tdc.reg1.pos1;
  pos2 = (tdc.reg2.strand == '+') ? tdc.reg2.pos2 : tdc.reg2.pos1;
  refID1 = tdc.reg1.chr;
  refID2 = tdc.reg2.chr;
  cname = tdc.cluster;
  strand1 = tdc.reg1.strand;
  strand2 = tdc.reg2.strand;

  mapq1 = tdc.getMeanMapq();
  mapq2 = mapq1;

  local1 = true;
  local2 = true; // by definition, for discordant alignment came from alignment location...
  
  if (refID1 != refID2)
    span = -1;
  else
    span = abs((int)pos1 - (int)pos2);
  
}

string BreakPoint::BreakPointHeader() {
  string sep = ",";
  stringstream header;
  header << "evidence" << sep << "chr1" << sep << "pos1" << sep << "strand1" << sep
         << "chr2" << sep << "pos2" << sep << "strand2" << sep
         << "mapq1" << sep << "mapq2" << sep 
	 << "nsplit" << sep << "tsplit" << sep 
	 << "ndisc" << sep << "tdisc" << sep 
    //<< "tumcount" << sep << "numcount" << sep 
	 << "homology" << sep << "insertion" << sep
         << "cname" << sep << "span" << sep << "num.dups" << sep << "num.parts" 
	 << sep << "confidence" << sep << "supporting.reads";
  return header.str();
}

bool BreakPoint::hasDiscordant() const {
  return !dc.reg1.isEmpty();
}
