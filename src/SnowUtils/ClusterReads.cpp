#include "ClusterReads.h"

// perform clustering of discordant reads, with supplied orientations
// orientation is one of FF, FR, RR, RF, AA (for any)
void ClusterReads::clusterReads(BamAlignmentVector &bav, GenomicRegionVector &grv, RMap &rmap, string orientation, int isize, int cluster_buffer) {

  if (bav.size() < 3)
    return;

  char anchor_strand  = (orientation.at(0) == 'F') ? '+' : '-';
  char partner_strand = (orientation.at(1) == 'F') ? '+' : '-';
  if (orientation == "BB") {
    anchor_strand  = '*';
    partner_strand = '*';
  }

  int dref = -1, dpos1 = -1, dpos2 = -1;
  int dtcount = 0, dncount = 0; 

  BamAlignmentVector::const_iterator it = bav.begin();

  // find the initial read with correct orientation
  for (; it != bav.end(); it++) {

    bool isdisc = validDiscordant(*it, isize, orientation);
    if (isdisc) {
      dref  = it->MateRefID;
      dpos1 = it->MatePosition;
      dpos2 = it->MatePosition;
      if (SVBamReader::IsTumorRead(*it))
	dtcount = 1;
      else
	dncount = 1;
      break;
    }
  }

  // no clusters found
  if (dref < 0)
    return;

  // start the first cluster
  grv.push_back(GenomicRegion(dref, dpos1, dpos2)); 

  // keep track of the anchor breakpoint
  GenomicRegion anc(it->RefID, it->Position, it->Position);
  anc.strand = anchor_strand;
  GenomicRegion par(dref, it->MatePosition, it->MatePosition);
  par.strand = partner_strand; 

  stringstream clusttag;
  for (it = it + 1; it != bav.end(); it++) {

    bool isdisc = validDiscordant(*it, isize, orientation);

    if (isdisc) {

      int diff = it->MatePosition - dpos2;
      bool different_chr = dref != it->MateRefID;
      assert(diff >= 0 || different_chr); // make sure ordering of mate position is true

      if ( (diff >= cluster_buffer) || different_chr) { // set new cluster

	dref = it->MateRefID;
	
	// finish off the old cluster
	finalizeCluster(grv, rmap, dpos2, anc, par, dtcount, dncount);
	
        // start the new one
	GenomicRegion gr(it->MateRefID, it->MatePosition); //leaves pos2 undefined
	gr.strand = anchor_strand;
	gr.mapq.push_back(it->MapQuality);
	grv.push_back(gr);
	dtcount = 0;
	dncount = 0;
	anc.chr = it->RefID;
	anc.pos1 = it->Position;
	anc.pos2 = it->Position;
	par.chr = it->MateRefID;
	par.pos1 = it->MatePosition;
	par.pos2 = it->MatePosition;
      }
   
      if (SVBamReader::IsTumorRead(*it))
	dtcount++;
      else
	dncount++;
      
      // add the read name to the GenomicRegion
      grv.back().rname.insert(pair<string, size_t>(it->Name, 0));
      grv.back().mapq.push_back(it->MapQuality);

      dpos2 = it->MatePosition;
      
      // update the anchor and partner breakpoint position
      //int anc_tip_pos = it->IsReverseStrand() ? it->Position : it->GetEndPosition();
      //int par_tip_pos = it->IsMateReverseStrand() ? it->MatePosition : it->GetEndPosition();
      
      anc.pos1 = min(anc.pos1, it->Position);
      par.pos1 = min(par.pos1, it->MatePosition);
      anc.pos2 = max(anc.pos2, it->Position);
      par.pos2 = max(par.pos2, it->MatePosition);

      // discordant clusters shouldn't grow too big
      int lim = 20000;
      if (abs(anc.width()) >= lim || abs(par.width()) >= lim)
	cerr << anc << " " << par << endl;
      //assert(abs(anc.pos1 - anc.pos2) < lim);
      //assert(abs(par.pos1 - par.pos2) < lim);
    }
  }
    
  // finish the last one
  finalizeCluster(grv, rmap, dpos2, anc, par, dtcount, dncount);

  return;
}

void ClusterReads::finalizeCluster(GenomicRegionVector &grv, RMap &rmap, int pos, GenomicRegion anc, GenomicRegion par, int dtcount, int dncount) {

  if (grv.size() == 0)
    return;

  // if cluster has only 1-2 reads, remove it
  if ( (dtcount + dncount) < 3) {
    grv.pop_back();
    return;
  }

  // back sure the mapq is sufficient
  //double mean_mapq = accumulate(grv.back().mapq.begin(), grv.back().mapq.end(), 0.0) / grv.back().mapq.size();
  //if (mean_mapq < 5) {
    //cerr << "...discarding due to low mapq " << mean_mapq << endl;
  //  grv.pop_back();
  //  return;
  // }
    
  // regions that are too big are probably false positives, skip
  if (abs(grv.back().width()) > 3000) {
    grv.pop_back();
    return;
  }

  // finish the last one
  grv.back().pos2 = pos;
  grv.back().tcount = dtcount;					       
  grv.back().ncount = dncount;					       

  // update the string
  stringstream clusttag;
  if (anc < par)
    clusttag << anc << "_" << par;
  else
    clusttag << par << "_" << anc;     
  grv.back().cluster = clusttag.str();

  // write the read map
  string ctag = clusttag.str();
  GMap nams = grv.back().rname;
  for (GMap::iterator jt = nams.begin(); jt != nams.end(); jt++)
    rmap.insert(pair<string, string>(jt->first, ctag));
  grv.back().rname.clear();

}

// return whether a given read should even be considered to be discordant
bool ClusterReads::validDiscordant(const BamAlignment &it, int isize, string orientation) {

  bool ancrev = (orientation.at(0) == 'R');
  bool parrev = (orientation.at(1) == 'R');
  bool any_orientation = orientation.at(0) == 'B';

  bool isNotR2    = !it.HasTag("IR");
  bool mapped     = it.IsMapped() && it.IsMateMapped();
  bool oriented   = (it.IsReverseStrand() == ancrev && it.IsMateReverseStrand() == parrev) || (any_orientation);
  bool discordant = it.InsertSize >= isize || (it.RefID != it.MateRefID);
  return isNotR2 && mapped && oriented && discordant && oriented;

}
