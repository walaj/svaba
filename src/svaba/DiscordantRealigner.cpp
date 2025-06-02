#include "DiscordantRealigner.h"
#include "svabaUtils.h"
#include "SeqLib/BWAAligner.h"

#ifdef QNAME
#define DEBUG(msg, read)				\
  if (read.Qname() == QNAME && (read.AlignmentFlag() == QFLAG || QFLAG == -1)) { std::cerr << (msg) << " read " << r << std::endl; }
#else
#define DEBUG(msg, read)
#endif

//bool DiscordantRealigner::ShouldRealign(const SeqLib::BamRecord& r) const {
bool DiscordantRealigner::ShouldRealign(const svabaRead& r) const {

  //if (r.SR() == "t000_147_H01PEALXX140819:3:1219:25827:26572")
  //  std::cerr << " DD " << r.GetDD() << " fi " << r.FullInsertSize() << " lim " << min_isize_for_discordant_realignment << std::endl;

  // read was already realigned
  if (r.GetDD() != 0)
    return false;
  
  // read can't be discordant anyway, since not both mapped
  if (!r.MappedFlag() || !r.MateMappedFlag())
    return false;
  
  if (r.Interchromosomal())
    return true;

  // get the insert size
  int fi = r.FullInsertSize();

  // if too small, don't realign
  if (fi < min_isize_for_discordant_realignment && !r.Interchromosomal())
    return false;

  return true;
}

//bool DiscordantRealigner::RealignRead(SeqLib::BamRecord& r, const SeqLib::BWAWrapper* bwa) const {
bool DiscordantRealigner::RealignRead(svabaRead& r, const SeqLib::BWAAligner& bwa) const {

  DEBUG("Realigning original read", r);

  SeqLib::BamRecordVector als;
  bwa.alignSequence(r.Seq(), r.Qname(), als, false, 0.60, secondary_cap);

  // no alignments, so label as bad
  if (!als.size()) {
    r.SetDD(MAPS_NOWHERE);
    return false;
  }

  // discordant alignments
  SeqLib::GenomicRegion gr  = r.AsGenomicRegion();
  SeqLib::GenomicRegion grm = r.AsGenomicRegionMate();
  grm.Pad(discordant_realign_mate_pad);
  
  bool has_orig = false;
  bool maps_near_mate = false;
  for (auto& i : als) {

    DEBUG("   Realigned hit (one of) " + std::to_string(als.size()), i);

    // if there is another alignment that overlaps with the original
    // but has a lower MAPQ, take the new mapq
    if (gr.GetOverlap(i.AsGenomicRegion())) {
      has_orig = true;
      r.SetMapQuality(std::min(i.MapQuality(), r.MapQuality()));
    }
    
    DEBUG("   Realigned hit has overlap with mate of N bases: " + std::to_string(r.OverlappingCoverage), i);
    
    // if mate is mapped, and read maps to near mate region wo clips, it's not disc
    if (grm.GetOverlap(i.AsGenomicRegion())) {
      // if not clipped or overlaps same covered region as original, 
      // then it has a secondary mapping
      if (r.NumClip() < 20 || r.OverlappingCoverage(i) >= 20) 
	maps_near_mate = true;
      DEBUG("   Found realignment hit near mate", i);
    }
  } // end als loop
  
  // add tags
  if (!has_orig) 
    r.SetDD(MAPS_NOT_NEAR_ORIG);
  else if (maps_near_mate) 
    r.SetDD(MAPS_NEAR_MATE);
  else
    r.SetDD(als.size());
  
  return r.GetDD() < 0; // was it realigned?
}

/*void DiscordantRealigner::ReassignRead(svabaRead& r, const svabaRead& s) const {

  r.Reassign(s);
  r.SetDD(REASSIGNED_READ);
  return;

  }*/
