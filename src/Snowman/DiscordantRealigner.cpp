#include "DiscordantRealigner.h"

#ifdef QNAME
#define DEBUG(msg, read)				\
  if (read.Qname() == QNAME && (read.AlignmentFlag() == QFLAG || QFLAG == -1)) { std::cerr << (msg) << " read " << r << std::endl; }
#else
#define DEBUG(msg, read)
#endif

bool DiscordantRealigner::ShouldRealign(const SeqLib::BamRecord& r) const {

  // read was already realigned
  if (r.GetIntTag("DD") != 0)
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

bool DiscordantRealigner::RealignRead(SeqLib::BamRecord& r, const SeqLib::BWAWrapper* bwa) const {

  DEBUG("Realigning original read", r);

  SeqLib::BamRecordVector als;
  bwa->AlignSequence(r.Sequence(), r.Qname(), als, false, 0.60, secondary_cap);
  
  // no alignments, so label as bad
  if (!als.size()) {
    r.AddIntTag("DD", -3);
    return false;
  }

  // discordant alignments
  SeqLib::GenomicRegion gr  = r.asGenomicRegion();
  SeqLib::GenomicRegion grm = r.asGenomicRegionMate();
  grm.Pad(discordant_realign_mate_pad);
  
  bool has_orig = false;
  bool maps_near_mate = false;
  for (auto& i : als) {

    DEBUG("   Realigned hit (one of) " + std::to_string(als.size()), i);

    // if the realignment has a better mapq at differnt location
    // then take that, and say that it is too uncertain to be a 
    // reliable alignment
    if (i.MapQuality() > r.MapQuality() && gr.GetOverlap(i.asGenomicRegion())) {
      ReassignRead(r, i);
      continue;
    }
    
    // if there is another alignment that overlaps with the original
    // but has a lower MAPQ, take the new mapq
    if (gr.GetOverlap(i.asGenomicRegion())) {
      has_orig = true;
      r.SetMapQuality(std::min(i.MapQuality(), r.MapQuality()));
    }
    
    DEBUG("   Realigned hit has overlap with mate of N bases: " + std::to_string(r.OverlappingCoverage), i);
    
    // if mate is mapped, and read maps to near mate region wo clips, it's not disc
    if (grm.GetOverlap(i.asGenomicRegion())) {
      // if not clipped or overlaps same covered region as original, 
      // then it has a secondary mapping
      if (r.NumClip() < 20 || r.OverlappingCoverage(i) >= 20) 
	maps_near_mate = true;
      DEBUG("   Found realignment hit near mate", i);
    }
  } // end als loop
  
  // add tags
  if (!has_orig) 
    r.AddIntTag("DD", -1);
  else if (maps_near_mate) 
    r.AddIntTag("DD", -2);
  else
    r.AddIntTag("DD", als.size());
  
  if (r.GetIntTag("DD") < 0)
    return true; // was realigned
  else 
    return false;
  
}

void DiscordantRealigner::ReassignRead(SeqLib::BamRecord& r, SeqLib::BamRecord& s) const {

  s.SetChrIDMate(r.MateChrID());
  s.SetPositionMate(r.MatePosition());
  s.SetPairMappedFlag();

  if (r.MateReverseFlag())
    s.SetMateReverseFlag();

  s.AddZTag("SR", r.GetZTag("SR"));
  r = s; // reassign the read
  s.AddIntTag("DD", -5); // read is re-assigned, so too worrisome for discordant
  return;

}
