#pragma once

namespace SeqLib {
  class BWAAligner;
}
#include "svabaRead.h"

// object that holds functions for re-aligning a discordant read
// and checking if it is a true discordant read
class DiscordantRealigner {

 public:

  DiscordantRealigner() {}

  ~DiscordantRealigner() {}

  // check if a read should be realigned
  //bool ShouldRealign(const SeqLib::BamRecord& r) const;
  bool ShouldRealign(const svabaRead& r) const;

  // realign read
  bool RealignRead(svabaRead& r, const SeqLib::BWAAligner& bwa) const;

  // reassign a read if it has a better mapping
  // r is the read to be reassigned, s source alignment is the new alignment
  //void ReassignRead(SeqLib::BamRecord& r, const SeqLib::BamRecord& s) const;
  //void ReassignRead(svabaRead& r, const svabaRead& s) const;

  static const int MAPS_NOT_NEAR_ORIG = -1;  
  static const int MAPS_NEAR_MATE = -2;
  static const int MAPS_NOWHERE = -3;
  static const int MATE_BAD_DISC = -4;
  static const int REASSIGNED_READ = -5;

 private:
  
  // max number of secondary reads to consider
  int secondary_cap = 20;

  // minimum size to check
  int min_isize_for_discordant_realignment = 1000;

  // how much to search aroudn mate region for read alignent
  int discordant_realign_mate_pad = 100;
  
};
