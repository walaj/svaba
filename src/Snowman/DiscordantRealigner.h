#ifndef DISCORDANT_REALIGNER_H__
#define DISCORDANT_REALIGNER_H__

#include "SeqLib/BamRecord.h"
#include "SeqLib/BWAWrapper.h"

// object that holds functions for re-aligning a discordant read
// and checking if it is a true discordant read
class DiscordantRealigner {

 public:

  DiscordantRealigner() {}

  ~DiscordantRealigner() {}

  // check if a read should be realigned
  bool ShouldRealign(const SeqLib::BamRecord& r) const;

  // realign read
  bool RealignRead(SeqLib::BamRecord& r, const SeqLib::BWAWrapper* bwa) const;

  // reassign a read if it has a better mapping
  // r is the original read to be reassigned, s is the new alignment
  void ReassignRead(SeqLib::BamRecord& r, SeqLib::BamRecord& s) const;
  

 private:
  
  // max number of secondary reads to consider
  int secondary_cap = 20;

  // minimum size to check
  int min_isize_for_discordant_realignment = 1000;

  // how much to search aroudn mate region for read alignent
  int discordant_realign_mate_pad = 100;
  
};


#endif
