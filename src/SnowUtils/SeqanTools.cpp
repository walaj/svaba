#include "SeqanTools.h"
#include "SnowUtils.h"

double SeqanTools::SWalign(TSequence &ref,int32_t &pos, string &rseq, int32_t &score, bool revcomp /* false */) {
  
  int match = 4;
  int mismatch = -2;
  int gapopen = -4;
  int gapextend = -2;

  if (revcomp)
    SnowUtils::rcomplement(rseq);
  TSequence read = rseq; 
    
  TAlign align;
  resize(rows(align), 2); 
  assignSource(row(align,0), ref); 
  assignSource(row(align,1), read);
  score = localAlignment(align, Score<int,Simple>(match, mismatch, gapopen, gapextend));
  
  //r.sw_score = score;
  
  unsigned new_pos = clippedBeginPosition(row(align,0));
  //unsigned read_start_pos = clippedBeginPosition(row(align,1));
  //unsigned read_end_pos = clippedEndPosition(row(align, 1));
  pos = new_pos;
  //r.pos = new_pos;
  
  // trim the front of the read if it falls off the front
  //int end = read_end_pos - read_start_pos;
  //if (end > 20) 
  //  r.seq = r.seq.substr(read_start_pos, end);
  
  return score;
}
