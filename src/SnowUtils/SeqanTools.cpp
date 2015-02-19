#include "SeqanTools.h"
#include "SnowUtils.h"

double SeqanTools::SWalign(TSequence &ref,int32_t &pos, string &rseq, int32_t &score, int cutoff, bool revcomp /* false */) {
  
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

  //debug
  //cout << "score for "  << rseq << " is " << score << " revcomp " << revcomp << endl;

  if (score < cutoff)
    return false;

  //r.sw_score = score;
  //  cout << align << " " << read << endl;
  unsigned new_pos = clippedBeginPosition(row(align,0));
  unsigned clip_end_pos = clippedEndPosition(row(align,0));
  unsigned read_start_pos = clippedBeginPosition(row(align,1));
  unsigned read_end_pos = clippedEndPosition(row(align, 1));
  //cout << align << " " << read << " read_start " << read_start_pos << 
  //  " read_end_pos " << read_end_pos << " contig_start_pos " << new_pos << 
  //     " clip_end_pos " << clip_end_pos << endl;

    pos = new_pos;

  // trim the front of the read if it falls off the front
  int align_span = read_end_pos - read_start_pos + 1;
  //if (end > 20) 
  if (align_span != rseq.length()) {
    rseq = rseq.substr(read_start_pos, align_span);
    pos = new_pos -read_start_pos;
    //assert(pos >= 0);
  }
  
  return score;
}
