#ifdef HAVE_SEQAN_BASIC_H

#include "SeqanTools.h"
#include "SnowUtils.h"

double SeqanTools::SWalign(TSequence &ref,int32_t &pos, string &rseq, int32_t &score, int cutoff, bool revcomp, bool indel /* false */) {

  int match = 4;
  int mismatch = -2;
  int gapopen, gapextend;
  //if (indel) {
  //  gapopen = -100;
  //  gapextend = -100;
  //} else {
    gapopen = -16;
    gapextend = -16;
    //}

  Score<int, Simple> sc(match, mismatch, gapopen, gapextend);

  if (revcomp) 
    SnowUtils::rcomplement(rseq);
  TSequence read = rseq; 

  // try the MeyersBitVector first
  //TAlign malign;
  //resize(rows(malign),2);
  //assignSource(row(malign,0),ref);
  //assignSource(row(malign,1),read);
  //int mscore = globalAlignmentScore(malign, MyersBitVector());
    
  TAlign align;
  resize(rows(align), 2); 
  assignSource(row(align,0), ref); 
  assignSource(row(align,1), read);
  score = localAlignment(align, sc);

  //debug 
  //if (rseq == "GAGCTGGGGGTTAGGTGAAGGAAATTGTGTGTGTGTGTGTGTGTGTGTCTGTGCATGCACATGCGTGTGTGCACGC")
  //  cout << align << " score " << score << endl;

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

  // left end of read hangs off
  if (read_start_pos > 0 && pos == 0 && rseq.length() > 5)
    rseq = rseq.substr(read_start_pos, rseq.length() - read_start_pos);
  // right end of read hangs off
  else if (read_end_pos > clip_end_pos) 
    rseq = rseq.substr(0, read_end_pos);

  //debug 
  //score = 0;

  return score;
}

#endif
