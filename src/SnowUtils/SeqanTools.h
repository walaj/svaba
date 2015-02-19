#ifndef SEQAN_TOOLS_H
#define SEQAN_TOOLS_H

#include <string>
#include <time.h>
#include <ctime>
#include <vector>

#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace std;
using namespace seqan;

typedef String<Dna5> TSequence;
typedef Align<TSequence,ArrayGaps> TAlign; 

namespace SeqanTools {
  
  double SWalign(TSequence &ref,int32_t &pos, string &rseq, int32_t &score, int cutoff, bool revcomp = false);

}

#endif
