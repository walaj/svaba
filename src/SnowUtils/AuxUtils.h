#ifndef AUXUTILS_H
#define AUXUTILS_H

#include <string>
#include <vector>
#include <cstdlib>
#include <sstream>

using namespace std;

struct SBlat {

  SBlat(string input);
  
  int match, mismatch, repmatch, Ns, Qgap_count, Qgap_bases, Tgap_count, Tgap_bases;
  char strand;
  string query_name;
  int query_size, query_start, query_end;
  string hit_name;
  int hit_size, hit_start, hit_end;
  int blockCount, blockSizes; 

  // define how these are to be sorted. Sort by biggest match first
  bool operator < (const SBlat& b) const { return (match > b.match); }

};

typedef vector<SBlat> SBlatVec;

struct RepeatMasker {

  RepeatMasker(string input);

  //RepeatMasker(int tsw, string tcontig, string trep_name, string trep_class, 
  //	       double tdiv, double tdel, double tins, int beg, int end) :
  // sw(tsw), contig(tcontig), repeat_name(trep_name), repeat_class(trep_class), perc_div(tdiv), perc_del(tdel), perc_ins(tins), 
  // query_begin(beg), query_end(end) {}

  int sw;
  string contig;
  string repeat_name;
  string repeat_class;
  
  double perc_div;
  double perc_del;
  double perc_ins;
  
  int query_begin;
  int query_end;

  // define how these are to be sorted. Sort by biggest sw first
  bool operator < (const RepeatMasker& rep) const { 
    return (sw > rep.sw); 
  }

  // print it out
  string toString() {
    stringstream ss;
    ss << "SW: " << sw << " Contig: " << contig << " RepeatName: " << repeat_name << 
      " RepeatClass: " << repeat_class << " %Div: " << perc_div << " %Del: " << perc_del << 
      " %Ins: " << perc_ins << " QueryBegin: " << query_begin << " QueryEnd: " << query_end;
    return ss.str();
  }

};

typedef vector<RepeatMasker> RepeatMaskerVec;

#endif
