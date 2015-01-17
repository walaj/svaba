#include "AuxUtils.h"

#include <iostream>
#include <sstream>

using namespace std;

// constructor that takes a line from a BLAT hit
SBlat::SBlat(string input) {
  
  string s_match, s_mismatch, s_repmatch, s_Ns, s_Qgap_count, s_Qgap_bases, 
    s_Tgap_count, s_Tgap_bases, s_strand, query_name, s_query_size, s_query_start,
    s_query_end, hit_name, s_hit_size, s_hit_start, s_hit_end, s_blockCount, s_blockSizes;
  
  istringstream iss(input);
  
  if (!(iss >> s_match >> s_mismatch >> s_repmatch >> s_Ns >> s_Qgap_count >>
	s_Qgap_bases >> s_Tgap_count >> s_Tgap_bases >> s_strand >> query_name >>
	s_query_size >> s_query_start >> s_query_end >> hit_name >> s_hit_size >> 
	s_hit_start >> s_hit_end >> s_blockCount >> s_blockSizes))
    cerr << "Error in making SBlat object" << endl; 
  
  match = stoi(s_match);
  mismatch = stoi(s_mismatch);
  repmatch = stoi(s_repmatch);
  Ns = stoi(s_Ns);
  Qgap_count = stoi(s_Qgap_count);	  
  Qgap_bases = stoi(s_Qgap_bases);
  Tgap_count = stoi(s_Tgap_count);
  Tgap_bases = stoi(s_Tgap_bases);
  query_size = stoi(s_query_size);
  query_start = stoi(s_query_start);
  query_end   = stoi(s_query_end);
  hit_size = stoi(s_hit_size);
  hit_start = stoi(s_hit_start);
  hit_end = stoi(s_hit_end);
  blockCount = stoi(s_blockCount);
  blockSizes = stoi(s_blockSizes);
}

RepeatMasker::RepeatMasker(string input) {
  
  string s_sw, s_perc_div, s_perc_del, s_perc_ins, s_contig, 
    s_query_begin, s_query_end, s_query_left, dum, s_repeat_class, s_repeat_name, 
    s_repeat_begin, s_repeat_end, s_repeat_left, s_id, dum2;
  
  istringstream iss(input);
  if (!(iss >> s_sw >> s_perc_div >> s_perc_del >> s_perc_ins >> s_contig >> s_query_begin 
	>> s_query_end >> s_query_left >> dum >> s_repeat_name >> s_repeat_class >> 
	s_repeat_begin >> s_repeat_end >> s_repeat_left >> s_id)) 
    cerr << "Error in making RepeatMasker object" << endl; 
  
  sw = stoi(s_sw);
  perc_div = stod(s_perc_div);
  perc_del = stod(s_perc_del);
  perc_ins = stod(s_perc_ins);
  query_begin = stoi(s_query_begin);
  query_end = stoi(s_query_end);
  contig = s_contig;
  repeat_class = s_repeat_class;
  repeat_name = s_repeat_name;
  
}
