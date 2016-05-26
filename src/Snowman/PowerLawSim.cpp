#include "PowerLawSim.h"

#include <cmath>
#include <string>
#include <cassert>
#include <algorithm>
#include <unordered_set>

#define MAX_DUP_DEL_INV_SIZE 2000
#define MAX_RAR_SIZE 35e6
#define EVENT_BUFFER 4000

void PowerLawSim(faidx_t* findex, int num_breaks, double power_law, SnowTools::GRC& grc, std::ofstream& outfasta, std::ofstream& events) {

  double frac_inter = 0.05;

  bool verbose = false;

  // generate the random numbers
  std::vector<int> rpower = drawFromPower(1, MAX_RAR_SIZE, power_law, num_breaks);
  // generate again, but for small indels
  std::vector<int> small_indels2 = drawFromPower(1, 10, -3, num_breaks / 4);
  // generate inter-chr
  size_t num_intra = rpower.size() + small_indels2.size();
  std::vector<int> inter_chr((double)num_intra / (1-frac_inter) - num_intra, -1);
  rpower.insert(rpower.end(), small_indels2.begin(), small_indels2.end());
  rpower.insert(rpower.end(), inter_chr.begin(), inter_chr.end());
  std::random_shuffle(rpower.begin(), rpower.end());

  // check
  for (auto& i : rpower)
    assert(i <= MAX_RAR_SIZE);

  const char TCGA[5] = "TCGA";

  // get random inserts
  std::vector<std::string> inserts;
  std::vector<int> rsize = drawFromPower(1, 20, -2, num_breaks);
  for (int i = 0; i < num_breaks; ++i) {
    if (rand() % 2 == 0) {
      inserts.push_back(std::string());
      continue;
    }
    std::string is(rsize[i], 'N');
    for (size_t j = 0; j < is.length(); ++j)
      is[j] = TCGA[rand() % 3];      
    inserts.push_back(is);
  }

  SVEvent e;
  std::string ename, etype;
  for (int i = 0; i < num_breaks; ++i) {

    e.number = i;
    if (i % 500 == 0)
      std::cerr << "working on " << i << " of " << num_breaks << std::endl;
    
    int rval = rand() % 3;
    std::string frag, outstring;

    if (rpower[i] == -1) {

      etype = "INT";
      ename = "tumor-" + std::to_string(i) + "-" + etype;
      if (verbose) std::cerr << "...generating inter-chr " << std::endl;
      genRandomSequence(frag, e.reg1, EVENT_BUFFER, findex, grc);
      std::string frag_ic;
      do { 
	genRandomSequence(frag_ic, e.reg2, EVENT_BUFFER, findex, grc);
      } while (e.reg1.chr == e.reg2.chr);

      outstring = frag + inserts[i] + frag_ic;

      events << SnowTools::CHR_NAME[e.reg1.chr] << "\t" << e.reg1.pos2 << "\t" << SnowTools::CHR_NAME[e.reg2.chr] << "\t" << e.reg2.pos1 << "\t+\t-\t" << "N\t-1\tINT\t" << (inserts[i].empty() ? "N" : inserts[i]) << "\t" << ename << std::endl;

    // deletion
    } else if (rpower[i] <= 50 && rand() % 2) {
      etype = "del";
      ename = "tumor-" + std::to_string(i) + "-" + etype;
      if (verbose) std::cerr << "...generating deletion of length " << rpower[i] << std::endl;      
      genRandomSequence(frag, e.reg1, EVENT_BUFFER * 2 + rpower[i], findex, grc);
      try {outstring = frag.substr(0, EVENT_BUFFER) + frag.substr(EVENT_BUFFER + rpower[i], frag.length() - rpower[i] - EVENT_BUFFER);} catch (...) { std::cerr << " len " << (frag.length() - rpower[i] - EVENT_BUFFER)  << std::endl; } 
      assert(outstring.length() < EVENT_BUFFER * 10);
      events << SnowTools::CHR_NAME[e.reg1.chr] << "\t" << (e.reg1.pos1+EVENT_BUFFER) << "\t" << SnowTools::CHR_NAME[e.reg1.chr] << "\t" << (e.reg1.pos1 + EVENT_BUFFER + rpower[i]) << "\t+\t-\tN\t" << rpower[i] << "\tdel\tN\t" << ename << std::endl;
    // insertion
    } else if (rpower[i] <= 50) {
      etype = "ins";
      ename = "tumor-" + std::to_string(i) + "-" + etype;
      if (verbose) std::cerr << "...generating insertion of length " << rpower[i] << std::endl;      
      genRandomSequence(frag, e.reg1, EVENT_BUFFER * 2, findex, grc);
      std::string ins(rpower[i], 'N');
      for (int y = 0; y < rpower[i]; ++y) 
	ins[y] = TCGA[rand() % 3];
      outstring = frag.substr(0, EVENT_BUFFER) + ins + frag.substr(EVENT_BUFFER, frag.length() - EVENT_BUFFER);
      assert(outstring.length() < EVENT_BUFFER * 10);
      events << SnowTools::CHR_NAME[e.reg1.chr] << "\t" << (e.reg1.pos1+EVENT_BUFFER) << "\t" << SnowTools::CHR_NAME[e.reg1.chr] << "\t" << (e.reg1.pos1 + EVENT_BUFFER + 1) << "\t+\t-\t" << ins << "\t" << rpower[i] << "\tins\tN\t" << ename << std::endl;
    // TANDEM DUPLICATION
    } else if (rpower[i] < MAX_DUP_DEL_INV_SIZE && rval == 0) {
      etype = "DUP";
      ename = "tumor-" + std::to_string(i) + "-" + etype;
      if (verbose) std::cerr << "...generating DUP of length " << rpower[i] << std::endl;      
      genRandomSequence(frag, e.reg1, EVENT_BUFFER * 2 + rpower[i], findex, grc);
      outstring = frag.substr(0, EVENT_BUFFER + rpower[i]) + inserts[i] + frag.substr(EVENT_BUFFER, frag.length() - EVENT_BUFFER); //rpower[i]) + frag.substr(EVENT_BUFFER, frag.length() - EVENT_BUFFER);
      assert(outstring.length() < EVENT_BUFFER * 10);
      events << SnowTools::CHR_NAME[e.reg1.chr] << "\t" << (e.reg1.pos1+EVENT_BUFFER+rpower[i]) << "\t" << SnowTools::CHR_NAME[e.reg1.chr] << "\t" << (e.reg1.pos1 + EVENT_BUFFER) << "\t+\t-\t" << "N" << "\t" << rpower[i] << "\tDUP\t" << (inserts[i].empty() ? "N" : inserts[i]) << "\t" << ename << std::endl;
    // DELETION
    } else if (rpower[i] < MAX_DUP_DEL_INV_SIZE && rval == 1) {
      etype = "DEL";
      ename = "tumor-" + std::to_string(i) + "-" + etype;
      if (verbose) std::cerr << "...generating DEL of length " << rpower[i] << std::endl;      
      genRandomSequence(frag, e.reg1, rpower[i] + EVENT_BUFFER * 2, findex, grc);
      outstring = frag.substr(0, EVENT_BUFFER) + inserts[i] + frag.substr(EVENT_BUFFER + rpower[i], frag.length() - EVENT_BUFFER - rpower[i]);
      assert(outstring.length() < EVENT_BUFFER * 10);
      events << SnowTools::CHR_NAME[e.reg1.chr] << "\t" << (e.reg1.pos1+EVENT_BUFFER) << "\t" << SnowTools::CHR_NAME[e.reg1.chr] << "\t" << (e.reg1.pos1 + EVENT_BUFFER + rpower[i]) << "\t+\t-\tN\t" << rpower[i] << "\tDEL\t" << (inserts[i].empty() ? "N" : inserts[i]) << "\t" << ename << std::endl;
    // INV
    } else if (rpower[i] < MAX_DUP_DEL_INV_SIZE && rval == 2 && false) {
      etype = "INV";
      ename = "tumor-" + std::to_string(i) + "-" + etype;
      if (verbose) std::cerr << "...generating INV of length " << rpower[i] << std::endl;      
      genRandomSequence(frag, e.reg1, rpower[i] + EVENT_BUFFER * 2, findex, grc);
      std::string inv_frag = frag.substr(EVENT_BUFFER, rpower[i]);
      std::reverse(inv_frag.begin(), inv_frag.end());
      outstring = frag.substr(0, EVENT_BUFFER) + inv_frag + frag.substr(EVENT_BUFFER, frag.length() - EVENT_BUFFER - rpower[i]);
      assert(outstring.length() < EVENT_BUFFER * 10);
      events << SnowTools::CHR_NAME[e.reg1.chr] << "\t" << (e.reg1.pos1+EVENT_BUFFER) << "\t" << SnowTools::CHR_NAME[e.reg1.chr] << "\t" << (e.reg1.pos1 + EVENT_BUFFER + rpower[i]) << "\t+\t-\tN\t" << rpower[i] << "\tINV\t" << (inserts[i].empty() ? "N" : inserts[i]) << "\t" << ename << std::endl;
    } else  {
      etype = "RAR";
      ename = "tumor-" + std::to_string(i) + "-" + etype;
      if (verbose) std::cerr << "...generating RAR of length " << rpower[i] << std::endl;      
      genRandomSequence(frag, e.reg1, EVENT_BUFFER, findex, grc);
      bool leftfit  = e.reg1.pos1 - rpower[i] > 1e6;
      bool rightfit = e.reg1.pos2 + rpower[i] < SnowTools::CHR_LEN[e.reg1.chr];
      int rrr = -1;
      if (leftfit && rightfit)
	rrr = rand() % 2;
      else if (leftfit)
	rrr = 0;
      else if (rightfit)
	rrr = 1;
      
      // rearrange to left
      std::string seq2;
      if (rrr == 0) {

	if (verbose) 
	  std::cerr << "...generating LEFT RAR of length " << rpower[i] << std::endl;      

	SnowTools::GenomicRegion gr(e.reg1.chr, e.reg1.pos1 - rpower[i] - EVENT_BUFFER, e.reg1.pos1 - rpower[i]);
	std::string chrstring = SnowTools::CHR_NAME[e.reg1.chr];
	int len;
	char * seq = faidx_fetch_seq(findex, const_cast<char*>(chrstring.c_str()), gr.pos1, gr.pos2 - 1, &len);
	if (!seq)
	  break;
	else
	  seq2 = std::string(seq);

        outstring = seq2 + inserts[i] + frag;
	assert(outstring.length() < EVENT_BUFFER * 10);
	events << SnowTools::CHR_NAME[gr.chr] << "\t" << gr.pos2 << "\t" << SnowTools::CHR_NAME[e.reg1.chr] << "\t" << e.reg1.pos1 << "\t+\t-\tN\t" << rpower[i] << "\tRAR\t" << (inserts[i].empty() ? "N" : inserts[i]) << "\t" << ename << std::endl;
	//std::cerr << rpower[i] << " span " << (std::abs(gr.pos2 - e.reg1.pos1) - rpower[i]) << std::endl;
	assert(rpower[i] == (e.reg1.pos1 - gr.pos2));


      // rearrange to right
      } else if (rrr == 1) {

	if (verbose) 
	  std::cerr << "...generating RIGHT RAR of length " << rpower[i] << std::endl;      
	SnowTools::GenomicRegion gr(e.reg1.chr, e.reg1.pos2 + rpower[i], e.reg1.pos2 + rpower[i] + EVENT_BUFFER);
	std::string chrstring = SnowTools::CHR_NAME[gr.chr];
	int len;
	char * seq = faidx_fetch_seq(findex, const_cast<char*>(chrstring.c_str()), gr.pos1, gr.pos2 - 1, &len);
	if (!seq)
	  break;
	else
	  seq2 = std::string(seq);

	events << SnowTools::CHR_NAME[e.reg1.chr] << "\t" << e.reg1.pos2 << "\t" << SnowTools::CHR_NAME[gr.chr] << "\t" << gr.pos1 << "\t+\t-\tN\t" << rpower[i] << "\tRAR\t" << (inserts[i].empty() ? "N" : inserts[i]) << "\t" << ename << std::endl;	
        outstring = frag + inserts[i] + seq2;
	assert(outstring.length() < EVENT_BUFFER * 10);
	//std::cerr << rpower[i] << " span " << (std::abs(gr.pos1 - e.reg1.pos2) - rpower[i]) << std::endl;
	assert(rpower[i] == gr.pos1 - e.reg1.pos2);

      // rearrangement doesn't fit
      } else {
	if (verbose) std::cerr << "FAILED TO FIT RAR of length " << rpower[i] << " on break " << e.reg1 << std::endl;
      }

      
    }

    assert(outstring.length() < 50000);
    if (!outstring.empty())
      outfasta << ">" << ename << std::endl << outstring << std::endl;
    else 
      std::cerr << "...empty string on rar " << i << std::endl;
    
  }

}

std::vector<int> drawFromPower(double x0, double x1, double power, int n_draws) {

  assert(power != -1);

  const int PRECISION = 1e6;

  std::vector<int> rpower(n_draws, 0);
  
  for (int i = 0; i < n_draws; ++i) {
    double r = (double)(rand() % PRECISION)/(double)PRECISION;
    double t1 = std::pow(x1, power+1);
    double t2 = std::pow(x0, power+1);
    double tsum = (t1-t2) * r + t2;
    rpower[i] = std::floor(std::pow(tsum, 1 / (power + 1)));
  }

  return rpower;
}

void genRandomSequence(std::string& s, SnowTools::GenomicRegion& gr, int width, faidx_t * findex, SnowTools::GRC& grc) {

  gr.chr = SnowTools::weightedRandom(SnowTools::CHR_CUMSUM_WEIGHT_X);
  std::string chrstring = SnowTools::CHR_NAME[gr.chr];
  char * seq = nullptr;
  int len;
  // get the first sequence
  do {
    if (seq)
      free(seq);
    seq = nullptr;
    s = std::string();
    gr.pos1 = 1e6 + rand() % (int)(SnowTools::CHR_LEN_VEC[gr.chr] - 2e6);
    gr.pos2 = gr.pos1 + width;
    seq = faidx_fetch_seq(findex, const_cast<char*>(chrstring.c_str()), gr.pos1, gr.pos2 - 1, &len);
    if (seq)
      s = std::string(seq);
  } while (!seq || s.find("N") != std::string::npos || grc.findOverlapping(gr));

}
