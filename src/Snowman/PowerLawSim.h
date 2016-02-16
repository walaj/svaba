#ifndef SNOWMAN_POWER_LAW_SIM_HH
#define SNOWMAN_POWER_LAW_SIM_HH

#include <fstream>

#include "SnowTools/GenomicRegionCollection.h"
#include "SnowTools/BWAWrapper.h"

void PowerLawSim(faidx_t* findex, int num_breaks, double power_law, SnowTools::GRC& grc, 
		 std::ofstream& outfasta, std::ofstream& events);
std::vector<int> drawFromPower(double x0, double x1, double power, int n_draws);
void genRandomSequence(std::string& s, SnowTools::GenomicRegion& gr, int width, faidx_t * findex, SnowTools::GRC& grc);

struct SVEvent {
  
  SnowTools::GenomicRegion reg1;
  SnowTools::GenomicRegion reg2;
  
  int break1;
  int break2;

  std::string ins;
  std::string r_ins;

  std::string etype;

  int number;
};

#endif
