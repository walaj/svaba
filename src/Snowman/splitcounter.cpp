#include "splitcounter.h"

#include <getopt.h>
#include <iostream>
#include <sstream>

#include "SnowmanUtils.h"
#include "SnowTools/BamWalker.h"

SnowTools::GRC bks;

namespace opt {
  static std::string bam;
  static std::string breaks_file;
  static std::string header;
}

static const char* shortopts = "hb:B:H:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "bam",               required_argument, NULL, 'b' },
  { "breaks",               required_argument, NULL, 'B' },
  { "breaks-header",               required_argument, NULL, 'H' },
  { NULL, 0, NULL, 0 }
};

static const char *RUN_SPLITCOUNTER_MESSAGE =
"Usage: snowman splitcounter [OPTION] -b BAM\n\n"
"  Description: Simple utility to count locations of split-alignment breakpoints (eg for PacBio)\n"
"\n"
"  General options\n"
"  -h, --help                          Display this help and exit\n"
"  Required input\n"
"  -b, --bam                           BAM file containing split alignments to count\n"
"  Optional input\n"
"  -B, --breaks                        BED file containing regions to limit analysis to\n"
"  -H, --breaks-header                 BAM file to pair with the BED file (to get header from)\n"
"\n";

void runSplitCounter(int argc, char** argv) {

  parseSplitCounterOptions(argc, argv);

  // open the header BAM
  std::cerr << "...opening the header BAM " << opt::header << std::endl;
  SnowTools::BamWalker bwh(opt::header);

  // open the breaks
  std::cerr << "...reading " << opt::breaks_file << std::endl;
  SnowmanUtils::__open_bed(opt::breaks_file, bks, bwh.header());
  std::cerr << "...read in " << SnowTools::AddCommas(bks.size()) << " break regions " << std::endl;

  // open the BAM
  SnowTools::BamWalker walk(opt::bam);

  SnowTools::BamRead r;
  bool rule;
									
  std::string rules_string = "{\"\" : {\"rules\" : [{\"clip\" : 20}]}}";
  //std::cerr << "Rules: " << rules_string << std::endl;
  SnowTools::MiniRulesCollection mr(rules_string, walk.header());
  walk.SetMiniRulesCollection(mr);

  size_t countr = 0;
  while(walk.GetNextRead(r, rule)) {

    if (!rule)
      continue;
    
    ++countr;

    SnowTools::Cigar cig = r.GetCigar();
    
    assert(cig.size() > 1);

    //if (countr % 100000 == 0)
    //  std::cerr << "...at split read " << SnowTools::AddCommas(countr) << " at position " << r.Brief(walk.header()) << std::endl;

    // if (cig[0].Type() == 'H' || cig[0].Type() == 'S')  // split on left
	 // std::cout << "L\t" << r.ChrName(walk.header()) << "\t" << r.Position() << "\t" << (r.ReverseFlag() ? "-" : "+") << "\t" << r.Qname() << "\t" << r.MapQuality() << "\t" << (r.SecondaryFlag() ? "SEC" : "PRI") << "\t" << (bks.findOverlapping(r.asGenomicRegion()) ? "IN_BK" : "NO_BK") << "\t0" << std::endl;

    //if (cig.back().Type() == 'H' || cig.back().Type() == 'S')  // split on left
    //  std::cout << "R\t" << r.ChrName(walk.header()) << "\t" << r.PositionEnd() << "\t" << (r.ReverseFlag() ? "-" : "+") << "\t" << r.Qname() << "\t" << r.MapQuality() << "\t" << (r.SecondaryFlag() ? "SEC" : "PRI") << "\t" << (bks.findOverlapping(r.asGenomicRegion()) ? "IN_BK" : "NO_BK") << "\t0" << std::endl;

    // look for indels
    int pos = r.Position();
    for (auto& c : cig) {
      //if (c.Type() == 'D' && c.Length() >= 30) 
      if (c.Type() == 'D' || c.Type() == 'I')
	std::cout << c.Type() << "\t" << r.ChrName(walk.header()) << "\t" << pos << "\t" << "+" << "\t" << r.Qname() << "\t" << r.MapQuality() << "\t" << (r.SecondaryFlag() ? "SEC" : "PRI") << "\t" << (bks.findOverlapping(r.asGenomicRegion()) ? "IN_BK" : "NO_BK") << "\t" << c.Length() << std::endl;
      if (c.ConsumesReference())
	pos += c.Length();
    }
    
  }

  std::cerr << "...done" << std::endl;
}

void parseSplitCounterOptions(int argc, char** argv) {

  bool die = false;

  if (argc <= 2) 
    die = true;

  bool help = false;
  std::stringstream ss;

  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'b': arg >> opt::bam; break;
    case 'B': arg >> opt::breaks_file; break;
    case 'H': arg >> opt::header; break;
    case 'h': help = true; break;
    }
  }

  if (die || help) 
    {
      std::cerr << "\n" << RUN_SPLITCOUNTER_MESSAGE;
      die ? exit(EXIT_FAILURE) : exit(EXIT_SUCCESS);
    }


}
