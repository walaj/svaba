#include "splitcounter.h"

#include <getopt.h>
#include <iostream>
#include <sstream>
#include <cassert>
#include "gzstream.h"

#include "SnowmanUtils.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/ReadFilter.h"


#define MIN_CLIP_FACTOR 2
#define MIN_MATCH_FACTOR 4
#define MIN_MAPQ 10
#define READ_LEN 151

SeqLib::GRC bks;

namespace opt {
  static std::string bam;
  static std::string breaks_file;
  static std::string header;
  static std::string fasta;
  static int split_num = 1000;
}

static const char* shortopts = "hb:B:H:f:n:o:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "bam",               required_argument, NULL, 'b' },
  { "fasta",               required_argument, NULL, 'f' },
  { "breaks",               required_argument, NULL, 'B' },
  { "breaks-header",               required_argument, NULL, 'H' },
  { "num-bases",               required_argument, NULL, 'n' },
  { "out-fasta",               required_argument, NULL, 'o' },
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

static const char *RUN_FASTASPLIT_MESSAGE =
"Usage: snowman fastasplit -f fasta -n num_bases > outfile.fa\n\n"
"  Description: Divide a fasta into smaller sub-sequences, splitting some seqs in the middle\n"
"\n"
"  General options\n"
"  -h, --help                          Display this help and exit\n"
"  Required input\n"
"  -f, --fasta                         BAM file containing split alignments to count\n"
"  -n, --num-bases                     Number of bases to split into\n"
"\n";


void parseFastaSplitOptions(int argc, char** argv);

void runSplitFasta(int argc, char** argv) {

  parseFastaSplitOptions(argc, argv);

  std::cerr << " SPLITTING: " << opt::fasta << "\n" <<
    " NUM: " << SeqLib::AddCommas(opt::split_num) << std::endl;
  
  // open the fasta
  igzstream infile(opt::fasta.c_str(), std::ios::in);
  
  // confirm that it is open
  if (!infile) {
    std::cerr << "Can't read file " << opt::fasta << " for parsing fasta" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  int line_count = 0;

  //
  std::string line;
  std::string curr_line;
  size_t this_line_count = 0; // new line count
  while (std::getline(infile, line, '\n')) {
    
    ++this_line_count;
    if (this_line_count % 20000)
      std::cerr << "...on fasta file line (new-line separated) " << SeqLib::AddCommas(this_line_count) << std::endl;

    if (line.empty() || line.find(">") != std::string::npos)
      continue;
    
    // doesn't hit limit
    if ((int)(line.length() + curr_line.length()) < opt::split_num)
      curr_line = curr_line + line;
    
    // hits limit, but only if adding the new line
    else if ((int)line.length() <= opt::split_num) {
      std::cout << ">" << line_count << std::endl;
      ++line_count;
      std::cout << curr_line << std::endl;
      curr_line = line;
    }
      
    // the line itself is too big, split up
    else if ((int)line.length() > opt::split_num) {
      int lsize = line.size();
      int lstart = 0;
      while (lstart < (int)line.size()) {
	std::cout << ">" << line_count << std::endl; 
	++line_count;
	std::cout << line.substr(lstart, std::min(opt::split_num, lsize - lstart)) << std::endl;
	lstart += opt::split_num;
      }

    }
    
    
  }

}

void runSplitCounter(int argc, char** argv) {

  parseSplitCounterOptions(argc, argv);

  // open the header BAM
  std::cerr << "...opening the header BAM " << opt::header << std::endl;
  SeqLib::BamReader bwh;
  assert(bwh.Open(opt::header));

  // open the breaks
  std::cerr << "...reading " << opt::breaks_file << std::endl;
  SnowmanUtils::__open_bed(opt::breaks_file, bks, bwh.Header());
  std::cerr << "...read in " << SeqLib::AddCommas(bks.size()) << " break regions " << std::endl;

  // open the BAM
  SeqLib::BamReader walk;
  assert(walk.Open(opt::bam));

  SeqLib::BamRecord r;
	
  std::string rules_string = "{\"\" : {\"rules\" : [{\"mapq\" : 10}]}}";

  std::cerr << "Rules: " << rules_string << std::endl;
  SeqLib::Filter::ReadFilterCollection mr(rules_string, walk.Header());

  // convert qname (big) to in (small)
  std::unordered_map<std::string, int> qname_set;

  size_t countr = 0;
  while(walk.GetNextRecord(r)) {

    if (r.ChrID() > 23)
      break;

    if (!mr.isValid(r))
      continue;
    
    ++countr;

    SeqLib::Cigar cig = r.GetCigar();
    
    if (cig.size() <= 1)
      continue;
    
    int matchlen = 0;
    for (auto& c : cig)
      if (c.Type() == 'M')
	matchlen += c.Length();
    if (matchlen < MIN_MATCH_FACTOR * READ_LEN)
      continue;

    if (!qname_set.count(r.Qname()))
      qname_set[r.Qname()] = qname_set.size();
    std::string QN = std::to_string(qname_set[r.Qname()]);

    if (countr % 100000 == 0 )
      std::cerr << "...at split read " << SeqLib::AddCommas(countr) << " at position " << r.Brief() << std::endl;

    if ( (cig[0].Type() == 'H' || cig[0].Type() == 'S') && cig[0].Length() >= MIN_CLIP_FACTOR * READ_LEN)  // split on left
      std::cout << "L\t" << r.ChrName(walk.Header()) << "\t" << r.Position() << 
	"\t" << (r.ReverseFlag() ? "-" : "+") << "\t" << QN << "\t" << r.MapQuality() << "\t" << (r.SecondaryFlag() ? "SEC" : "PRI") << "\t" << (bks.CountOverlaps(r.asGenomicRegion()) ? "BK" : "NBK") << "\t0" << std::endl;
    
    if ( (cig.back().Type() == 'H' || cig.back().Type() == 'S') && cig.back().Length() >= MIN_CLIP_FACTOR * READ_LEN)  // split on right
      std::cout << "R\t" << r.ChrName(walk.Header()) << "\t" << r.PositionEnd() 
		<<  "\t" << (r.ReverseFlag() ? "-" : "+") << "\t" << QN << "\t" << r.MapQuality() << "\t" << (r.SecondaryFlag() ? "SEC" : "PRI") << "\t" << (bks.CountOverlaps(r.asGenomicRegion()) ? "BK" : "NBK") << "\t0" << std::endl;

    // look for indels
    int pos = r.Position();
    for (auto& c : cig) {
      if ((c.Type() == 'D' || c.Type() == 'I') && c.Length() >= 30) 
	std::cout << c.Type() << "\t" << r.ChrName(walk.Header()) << "\t" << pos << "\t+\t" 
		  << QN << "\t" << r.MapQuality() << "\t" << (r.SecondaryFlag() ? "SEC" : "PRI") << "\t" 
		  << (bks.CountOverlaps(r.asGenomicRegion()) ? "BK" : "NBK") << "\t" << c.Length() << std::endl;
      if (c.ConsumesReference())
	pos += c.Length();
    }

    // try BLAT realignment
    //SeqLib::BamReadVector brv;
    //blat.querySequence(r.Qname(), r.Sequence(), brv);
    
    //for(auto& a : brv)
    // w.writeAlignment(a);
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


void parseFastaSplitOptions(int argc, char** argv) {

  bool die = false;

  if (argc <= 2) 
    die = true;

  bool help = false;
  std::stringstream ss;

  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'f': arg >> opt::fasta; break;
    case 'n': arg >> opt::split_num; break;
    case 'h': help = true; break;
    }
  }

  if (die || help) 
    {
      std::cerr << "\n" << RUN_FASTASPLIT_MESSAGE;
      die ? exit(EXIT_FAILURE) : exit(EXIT_SUCCESS);
    }


}
