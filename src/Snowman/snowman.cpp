//-----------------------------------------------
// Copyright 2014 Broad Institute 
// Modified by Jeremiah Wala (jwala@broadinstitute.org)
// Released under the GPL
// Incorportes the core of String Graph Assembler, 
// -- Copyright 2009 Wellcome Trust Sanger Institute
// -- Written by Jared Simpson
// -- Released under the GPL
//-----------------------------------------------

/* CONFIGURE SNOWMAN

sparse=/xchip/gistic/Jeremiah/sparsehash-2.0.2
bamt=/broad/software/free/Linux/redhat_5_x86_64/pkgs/pezmaster31_bamtools-6708a21
seqan=/xchip/gistic/Jeremiah/seqan-trunk/core
htslib=/xchip/gistic/Jeremiah/htslib-1.1
shome=/home/unix/jwala/GIT/isva/Snowman
./configure --with-sparsehash=$sparse --with-bamtools=$bamt --prefix=$shome --with-seqan=$seqan --with-htslib=$htslib

*/

//in core/include/seqan/bam_io/cigar.h had to add parenthesis around lhs.count and rhs.count

#include <string>
#include <iostream>
#include "run.h"
//#include "gather.h"
#include "vcf.h"
//#include "prep_pon.h"
#include <vector>
#include "BreakPoint.h"
#include "genpon.h"
//#include "clean.h"
//#include "nozzle.h"
//#include "bamprocess.h"
//#include "cluster.h"

//debug
#include "SnowUtils.h"
#include "SeqanTools.h"

#define PROGRAM_BIN "snowman"
#define AUTHOR "Jeremiah Wala <jwala@broadinstitute.org>"

//Jeremiah
using namespace std;

static const char *SNOWMAN_USAGE_MESSAGE =
"Program: " PACKAGE_NAME "\n"
"Version: " PACKAGE_VERSION "\n"
"Contact: " AUTHOR " [" PACKAGE_BUGREPORT "]\n"
"Usage: " PROGRAM_BIN " <command> [options]\n\n"
"Commands:\n"
"           run            Run the Snowman SV algorithm on an SV bam pair\n"
  //"           gather         Gather the contigs and alignments, send to BWA-MEM and produce outputs\n"
"           vcf            Tools for handling VCF files within Snowman and across centers\n"
"           refilter       Specify a new set of confidence / evidence criteria for a bp file. Can re-run snowman vcf to make new calls\n"
"           pon            Generate a panel of normal files from a list of germline vcf files. For use with <snowman refilter>\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

int main(int argc, char** argv) {

  /*
  //debug
  TSequence read =     "CGGAGAAAGGGGGATCTGGTGATCCTAATTTTGTCACCAGCTTTGAAATCTTTAGCAACGGATTCTCCATGGAGAG";
  TSequence ref  = "AGGAAATTTCATGTCATTGTACAAAGATGCTACACCGGAGAAAGGGGGATCTGGTGATCCTAATTTTGTCACCAGCTTTGAAATCTTTAGCAACGGATTCTCCATGGAGAGCAGTAGAGGGGTGGTGATGGAACTGCTGCACCGAGTGACGCAAG";
  string rr = "CGGAGAAAGGGGGATCTGGTGATCCTAATTTTGTCACCAGCTTTGAAATCTTTAGCAACGGATTCTCCATGGAGAG";
  string cc = "AGGAAATTTCATGTCATTGTACAAAGATGCTACACCGGAGAAAGGGGGATCTGGTGATCCTAATTTTGTCACCAGCTTTGAAATCTTTAGCAACGGATTCTCCATGGAGAGCAGTAGAGGGGTGGTGATGGAACTGCTGCACCGAGTGACGCAAG";
  TAlign align;
  resize(rows(align), 2); 
  assignSource(row(align,0), ref); 
  assignSource(row(align,1), read);
  Score<int, Simple> sc(4, -2, -4, -2);
  int32_t score = localAlignment(align, sc, 0,0);

  int cl = clock();
  for (int i = 0; i < 10000; i++) {
    cc.find(rr) != string::npos;
    reverse(rr);
    //cc.find(rr) != string::npos;
  }
  
  cout << "score " << score << endl;
  cout << clock() - cl << endl;

  int rll = rr.length();
  int cll = cc.length();
  for (int i = 0; i < 10000; i++) {
    score = localAlignment(align, sc, rll-cll-2, cll-rll+2);
    //score = localAlignment(align, sc, -5,5);
  }

  cout << "score " << score << endl;
  cout << clock() - cl << endl;

  exit(0);

*/
  if (argc <= 1) {
    std::cerr << SNOWMAN_USAGE_MESSAGE;
    return 0;
  } else {
    std::string command(argv[1]);
    if (command == "help" || command == "--help") {
      std::cerr << SNOWMAN_USAGE_MESSAGE;
      return 0;
      //} else if (command == "preprocess") {
      //runPrep(argc-1, argv+1);
    } else if (command == "run") {
      runSnowman(argc -1, argv + 1);
      //} else if (command == "gather") {
      //  runConcat(argc-1, argv+1);
    } else if (command == "refilter") {
      runRefilterBreakpoints(argc-1, argv+1);
    } else if (command == "pon") {
      runGeneratePON(argc-1, argv+1);
    }
  else if (command == "vcf") {
      runVCF(argc-1, argv+1);
    } else {
      std::cerr << SNOWMAN_USAGE_MESSAGE;
      return 0;
    }
    
//else if (command == "cluster") {
      //runCluster(argc-1, argv+1);
    //}

  }
  
  std::cerr << "Done with snowman" << std::endl;
  return 0;

}
