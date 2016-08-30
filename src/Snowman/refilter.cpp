#include "refilter.h"

#include <getopt.h>
#include <sstream>
#include <iostream>

#include "gzstream.h"
#include "SeqLib/BamReader.h"

#include "vcf.h"
#include "BreakPoint2.h"
#include "SnowmanUtils.h"


static ogzstream os_allbps_r;

namespace opt {

  static std::string input_file;
  static std::string output_file;
  static std::string pon = "";
  static std::string analysis_id = "refilter";
  static bool noreads = false;
  static std::string indel_mask;

  static std::string ref_index = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";

  static std::string normal_bam;
  static std::string tumor_bam;

  static std::string bam;

  static int verbose = 1;

  static double lod = 8;
  static double lod_db = 7;
  static double lod_no_db = 2.5;
  static double lod_germ = 3;
}

enum { 
  OPT_LOD,
  OPT_DLOD,
  OPT_NDLOD,
  OPT_LODGERM
};


static const char* shortopts = "hxi:a:v:q:g:m:b:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "input-bps",               required_argument, NULL, 'i'},
  { "bam",                     required_argument, NULL, 'b'},
  { "panel-of-normals",        required_argument, NULL, 'q'},
  { "tumor-bam",               required_argument, NULL, 't' },
  { "normal-bam",               required_argument, NULL, 'n' },
  { "indel-mask",              required_argument, NULL, 'm'},
  { "reference-genome",        required_argument, NULL, 'g'},
  { "analysis-id",             required_argument, NULL, 'a'},
  { "no-reads",                no_argument, NULL, 'x'},
  { "verbose",                 required_argument, NULL, 'v' },
  { "lod",                   required_argument, NULL, OPT_LOD },
  { "lod-dbsnp",             required_argument, NULL, OPT_DLOD },
  { "lod-no-dbsnp",          required_argument, NULL, OPT_NDLOD },
  { "lod-germ",          required_argument, NULL, OPT_LODGERM },

  { NULL, 0, NULL, 0 }
};


static const char *BP_USAGE_MESSAGE =
"Usage: snowman refilter [OPTION] -i bps.txt.gz -o bps.new.txt.gz\n\n"
"  Description: \n"
"\n"
"  General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h, --help                           Display this help and exit\n"
"  -g, --reference-genome               Path to indexed reference genome to be used by BWA-MEM. Default is Broad hg19 (/seq/reference/...)\n"
"  -b, --opt-bam                        Input BAM file to get header from\n"
"  Required input\n"
"  -i, --input-bps                      Original bps.txt.gz file\n"
"  -b, --bam                            BAM file used to grab header from\n"
"  Variant filtering and classification\n"
"      --lod                            LOD cutoff to classify indel as real / artifact [8]\n"
"      --lod-dbsnp                      LOD cutoff to classify indel as somatic (at DBsnp site) [7]\n"
"      --lod-no-dbsnp                   LOD cutoff to classify indel as somatic (not at DBsnp site) [2.5]\n"
"  Optional input\n"                       
"  -q, --panel-of-normals               Panel of normals files generated from snowman pon\n"                       
"  -a, --id-string                      String specifying the analysis ID to be used as part of ID common.\n"
"  -x, --no-reads                       Flag to turn off recording of supporting reads. Setting this flag greatly reduces file size.\n"
"  -m, --indel-mask                     BED-file with blacklisted regions for indel calling. Default none\n"
"\n";

// parse the command line options
void parseBreakOptions(int argc, char** argv) {
  bool die = false;
  
  if (argc <= 2) 
    die = true;
  
  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'h': die = true; break;
    case 'g': arg >> opt::ref_index; break;
    case 'm': arg >> opt::indel_mask; break;
    case 'i': arg >> opt::input_file; break;
    case 'v': arg >> opt::verbose; break;
    case 'q': arg >> opt::pon; break;
    case 'a': arg >> opt::analysis_id; break;
    case 'x': opt::noreads = true; break;
    case OPT_LOD: arg >> opt::lod; break;
    case OPT_DLOD: arg >> opt::lod_db; break;
    case OPT_NDLOD: arg >> opt::lod_no_db; break;
    case OPT_LODGERM: arg >> opt::lod_germ; break;
    case 'b': arg >> opt::bam; break; 
    }
  }
  
  if (opt::input_file.length() == 0)
    die = true;
  
  if (opt::bam.length() == 0) {
    std::cerr << "BAM is required (for the header)" << std::endl;
    die = true;
  }

  if (die) {
    std::cerr << "\n" << BP_USAGE_MESSAGE;
    exit(1);
  }
}

void runRefilterBreakpoints(int argc, char** argv) {
  
  parseBreakOptions(argc, argv);
  
  opt::output_file = opt::analysis_id + ".filtered.bps.txt.gz";
  if (opt::verbose > 0) {
    std::cerr << "Input bps file:  " << opt::input_file << std::endl;
    std::cerr << "Output bps file: " << opt::output_file << std::endl;
    std::cerr << "Panel of normals file: " << opt::pon << std::endl;
    std::cerr << "Indel mask BED:      " << opt::indel_mask << std::endl;
    std::cerr << "Analysis id: " << opt::analysis_id << std::endl;
    std::cerr << "    LOD cutoff (artifact vs real) :  " << opt::lod << std::endl << 
    "    LOD cutoff (somatic, at DBSNP):  " << opt::lod_db << std::endl << 
    "    LOD cutoff (somatic, no DBSNP):  " << opt::lod_no_db << std::endl <<
    "    LOD cutoff (germline, AF>=0.5 vs AF=0):  " << opt::lod_germ << std::endl;

  }

    
  if (!SeqLib::read_access_test(opt::input_file)) {
    std::cerr << "ERROR: Cannot read file " << opt::input_file  << std::endl;
    exit(EXIT_FAILURE);
  }
  
  SeqLib::BamReader bwalker;
  assert(bwalker.Open(opt::bam));

  // load the reference index
  //if (opt::verbose > 0)
  //  std::cerr << "attempting to load: " << opt::ref_index << std::endl;
  //faidx_t * findex = fai_load(opt::ref_index.c_str());  // load the reference
  
  // open the output file
  /*igzstream iz(opt::input_file.c_str());
  if (!iz) {
    std::cerr << "Can't read file " << opt::input_file << std::endl;
    exit(EXIT_FAILURE);
  }
  */
  // read in the PON
  //std::unique_ptr<SnowTools::PON> pmap;
  //SnowTools::BreakPoint::readPON(opt::pon, pmap);
  
  // read the indel mask
  //GenomicIntervalTreeMap grm_mask;
  if (!SeqLib::read_access_test(opt::indel_mask) && opt::indel_mask.length()) {
    std::cerr << "indel mask " << opt::indel_mask << " does not exist / is not readable. Skipping indel masking."  << std::endl;
    opt::indel_mask = "";
  }
  
  SeqLib::GRC grv_mask;
  if (opt::indel_mask.length()) 
    grv_mask = SeqLib::GRC(opt::indel_mask, bwalker.Header());
  //grv_mask.regionFileToGRV(opt::indel_mask, 0, NULL);
 
  VCFHeader header;
  header.filedate = SnowmanUtils::fileDateString();
  header.source = "";//opt::args;
  header.reference = "";//opt::refgenome;

  // open bps file
  std::string new_bps_file = opt::analysis_id + ".bps.txt.gz";
  SnowmanUtils::fopen(new_bps_file, os_allbps_r);

  // read in the BPS
  std::vector<std::string> allele_names; // store with real name
  std::map<std::string, SampleInfo> tmp_alleles;
  std::string line, line2, val;
  igzstream infile(opt::input_file.c_str(), std::ios::in);
  size_t line_count = 0;
  while (getline(infile, line, '\n')) {
    if (line_count == 0) {
      os_allbps_r << line << std::endl;
      std::istringstream f(line);
      size_t scount = 0;
      while (std::getline(f, val, '\t')) {
	++scount;
	if (scount > 33) { // 34th column should be first sample ID
	  assert(val.at(0) == 't' || val.at(0) == 'n');
	    allele_names.push_back(val);
	}
      }

    } else {
	BreakPoint * bp = new BreakPoint(line, bwalker.Header());

	// fill in with the correct names from the header of bps.txt
	std::string id = "";
	for (auto& i : allele_names) {
	  id += "A";
	  tmp_alleles[i] = bp->allele[id];
	}
	bp->allele = tmp_alleles;

	// fill in discordant info
	for (auto& i : bp->allele) {
	  if (i.first.at(0) == 't')
	    bp->dc.tcount += i.second.disc;
	  else
	    bp->dc.tcount += i.second.disc;
	}

	// score them
	bp->scoreBreakpoint(opt::lod, opt::lod_db, opt::lod_no_db, opt::lod_germ, 0);
	os_allbps_r << bp->toFileString(/*opt::no_reads*/false) << std::endl;
	delete bp;
      }
    ++line_count;
  }

  os_allbps_r.close();
  
  // primary VCFs
  std::cerr << " input file " << opt::input_file << std::endl;
  if (SeqLib::read_access_test(new_bps_file)) {
    if (opt::verbose)
      std::cerr << "...making the primary VCFs (unfiltered and filtered) from file " << new_bps_file << std::endl;
    VCFFile snowvcf(new_bps_file, opt::analysis_id, bwalker.Header(), header);
    std::string basename = opt::analysis_id + ".snowman.unfiltered.";
    snowvcf.include_nonpass = true;
    snowvcf.writeIndels(basename, false, false);
    snowvcf.writeSVs(basename, false, false);

    basename = opt::analysis_id + ".snowman.";
    snowvcf.include_nonpass = false;
    snowvcf.writeIndels(basename, false, false);
    snowvcf.writeSVs(basename, false, false);

  } else {
    std::cerr << "Failed to make VCF. Could not file bps file " << opt::input_file << std::endl;
  }


  /*
 
  ogzstream oz(opt::output_file.c_str(), std::ios::out);
  if (!oz) {
    std::cerr << "Can't write to output file " << opt::output_file << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (opt::verbose)
    std::cerr << "...refiltering variants" << std::endl;
  
  // set the header
  oz << SnowTools::BreakPoint::header() << std::endl;
  
  std::string line;

  //skip the header
  std::getline(iz, line, '\n');
  
  size_t count = 0;
  while (std::getline(iz, line, '\n')) {
    
    if (++count % 10000 == 1 && opt::verbose > 0)
      std::cerr << "filtering breakpoint " << SnowTools::AddCommas(count) << std::endl;
    
    SnowTools::BreakPoint bp(line, bwalker.header());
    
    // check if in panel of normals
    //bp.checkPon(pmap);
    
    // check for blacklist
    bp.checkBlacklist(grv_mask);
    
    // check for repeat sequence
    if (bp.b1.gr.chr < 24) {
      SnowTools::GenomicRegion gr = bp.b1.gr;
      gr.pad(20);
      
      // get the reference seqeuence for this piece
      //int len;
      //std::string chrstring = SnowTools::GenomicRegion::chrToString(gr.chr);
      //char * seq = faidx_fetch_seq(findex, const_cast<char*>(chrstring.c_str()), gr.pos1-1, gr.pos2-1, &len);
      //if (!seq) {
      //std::cerr <<  "faidx_fetch_seq fail" << std::endl;
      //}
      //std::string seqr = std::string(seq);

      //for (auto& i : repr)
      //if (seqr.find(i) != std::string::npos)
      //	  bp.repeat_seq = i;
      
    }
    
    //if (bp.hasMinimal() && bp.b1.gr.chr < 24 && bp.gr2.chr < 24)
    oz << bp.toFileString(opt::noreads) << std::endl;
  }
  
  oz.close();
  */
  // make the VCF files
  /*
    if (opt::verbose)
    std::cerr << "...converting filtered.bps.txt.gz to vcf files" << std::endl;
    VCFFile filtered_vcf(opt::output_file, opt::ref_index.c_str(), '\t', opt::analysis_id);
    
    // output the vcfs
      string basename = opt::analysis_id + ".broad-snowman.DATECODE.filtered.";
      filtered_vcf.writeIndels(basename, true); // zip them -> true
      filtered_vcf.writeSVs(basename, true);
  */
}
