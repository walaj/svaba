#include "AssemblyBamWalker.h"
#include "SeqLib/UnalignedSequence.h"
#include "SeqLib/ReadFilter.h"

#include <chrono>
#include <thread>

#include "SnowmanUtils.h"
#include "run_snowman.h"
#include <fstream>

//#define QNAME "5427150"
#define MIN_ISIZE_FOR_DISC 700

//static std::vector<ContigElement*> contig_elems;
static int num_to_run;
static pthread_mutex_t snow_lock;
static ogzstream all_align, os_allbps;
//static std::ofstream all_align, os_allbps;
static std::string tt, nn; // so hacky
static std::shared_ptr<hts_idx_t> ttindex, nnindex;
static faidx_t* f;
static std::set<std::string> prefixes;

static SeqLib::Filter::ReadFilterCollection * mr;

//static ofstream os_allbps;
static struct timespec start;

static SeqLib::RefGenome * ref_genomeAW, * ref_genome_viral_dummy;

bool __good_contig(const SeqLib::BamRecordVector& brv, const SeqLib::GenomicRegionVector& regions, int max_len, int max_mapq) {
  // NO INDELS
  /*  return (brv.size() && regions.size() && brv.size() < 20  && (brv.size() > 1) && 
	  brv[0].Length() < 20000 && brv[0].CigarSize() < 50 &&
	  max_len > 250 && max_mapq >= 0);
  */
  
  // all hard clip?
  size_t hc = 0;
  for (auto& i : brv)
    if (i.NumHardClip())
      ++hc;
  
  return (brv.size() && hc != brv.size() && regions.size() && brv.size() < 20  && (brv.size() > 1 || brv[0].CigarSize() > 1) && 
	  brv[0].Length() < 20000 && brv[0].CigarSize() < 20 &&
	  max_len > 101 && max_mapq >= 0);
}

//bool runAC(SeqLib::BamRecordVector& brv, faidx_t * f, std::shared_ptr<hts_idx_t> pt, std::shared_ptr<hts_idx_t> pn,
//	   const std::string& t, const std::string& n, const SeqLib::GenomicRegionVector& regions) {
bool runAC(const ContigElement * c) {
  
  SeqLib::GRC trimmed_regions;
  for (auto& r : c->regions)
    if (r.chr < 24)
      trimmed_regions.add(r);

  if (!trimmed_regions.size())
    return true;

  SnowmanBamWalker twalk, nwalk;
  
  if (!tt.empty()) {
    twalk = SnowmanBamWalker(tt);
    twalk.prefix = "t000";
    twalk.max_cov = 200;
    twalk.get_mate_regions = false;
    twalk.SetPreloadedIndex(ttindex);
    twalk.SetMultipleRegions(trimmed_regions);
    twalk.readBam();
  }
  if (!nn.empty()) {
    nwalk = SnowmanBamWalker(nn);
    nwalk.prefix = "n000";
    nwalk.max_cov = 200;
    nwalk.get_mate_regions = false;
    nwalk.SetPreloadedIndex(nnindex);
    nwalk.SetMultipleRegions(trimmed_regions);
    nwalk.readBam();
  }
  SeqLib::BamRecordVector bav_this;
  for (auto& q : twalk.reads)
    bav_this.push_back(q);
  for (auto& q : nwalk.reads)
    bav_this.push_back(q);

  // cluster the reads
  // set region to empty, just won't double-check that cluster overlaps regino. No big deal
  DiscordantClusterMap dmap = DiscordantCluster::clusterReads(bav_this, SeqLib::GenomicRegion(), 37, nullptr);

  std::vector<AlignedContig> this_alc;

  // contruct the index
  SeqLib::UnalignedSequenceVector usv;
  assert(c->brv[0].Qname().length());
  assert(c->brv[0].Sequence().find("N") == std::string::npos);
  
  // here we have to flip if it has (-) alignment to reference.
  // this is because the pipeline assumes it came from a 
  // de novo assembly, which is PRE alignment. Thus, we don't want 
  // to take the i.Sequeunce directly, as this has been reverse-complemented
  // BY THE ALIGNER.
  std::string sss = c->brv[0].Sequence();
  if (c->brv[0].ReverseFlag())
    SeqLib::rcomplement(sss);
  usv.push_back({c->brv[0].Qname(), sss, std::string()});
  
  SeqLib::BWAWrapper bw;
  bw.constructIndex(usv);
  
  if (!prefixes.size()) {
    prefixes.insert("t000"); 
    prefixes.insert("n000");
  }

  this_alc.push_back(AlignedContig(c->brv, prefixes));
  AlignedContig * ac = &this_alc.back();
  for (auto& kk : ac->m_frag_v)
    kk.m_max_indel = 20;

  // align the reads
  alignReadsToContigs(bw, usv, bav_this, this_alc, ref_genomeAW);
  ac->assignSupportCoverage(); // dummies
  ac->splitCoverage(); 
  
  ac->addDiscordantCluster(dmap);

  // get teh coverages
  std::unordered_map<std::string, STCoverage*> covs;
  std::unordered_map<std::string, STCoverage*> clip_covs;
  covs["t000"] = &twalk.cov;
  covs["n000"] = &nwalk.cov;
  clip_covs["t000"] = &twalk.cov;
  clip_covs["n000"] = &nwalk.cov;
  
  std::vector<BreakPoint> allbreaks = ac->getAllBreakPoints(false); // false says dont worry about "local"
  for (auto& i : allbreaks)
    i.repeatFilter();
  for (auto& i : allbreaks)
    i.addCovs(covs, clip_covs);
  for (auto& i : allbreaks)
    i.scoreBreakpoint(8, 2.5, 7, 3, 0);
  for (auto& i : allbreaks)
    i.setRefAlt(ref_genomeAW, ref_genome_viral_dummy);

  double region_width = 0;
  for (auto& i : trimmed_regions)
    region_width += i.width();
  double cov = (double)bav_this.size() / region_width * 250;

  // print message
  if (cov > 100 || num_to_run % 1000 == 0) {
    std::stringstream ss;
#ifndef __APPLE__
    ss << SeqLib::displayRuntime(start);
#endif
    std::cerr << "..read " << SeqLib::AddCommas(bav_this.size()) << " reads from " << 
      trimmed_regions.size() << " region starting at " << trimmed_regions.at(0) << " Queue-left " << 
      SeqLib::AddCommas(num_to_run) << " " << ss.str() << std::endl;
  }

  std::stringstream outr;
  bool no_reads = true;
  for (auto& i : allbreaks)
    outr << i.toFileString(no_reads) << std::endl;

  // MUTEX LOCKED
  ////////////////////////////////////
  pthread_mutex_lock(&snow_lock);  
  --num_to_run;
  //all_align << (*ac) << std::endl;
  os_allbps << outr.str();
  ////////////////////////////////////
  // MUTEX UNLOCKED
  ////////////////////////////////////
  pthread_mutex_unlock(&snow_lock);

  return true;

}


void AssemblyBamWalker::walkDiscovar()
{

  //std::string rules = "global@!hardclip;!duplicate;!qcfail;phred[4,100];length[25,1000]%region@WG%!isize[0,1200];mapq[0,1000]%clip[5,1000]%ins[1,1000];mapq[0,100]%del[1,1000];mapq[1,1000]%mapped;!mate_mapped;mapq[1,1000]%mate_mapped;!mapped";
  std::string rules = "{\"global\" : {\"duplicate\" : false, \"qcfail\" : false}, \"\" : { \"rules\" : [{\"isize\" : [MINS,0]},{\"rr\" : true},{\"ff\" : true}, {\"rf\" : true}, {\"ic\" : true}, {\"clip\" : 5, \"phred\" : 4}, {\"ins\" : true}, {\"del\" : true}, {\"mapped\": true , \"mate_mapped\" : false}, {\"mate_mapped\" : true, \"mapped\" : false}]}}";

  rules.replace(rules.find("MINS"), 4, std::to_string(MIN_ISIZE_FOR_DISC));

  mr = new SeqLib::Filter::ReadFilterCollection(rules, SeqLib::BamHeader());

  f = findex;
  tt = tbam;
  nn = nbam;
  nnindex = nindex;
  ttindex = tindex;

  SnowmanUtils::fopen(id + ".alignments.txt.gz", all_align);
  SnowmanUtils::fopen(id + ".bps.txt.gz", os_allbps);
  os_allbps << BreakPoint::header() << std::endl;

  SeqLib::BamRecord r;
  std::cerr << "...starting to walk assembly BAM" << std::endl;

  std::vector<AlignedContig> ac_vec;
  
  std::unordered_map<std::string, SeqLib::BamRecordVector> map;

  SeqLib::BamRecordVector brv;

  SeqLib::GenomicRegionVector regions;

  // start the timer
#ifndef __APPLE__
  clock_gettime(CLOCK_MONOTONIC, &start);
#endif

  // open the mutex
  if (pthread_mutex_init(&snow_lock, NULL) != 0) {
      printf("\n mutex init failed\n");
      return;
  }

  // Create the queue and consumer (worker) threads
  wqueue<AssemblyWalkerWorkItem*>  queue;
  std::vector<ConsumerThread<AssemblyWalkerWorkItem>*> threadqueue;
  for (int i = 0; i < numThreads; i++) {
    ConsumerThread<AssemblyWalkerWorkItem>* threadr = new ConsumerThread<AssemblyWalkerWorkItem>(queue, true);
    threadr->start();
    threadqueue.push_back(threadr);
  }
  
  int count = 0, acc_count = 0;
  std::string curr_name;

  std::vector<AssemblyWalkerWorkItem*> tmp_queue;

  // open ref genome for extracting sequences
  ref_genomeAW = new SeqLib::RefGenome(refGenome);
  
  int max_mapq = 0, max_len = 0;
  while(GetNextRecord(r)) {

    if (++count % 200000 == 0) 
      std::cerr << "...checking contig " << SeqLib::AddCommas(count) << " to see if we should parse" << std::endl;

    // give each contig a unique prefix
    //r.SetQname("contig_" + r.Qname());

    // first one
    if (curr_name.empty())
      curr_name = r.Qname();

    // add the MC tag
    if (r.ChrID() < br->n_targets && r.ChrID() >= 0) {
      r.AddZTag("MC", std::string(br->target_name[r.ChrID()]));
    }
    else 
      continue;

    // add to existing
    if (curr_name == r.Qname()) {

      SeqLib::GenomicRegion gr = r.asGenomicRegion();
      gr.Pad(100);
      regions.push_back(gr);
      brv.push_back(r);
      max_mapq = std::max(r.MapQuality(), max_mapq);
      max_len = std::max(r.Length(), max_len);
      
    } else {

#ifdef QNAME      
      if (brv.size())
	if (brv[0].Qname() == QNAME || true) {
	  std::cerr << " found it here -- checking " << std::endl;
	  std::cerr << " brv.size() " << brv.size() << " regions.size() " << regions.size() << " brv.CigarSize() " << brv[0].CigarSize() << 
	    " brv.Length() " << brv[0].Length() << " max_len " << max_len << " max_mapq " << max_mapq << std::endl;
	}
#endif

      if (true || __good_contig(brv, regions, max_len, max_mapq)) { 

	++acc_count;
	//AssemblyWalkerWorkItem * item = new AssemblyWalkerWorkItem(brv, findex, tindex, nindex, regions, tbam, nbam);
	AssemblyWalkerWorkItem * item = new AssemblyWalkerWorkItem(new ContigElement(brv, regions)); 
	
#ifdef QNAME      
	if (brv[0].Qname() == QNAME)
	  std::cerr << " found it here -- adding " << std::endl;
#endif

	tmp_queue.push_back(item);
	regions.clear();
      }

      // prepare the next one
      brv.clear();
      brv.push_back(r);
      curr_name = r.Qname();
      regions.clear();
      SeqLib::GenomicRegion gr = r.asGenomicRegion();
      gr.Pad(1000);
      regions.push_back(gr);
      max_mapq = r.MapQuality();
      max_len  = r.Length();
    }
  }


  // add the last one
  if (__good_contig(brv, regions, max_len, max_mapq)) {
    AssemblyWalkerWorkItem * item = new AssemblyWalkerWorkItem(new ContigElement(brv, regions)); 
    tmp_queue.push_back(item);
  }

  std::cerr << "...done adding " << SeqLib::AddCommas(tmp_queue.size()) << " contigs. Firing up " << 
    numThreads << " re-alignment threads" << std::endl;
  num_to_run = tmp_queue.size();
  
  if (!num_to_run)
    return;

  for (auto& i : tmp_queue) 
    queue.add(i);
  
  // wait for the threads to finish
  for (int i = 0; i < numThreads; i++) 
    threadqueue[i]->join();  

  all_align.close();
  os_allbps.close();

}

