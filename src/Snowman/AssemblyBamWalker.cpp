#include "AssemblyBamWalker.h"
#include "SnowmanUtils.h"
#include "run_snowman.h"

#include <chrono>
#include <thread>

//#define QNAME "5427150"

//static std::vector<ContigElement*> contig_elems;
static int num_to_run;
static pthread_mutex_t snow_lock;
static ogzstream all_align, os_allbps;
static std::string tt, nn; // so hacky
static std::shared_ptr<hts_idx_t> ttindex, nnindex;
static faidx_t* f;

//static ofstream os_allbps;
static struct timespec start;

static SnowTools::MiniRulesCollection * mr;

bool __good_contig(const SnowTools::BamReadVector& brv, const SnowTools::GenomicRegionVector& regions, int max_len, int max_mapq) {
  // NO INDELS
  return (brv.size() && regions.size() && brv.size() < 20  && (brv.size() > 1) && 
	  brv[0].Length() < 20000 && brv[0].CigarSize() < 50 &&
	  max_len > 250 && max_mapq >= 0);

  return (brv.size() && regions.size() && brv.size() < 20  && (brv.size() > 1 || brv[0].CigarSize() > 1) && 
	  brv[0].Length() < 20000 && brv[0].CigarSize() < 20 &&
	  max_len > 250 && max_mapq >= 0);
}

//bool runAC(SnowTools::BamReadVector& brv, faidx_t * f, std::shared_ptr<hts_idx_t> pt, std::shared_ptr<hts_idx_t> pn,
//	   const std::string& t, const std::string& n, const SnowTools::GenomicRegionVector& regions) {
bool runAC(const ContigElement * c) {
	   
  
  SnowmanBamWalker twalk(tt);
  SnowmanBamWalker nwalk(nn);
  
  twalk.prefix = "t000";
  nwalk.prefix = "n000";
  twalk.max_cov = 200;
  nwalk.max_cov = 200;
  twalk.get_mate_regions = false;
  nwalk.get_mate_regions = false;
  twalk.setBamWalkerRegions(c->regions, ttindex);
  nwalk.setBamWalkerRegions(c->regions, nnindex);

  twalk.SetMiniRulesCollection(*mr);
  nwalk.SetMiniRulesCollection(*mr);

  twalk.readBam();
  nwalk.readBam();

  //if (c->brv[0].Qname() == "9689804")
  //  for (auto& i : nwalk.reads)
  //    std::cerr << i << std::endl;

  SnowTools::BamReadVector bav_this;
  for (auto& q : twalk.reads)
    bav_this.push_back(q);
  for (auto& q : nwalk.reads)
    bav_this.push_back(q);

  // cluster the reads
  // set region to empty, just won't double-check that cluster overlaps regino. No big deal
  SnowTools::DiscordantClusterMap dmap = SnowTools::DiscordantCluster::clusterReads(bav_this, SnowTools::GenomicRegion());

  std::vector<SnowTools::AlignedContig> this_alc;

  // contruct the index
  SnowTools::USeqVector usv;
  assert(c->brv[0].Qname().length());
  assert(c->brv[0].Sequence().find("N") == std::string::npos);
  
  // here we have to flip if it has (-) alignment to reference.
  // this is because the pipeline assumes it came from a 
  // de novo assembly, which is PRE alignment. Thus, we don't want 
  // to take the i.Sequeunce directly, as this has been reverse-complemented
  // BY THE ALIGNER.
  std::string sss = c->brv[0].Sequence();
  if (c->brv[0].ReverseFlag())
    SnowTools::rcomplement(sss);
  usv.push_back({c->brv[0].Qname(), sss});
  
  SnowTools::BWAWrapper bw;
  bw.constructIndex(usv);
  
  this_alc.push_back(SnowTools::AlignedContig(c->brv));
  SnowTools::AlignedContig * ac = &this_alc.back();
  for (auto& kk : ac->m_frag_v)
    kk.m_max_indel = 20;

  // align the reads
  alignReadsToContigs(bw, usv, bav_this, this_alc);
  ac->assignSupportCoverage(); // dummies
  ac->splitCoverage(); 
  
  ac->addDiscordantCluster(dmap);

  // get teh coverages
  std::unordered_map<std::string, SnowTools::STCoverage*> covs;
  std::unordered_map<std::string, SnowTools::STCoverage*> clip_covs;
  covs["t000"] = &twalk.cov;
  covs["n000"] = &nwalk.cov;
  clip_covs["t000"] = &twalk.cov;
  clip_covs["n000"] = &nwalk.cov;
  
  std::vector<SnowTools::BreakPoint> allbreaks = ac->getAllBreakPoints(false); // false says dont worry about "local"
  for (auto& i : allbreaks)
    i.repeatFilter();
  for (auto& i : allbreaks)
    i.addCovs(covs, clip_covs);
  for (auto& i : allbreaks)
    i.scoreBreakpoint(8, 2.5, 7);
  for (auto& i : allbreaks)
    i.setRefAlt(f, nullptr);

  double region_width = 0;
  for (auto& i : c->regions)
    region_width += i.width();
  double cov = (double)bav_this.size() / region_width * 250;

  // print message
  if (cov > 100 || num_to_run % 1000 == 0) {
    std::stringstream ss;
    ss << SnowTools::displayRuntime(start);
    std::cerr << "..read " << SnowTools::AddCommas(bav_this.size()) << " reads from " << 
      c->regions.size() << " region starting at " << c->regions.at(0) << " Queue-left " << 
      SnowTools::AddCommas(num_to_run) << " " << ss.str() << std::endl;
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

  std::string rules = "global@!hardclip;!duplicate;!qcfail;phred[4,100];length[25,1000]%region@WG%!isize[0,1200];mapq[0,1000]%clip[5,1000]%ins[1,1000];mapq[0,100]%del[1,1000];mapq[1,1000]%mapped;!mate_mapped;mapq[1,1000]%mate_mapped;!mapped";
  mr = new SnowTools::MiniRulesCollection(rules);

  f = findex;
  tt = tbam;
  nn = nbam;
  nnindex = nindex;
  ttindex = tindex;

  //SnowmanUtils::fopen("assembly.alignments.txt.gz", all_align);
  SnowmanUtils::fopen("assembly.bps.txt.gz", os_allbps);
  os_allbps << SnowTools::BreakPoint::header() << endl;

  SnowTools::BamRead r;
  std::cerr << "...starting to walk assembly BAM" << std::endl;

  std::vector<SnowTools::AlignedContig> ac_vec;
  
  std::unordered_map<std::string, SnowTools::BamReadVector> map;

  SnowTools::BamReadVector brv;

  SnowTools::GenomicRegionVector regions;

  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);

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
  
  bool rule;
  int count = 0, acc_count = 0;
  std::string curr_name;

  std::vector<AssemblyWalkerWorkItem*> tmp_queue;
  
  int max_mapq = 0, max_len = 0;
  while(GetNextRead(r, rule)) {
    if (++count % 200000 == 0) 
      std::cerr << "...checking contig " << SnowTools::AddCommas(count) << " to see if we should parse" << std::endl;

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

      SnowTools::GenomicRegion gr = r.asGenomicRegion();
      gr.pad(100);
      regions.push_back(gr);
      brv.push_back(r);
      max_mapq = std::max(r.MapQuality(), max_mapq);
      max_len = std::max(r.Length(), max_len);
      
    } else {

#ifdef QNAME      
      if (brv.size())
	if (brv[0].Qname() == QNAME) {
	  std::cerr << " found it here -- checking " << std::endl;
	  std::cerr << " brv.size() " << brv.size() << " regions.size() " << regions.size() << " brv.CigarSize() " << brv[0].CigarSize() << 
	    " brv.Length() " << brv[0].Length() << " max_len " << max_len << " max_mapq " << max_mapq << std::endl;
	}
#endif

      if (__good_contig(brv, regions, max_len, max_mapq)) { 

	++acc_count;
	//AssemblyWalkerWorkItem * item = new AssemblyWalkerWorkItem(brv, findex, tindex, nindex, regions, tbam, nbam);
	AssemblyWalkerWorkItem * item = new AssemblyWalkerWorkItem(new ContigElement(brv, regions)); 
	
#ifdef QNAME      
	if (brv[0].Qname() == QNAME)
	  std::cerr << " found it here -- adding " << std::endl;
#endif

	tmp_queue.push_back(item);
	//queue.add(item);
	//std::cerr << " qua " << queue.size() <<std::endl;

	regions.clear();
	/// hack
	//while (queue.size() > 1000) {
	//  std::cerr << queue.size() << std::endl;
	//  std::this_thread::sleep_for(std::chrono::milliseconds(1));
	//	}

	/*if (acc_count % 10000 == 0)
	  std::cerr << " adding count " << SnowTools::AddCommas(count) << " Nmap: " << brv.size() << " Name: " << 
	    brv[0].Qname() << " " << brv[0].asGenomicRegion() << " qsize " << tmp_queue.size() << std::endl;	
	*/
	/*if (count > 5e6) {
	  all_align.close();
	  os_allbps.close();
	  return; //exit(0);
	  break;
	  }*/


	/*
	std::vector<SnowTools::AlignedContig> this_alc;

	int rsize = regions.size();

	// grab the reads
	assert(regions.size());
	twalk.setBamWalkerRegions(regions, tindex);
	nwalk.setBamWalkerRegions(regions, nindex);
	twalk.readBam();
	nwalk.readBam();
	SnowTools::BamReadVector bav_this;
	for (auto& q : twalk.reads)
	  bav_this.push_back(q);
	for (auto& q : nwalk.reads)
	  bav_this.push_back(q);
	twalk.resetAll(); twalk.reads.clear();
	nwalk.resetAll(); nwalk.reads.clear();

	// contruct the index
	SnowTools::USeqVector usv;
	for (auto i : brv) {
	  assert(i.Qname().length());
	  assert(i.Sequence().find("N") == std::string::npos);
	  usv.push_back({i.Qname(), i.Sequence()});
	}
	SnowTools::BWAWrapper bw;
	bw.constructIndex(usv);

	this_alc.push_back(SnowTools::AlignedContig(brv));
	SnowTools::AlignedContig * ac = &this_alc.back();

	// align the reads
	alignReadsToContigs(bw, usv, bav_this, this_alc);
	ac->assignSupportCoverage(); // dummies
	ac->splitCoverage(); 

	// get teh coverages
	std::unordered_map<std::string, SnowTools::STCoverage*> covs;
	std::unordered_map<std::string, SnowTools::STCoverage*> clip_covs;
	covs["t000"] = &twalk.cov;
	covs["n000"] = &nwalk.cov;
	clip_covs["t000"] = &twalk.cov;
	clip_covs["n000"] = &nwalk.cov;

	std::vector<SnowTools::BreakPoint> allbreaks = ac->getAllBreakPoints(false); // false says dont worry about "local"
	for (auto& i : allbreaks)
	  i.addCovs(covs, clip_covs);
	for (auto& i : allbreaks)
	  i.scoreBreakpoint();
	for (auto& i : allbreaks)
	  i.setRefAlt(findex, nullptr);

	all_align << (*ac) << std::endl;
	bool no_reads = true;
	for (auto& i : allbreaks)
	  os_allbps << i.toFileString(no_reads) << std::endl;

	std::stringstream ss;
	ss << SnowTools::displayRuntime(start);
	std::cerr << "..read " << SnowTools::AddCommas(bav_this.size()) << " reads from " << rsize << " regions for contig " << SnowTools::AddCommas(count) << " " << ss.str() << std::endl;
	*/
	
      }

      // prepare the next one
      brv.clear();
      brv.push_back(r);
      curr_name = r.Qname();
      regions.clear();
      SnowTools::GenomicRegion gr = r.asGenomicRegion();
      gr.pad(1000);
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

  std::cerr << "...done adding " << SnowTools::AddCommas(tmp_queue.size()) << " contigs. Firing up " << 
    numThreads << " re-alignment threads" << std::endl;
  num_to_run = tmp_queue.size();
  for (auto& i : tmp_queue) 
    queue.add(i);
  
  // wait for the threads to finish
  for (int i = 0; i < numThreads; i++) 
    threadqueue[i]->join();  

  all_align.close();
  os_allbps.close();

}

