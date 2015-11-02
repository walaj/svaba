#include "AssemblyBamWalker.h"
#include "SnowmanUtils.h"
#include "run_snowman.h"

#include <chrono>
#include <thread>

static pthread_mutex_t snow_lock;
static ogzstream all_align, os_allbps;
static struct timespec start;

static int queue_size = 0;

bool runAC(SnowTools::BamReadVector& brv, faidx_t * f, std::shared_ptr<hts_idx_t> pt, std::shared_ptr<hts_idx_t> pn,
	   const std::string& t, const std::string& n, const SnowTools::GenomicRegionVector& regions) {
  
  SnowmanBamWalker twalk(t);
  SnowmanBamWalker nwalk(n);
  twalk.prefix = "t000";
  nwalk.prefix = "n000";
  twalk.max_cov = 500;
  nwalk.max_cov = 500;
  twalk.get_mate_regions = false;
  nwalk.get_mate_regions = false;
  twalk.setBamWalkerRegions(regions, pt);
  nwalk.setBamWalkerRegions(regions, pn);

  twalk.readBam();
  nwalk.readBam();

  SnowTools::BamReadVector bav_this;
  for (auto& q : twalk.reads)
    bav_this.push_back(q);
  for (auto& q : nwalk.reads)
    bav_this.push_back(q);

  std::vector<SnowTools::AlignedContig> this_alc;


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
    i.setRefAlt(f, nullptr);

  // MUTEX LOCKED
  ////////////////////////////////////
  pthread_mutex_lock(&snow_lock);  

  all_align << (*ac) << std::endl;
  bool no_reads = true;
  for (auto& i : allbreaks)
    os_allbps << i.toFileString(no_reads) << std::endl;
  
  std::stringstream ss;
  ss << SnowTools::displayRuntime(start);
  std::cerr << "..read " << SnowTools::AddCommas(bav_this.size()) << " reads from " << regions.size() << " regions " << ss.str() << std::endl;

  ////////////////////////////////////
  // MUTEX UNLOCKED
  ////////////////////////////////////
  pthread_mutex_unlock(&snow_lock);
  
  return true;

}


void AssemblyBamWalker::walkDiscovar()
{

  SnowmanUtils::fopen("assembly.alignments.txt.gz", all_align);
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
  

  std::vector<AssemblyWalkerWorkItem*> tmp_queue;

  bool rule;
  int count = 0;
  std::string curr_name;

  int max_mapq = 0;
  while(GetNextRead(r, rule)) {
    if (++count % 50000 == 0) 
      std::cerr << "...working on contig " << SnowTools::AddCommas(count) << std::endl;

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
      brv.push_back(r);
      SnowTools::GenomicRegion gr = r.asGenomicRegion();
      gr.pad(400);
      regions.push_back(gr);
      max_mapq = std::max(r.MapQuality(), max_mapq);
    } else {
      if (/*brv[0].CigarSize() > 1 || */brv.size() > 1 && brv.size() < 20 && regions.size() < 20 && max_mapq >= 50) {

	AssemblyWalkerWorkItem * item     = new AssemblyWalkerWorkItem(brv, count, findex, tindex, nindex, regions, 
								       tbam, nbam);
	queue.add(item);

	regions.clear();
	/// hack
	while (queue.size() > 50) {
	  //std::cerr << queue.size() << std::endl;
	  std::this_thread::sleep_for(std::chrono::milliseconds(1));
	}

	std::cerr << " count " << SnowTools::AddCommas(count) << std::endl;	
	if (count > 5e6) {
	  all_align.close();
	  os_allbps.close();
	  return; //exit(0);
	  break;
	}


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
      max_mapq = 0;
    }
  }

  // wait for the threads to finish
  for (int i = 0; i < numThreads; i++) 
    threadqueue[i]->join();  

  all_align.close();
  os_allbps.close();

}

