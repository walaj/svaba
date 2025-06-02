#include "DiscordantCluster.h"

#include <set>
#include <cassert>
#include <numeric>
#include <unordered_set>

#include "SeqLib/BamWalker.h"
#include "svabaUtils.h"

#define DISC_PAD 150
#define MIN_PER_CLUST 2
#define DEFAULT_ISIZE_THRESHOLD 2000 // shouldn't be hit if isize was learned

//#define DEBUG_CLUSTER 1

namespace SeqLib {
  class BamHeader;
  class GenomicRegion;
}

using SeqLib::GenomicRegion;

DiscordantClusterMap DiscordantCluster::clusterReads(
						     svabaReadVector& bav,
						     const SeqLib::GenomicRegion& interval,
						     int max_mapq_possible) {
  
  
#ifdef DEBUG_CLUSTER    
  //for (auto& i : bav)
  //  std::cerr << " PRE DEDUPED CLUSTER " << i << std::endl;
#endif
  
  // sort by position
  std::sort(bav.begin(), bav.end(), SeqLib::BamRecordSort::ByReadPosition());
  
#ifdef DEBUG_CLUSTER    
  //for (auto& i : bav_dd)
  //  std::cerr << " DEDUPED CLUSTER " << i << std::endl;
#endif
  
  // clear the tmp map. Now we want to use it to store if we already clustered read
  //tmp_map.clear();
  
  svabaReadClusterVector fwd, rev, fwdfwd, revrev, fwdrev, revfwd;
  
  // make the fwd and reverse READ clusters. dont consider mate yet
  __cluster_reads(bav, fwd, rev, FRORIENTATION);
  __cluster_reads(bav, fwd, rev, FFORIENTATION);
  __cluster_reads(bav, fwd, rev, RFORIENTATION);
  __cluster_reads(bav, fwd, rev, RRORIENTATION);
  
  // remove singletons
  __remove_singletons(fwd);
  __remove_singletons(rev);
  
#ifdef DEBUG_CLUSTER
  for (auto& i : fwd) {
    std::cerr << "fwd cluster " << std::endl;
    for (auto& j : i)
      std::cerr << "fwd " << j << std::endl;
  }
  for (auto& i : rev) {
    std::cerr << "rev cluster " << std::endl;
    for (auto& j : i)
      std::cerr << "rev " << j << std::endl;
  }
#endif
  
  // within the forward read clusters, cluster mates on fwd and rev
  __cluster_mate_reads(fwd, fwdfwd, fwdrev);
  
  // within the reverse read clusters, cluster mates on fwd and rev
  __cluster_mate_reads(rev, revfwd, revrev); 
  
  // remove singletons
  __remove_singletons(fwdfwd);
  __remove_singletons(revfwd);
  __remove_singletons(fwdrev);
  __remove_singletons(revrev);
  
  // we have the reads in their clusters. Just convert to discordant reads clusters
  DiscordantClusterMap dcm;
  __convertToDiscordantCluster(dcm, fwdfwd, bav, max_mapq_possible);
  __convertToDiscordantCluster(dcm, fwdrev, bav, max_mapq_possible);
  __convertToDiscordantCluster(dcm, revfwd, bav, max_mapq_possible);
  __convertToDiscordantCluster(dcm, revrev, bav, max_mapq_possible);
  
#ifdef DEBUG_CLUSTER
  std::cerr << "----fwd cluster count: " << fwd.size() << std::endl;
  std::cerr << "----rev cluster count: " << rev.size() << std::endl;
    std::cerr << "----fwdfwd cluster count: " << fwdfwd.size() << std::endl;
    std::cerr << "----fwdrev cluster count: " << fwdrev.size() << std::endl;
    std::cerr << "----revfwd cluster count: " << revfwd.size() << std::endl;
    std::cerr << "----revrev cluster count: " << revrev.size() << std::endl;

    for (auto& ii : fwdrev) {
      std::cerr << " ____________ CLUSTER ______________" << std::endl;
      for (auto& jj : ii)
	std::cerr << "FWDREV _____ " << jj << std::endl;
    }
    for (auto& ii : revfwd) {
      std::cerr << " ____________ CLUSTER ______________" << std::endl;
      for (auto& jj : ii)
	std::cerr << "FWDREV _____ " << jj << std::endl;
    }

#endif

    // // score by number of maps
    // for (auto d : dd_clean) {
    //   for (auto& r : d.second.reads) {
    // 	double rr = r.second.GetDD();
    // 	d.second.read_score += (rr > 0) ? 1/rr : 1;
    //   }
    //   for (auto& r : d.second.mates) {
    // 	double rr = r.second.GetDD();
    // 	d.second.mate_score += (rr > 0) ? 1/rr : 1;
    //   }
    // }

    return dcm;
    
  }
  
  // this reads is reads in the cluster. all_reads is big pile where all the clusters came from
  DiscordantCluster::DiscordantCluster(const svabaReadVector& this_reads,
				       const svabaReadVector& all_reads,
				       int max_mapq_possible) {
    
    if (this_reads.size() == 0)
      return;
    
    // check the orientations, fill the reads
    bool rev = this_reads[0].ReverseFlag();
    bool mrev = this_reads[0].MateReverseFlag();
    
    // the ID of the discordant cluster is just the first reads Qname
    m_id = this_reads[0].Qname();
    assert(m_id.length());
    
    assert(this_reads.back().MatePosition() - this_reads[0].MatePosition() < 10000);
    assert(this_reads.back().Position() - this_reads[0].Position() < 10000);
    std::vector<int> isizer;
    
    // // get distribution of isizes. Reject outliers

    // for (auto& i : this_reads)
    //   isizer.push_back(std::abs(i.FullInsertSize()));
    // double sd = 0, mm = 0, medr = 0;
    // std::sort(isizer.begin(), isizer.end());
    // if (isizer.size() >= 5 && isizer.back() - isizer[0] > 400) {

    //   medr = svabaUtils::CalcMHWScore(isizer);

    //   // get the mean
    //   double sum = std::accumulate(isizer.begin(), isizer.end(), 0.0);
    //   mm = isizer.size() > 0 ? sum / isizer.size() : 0;

    //   // get isize stdev
    //   std::vector<double> diff(isizer.size());
    //   std::transform(isizer.begin(), isizer.end(), diff.begin(), [mm](double x) { return x - mm; });
    //   double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    //   sd = std::sqrt(sq_sum / isizer.size());      
    // }

    // double min_isize = 0;
    // if (sd > 0) { // ? (medr / sd) < 4 : false) {
    //   min_isize = medr - 2 * sd;

    //   //std::cout << (*this ) << std::endl;
    //   //std::cout << " MED " << medr << " SD " << sd << " MEAN " << mm << " CUTOFF " << min_isize << " medr/sd " << (medr/sd) << std::endl;;
    // }
    // min_isize = min_isize < 1000 ? min_isize : 0;
    
    for (auto& i : this_reads) {
	// double check that we did the clustering correctly. All read orientations should be same
	assert(rev == i.ReverseFlag() && mrev == i.MateReverseFlag()); 

	// add the read to the read map
	std::string tmp = i.SR(); // this name like t001_165_qname
	assert(tmp.length());
	reads[tmp] = i;

	// count number of reads per BAM
	++counts[i.Prefix()];
	
	// the ID is the lexographically lowest qname
	std::string qn = i.Qname();
	if (qn < m_id)
	  m_id = qn;

	if (tmp.at(0) == 't') {
	  ++tcount;
	} else {
	  ++ncount;
	}

	isizer.push_back(i.FullInsertSize());

      }
       
    // loop through the big stack of reads and find the mates
    addMateReads(all_reads);
    assert(reads.size());

    // set the regions
    m_reg1 = GenomicRegion(-1, -1, -1);
    m_reg1.pos1 = INT_MAX;
    m_reg2 = GenomicRegion(-1, -1, -1); // mate region
    m_reg2.pos1 = INT_MAX;
    for (auto& i : reads) 
      {
	m_reg1.strand = i.second.ReverseFlag() ? '-' : '+'; //r_strand(i.second) == '+'; //(!i.second->IsReverseStrand()) ? '+' : '-';
	m_reg1.chr = i.second.ChrID(); //r_id(i.second); //i.second->RefID;
	if (i.second.Position() < m_reg1.pos1)
	  m_reg1.pos1 = i.second.Position(); //r_pos(i.second); //i.second->Position;
	int endpos = i.second.PositionEnd(); //r_endpos(i.second);
	if (endpos > m_reg1.pos2)
	  m_reg1.pos2 = endpos;
	assert(m_reg1.Width() < 5000);
      }
    
    for (auto& i : mates) 
      {
	m_reg2.strand = i.second.ReverseFlag() ? '-' : '+'; //r_strand(i.second) == '+'; //(!i.second->IsReverseStrand()) ? '+' : '-';
	m_reg2.chr = i.second.ChrID(); //i.second->RefID;
	if (i.second.Position() < m_reg2.pos1)
	  m_reg2.pos1 = i.second.Position();
	int endpos = i.second.PositionEnd();
	if (endpos > m_reg2.pos2)
	  m_reg2.pos2 = endpos;
	assert(m_reg2.Width() < 5000);
      }
    
    // if no mates (we didn't look up), still make the mate region
    if (mates.size() == 0) {
      for (auto& i : reads) {
	m_reg2.strand = i.second.MateReverseFlag() ? '-' : '+'; 
	m_reg2.chr = i.second.MateChrID();
	if (i.second.MatePosition() < m_reg2.pos1)
	  m_reg2.pos1 = i.second.MatePosition();
	int endpos = i.second.MatePosition() + i.second.Length();
	if (endpos > m_reg2.pos2)
	  m_reg2.pos2 = endpos;
	assert(m_reg2.Width() < 5000);
      }
    }

    int HQMAPQ = max_mapq_possible == 60 ? 25 : std::floor(0.70 * (double)max_mapq_possible) ;

    std::unordered_set<std::string> hqq;
    // set which reads are HQ
    for (auto& i : mates) {
      int nm=0;
      i.second.GetIntTag("NM", nm);
      if (i.second.MapQuality() >= HQMAPQ && nm < 3) {
	hqq.insert(i.second.Qname());
      }
    }
    for (auto& i : reads) {
      int nm=0;
      i.second.GetIntTag("NM", nm);
      if (i.second.MapQuality() >= HQMAPQ && hqq.count(i.second.Qname()) && nm < 3) {
	//if(i.second.GetZTag("SR").at(0) == 't')
	if(i.second.Tumor())
	  ++tcount_hq;
	else
	  ++ncount_hq;
      }
    }

    mapq1 = __getMeanMapq(false);
    mapq2 = __getMeanMapq(true);
    assert(mapq1 >= 0);
    assert(mapq2 >= -1); // can have -1 as placeholder if no mate reads (bc didnt do lookup)

    // orient them correctly so that left end is first
    if (m_reg2 < m_reg1) {
      std::swap(m_reg1, m_reg2);
      std::swap(reads, mates);
      std::swap(mapq1, mapq2);
    }

  }
  
  void DiscordantCluster::addMateReads(const svabaReadVector& bav) 
  { 
    
    if (!reads.size())
      return;
    
    // log the qnames
    std::unordered_set<std::string> qnames;
    for (auto& i : reads) 
      qnames.insert(i.second.Qname());

    // get region around one of the reads.
    // OK, so this is necessary because...
    // we are looping through and trying to fill mate reads by comparing
    // against qname of reads. BUT this can be problematic if there are secondary 
    // or supplementarty alignments. So use this region to check that mate ALSO agrees
    // with mate regions
    GenomicRegion g(reads.begin()->second.MateChrID(), reads.begin()->second.MatePosition(), reads.begin()->second.MatePosition()); 
    bool st = reads.begin()->second.MateReverseFlag();
    g.Pad(DISC_PAD + 1000);
    
    for (auto& i : bav) {
      std::string sr;
      if (qnames.count(i.Qname())) {
	std::string tmp = i.SR();
	  if (reads.count(tmp) == 0)  {// only add if this is a mate read
	    if (i.ReverseFlag() == st && g.GetOverlap(i.AsGenomicRegion()) > 0) // agrees with intiial mate orientation and position
	      mates[tmp] = i;
	  }
	}
    }
  }
  
  double DiscordantCluster::__getMeanMapq(bool mate) const 
  {
    double mean = 0;
    std::vector<int> tmapq;
    if (mate) {
      for (auto& i : mates)
	tmapq.push_back(i.second.MapQuality());
    } else {
      for (auto& i : reads)
	tmapq.push_back(i.second.MapQuality());
    }
    
    // if no mates, set to -1
    if (!mates.size() && mate)
      return -1;

    if (tmapq.size() > 0)
      mean = ((double)std::accumulate(tmapq.begin(), tmapq.end(), 0.0)) / ((double)tmapq.size());

    // actually, get the median
    //if (tmapq.size() > 0)
    // mean = CalcMHWScore(tmapq); //std::accumulate(tmapq.begin(), tmapq.end(), 0.0) / tmapq.size();

    return mean;
  }
  
std::string DiscordantCluster::toRegionString(const SeqLib::BamHeader& h) const 
  {
    int pos1 = (m_reg1.strand == '+') ? m_reg1.pos2 : m_reg1.pos1;
    int pos2 = (m_reg2.strand == '+') ? m_reg2.pos2 : m_reg2.pos1;
    
    std::stringstream ss;
    ss << h.IDtoName(m_reg1.chr) << ":" << pos1 << "(" << m_reg1.strand << ")" << "-" << 
      h.IDtoName(m_reg2.chr) << ":" << pos2 << "(" << m_reg2.strand << ")";
    return ss.str();
    
  }
  
// define how to print this
std::string DiscordantCluster::print(const SeqLib::BamHeader& h) const {
    std::stringstream ss;
    ss << toRegionString(h) << " Tcount: " << tcount << 
      " Ncount: "  << ncount << " Mean MAPQ: " 
       << mapq1 << " Mean Mate MAPQ: " << mapq2 << " Valid: " << (valid() ? "TRUE" : "FALSE");
    return ss.str();
  }
  
  bool DiscordantCluster::valid() const {

    // it's OK if the clusters are inter-chromosomal
    if (m_reg1.chr != m_reg2.chr)
      return true;

    // if its RF orientation (<----->) then bc bps are on outside, regions can overlaps
    if (m_reg1.strand == '-' && m_reg2.strand == '+')
      return true;

    // the clusters overlap, doesn't make sense for del and inversion type
    if (m_reg1.pos2 > m_reg2.pos1) 
      return false;

    return true;

  }

  
  // define how to print to file
  std::string DiscordantCluster::toFileString(const SeqLib::BamHeader& h, bool with_read_names) const { 
    
    std::string sep = "\t";
    
    // add the reads names (currently off)
    std::string reads_string;
    if (with_read_names) {
	
	std::unordered_set<std::string> qnset;
	for (auto& i : reads) 
	  {
	    if (qnset.count(i.second.Qname()))
	      continue;
	    std::string tmp = i.second.SR();
	    qnset.insert(i.second.Qname());
	    reads_string += tmp + ",";
	  }
	for (auto& i : mates) {
	    if (qnset.count(i.second.Qname()))
	      continue;
	    std::string tmp = i.second.SR();
	    qnset.insert(i.second.Qname());
	    reads_string += tmp + ",";
	  }

	
	if (reads_string.empty())
	  reads_string = "x";
	else
	  reads_string.pop_back(); // delete last comma
      }
    
    int pos1 = m_reg1.strand == '+' ? m_reg1.pos2 : m_reg1.pos1; // get the edge of the cluster
    int pos2 = m_reg2.strand == '+' ? m_reg2.pos2 : m_reg2.pos1;

    std::stringstream out;
    out << h.IDtoName(m_reg1.chr) << sep << pos1 << sep << m_reg1.strand << sep 
	<< h.IDtoName(m_reg2.chr) << sep << pos2 << sep << m_reg2.strand << sep 
	<< tcount << sep << ncount << sep << tcount_hq << sep << ncount_hq
	<< sep << mapq1 << sep 
	<< mapq2 << sep << (m_contig.length() ? m_contig : "x") << sep << toRegionString(h)
	<< sep << (reads_string.length() ? reads_string : "x");

    return (out.str());
    
  }
  
  // define how to sort theses
  bool DiscordantCluster::operator<(const DiscordantCluster &b) const 
  {
    if (m_reg1.chr < b.m_reg1.chr)
      return true;
    if (m_reg1.pos1 < b.m_reg1.pos1)
      return true;
    return false;
  }
  
  /**
   * Cluster reads by alignment position 
   * 
   * Checks whether a read belongs to a cluster. If so, adds it. If not, ends
   * and stores cluster, adds a new one.
   *
   * @param cvec Stores the vector of clusters, which themselves are vectors of read pointers
   * @param clust The current cluster that is being added to
   * @param a Read to add to cluster
   * @param mate Flag to specify if we should cluster on mate position instead of read position
   * @return Description of the return value
   */
  bool DiscordantCluster::__add_read_to_cluster(svabaReadClusterVector &cvec,
						svabaReadVector &clust,
						const svabaRead &a,
						bool mate) {

    // get the position of the previous read. If none, we're starting a new one so make a dummy
    std::pair<int,int> last_info;
    if (clust.size() == 0)
      last_info = {-1, -1};
    else if (mate)
      last_info = {clust.back().MateChrID(), clust.back().MatePosition()};
    else
      last_info = {clust.back().ChrID(), clust.back().Position()};
    
    // get the position of the current read
    std::pair<int,int> this_info;
    if (mate)
      this_info = {a.MateChrID(), a.MatePosition()};
    else 
      this_info = {a.ChrID(), a.Position()};

    // is this cluster too big? happens if too many discordant reads. Enforce a hard cutoff
    bool too_big;
    if (mate)
      too_big = (clust.size() > 1 && (clust.back().MatePosition() - clust[0].MatePosition()) > 3000);
    else
      too_big = (clust.size() > 1 && (clust.back().Position() - clust[0].Position()) > 3000);      
    
    // note to self. this rarely gets hit?
    //if (too_big) 
    //  std::cerr << "Cluster too big at " << clust[0].Brief() << " to " << clust.back().Brief() << ". Breaking off after 3000bp" << std::endl;

    // check if this read is close enough to the last
    if (!too_big &&  (this_info.first == last_info.first) && (this_info.second - last_info.second) <= DISC_PAD) {

      // read belongs in current cluster, so add
      clust.push_back(a);
      last_info = this_info;
      return true;
      
    // read does not belong to cluster. close this cluster and add to cvec
    } else {
      
      // if enough supporting reads, add as a cluster
      if (clust.size() >= MIN_PER_CLUST) {
	cvec.push_back(clust);
      }
      
      // clear this cluster and start a new one
      clust.clear();
      clust.push_back(a);
      
      return false;
    }
  }

  void DiscordantCluster::__cluster_mate_reads(svabaReadClusterVector& brcv, svabaReadClusterVector& fwd, svabaReadClusterVector& rev)
  {
    // loop through the clusters, and cluster within clusters based on mate read
    for (auto& v : brcv) 
      {
	
	svabaReadVector this_fwd, this_rev;
	std::sort(v.begin(), v.end(), SeqLib::BamRecordSort::ByMatePosition());

	for (auto& r : v) 
	  {
	    
	    // not a discordant read, ignore it
	    if (r.dd <= 0)
	      continue;

	    // forward clustering
	    if (!r.MateReverseFlag())
	      __add_read_to_cluster(fwd, this_fwd, r, true);
	    // reverse clustering 
	    else if (r.MateReverseFlag()) 
	      __add_read_to_cluster(rev, this_rev, r, true);
	    
	  }
	// finish the last clusters
	if (this_fwd.size() > 0)
	  fwd.push_back(this_fwd);
	if (this_rev.size() > 0)
	  rev.push_back(this_rev);
      } // finish main cluster loop
  }
  
  void DiscordantCluster::__cluster_reads(const svabaReadVector& brv,
					  svabaReadClusterVector& fwd,
					  svabaReadClusterVector& rev,
					  int orientation) 
  {

    // hold the current cluster
    svabaReadVector this_fwd, this_rev;

    std::unordered_set<std::string> tmp_set;

    // cluster in the READ direction, separately for fwd and rev
    for (auto& i : brv) {

      // not a discordant read, ignore it
      if (i.dd <= 0)
	continue;
      
      // only cluster FR reads together, RF reads together, FF together and RR together
      if (i.PairOrientation() != orientation) {
	continue;
      }

      std::string qq = i.Qname();

      // only cluster if not seen before (e.g. left-most is READ, right most is MATE)
      if (i.PairMappedFlag() && tmp_set.count(qq) == 0) {

	tmp_set.insert(qq);

	// forward clustering
	if (!i.ReverseFlag()) 
	  __add_read_to_cluster(fwd, this_fwd, i, false);
	// reverse clustering 
	else 
	  __add_read_to_cluster(rev, this_rev, i, false);
      }
    }

    // finish the last clusters
    if (this_fwd.size() > 0) 
      fwd.push_back(this_fwd);
    if (this_rev.size() > 0)
      rev.push_back(this_rev);

  }

void DiscordantCluster::__convertToDiscordantCluster(DiscordantClusterMap &dd,
						     const svabaReadClusterVector& cvec,
						     const svabaReadVector& bav,
						     int max_mapq_possible) {
  
  for (auto& v : cvec) {
    if (v.size() > 1) { // no clusters with just one read
      DiscordantCluster d(v, bav, max_mapq_possible); /// slow but works (erm, not really slow)
      dd[d.m_id] = d;
    }
  }
}

GenomicRegion DiscordantCluster::GetMateRegionOfOverlap(const GenomicRegion& gr) const {
  
    if (gr.GetOverlap(m_reg1))
      return m_reg2;
    if (gr.GetOverlap(m_reg2))
      return m_reg1;
    return GenomicRegion();

  }

  bool DiscordantCluster::isEmpty() const {
    return m_reg1.IsEmpty() || m_reg2.IsEmpty() || m_reg1.chr == -1 || m_reg2.chr == -1;
  }

  void DiscordantCluster::__remove_singletons(svabaReadClusterVector& b)  {

    // remove cluster with only one read
    b.erase(
	    std::remove_if(b.begin(), b.end(),
			   [](const auto& subvec) { return subvec.size() <= 1; }),
	    b.end());
  }

