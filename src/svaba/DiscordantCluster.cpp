#include "DiscordantCluster.h"

#include <set>
#include <cassert>
#include <numeric>
#include <unordered_set>

#include "SeqLib/BamWalker.h"
#include "svabaUtils.h"

#define DISC_PAD 150
#define MIN_PER_CLUST 2
#define MAX_CLUSTER_WIDTH 3000
#define MAX_REGION_WIDTH 5000 
#define DEFAULT_ISIZE_THRESHOLD 2000 // shouldn't be hit if isize was learned

//#define DEBUG_CLUSTER 1

namespace SeqLib {
  class BamHeader;
  class GenomicRegion;
}

using SeqLib::GenomicRegion;

DiscordantClusterMap DiscordantCluster::clusterReads(svabaReadPtrVector& bav,
						     const SeqLib::GenomicRegion& interval,
						     const SeqLib::BamHeader& header) {
  
  
  // sort by position
  std::sort(bav.begin(), bav.end(), SeqLib::BamRecordSort::ByReadPositionSharedPtr());
  
  svabaReadClusterVector fwd, rev, fwdfwd, revrev, fwdrev, revfwd;

  //debug
  // for (const auto& b : bav) {
  //   if (b->Qname() == "LH00306:129:227V5CLT4:6:2114:6074:25724")
  // 	std::cerr << " CLUSTER " << *b << std::endl;
  // }
  
  // make the fwd and reverse READ clusters. dont consider mate yet
  __cluster_reads(bav, fwd, rev, SeqLib::Orientation::FR);
  __cluster_reads(bav, fwd, rev, SeqLib::Orientation::FF);
  __cluster_reads(bav, fwd, rev, SeqLib::Orientation::RF);
  __cluster_reads(bav, fwd, rev, SeqLib::Orientation::RR);
  
#ifdef DEBUG_CLUSTER
  for (auto& i : fwd) {
    std::cerr << "fwd cluster " << std::endl;
    for (auto& j : i)
      std::cerr << "fwd " << *j << std::endl;
  }
  for (auto& i : rev) {
    std::cerr << "rev cluster " << std::endl;
    for (auto& j : i)
      std::cerr << "rev " << *j << std::endl;
  }
#endif
  
  // within the forward read clusters, cluster mates on fwd and rev
  __cluster_mate_reads(fwd, fwdfwd, fwdrev);
  
  // within the reverse read clusters, cluster mates on fwd and rev
  __cluster_mate_reads(rev, revfwd, revrev); 

  // we have the reads in their clusters. Just convert to discordant reads clusters
  DiscordantClusterMap dcm;
  __convertToDiscordantCluster(dcm, fwdfwd, bav, header);
  __convertToDiscordantCluster(dcm, fwdrev, bav, header);
  __convertToDiscordantCluster(dcm, revfwd, bav, header);
  __convertToDiscordantCluster(dcm, revrev, bav, header);

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
	std::cerr << "FWDREV _____ " << *jj << std::endl;
    }
    for (auto& ii : revfwd) {
      std::cerr << " ____________ CLUSTER ______________" << std::endl;
      for (auto& jj : ii)
	std::cerr << "FWDREV _____ " << *jj << std::endl;
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
  DiscordantCluster::DiscordantCluster(const svabaReadPtrVector& this_reads,
				       const svabaReadPtrVector& all_reads,
				       const SeqLib::BamHeader& header) {

    assert(this_reads.size());
    assert(all_reads.size());

    // check the orientations, fill the reads
    auto first_read = this_reads.at(0);
    auto last_read = this_reads.back();
    bool rev = first_read->ReverseFlag();
    bool mrev = first_read->MateReverseFlag();

    // loop the reads in this cluster
    for (auto& i : this_reads) {
      
      // double check that we did the clustering correctly. All read orientations should be same
      assert(rev == i->ReverseFlag() && mrev == i->MateReverseFlag()); 
      
      // add the read to the read map
      std::string tmp = i->UniqueName(); // this name like t001_165_qname
      assert(tmp.length());
      reads[tmp] = i; 
      
      // count number of reads per BAM
      ++counts[i->Prefix()];
      
      if (tmp.at(0) == 't') {
	++tcount;
      } else {
	++ncount;
      }
      
    }

    // loop through the big stack of reads and find the mates
    addMateReads(all_reads);
    assert(reads.size());

    // set the regions
    m_reg1 = GenomicRegion(-1, -1, -1);
    m_reg1.pos1 = INT_MAX;
    m_reg2 = GenomicRegion(-1, -1, -1); // mate region
    m_reg2.pos1 = INT_MAX;
    SeqLib::Orientation por;
    for (const auto& [_, read] : reads) 
      {
        por = read->PairOrientation();
        m_reg1.strand = read->ReverseFlag()     ? '-' : '+'; 
        m_reg2.strand = read->MateReverseFlag() ? '-' : '+'; 	
        m_reg1.chr = read->ChrID(); 
        m_reg2.chr = read->MateChrID();

        // get left side 
        if (read->Position() < m_reg1.pos1)
          m_reg1.pos1 = read->Position();
        if (read->MatePosition() < m_reg2.pos1)
          m_reg2.pos1 = read->MatePosition();
        
        // get right side
        if (read->PositionEnd() > m_reg1.pos2)
          m_reg1.pos2 = read->PositionEnd();
        if (read->MatePosition() > m_reg2.pos2) 
          m_reg2.pos2 = read->MatePosition() + read->Length(); // since don't have mate end
        
        // Check region width constraints after building regions
        if (m_reg1.Width() >= MAX_REGION_WIDTH) {
          std::cerr << "Warning: Region 1 width (" << m_reg1.Width() << "bp) exceeds " << MAX_REGION_WIDTH << "bp limit. Cluster is too large." << std::endl;
        }
        if (m_reg2.Width() >= MAX_REGION_WIDTH) {
          std::cerr << "Warning: Region 2 width (" << m_reg2.Width() << "bp) exceeds " << MAX_REGION_WIDTH << "bp limit. Cluster is too large." << std::endl;
        }
        assert(m_reg1.Width() < 5000);
        assert(m_reg2.Width() < 5000);	
      }
  
    mapq1 = __getMeanMapq(reads);
    mapq2 = __getMeanMapq(mates);
    assert(mapq1 >= 0);
    assert(mapq2 >= -1); // can have -1 as placeholder if no mate reads (bc didnt do lookup)
    
    // orient them correctly so that left end is first
    if (m_reg2 < m_reg1) {
      std::swap(m_reg1, m_reg2);
      std::swap(reads, mates);
      std::swap(mapq1, mapq2);
    }
    
    // set the ID
    std::string ortid;

    switch (por) {
    case SeqLib::Orientation::FR: ortid = "FR"; break;
    case SeqLib::Orientation::FF: ortid = "FF"; break;
    case SeqLib::Orientation::RF: ortid = "RF"; break;
    case SeqLib::Orientation::RR: ortid = "RR"; break;
    default: ortid = "UD"; break;  // fallback for UD or anything unexpected
    }    
    
    m_id = ortid + "_" + m_reg1.ChrName(header) + "_" + std::to_string(m_reg1.pos1) +
      "___" + m_reg2.ChrName(header) + "_" + std::to_string(m_reg2.pos1);
    
  }

void DiscordantCluster::addMateReads(const svabaReadPtrVector& bav) 
  { 
    
    if (!reads.size())
      return;
    
    // log the qnames
    // so qnames is just the same map as reads except:
    // -- qnames: qname : svabaRead
    // -- reads: uniquesame : svabaRead
    // the qnames structure is useful here because when we search the pile
    // of reads for the mate, the read and mate should have same qname
    // but then it's good to have the read when we are looking at the mate
    // to check some things (e.g. read and mate are different members of pair)
    std::unordered_map<std::string, svabaReadPtr> qnames;
    for (auto& [_,read] : reads) 
      qnames.insert({read->Qname(), read});

    // OK, so this is necessary because...
    // we are looping through and trying to fill mate reads by comparing
    // against qname of reads. BUT this can be problematic if there are secondary 
    // or supplementarty alignments. So use this region to check that mate ALSO agrees
    // with mate regions
    auto& first_read = reads.begin()->second;
    GenomicRegion g(first_read->MateChrID(), first_read->MatePosition(),
		    first_read->MatePosition()); 
    g.Pad(DISC_PAD + 1000);

    // loop all of the reads and find the mates,
    // if we have them (don't have to)

    for (const auto& i : bav) {
      
      // first check if this read has same name as one already stored in reads
      auto it_read = qnames.find(i->Qname());
      if (it_read == qnames.end()) // qname not found, so this read can't be a matching mate
	continue;
      
      // now, check if this is the same read as we already now about (in reads)
      // or if it is not (the its the mate)
      std::string tmp = i->UniqueName();
      auto it = reads.find(tmp);
      if (it != reads.end()) // shoot, this it just the read as "read", so we're done
	continue;
      
      // now check that the read and this one (i = candidate mate) are opposite in pair
      if (i->FirstFlag() == it_read->second->FirstFlag())
	continue;

      // check that the orientation 
      if (i->PairOrientation() != it_read->second->PairOrientation())
	continue;

      // don't want to deal with these
      if (i->SecondaryFlag() || it_read->second->SecondaryFlag())
	continue;
      
      // agrees with intiial mate orientation and position      
      if (g.GetOverlap(i->AsGenomicRegion()) > 0)
	mates[tmp] = i;
      
    }
  }

double DiscordantCluster::__getMeanMapq(const DiscordantReadMap& m) const {
  if (m.empty()) return -1.0;
  // sum them up
  double total = 0.0;
  for (auto const& kv : m) {
    total += kv.second->MapQuality();
  }
  return total / m.size();
}

void DiscordantCluster::labelReads() {

  for (auto& [_, r] : reads) {
    std::string dcstring;
    r->GetZTag("DC", dcstring);
    if (!dcstring.empty())
      dcstring.append("d");
    dcstring.append(m_id);
    r->AddZTag("DC", dcstring);
  }

  for (auto& [_, r] : mates) {
    std::string dcstring;
    r->GetZTag("DC", dcstring);
    if (!dcstring.empty())
      dcstring.append("d");
    dcstring.append(m_id);
    r->AddZTag("DC", dcstring);
  }
  
  
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
  
  int pos1 = m_reg1.strand == '+' ? m_reg1.pos2 : m_reg1.pos1; // get the edge of the cluster
  int pos2 = m_reg2.strand == '+' ? m_reg2.pos2 : m_reg2.pos1;
  
  std::stringstream out;
  out << h.IDtoName(m_reg1.chr) << sep << pos1 << sep << m_reg1.strand << sep 
      << h.IDtoName(m_reg2.chr) << sep << pos2 << sep << m_reg2.strand << sep 
      << tcount << sep << ncount
      << sep << mapq1 << sep 
      << mapq2 << sep << (m_contig.length() ? m_contig : "x") << sep << m_id;
  
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
 * @param cvec Stores the vector of clusters, which grows here
 * @param clust The current cluster that is being added to
 * @param a Read to add to cluster
 * @param mate Flag to specify if we should cluster on mate position instead of read position
 * @return Description of the return value
 */
bool DiscordantCluster::__add_read_to_cluster(svabaReadClusterVector& cvec, 
                                              svabaReadPtrVector& clust,
                                              svabaReadPtr& a,
                                              bool ismate)
{

  // when ismate is false - just look at this read, ignore the mate
  // when ismate is true - just consider the mate of this read, ignore this pairmate
  
  // last info stores the latest position of the "end" of this cluster
  std::pair<int, int> last_info;
  if (clust.empty()) {
    last_info = {-1, -1};
  } else if (ismate) {
    last_info = {clust.back()->MateChrID(), clust.back()->MatePosition()};
  } else {
    last_info = {clust.back()->ChrID(), clust.back()->Position()};
  }

  // what is the position of the input read
  std::pair<int, int> this_info;
  if (ismate)
    this_info = {a->MateChrID(), a->MatePosition()};
  else
    this_info = {a->ChrID(), a->Position()};

  // check if the cluster is getting too wide
  bool too_wide = !__valid_cluster(clust, ismate);

  // if cluster not too wide, and read is on the same chromosome
  // and the difference in position is less than the discordant "coupling"
  // parameter, then add it to this cluster
  if (!too_wide && this_info.first == last_info.first &&
      (this_info.second - last_info.second) <= DISC_PAD) {

    clust.push_back(a);
    return true;

  // if not, then we're done with the cluster. The key is that the
  // input reads are position sorted, so that if this didn't make it, then
  // the rest won't, and we whould start the next cluster with the given read
  } else {
    if (clust.size() >= MIN_PER_CLUST)
      cvec.push_back(clust);
    
    clust.clear();
    clust.push_back(a);
    return false;
  }
}

void DiscordantCluster::__cluster_mate_reads(svabaReadClusterVector& brcv,
                                             svabaReadClusterVector& fwd,
                                             svabaReadClusterVector& rev)
{
  // loop the nascent "clustesr" which are just clustered based on
  // left read, not on the right (mate) read yet
  for (auto& cluster : brcv) {
    
    svabaReadPtrVector this_fwd, this_rev;

    // sort the reads in the 
    std::sort(cluster.begin(), cluster.end(),
	      [](const svabaReadPtr& a, const svabaReadPtr& b) {
		return (a->MateChrID() < b->MateChrID()) ||
		  (a->MateChrID() == b->MateChrID() &&
		   a->MatePosition() < b->MatePosition());
	      });

    // loop the reads in the nascent cluster
    for (svabaReadPtr& r : cluster) {

      assert(r->dd > 0);

      // see __cluster_reads for logic. Exactly the same, just now for mates
      if (!r->MateReverseFlag())
        __add_read_to_cluster(fwd, this_fwd, r, true);
      else
        __add_read_to_cluster(rev, this_rev, r, true);
    }

    // finish the last clusters
    if (__valid_cluster(this_fwd, true) && this_fwd.size() >= MIN_PER_CLUST)
      fwd.push_back(this_fwd);
    if (__valid_cluster(this_rev, true) && this_rev.size() >= MIN_PER_CLUST)
      rev.push_back(this_rev);
  }
}

void DiscordantCluster::__cluster_reads(svabaReadPtrVector& brv,
					svabaReadClusterVector& fwd,
					svabaReadClusterVector& rev,
					SeqLib::Orientation orientation) 
{
  
  // hold the current cluster
  svabaReadPtrVector this_fwd, this_rev;
  
  std::unordered_set<std::string> tmp_set;

  // int pair_orientation = 0;
  // if (brv.size())
  //   pair_orientation = brv.front()->PairOrientation();
  
  // cluster in the READ direction, separately for fwd and rev
  for (svabaReadPtr& i : brv) {

    // don't want to deal with secondary alignments (multiple mappings, but not split)
    if (i->SecondaryFlag())
      continue;
    
    // if (pair_orientation != i->PairOrientation()) {
    //   std::cerr << *i << std::endl;
    //   std::cerr << i->PairOrientation() << std::endl;
    //   std::cerr << pair_orientation << std::endl;
    // }

    // not a discordant read, ignore it
    if (i->dd <= 0)
      continue;
    
    // only cluster FR reads together, RF reads together, FF together and RR together
    if (i->PairOrientation() != orientation) {
      continue;
    }
      
    std::string qq = i->Qname();
    
    // only cluster if not seen before (e.g. left-most is READ, right most is MATE)
    if (i->PairMappedFlag() && tmp_set.count(qq) == 0) {
      
      tmp_set.insert(qq);
      
      // forward clustering -- this_fwd is just the most current cluster and
      // fwd (or rev) is the collection of clusters that grows here
      if (!i->ReverseFlag()) 
	__add_read_to_cluster(fwd, this_fwd, i, false);
      // reverse clustering 
      else 
	__add_read_to_cluster(rev, this_rev, i, false);
    }
  }
  
  // finish the last clusters
  if (__valid_cluster(this_fwd, false) && this_fwd.size() >= MIN_PER_CLUST)
    fwd.push_back(this_fwd);
  if (__valid_cluster(this_rev, false) && this_rev.size() >= MIN_PER_CLUST)
    rev.push_back(this_rev);
  
}

bool DiscordantCluster::__valid_cluster(svabaReadPtrVector& clust, bool ismate) {

  if (clust.size() < 2)
    return true;
  
  // check if getting too wide
  bool too_wide;
  if (ismate)
    too_wide = (clust.size() > 1 && (clust.back()->MatePosition() - clust[0]->MatePosition()) > MAX_CLUSTER_WIDTH);
  else
    too_wide = (clust.size() > 1 && (clust.back()->Position() - clust[0]->Position()) > MAX_CLUSTER_WIDTH); 

  return !too_wide;
}

void DiscordantCluster::__convertToDiscordantCluster(DiscordantClusterMap &dd,
						     const svabaReadClusterVector& cvec,
						     const svabaReadPtrVector& bav,
						     const SeqLib::BamHeader& header) {
  
  // nothign to do if no clusters
  if (cvec.size() == 0)
    return; 

  // always have to be more input reads than clusters
  assert(bav.size() > cvec.size());
  
  // no reads, nothign to do
  if (!bav.size())
    return;

  // loop through the clusters
  for (auto& v : cvec) {
    DiscordantCluster d(v, bav, header);
    dd[d.m_id] = d;
    
    // Warn if widths are too large
    if (d.m_reg1.Width() >= MAX_REGION_WIDTH || d.m_reg2.Width() >= MAX_REGION_WIDTH) {
      std::cerr << "Warning: Wide cluster " << d.m_id << " with region widths "
                << d.m_reg1.Width() << "bp and " << d.m_reg2.Width() << "bp (>= " << MAX_REGION_WIDTH << "bp limit)" << std::endl;
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
 
