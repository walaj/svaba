#include "BamStats.h"

#include <cmath>

using namespace SeqLib;

//#define DEBUG_STATS 1

BamReadGroup::BamReadGroup(const std::string& name) : reads(0), supp(0), unmap(0), qcfail(0), 
						      duplicate(0), mate_unmap(0), m_name(name)
{

  mapq = Histogram(0,100,1);
  nm = Histogram(0,100,1);
  isize = Histogram(-2,2000,10);
  clip = Histogram(0,100,5);
  phred = Histogram(0,100,1);
  len = Histogram(0,250,1);

}

  std::ostream& operator<<(std::ostream& out, const BamStats& qc) {
    out << "ReadGroup\tReadCount\tSupplementary\tUnmapped\tMateUnmapped\tQCFailed\tDuplicate\tMappingQuality\tNM\tInsertSize\tClippedBases\tMeanPhredScore\tReadLength" << std::endl;
    for (auto& i : qc.m_group_map)
      out << i.second << std::endl;
    return out;
  }

  std::ostream& operator<<(std::ostream& out, const BamReadGroup& qc) {
    std::string sep = "\t";
    out << qc.m_name << sep << qc.reads << sep <<
      qc.supp << sep << 
      qc.unmap << sep <<
      qc.mate_unmap << sep <<
      qc.qcfail << sep <<
      qc.duplicate << sep << 
      qc.mapq.toFileString() << sep << 
      qc.nm.toFileString() << sep <<
      qc.isize.toFileString() << sep <<
      qc.clip.toFileString() <<  sep << 
      qc.phred.toFileString() << sep <<
      qc.len.toFileString();
    return out;
  }

void BamReadGroup::addRead(BamRecord &r)
{
  
  ++reads;
  if (r.SecondaryFlag())
    ++supp;
  if (r.QCFailFlag())
    ++qcfail;
  if (r.DuplicateFlag())
    ++duplicate;
  if (!r.MappedFlag())
    ++unmap;
  if (!r.MateMappedFlag())
    ++mate_unmap;

  int mapqr = r.MapQuality();
  if (mapqr >=0 && mapqr <= 100)
    mapq.addElem(mapqr);
  
  int32_t this_nm = r.GetIntTag("NM");;
  //r_get_int32_tag(r, "NM", this_nm);
  if (this_nm <= 100)
    nm.addElem(this_nm);
  
  int32_t isizer = -1;
  if (!r.PairMappedFlag())
    isizer = -2;
  else if (!r.Interchromosomal())
    isizer = std::abs(r.InsertSize());
  isize.addElem(isizer);

  int32_t c = r.NumClip();
  //r_get_clip(r,c);
  clip.addElem(c);
  
  len.addElem(r.Length());

  phred.addElem((int)r.MeanPhred());

}

void BamStats::addRead(BamRecord &r)
{

  // get the read group
  std::string rg = r.GetZTag("RG"); 
  if (rg.empty()) // try grabbing from QNAME
    rg = "QNAMED_" + r.ParseReadGroup();

#ifdef DEBUG_STATS
  std::cout << "got read group tag " << rg << std::endl;
#endif
  
  std::unordered_map<std::string, BamReadGroup>::iterator ff = m_group_map.find(rg);
  
  if (ff == m_group_map.end())
    {
      m_group_map[rg] = BamReadGroup(rg);
      m_group_map[rg].addRead(r);
    } 
  else
    {
      ff->second.addRead(r);
    }
}

