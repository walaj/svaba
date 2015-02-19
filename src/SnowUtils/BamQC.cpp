#include "BamQC.h"
#include "VariantBamReader.h"

template <typename T> void printQCVec(ostream &out, const vector<T> &vec) {
  for (auto it = vec.begin(); it != vec.end(); it++)
    out << "," << *it;
  return;
}

// instantiate a new read group
BamQCReadGroup::BamQCReadGroup() {

  mapq = vector<size_t>(61,0); // 60 elems of 0
  nm   = vector<size_t>(102,0); // 101 elems of 0
  isize= vector<size_t>(2001,0); // (everything above 2000 is inter)
  clip = vector<size_t>(102,0); 
  as   = vector<size_t>(102,0);
  xp   = vector<size_t>(102,0);
  len  = vector<size_t>(102,0);
  phred= vector<size_t>(61, 0);
}

// make the output
std::ostream& operator<<(std::ostream& out, const BamQC& qc) {

  for (auto it = qc.map.begin(); it != qc.map.end(); it++)
    out << "READGROUP:" << it->first << endl << it->second << endl;
  return out;
}

// make the output
std::ostream& operator<<(std::ostream& out, const BamQCReadGroup& rg) {

  out << "total," << rg.num_reads << endl;
  out << "unmap," << rg.unmap << endl;
  out << "qcfail," << rg.qcfail << endl;
  out << "duplicate," << rg.duplicate << endl;
  out << "supplementary," << rg.supp << endl;

  out << "mapq";
  printQCVec<size_t>(out, rg.mapq);
  out << endl;

  out << "nm";
  printQCVec<size_t>(out, rg.nm);
  out << endl;

  out << "isize";
  printQCVec<size_t>(out, rg.isize);
  out << endl;

  out << "as";
  printQCVec<size_t>(out, rg.as);
  out << endl;
  
  out << "xp";
  printQCVec<size_t>(out, rg.xp);
  out << endl;

  out << "clip";
  printQCVec<size_t>(out, rg.clip);
  out << endl;

  out << "len";
  printQCVec<size_t>(out, rg.len);
  out << endl;

  out << "phred";
  printQCVec<size_t>(out, rg.phred);

  return out;
}

// add an addional alignment
void BamQC::addRead(BamTools::BamAlignment &a) {

  try {

    if (a.Name == "")
      	a.BuildCharData();
      string rgroup;
      if (!a.GetTag("RG",rgroup))
	cerr << "Failed to read rgroup" << endl;
      
      int this_isize = a.InsertSize;
      this_isize = (a.MateRefID != a.RefID || this_isize > 2000) ? 2000 : this_isize;
      
      // get clip num
      unsigned clipnum = VariantBamReader::getClipCount(a);
      
      // get the mean phred quality
      size_t i = 0;
      int phred = 0;
      while(i < a.Qualities.length()) {
        phred += char2phred(a.Qualities[i]);
	i++;
      }
      if (a.Qualities.length() > 0)
	phred = static_cast<int>(floor(static_cast<float>(phred) / a.Qualities.length()));
      
      // get the NM tag
      uint32_t nm;
      if (a.GetTag("NM", nm)) {} else { nm = 0; }
      
      assert(a.MapQuality <= 60);
      assert(clipnum <= 101);
      //assert(as <= 101);
      //assert(xp <= 101);
      assert(a.Length <= 101 && a.Length  >= 0);
      assert(phred <= 60 && phred  >= 0);
      assert(nm <= 101);
      
      // discordant
      bool FR_f = !a.IsReverseStrand() && (a.Position < a.MatePosition) && (a.RefID == a.MateRefID) &&  a.IsMateReverseStrand();
      bool FR_r =  a.IsReverseStrand() && (a.Position > a.MatePosition) && (a.RefID == a.MateRefID) && !a.IsMateReverseStrand();
      bool FR = FR_f || FR_r;
      if (a.InsertSize > 0 && this_isize <= 2000 && a.IsPaired() && FR ) // only count "proper" reads 
	map[rgroup].isize[this_isize]++;
      
      // all the rest
      //map[rgroup].xp[xp]++;
      map[rgroup].len[a.Length]++;
      //map[rgroup].as[as]++;
      map[rgroup].clip[clipnum]++;
      map[rgroup].phred[phred]++;
      map[rgroup].num_reads++;
      map[rgroup].mapq[a.MapQuality]++;
      map[rgroup].nm[nm]++;
      
      if (!a.IsMapped())
	map[rgroup].unmap++;
      if (a.IsFailedQC()) 
	map[rgroup].qcfail++;
      if (a.IsDuplicate())
	map[rgroup].duplicate++;
      if (!a.IsPrimaryAlignment())
	map[rgroup].supp++;
      
  } catch (...) {
    cerr << "Failed at adding to QC" << endl;
    //cerr << "Readgroup " << "NM " << nm << " mapq " << a.MapQuality << 
    //" xp " << xp << 
    //" len " << a.Length <<
    //" as " << as << 
    //" phred " << phred << endl;
  }
  
}
