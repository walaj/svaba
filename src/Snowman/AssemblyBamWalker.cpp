#include "AssemblyBamWalker.h"

#include "BamToolsUtils.h"

void AssemblyBamWalker::walkDiscovar()
{
  Read r;
  std::cout << "...starting to walk Discovar BAM" << std::endl;
  //std::string rule;

  std::vector<AlignedContig> ac_vec;
  
  bool rule;
  while(GetNextRead(r, rule))
    {
      BamTools::BamAlignment a;

      a.Name = r_qname(r);
      a.AlignmentFlag = r_flag(r);
      a.Position = r_pos(r);
      //a.MatePosition = r_mpos(r);
      a.RefID = r_id(r);
      //a.MateRefID = r_mid(r);
      r_seq(r, a.QueryBases);
      a.MapQuality = r_mapq(r);
      
      // parse the cigar
      for (size_t i = 0; i < r_cig_size(r); ++i) 
	{
	  BamTools::CigarOp op;
	  op.Type = r_cig_type(r, i);
	  op.Length = r_cig_len(r, i);
	  a.CigarData.push_back(op);
	}

      // add the primary alignment
      AlignedContig ac(a);
      ac_vec.push_back(ac);

      // parse the XP field
      std::string xp;
      r_get_Z_tag(r, "XP", xp);
      std::vector<BamTools::BamAlignment> bv;
      BamToolsUtils::parseXP2BamAlignment(xp, bv);
      for (auto& i : bv) {
	i.Name = a.Name;
	if (a.IsReverseStrand() == i.IsReverseStrand()) {
	  i.QueryBases = a.QueryBases;
	} else {
	  std::string qb = a.QueryBases;
	  SnowTools::rcomplement(qb);
	  i.QueryBases = qb;
	}
	
	ac_vec.back().addAlignment(i);
      }
      ac_vec.back().setMultiMapBreakPairs();
      if (ac_vec.back().hasVariant() && ac_vec.back().getMaxMapq() >= 55)
	  std::cout << ac_vec.back() << std::endl;
    }
}
