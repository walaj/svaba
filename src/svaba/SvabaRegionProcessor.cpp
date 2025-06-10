#include "SvabaRegionProcessor.h"

#include "svabaUtils.h"
#include "svabaOutputWriter.h"

#include "SeqLib/BWAAligner.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/BFC.h"
#include "SeqLib/BamRecord.h"

#include "svabaRead.h"
#include "svabaAssemblerEngine.h"

using SeqLib::GenomicRegion;
using SeqLib::GRC;
using SeqLib::AddCommas;
using SeqLib::BFC;
using SeqLib::BamRecordVector;
using SeqLib::CigarMap;
using SeqLib::UnalignedSequenceVector;
using SeqLib::BWAIndexPtr;
using SeqLib::BWAAligner;
using SeqLib::BWAIndex;

using std::vector;
using std::unordered_map;
using std::unordered_set;
using std::string;

using CountPair    = std::pair<int, int>;

SvabaRegionProcessor::SvabaRegionProcessor(SvabaSharedConfig& sh_cf) : sc(sh_cf)
{ }

void SvabaRegionProcessor::runMateCollectionLoop(const GenomicRegion& region,
						 svabaThreadUnit& stu)
{
  
  // store the newly found bad mate regions
  GRC this_bad_mate_regions; 
  
  CountPair counts = {0,0};
  
  MateRegionVector all_somatic_mate_regions;
  // add the origional, don't want to double back    
  all_somatic_mate_regions.add(MateRegion(region.chr, region.pos1, region.pos2));
  
  // collect together mate regions from control
  MateRegionVector normal_mate_regions;
  for (const auto& w : stu.walkers)
    if (w.first.at(0) == 'n')
      normal_mate_regions.Concat(w.second.mate_regions);
  normal_mate_regions.MergeOverlappingIntervals();
  normal_mate_regions.CreateTreeMap();
  
  // get the mates from somatic 3+ mate regions
  // that don't overlap with normal mate region
  MateRegionVector somatic_mate_regions;
    for (const auto& w : stu.walkers) {
      if (w.first.at(0) == 't')
	for (const auto& i : w.second.mate_regions) {
	  if (i.count >= sc.opts.mateLookupMin &&
	      !normal_mate_regions.CountOverlaps(i) &&
	      !stu.badd.CountOverlaps(i) &&
	      (!sc.blacklist.size() || !sc.blacklist.CountOverlaps(i)))
	    somatic_mate_regions.add(i); 
	}
    }
    somatic_mate_regions.MergeOverlappingIntervals();

    // if none, then we're done
    if (somatic_mate_regions.size() == 0)
      return;
    
    // print out to log
    for (auto& i : somatic_mate_regions)
      sc.logger.log(sc.opts.verbose > 1, sc.opts.verbose_log, 
		    "......mate region ",
		    i.ToString(sc.header) ,
		    " case read count that triggered lookup: ",
		    i.count);

    // convert MateRegionVector to GRC
    GRC somatic_mate_region_collection;
    for (const auto& s : somatic_mate_regions) 
      somatic_mate_region_collection.add(
	     GenomicRegion(s.chr, s.pos1, s.pos2, s.strand));
    somatic_mate_region_collection.MergeOverlappingIntervals();
    somatic_mate_region_collection.CreateTreeMap();

    // convert MateRegionVector to GRC
    GRC gg;
    for (auto& s : somatic_mate_regions) 
      gg.add(GenomicRegion(s.chr, s.pos1, s.pos2, s.strand));
    
    // collect the reads for this round
    CountPair mate_read_counts = {0,0};
    for (auto& w : stu.walkers) {

      svabaBamWalker& sbw = w.second;
      int oreads = sbw.reads.size();
      sbw.m_limit = 5000; //opts.sc.mate_region_lookup_limit;

      if (!sbw.SetRegions(gg))
	continue; // if gg is empty
      
      sbw.get_coverage = false;
      sbw.get_mate_regions = false; 
      
      // clear out the current store of mate regions, since we 
      // already added these to the to-do pile
      sbw.mate_regions.clear();
      
      // read the mate regions
      GRC bad_mate_regions = sbw.readBam();

      // add the bad regions to the total bad region tracker
      stu.badd.Concat(bad_mate_regions);
      stu.badd.MergeOverlappingIntervals();
      stu.badd.CreateTreeMap();

      // update the counts
      if (w.first.at(0) == 't') 
	mate_read_counts.first += (sbw.reads.size() - oreads);
      else
	mate_read_counts.second += (sbw.reads.size() - oreads);
      
    } // end walker read colletion
    
    sc.logger.log(sc.opts.verbose > 1, sc.opts.verbose_log, 
		  "......Mate region reads: <case, control>: <",
	       AddCommas(counts.first), ",",
		  AddCommas(counts.second), ">");
}


bool SvabaRegionProcessor::process(const SeqLib::GenomicRegion& region,
                                   svabaThreadUnit&             unit,
                                   size_t                       threadId)
{

  // count for this unit
  unit.processed_count++;
  unit.processed_since_memory_dump++;
  /*  if (unit.processed_count % 25 == 0) {
    sc.logger.log(true, true, "...processing ", SeqLib::AddCommas(unit.processed_count),
	       " of ", SeqLib::AddCommas(unit.total_count), " for thread ", unit.threadId);
  }
  */
  sc.logger.log(sc.opts.verbose > 1, sc.opts.verbose_log,
		"===Running region ", region.ToString(sc.header),
	      " on thread ", unit.threadId);
  
  // create a new BFC read error corrector for this
  std::shared_ptr<BFC> bfc;
  if (sc.opts.ecCorrectType == "f") {
    bfc = std::make_shared<BFC>();
    for (auto& w : unit.walkers)
      w.second.bfc = bfc;
  }
 
  // setup structures to store the final data for this region
  BamRecordVector all_contigs;
  
  // start a timer
  svabaUtils::svabaTimer st;
  st.start();

  // setup for the BAM walkers
  CountPair read_counts = {0,0}; //tumor, normal

  // holder for Bam: (read : cigar)
  unordered_map<string, CigarMap> cigmap;

  // create a new BFC read error corrector for this
  // its one BFC since correction is across BAMs
  for (auto& w : unit.walkers)
    w.second.bfc = bfc;
  
  // loop all of the BAMs (walkers)
  for (auto& w : unit.walkers) {
    
    svabaBamWalker& sbw = w.second;
    
    // set the region to jump to
    if (!region.IsEmpty()) {
      sbw.SetRegion(region);
    } else { // whole BAM analysis. If region file set, then set regions
      if (sc.file_regions.size()) {
	sbw.SetRegions(sc.file_regions);
      } else { // no regions, literally take every read in bam
	// default is walk all regions, so leave empty to cruise entire BAM
      }
    }

    sc.logger.log(sc.opts.verbose > 1, sc.opts.verbose_log, "---running svabaBamWalker",
		  sbw);
    
    // do the BAM reading, and store the bad mate regions
    SeqLib::GRC bad_regions = sbw.readBam();
    unit.badd.Concat(bad_regions);

    // merge overlapping intervals and remap
    unit.badd.MergeOverlappingIntervals();
    unit.badd.CreateTreeMap();

    // adjust the counts
    if (w.first.at(0) == 't') {
      read_counts.first += w.second.reads.size();
    } else {
      read_counts.second += w.second.reads.size();
    }

    // print if verbose
    sc.logger.log(sc.opts.verbose > 1, sc.opts.verbose_log,
		  "...main region reads <case,control> ",
		  "<",read_counts.first, ",", read_counts.second,">");
    
    // collect cigar strings
    cigmap[w.first] = w.second.cigmap;
    
  }// end the BAM loop
    
  // adjust counts and timer
  st.stop("r");
  
  // get the mate reads, if this is local assembly and has insert-size distro
  if (!region.IsEmpty()) {
    runMateCollectionLoop(region, unit);
    st.stop("m");
  }
  
  // do the discordant read clustering
  sc.logger.log(sc.opts.verbose > 1, false, 
		"...discordant read clustering");

  // tag the reads by discordant status
  svabaReadVector all_discordant_reads;
  for (auto& w : unit.walkers) {
    w.second.TagDiscordantReads();
    for (const auto& r : w.second.reads) {
      if (r.dd > 0)
	all_discordant_reads.push_back(r);
    }
  }

  // for dumping all discordant reads
  if (sc.opts.dump_discordant_reads) {
    unit.all_discordant_reads.insert(unit.all_discordant_reads.end(),
				    all_discordant_reads.begin(),
				    all_discordant_reads.end());
  }
  
  // do the discordant read clustering across BAMs
  DiscordantClusterMap dmap = DiscordantCluster::clusterReads(all_discordant_reads,
							      region,
							      60); // todo, max mapq
  all_discordant_reads.clear();



  // compile all of the raeds together for correction and assembly
  svabaReadVector all_reads_for_assembly;
  for (const auto& w : unit.walkers) {
    all_reads_for_assembly.insert(all_reads_for_assembly.end(),
				  w.second.reads.begin(),
				  w.second.reads.end());
  }

  // for dumping all reads
  if (sc.opts.dump_weird_reads) {
    unit.all_weird_reads.insert(unit.all_weird_reads.end(),
				all_reads_for_assembly.begin(),
				all_reads_for_assembly.end());
  }
  
  // do kmer correction
  //if (sc.opts.ec_correct_type == "s") {
  //  correctReads(all_seqs, input_reads);
  if (sc.opts.ecCorrectType == "f") {
    assert(bfc);
    int learn_reads_count = bfc->NumSequences();
    
    bfc->Train(); // training reads were added during initial retreival
    bfc->clear();  // clear memory and reads. Keeps training data
    
    st.stop("t");

    // reload with the reads to be corrected
    // NB s.CorrectedSeq() doesn't mean already corrected, just
    // this is the quality trimmed version
    for (auto& w : unit.walkers)
      for (const auto& s : w.second.reads)
	bfc->AddSequence(s.CorrectedSeq().c_str(), "", ""); 

    // error correct
    bfc->ErrorCorrect();

    // retrieve the corrected sequences
    string s, name_dum;
    size_t num_reads_corrected = 0;
    for (auto& r : all_reads_for_assembly) {
      assert(bfc->GetSequence(s, name_dum));
      r.SetCorrectedSeq(s);
      ++num_reads_corrected;
    }
    
    double kcov = bfc->GetKCov();
    int kmer    = bfc->GetKMer();

    sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log, "...BFC attempted correct ", num_reads_corrected,
		  " kmer: ", kmer);

    // clear it out, not needed anymore
    if (bfc)
      bfc->clear();
    
    st.stop("k");
  }

  // re-align the corrected reads and dump, just for debugging
  if (sc.opts.dump_corrected_reads) {

    // get the corrected sequences
    BamRecordVector corrected_alignments;
    for (const auto& rr : all_reads_for_assembly) {
      BamRecordVector this_corrected_alignments;
      unit.bwa_aligner->alignSequence(rr.CorrectedSeq(),
				      rr.UniqueName(),
				      this_corrected_alignments,
				      false,
				      SECONDARY_FRAC,
				      SECONDARY_CAP);
      corrected_alignments.insert(corrected_alignments.end(),
				  this_corrected_alignments.begin(),
				  this_corrected_alignments.end());
    }
    
    // add to svabaThreadUnit for output later
    unit.all_corrected_reads.insert(unit.all_corrected_reads.end(),
	        corrected_alignments.begin(),
		corrected_alignments.end());
  }

  vector<AlignedContig> alc; // where to put the contigs
  BamRecordVector all_aligned_contigs_this_region; // contig-to-genome alignments

  // get the reference sequence of the local region
  string lregion;
  try {
    lregion = unit.ref_genome->QueryRegion(sc.header.IDtoName(region.chr), region.pos1, region.pos2);
  } catch (...) {
    std::cerr << " Caught exception for lregion with reg " << region.ToString(sc.header);
    lregion = "";
  }
  
  // make a BWA index from the locally retrieved sequence
  UnalignedSequenceVector local_usv = {{"local", lregion, string()}};
  BWAIndexPtr local_bwa_index = std::make_shared<BWAIndex>();
  local_bwa_index->ConstructIndex(local_usv);
  BWAAligner local_bwa_aligner(local_bwa_index);

  sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log,
		"...running assemblies for region ",
		region); 
  
  // set the contig prefix
  string name = "c_" +
    std::to_string(region.chr+1) + "_" +
    std::to_string(region.pos1) + "_" +
    std::to_string(region.pos2);
  
  // where to store contigs
  vector<AlignedContig> all_AlignedContigs_this_region;
  
  // setup the engine
  svabaAssemblerEngine engine(name, sc.opts.sgaErrorRate,
			      sc.opts.sgaMinOverlap,
			      sc.readlen);
  engine.fillReadTable(all_reads_for_assembly);

  // do the actual assembly
  //DEBUG
  engine.performAssembly(1/*sc.opts.num_assembly_rounds*/);
  
  // retrieve contigs
  UnalignedSequenceVector all_unaligned_contigs_this_region =
    engine.getContigs();
  sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log, "...assembled ",
		all_unaligned_contigs_this_region.size(),
		" contigs for ",
		name);
  
  // loop and process unaligned contigs
  size_t count_contigs_of_size = 0;
  for (auto& i : all_unaligned_contigs_this_region) {
    
    // if too short, skip
    if ((int)i.Seq.length() < (sc.readlen * 1.2)) 
      continue;

    ++count_contigs_of_size;
    
    //// LOCAL REALIGNMENT
    // align to the local region
    BamRecordVector local_ct_alignments;
    local_bwa_aligner.alignSequence(i.Seq,
				    i.Name,
				    local_ct_alignments,
				    false,
				    SECONDARY_FRAC,
				    SECONDARY_CAP);
    
    // check if it has a non-local alignment
    /*bool valid_sv = true;
      for (auto& aa : local_ct_alignments) {
      if (aa.NumClip() < MIN_CLIP_FOR_LOCAL) 
      valid_sv = false; // has a non-clipped local alignment. can't be SV. Indel only
      }*/

    // do the main realignment of unaligned contigs to the reference genome
    BamRecordVector human_alignments;
    unit.bwa_aligner->alignSequence(i.Seq, i.Name,
				   human_alignments,
				   false, SECONDARY_FRAC, SECONDARY_CAP);

    // store all contig alignment
    all_aligned_contigs_this_region.insert(
	       all_aligned_contigs_this_region.end(),
	       human_alignments.begin(),
	       human_alignments.end());

    // sort the alignments by position
    std::sort(human_alignments.begin(), human_alignments.end());
    
    // add the chromosome string name
    for (auto& rr : human_alignments) {
      	rr.AddZTag("MC", sc.header.IDtoName(rr.ChrID()));
    }

    // add contig alignments to svabaThreadUnit for writing later
    unit.master_contigs.insert(unit.master_contigs.end(),
			       human_alignments.begin(),
			       human_alignments.end());

    // make the AlignedContig object for this contig
    AlignedContig ac(human_alignments, sc.prefixes);
    
    // assign the local variable to each
    ac.checkLocal(region);

    // add this
    all_AlignedContigs_this_region.push_back(ac);
  }

  
  sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log, "...assembled ",
		count_contigs_of_size,
		" contigs that meet size criteria for ",
		name);
  
  // didnt get any contigs that made it all the way through
  if (!all_AlignedContigs_this_region.size())
    return true;

  // Make a BWA mapper of the contigs themselves
  BWAIndexPtr contig_bwa_index = std::make_shared<BWAIndex>();
  contig_bwa_index->ConstructIndex(all_unaligned_contigs_this_region);
  BWAAligner contig_bwa_aligner(contig_bwa_index);

  // align the reads to the contigs
  for (auto& i : all_reads_for_assembly) {

    // do the read to contig alignment for a single read
    BamRecordVector reads2contigs_brv;
    contig_bwa_aligner.alignSequence(i.CorrectedSeq(),
				     i.Qname(),
				     reads2contigs_brv,
				     false,
				     0.60,
				     10000);
      
    // convert read to contig alignments to a svabaReadVector
    svabaReadVector reads2contigs_srv;
    for (auto& r : reads2contigs_brv)
      reads2contigs_srv.push_back(svabaRead(r, i.Prefix()));
    reads2contigs_brv.clear();
    
    // make sure we have only one alignment per contig
    unordered_set<string> cc;
    
    // check which ones pass
    svabaReadVector bpass;
    for (auto& r : reads2contigs_srv) {
      
      // make sure alignment score is OK
      int thisas = 0;
      r.GetIntTag("AS", thisas);
      if ((double)r.NumMatchBases() * 0.5 > thisas)
      	continue;
      
      //bool length_pass = (r.PositionEnd() - r.Position()) >=
      //	((double)seqr.length() * 0.75);

      string contig_name = all_unaligned_contigs_this_region[r.ChrID()].Name;
      if (/*length_pass && */!cc.count(contig_name)) {
	bpass.push_back(r);
	cc.insert(contig_name);
      }
    }

    // annotate the original read
    for (auto& r : bpass) {
      
      r2c this_r2c; // alignment of this read to this contig
      if (r.ReverseFlag())
	this_r2c.rc = true;
      
      this_r2c.AddAlignment(r);
      std::string contig_name = all_unaligned_contigs_this_region[r.ChrID()].Name;
      i.AddR2C(contig_name, this_r2c);
      
      // add the read to the right contig (loop to check for right contig)
      for (auto& a : all_AlignedContigs_this_region) {
	if (a.getContigName() != contig_name)
	  continue;
	a.AddAlignedRead(i);
      }
    } // read2contig alignment loop (per read)
    
  } // short read realignment loop (all reads)
  

  // Get contig coverage, discordant matching to contigs, etc
  for (auto& a : all_AlignedContigs_this_region) {
    
    // repeat sequence filter
    a.assessRepeats();
    
    a.splitCoverage();	
    // now that we have all the break support, check that the complex breaks are OK
    a.refilterComplex(); 
    // add discordant reads support to each of the breakpoints
    a.addDiscordantCluster(dmap);
    // add in the cigar matches
    a.checkAgainstCigarMatches(cigmap);
    // add to the final structure
    unit.master_alc.push_back(a);
    
  }
  
  st.stop("as");

  // get the breakpoints
  vector<BreakPoint> bp_glob;
  
  for (auto& i : alc) {
    vector<BreakPoint> allbreaks = i.getAllBreakPoints();
    bp_glob.insert(bp_glob.end(), allbreaks.begin(), allbreaks.end());
  }

  //if (dbsnp_filter && opt::dbsnp.length()) {
  //  WRITELOG("...DBSNP filtering", opt::verbose > 1, false);
  //  for (auto & i : bp_glob) 
  //   dbsnp_filter->queryBreakpoint(i);
  //}

  // filter against blacklist
  for (auto& i : bp_glob) 
    i.checkBlacklist(sc.blacklist);

  // add in the discordant clusters as breakpoints
  for (auto& i : dmap) {
    
    // dont send DSCRD if FR and below size
    /*bool below_size = 	i.second.m_reg1.strand == '+' && i.second.m_reg2.strand == '-' && 
      (i.second.m_reg2.pos1 - i.second.m_reg1.pos2) < min_dscrd_size_for_variant && 
      i.second.m_reg1.chr == i.second.m_reg2.chr;
    */
    
    // DiscordantCluster not associated with assembly BP and has 2+ read support
    if (!i.second.hasAssociatedAssemblyContig() && 
	(i.second.tcount + i.second.ncount) >= MIN_DSCRD_READS_DSCRD_ONLY &&
	i.second.valid()) {
      //	!below_size) {
      BreakPoint tmpbp(i.second, dmap,
		       region, sc.header, sc);
      bp_glob.push_back(tmpbp);
    }
  }

  // de duplicate the breakpoints
  std::sort(bp_glob.begin(), bp_glob.end());
  bp_glob.erase( unique( bp_glob.begin(), bp_glob.end() ), bp_glob.end() );

  // add the coverage data to breaks for allelic fraction computation
  unordered_map<string, STCoverage*> covs;
  for (auto& i : sc.opts.bams)
    covs[i.first] = &unit.walkers.at(i.first).cov;

  for (auto& i : bp_glob)
    i.addCovs(covs);
  
  for (auto& i : bp_glob) {
    i.readlen = sc.readlen; // set the readlength
    i.scoreBreakpoint(sc.opts.lod, sc.opts.lodDb, sc.opts.lodSomatic,
		      sc.opts.lodSomaticDb,
		      sc.opts.scaleError, /*min_dscrd_size_for_variant*/2000);
  }
  
  // label somatic breakpoints that intersect directly with normal as NOT somatic
  unordered_set<string> norm_hash;
  for (auto& i : bp_glob) // hash the normals
    if (!i.somatic_score && i.confidence == "PASS" && i.evidence == "INDEL") {
      norm_hash.insert(i.b1.hash());
      norm_hash.insert(i.b2.hash());
      norm_hash.insert(i.b1.hash(1));
      norm_hash.insert(i.b1.hash(-1));
    }

  // find somatic that intersect with norm. Set somatic = 0;
  for (auto& i : bp_glob)  
    if (i.somatic_score && i.evidence == "INDEL" && (norm_hash.count(i.b1.hash()) || norm_hash.count(i.b2.hash()))) {
      i.somatic_score = -3;
    }

  // remove indels at repeats that have multiple variants
  unordered_map<string, size_t> ccc;
  for (auto& i : bp_glob) {
    if (i.evidence == "INDEL" && i.repeat_seq.length() > 6) {
      ++ccc[i.b1.hash()];
    }
  }
  for (auto& i : bp_glob) {
    if (i.evidence == "INDEL" && ccc[i.b1.hash()] > 1)
      i.confidence = "REPVAR";
  }

  // remove somatic calls if they have a germline normal SV in them or indels with 
  // 2+germline normal in same contig
  unordered_set<string> bp_hash;
  for (auto& i : bp_glob) { // hash the normals
    if (!i.somatic_score && i.evidence != "INDEL" && i.confidence == "PASS") {
      bp_hash.insert(i.cname);
    }
  }
  for (auto& i : bp_glob)  // find somatic that intersect with norm. Set somatic = 0;
    if (i.somatic_score && i.num_align > 1 && bp_hash.count(i.cname)) {
      i.somatic_score = -2;
    }


  
  // remove somatic SVs that overlap with germline svs
  if (sc.germline_svs.size()) {
    for (auto& i : bp_glob) {
      if (i.somatic_score && i.b1.gr.chr == i.b2.gr.chr && i.evidence != "INDEL") {
	GenomicRegion gr1 = i.b1.gr;
	GenomicRegion gr2 = i.b2.gr;
	gr1.Pad(GERMLINE_CNV_PAD);
	gr2.Pad(GERMLINE_CNV_PAD);
	if (sc.germline_svs.OverlapSameInterval(gr1, gr2)) {
	  i.somatic_score = -1;
	}
      }
    }
      
  }

  // add the ref and alt tags
  // WHY IS THIS NOT THREAD SAFE?
  for (auto& i : bp_glob)
    i.setRefAlt(unit.ref_genome.get()); 

  // transfer local versions to thread store
  for (const auto& a : alc)
    if (a.hasVariant()) {
      unit.master_alc.push_back(a);
      //unit.m_bamreads_count += a.NumBamReads();
    }
  //for (const auto& d : dmap)
  //  unit.m_disc_reads += d.second.reads.size();

  // add the discordant clusters to the svabathreadunit for writing later
  unit.m_disc.insert(dmap.begin(), dmap.end());

  //for (const auto& a : alc)
  //  unit.m_bamreads_count += a.NumBamReads();
  for (auto& i : bp_glob) 
    if ( i.hasMinimal() && (i.confidence != "NOLOCAL" || i.complex_local ) ) 
      unit.m_bps.push_back(i);

  // clear out the reads and reset the walkers
  for (auto& w : unit.walkers) {
    w.second.clear(); 
  }

  // dump if getting to much memory
  if (unit.MemoryLimit(THREAD_READ_LIMIT, THREAD_CONTIG_LIMIT) && !sc.opts.hp) {
    sc.writer.writeUnit(unit, sc); // mutex and flush are inside this call
    unit.clear();
  }
   
  st.stop("pp");
  
  // display the run time
  /*  sc.logger.log(true, true,
		svabaUtils::runTimeString(read_counts.first,
					  read_counts.second,
					  alc.size(),
					  region,
					  sc.header,
					  st, start)); 
  */

  return true;
}
