#include "SvabaRegionProcessor.h"

#include "fml.h"
#include <iomanip>  // For std::setw and std::right
#include <sstream>  // For std::ostringstream
#include <chrono>

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
using SeqLib::BamRecordPtr;
using SeqLib::BamRecordPtrVector;

using std::vector;
using std::unordered_map;
using std::unordered_set;
using std::string;

SvabaRegionProcessor::SvabaRegionProcessor(SvabaSharedConfig& sh_cf) : sc(sh_cf)
{ }

void SvabaRegionProcessor::runMateCollectionLoop(const GenomicRegion& region,
						 svabaThreadUnit& stu)
{
  
  // store the newly found bad mate regions
  GRC this_bad_mate_regions; 

  MateRegionVector all_somatic_mate_regions;
  // add the origional, don't want to double back    
  all_somatic_mate_regions.add(MateRegion(region.chr, region.pos1, region.pos2));
  
  // collect together mate regions from control
  MateRegionVector normal_mate_regions;
  for (const auto& [key, walker] : stu.walkers) {
    if (!key.empty() && key[0] == 'n') {
      normal_mate_regions.Concat(walker->mate_regions);
    }
  }  
  normal_mate_regions.MergeOverlappingIntervals();
  normal_mate_regions.CreateTreeMap();
  
  // get the mates from somatic 3+ mate regions
  // that don't overlap with normal mate region
  MateRegionVector somatic_mate_regions;
  for (const auto& [key, walker] : stu.walkers) {
    if (!key.empty() && key[0] == 't') {
      for (const auto& i : walker->mate_regions) {
	if (i.count >= sc.opts.mateLookupMin &&
	    !normal_mate_regions.CountOverlaps(i) &&
	    //!stu.badd.CountOverlaps(i) && //TODO - don't use this anymore
	    (sc.blacklist.empty() || !sc.blacklist.CountOverlaps(i)))
	  {
	    somatic_mate_regions.add(i);
	  }
      }
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
  for (auto& [key, walker] : stu.walkers) {
    const int before = walker->reads.size();
    walker->m_limit = sc.opts.mate_region_lookup_limit;
    
    if (!walker->SetRegions(gg))
      continue;
    
    walker->get_coverage = false;
    walker->get_mate_regions = false;
    walker->mate_regions.clear();

    auto prior_end = walker->reads.end();
    GRC bad = walker->readBam(stu);

    // reset
    walker->m_limit = sc.opts.weird_read_limit;
    walker->get_coverage = true;
    walker->get_mate_regions = true;
    
    //stu.badd.Concat(bad);
    //stu.badd.MergeOverlappingIntervals();
    //stu.badd.CreateTreeMap();
    
    assert(!key.empty());
    auto& count = (key[0] == 't') ? stu.st.mate_read_count.first : stu.st.mate_read_count.second;    
    count += walker->reads.size() - before;
  }
  
  sc.logger.log(sc.opts.verbose > 1, sc.opts.verbose_log, 
		"......Mate region reads: <case, control>: <",
		stu.st.mate_read_count.first, ",",
		stu.st.mate_read_count.second, ">");

}


bool SvabaRegionProcessor::process(const SeqLib::GenomicRegion& region,
                                   svabaThreadUnit&             unit,
                                   size_t                       threadId)
{

  // count for this unit
  unit.processed_count++;
  sc.total_regions_done++; 
  if (sc.total_regions_done % 1000 == 0) {
    std::ostringstream msg;
    msg << "...processing "
	<< std::right << std::setw(5) << SeqLib::AddCommas(sc.total_regions_done)
	<< " of "
	<< std::right << std::setw(6) << SeqLib::AddCommas(sc.total_regions_to_process)
	<< " for thread "
	<< std::right << std::setw(2) << unit.threadId
	<< " for region "
	<< region;
    sc.logger.log(true, true, msg.str()); 
  }
  
  sc.logger.log(sc.opts.verbose > 1, sc.opts.verbose_log,
		"===Running region ", region.ToString(sc.header),
		" on thread ", unit.threadId);
  
  // start a new timer
  unit.st = svabaUtils::svabaTimer();
  unit.st.gr = region;
  unit.st.start();
  
  // holder for Bam: (read : cigar)
  unordered_map<string, CigarMap> cigmap;
  
  // loop all of the BAMs (walkers)
  CountPair read_counts = {0,0}; //tumor, normal
  for (auto& [key, walker] : unit.walkers) {
    
    // Set regions based on command-line region or file-defined intervals
    if (!region.IsEmpty())
      walker->SetRegion(region);
    else if (!sc.file_regions.empty())
      walker->SetRegions(sc.file_regions);
    // else: walk entire BAM (default)
    
    sc.logger.log(sc.opts.verbose > 1, sc.opts.verbose_log, "---running svabaBamWalker", walker);

    // reset the walker params
    walker->m_limit = sc.opts.weird_read_limit;
    walker->get_coverage = true;
    walker->get_mate_regions = true;
    SeqLib::GRC bad = walker->readBam(unit);
    //unit.badd.Concat(bad);

    sc.logger.log(sc.opts.verbose > 3, false, "...finished reading, found ", bad.size(), " bad regions");

    // TODO - right now we don't do anything with these bad regions
    //unit.badd.MergeOverlappingIntervals();
    //unit.badd.CreateTreeMap(); 
    
    // Update read count based on sample type
    auto& count = (key[0] == 't') ? unit.st.weird_read_count.first : unit.st.weird_read_count.second;
    count += walker->reads.size();
    
    // Store CIGARs
    cigmap[key] = walker->cigmap;
  } //end BAM loop

  sc.logger.log(sc.opts.verbose > 1, sc.opts.verbose_log,
		"...main region reads <case,control> ",
		"<", unit.st.weird_read_count.first, ",",
		unit.st.weird_read_count.second, ">");
  
  // adjust counts and timer
  unit.st.stop("r");

  // get the mate reads, if this is local assembly and has insert-size distro
  if (!region.IsEmpty()) { 
    runMateCollectionLoop(region, unit);
    unit.st.stop("m");
  }
  
  // do the discordant read clustering
  sc.logger.log(sc.opts.verbose > 1, false, 
		"...discordant read clustering");
  
  // tag the reads by discordant status
  svabaReadPtrVector all_discordant_reads;  

  for (auto& [_, walker] : unit.walkers) {
    //walker->RealignDiscordants(unit);
    for (auto& r : walker->reads) {
      if (r->dd > 0)
	all_discordant_reads.push_back(r);
    }
  }
  unit.st.dc_read_count = all_discordant_reads.size();

  // do the discordant read clustering across BAMs
  DiscordantClusterMap dmap = DiscordantCluster::clusterReads(all_discordant_reads, 
							      region,
							      sc.header);
  all_discordant_reads.clear();
  unit.st.dc_cluster_count = dmap.size();
  
  // add the discordant clusters to the svabathreadunit for writing later
  unit.m_disc.insert(dmap.begin(), dmap.end());
  
  // for dumping all reads
  if (sc.opts.dump_weird_reads) {
    for (auto& [_, dc] : dmap) {
      dc.labelReads(); // add the discordantcluster label to the read
    }
    for (const auto& [_, walker] : unit.walkers) {
      unit.all_weird_reads.insert(unit.all_weird_reads.end(),
				  walker->reads.begin(),
				  walker->reads.end());
    }
  }
  
  // do the discordant read clustering
  sc.logger.log(sc.opts.verbose > 1, false, "...error correcting");
  
  // do kmer correction
  if (sc.opts.ecCorrectType == "f") {
    
    // train and correct
    for (auto& [_, walker] : unit.walkers) {
      if (walker->reads.empty())
	continue;
      walker->Train();
      walker->ClearTraining();
      walker->AddBackReadsToCorrect();
      walker->ErrorCorrect();
    }
    
    // retrieve the corrected sequences
    string s, name_dum;
    size_t num_reads_corrected = 0;
    for (auto& [_, walker] : unit.walkers) {
      for (auto& r : walker->reads) {
	walker->bfc.GetSequence(s, name_dum);
	r->SetCorrectedSeq(s);
	++num_reads_corrected;
      }
      walker->bfc.clear(); // erase all reads and training dat      
    }
    
    sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log, "...BFC corrected ", num_reads_corrected);
    
    unit.st.stop("k");
  }
  
  // re-align the corrected reads and dump
  if (sc.opts.dump_corrected_reads) {

    // get the corrected sequences
    for (const auto& [_, walker] : unit.walkers) {
      for (const auto& r : walker->reads) {
	unit.bwa_aligner->alignSequence(r->CorrectedSeq(),
					r->UniqueName(),
					unit.all_corrected_reads,
					false,
					0.6,
					0);
      }
    }
    sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log,
		  "...realigned corrected reads to genome");
  }
  
  sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log,
		"...running assemblies for region ",
		region); 
  
  // set the contig prefix
  string ctg_prefix = "c_" +
    std::to_string(region.chr+1) + "_" +
    std::to_string(region.pos1) + "_" +
    std::to_string(region.pos2);
  
  // where to store the AlignedContigs for this region 
  std::unordered_map<std::string, AlignedContig> all_AlignedContigs_this_region;

  // setup the engine and peform assembly
  UnalignedSequenceVector all_unaligned_contigs_this_region;
  if (true)
  {
    
    svabaAssemblerEngine engine(ctg_prefix, sc.opts.sgaErrorRate,
				sc.opts.sgaMinOverlap,
				sc.readlen);

    size_t dbg = 0;
    for (const auto& [_, walker] : unit.walkers) {

      //debug
      if (region.pos1 == 4287501) {
	for (const auto& r : walker->reads)
	  std::cout << " RSEQ " << r->CorrectedSeq() << std::endl;
      }
      
      engine.fillReadTable(walker->reads);
      dbg += walker->reads.size();
    }
    
    // do the discordant read clustering
    sc.logger.log(sc.opts.verbose > 1, false, 
		  "...assembling reads");
    
    // do the actual assembly
    engine.performAssembly(1/*sc.opts.num_assembly_rounds*/);
    
    // retrieve contigs
    all_unaligned_contigs_this_region = engine.getContigs();

    if (region.pos1 == 4287501) 
      std::cerr << " contigs " <<
	all_unaligned_contigs_this_region.size() << std::endl;
  }
  if (false)
    {
      // build the reads structure
      size_t n_reads = 0;
      for (const auto& [_, walker] : unit.walkers) {
	n_reads += walker->reads.size();
      }
      fseq1_t *fseq = (fseq1_t*)calloc(n_reads, sizeof(fseq1_t));

      size_t i = 0;
      for (const auto& [_, walker] : unit.walkers) {
	for (const auto& r : walker->reads) {
	  fseq[i].seq     = strdup(r->CorrectedSeq().c_str());
	  fseq[i].l_seq   = r->CorrectedSeq().length();
	  
	  // we need a quality string for FASTQ mode give all high quality
	  std::string q(fseq[i].l_seq, 'I');
	  fseq[i].qual    = strdup(q.c_str());
	  ++i;
	}
      }
      fml_opt_t opt;
      fml_opt_init(&opt);

      // 3) Run the assembler
      int n_utgs = 0;
      fml_utg_t *utgs = nullptr;
      if (n_reads > 0)
	utgs = fml_assemble(&opt,n_reads,fseq,&n_utgs);
      
      // 4) Copy out contig names/sequences into a C++ vector
      all_unaligned_contigs_this_region.reserve(n_utgs);
      for (int i = 0; i < n_utgs; i++) {
	std::string name = ctg_prefix + "_" + std::to_string(i+1) + "C";
	std::string seq = utgs[i].seq;
	SeqLib::UnalignedSequence us(name, seq);
	all_unaligned_contigs_this_region.push_back(std::move(us));
      }

      // 5) Clean up
      fml_utg_destroy(n_utgs, utgs);
      for (int i = 0; i < n_reads; i++) {
        free(fseq[i].seq);
        free(fseq[i].qual);
      }
      free(fseq);
      
    }
  
  unit.st.contig_count = all_unaligned_contigs_this_region.size();
  
  sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log, "...assembled ",
		all_unaligned_contigs_this_region.size(),
		" contigs for ",
		ctg_prefix);
  
  /////////
  // SETUP FOR LOCAL ALIGNMENT
  // get the reference sequence of the local region
  /*  string lregion;
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
  */
  //////////////
  
  // loop and process unaligned contigs
  //vector<AlignedContig> alc; // where to put the contigs
  BamRecordPtrVector all_aligned_contigs_this_region; // contig-to-genome alignments
  size_t count_contigs_of_size = 0;
  for (auto& i : all_unaligned_contigs_this_region) {
    
    // if too short, skip
    if ((int)i.Seq.length() < (sc.readlen * 1.2)) 
      continue;
    
    ++count_contigs_of_size;
    
    //// LOCAL REALIGNMENT
    // align to the local region
    /*
      BamRecordVector local_ct_alignments;
      local_bwa_aligner.alignSequence(i.Seq,
      i.Name,
      local_ct_alignments,
      false,
      SECONDARY_FRAC,
      SECONDARY_CAP);
      
      // check if it has a non-local alignment
      bool valid_sv = true;
      for (auto& aa : local_ct_alignments) {
      if (aa.NumClip() < MIN_CLIP_FOR_LOCAL) 
      valid_sv = false; // has a non-clipped local alignment. can't be SV. Indel only
      }*/
    
    // do the main realignment of unaligned contigs to the reference genome
    BamRecordPtrVector human_alignments;
    unit.bwa_aligner->alignSequence(i.Seq, i.Name,
				    human_alignments,
				    false, SECONDARY_FRAC, SECONDARY_CAP);
    
    
    // sort the alignments by position
    std::sort(human_alignments.begin(),
	      human_alignments.end(),
	      SeqLib::BamRecordSort::ByReadPositionSharedPtr());
    
    // store all contig alignment
    all_aligned_contigs_this_region.insert(all_aligned_contigs_this_region.end(),
					   human_alignments.begin(),
					   human_alignments.end());
    
    // add contig alignments to svabaThreadUnit for writing later
    unit.master_contigs.insert(unit.master_contigs.end(),
			       human_alignments.begin(),
			       human_alignments.end());

    if (region.pos1 == 4287501)
      for (const auto& r : human_alignments)
     	std::cout << "HA " << *r << std::endl;
    
    // make the AlignedContig object for this contig
    AlignedContig ac(human_alignments, region, &sc);

    // add this
    all_AlignedContigs_this_region.insert_or_assign(ac.getContigName(), std::move(ac));
  }
  unit.st.aligned_contig_count = all_aligned_contigs_this_region.size(); 

  sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log, "...contigs that meet size criteria: ",
		count_contigs_of_size);

  if (region.pos1 == 4287501) {
    std::cerr <<" ALIGNED CONTIGS " << 
      all_aligned_contigs_this_region.size() << std::endl;
    for (const auto& r : all_unaligned_contigs_this_region)
      std::cout << "PASSED CONTIG " << r.Name << "\t" <<
	r.Seq << std::endl;
  }
	
  // didnt get any contigs that made it all the way through
  if (!all_AlignedContigs_this_region.size()) {
    unit.flush();
    return true;
  }

  // Make a BWA mapper of the contigs themselves
  BWAIndexPtr contig_bwa_index = std::make_shared<BWAIndex>();
  contig_bwa_index->ConstructIndex(all_unaligned_contigs_this_region);
  BWAAligner contig_bwa_aligner(contig_bwa_index);
  
  sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log, "...aligning reads to contigs");
  
  // align the reads to the contigs
  BamRecordPtrVector all_read2contigs;
  size_t seen = 0;

  for (auto& [_, walker] : unit.walkers) {
    for (auto& i : walker->reads) {
      ++seen;
      if (seen % 1000 == 0)
	sc.logger.log(sc.opts.verbose > 4, sc.opts.verbose_log, "- r2c ", SeqLib::AddCommas(seen));

      // do the read to contig alignment for a single read
      BamRecordPtrVector read2contigs;
      contig_bwa_aligner.alignSequence(i->CorrectedSeq(),
				       i->Qname(),
				       read2contigs,
				       false,
				       0.60,
				       1000);

      // make sure we have only one alignment per contig
      unordered_set<string> cc;
      
      // check which ones pass
      auto keep = [&](const BamRecordPtr& r) {
	int thisas = 0;
	r->GetIntTag("AS", thisas);
	if ((double)r->NumMatchBases() * 0.5 > thisas)
	  return false;
	
	std::string contig_name = all_unaligned_contigs_this_region[r->ChrID()].Name;
	if (cc.count(contig_name))
	  return false;
	
	cc.insert(contig_name);
	return true;
      };
	
      read2contigs.erase(
			 std::remove_if(read2contigs.begin(), read2contigs.end(),
					[&](const BamRecordPtr& r) { return !keep(r); }),
			 read2contigs.end());

      // annotate the original read
      for (auto& r : read2contigs) {
	all_read2contigs.push_back(r); // add to total	
	r2c this_r2c; // alignment of this read to this contig
	if (r->ReverseFlag())
	  this_r2c.rc = true;
	
	this_r2c.AddAlignment(r); // doesn't copy any memory
	std::string contig_name = all_unaligned_contigs_this_region[r->ChrID()].Name; //chromomomse here = contig
	i->AddR2C(contig_name, this_r2c); // i = AlignedContig

	if (region.pos1 == 4287501) {
	  std::cout << "R2C " << r->Qname() << " - " <<
	    r->Sequence() << " contig " << r->ChrName(sc.header) <<
	    std::endl;
	}
	
	// add the read to the right contig
	auto it = all_AlignedContigs_this_region.find(contig_name);
	if (it != all_AlignedContigs_this_region.end()) {
	  it->second.AddAlignedRead(i); // we are adding the svabaBamRead!! not the cread
	}
      } // read2contig alignment loop (per read)
    } // short read realignment loop (all reads)  
  } // short read realignment loop (all bam walkers)
  
  sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log, "...done aligning reads to contigs for ",
		all_read2contigs.size(), " reads. Processing variants");


  // BPS need a repeat filter
  // BPS can have addDiscordantFilter(dmap) applied directly
  
  // Get contig coverage, discordant matching to contigs, etc
  for (auto& [_, a] : all_AlignedContigs_this_region) {

    // right now all break points are stored in an AlignedContig
    // this routine calls all internally stored BreakPoint and
    // checks split read coverage, since all read2contig alignments
    // for this AlignedContig are stored in this AlignedContig
    a.splitCoverage();
    
    // now that we have all the break support,
    // check that the complex breaks are OK
    //a.refilterComplex();
    
    // add in the cigar matches
    //a.checkAgainstCigarMatches(cigmap);
    
    // add to the final structure
    //alc.push_back(std::move(a));
  }
  
  if (region.pos1 == 4287501) {
    for (const auto& [_, bbb] : all_AlignedContigs_this_region) {
      std::cout << " ALC " << bbb.printToAlignmentsFile(sc.header) << std::endl;
    }
  }
  
  unit.st.stop("as");

  // get the breakpoints
  BreakPointPtrVector bp_glob;
  for (auto& [_, i] : all_AlignedContigs_this_region) {
    BreakPointPtrVector tmp_allbreaks = i.getAllBreakPoints();
    bp_glob.insert(bp_glob.end(),
		   tmp_allbreaks.begin(),
		   tmp_allbreaks.end());
  }

  // set each interal BreakEnd as to whether it overlaps with the
  // "local" assembly window
  for (auto& i : bp_glob)
    i->setLocal(region);
  
  //if (dbsnp_filter && opt::dbsnp.length()) {
  //  WRITELOG("...DBSNP filtering", opt::verbose > 1, false);
  //  for (auto & i : bp_glob) 
  //   dbsnp_filter->queryBreakpoint(i);
  //}
  
  // filter against blacklist
  for (auto& i : bp_glob) 
    i->checkBlacklist(sc.blacklist);

  // add discordant reads
  for (auto& i : bp_glob) {
    i->CombineWithDiscordantClusterMap(dmap);
  }
  
  if (region.pos1 == 4287501) {
    for (const auto& bbb : bp_glob) {
      std::cout << " BPGLOB " << bbb->printSimple(sc.header) << std::endl;
    }
  }

  // add in the discordant clusters as breakpoints
  for (auto& [_, d] : dmap) {

    // dont send DSCRD if FR and below size
    /*bool below_size = 	i.second.m_reg1.strand == '+' && i.second.m_reg2.strand == '-' && 
      (i.second.m_reg2.pos1 - i.second.m_reg1.pos2) < min_dscrd_size_for_variant && 
      i.second.m_reg1.chr == i.second.m_reg2.chr;
    */

    // already associated with an assembly, so don't make a new variant
    // this information should be added to the DiscordantCluster
    // during "BreakPoint::CombineWithDiscordantClusterMap()"
    if (d.hasAssociatedAssemblyContig())
      continue;

    // too little read support, skip
    if ( (d.tcount + d.ncount) < MIN_DSCRD_READS_DSCRD_ONLY)
      continue;

    if (!d.valid())
      continue;
    
    // make the Discordant cluster
    BreakPointPtr tmpbp =
      std::make_shared<BreakPoint>(d, 
				   dmap,
				   region,
				   &sc);
    bp_glob.push_back(tmpbp);
  }
  
  // de duplicate the breakpoints
  std::sort(bp_glob.begin(), bp_glob.end(),
	    [](auto const& a, auto const& b){
	      return *a < *b;
	    });
  auto new_end = std::unique(
			     bp_glob.begin(), 
			     bp_glob.end(),
			     [](auto const& a, auto const& b){
			       // return true if *a and *b are "equal" in the sense of your operator<
			       return !(*a < *b) && !(*b < *a);
			     });
  
  if (region.pos1 == 4287501) {
    for (const auto& bbb : bp_glob) {
      std::cout << " AFTER BPGLOB " << bbb->printSimple(sc.header) << std::endl;
    }
  }
  
  // add the coverage data to breaks for allelic fraction computation
  unordered_map<string, STCoverage*> covs;
  for (auto& [key, walker] : unit.walkers) {
    covs[key] = &walker->cov;
  }
  for (auto& i : bp_glob)
    i->addCovs(covs);

  // score the breakpoints
  for (auto& i : bp_glob) {
    i->scoreBreakpoint(sc.opts.lod, sc.opts.lodDb, sc.opts.lodSomatic,
		      sc.opts.lodSomaticDb,
		      sc.opts.scaleError, /*min_dscrd_size_for_variant*/2000);
  }
  
  // label somatic breakpoints that intersect directly with normal as NOT somatic
  unordered_set<string> norm_hash;
  for (auto& i : bp_glob) { // hash the normals
    if (i->somatic_score <= 0 && i->confidence == "PASS") {
      for (const auto& h : i->getBreakEndHashes())
	norm_hash.insert(h);
    }
  }
  
  // find somatic that intersect with norm. Set somatic = 0;
  for (auto& i : bp_glob) {
    // don't need to check normal against normal
    if (i->somatic_score <= 0)
      continue;

    // if overlaps then turn off somatic
    for (const auto& h : i->getBreakEndHashes())
      if (norm_hash.count(h))
	i->somatic_score = -3; // not somatic
  }
  
  // remove indels at repeats that have multiple variants
  // unordered_map<string, size_t> ccc;
  // for (auto& i : bp_glob) {
  //   if (i->svtype == SVType::INDEL && i->repeat_seq.length() > 6) {
  //     ++ccc[i->b1.hash()];
  //   }
  // }
  // for (auto& i : bp_glob) {
  //   if (i->svtype == SVType::INDEL && ccc[i->b1.hash()] > 1)
  //     i->confidence = "REPVAR";
  // }

  // remove somatic calls if they have a germline normal SV in
  // them or indels with 
  // 2+ germline normal in same contig
  // unordered_set<string> bp_hash;
  // for (auto& i : bp_glob) { // hash the normals
  //   if (!i->somatic_score && i->svtype != SVType::INDEL && i->confidence == "PASS") {
  //     bp_hash.insert(i->cname);
  //   }
  // }
  // for (auto& i : bp_glob)  // find somatic that intersect with norm. Set somatic = 0;
  //   if (i->somatic_score && i->num_align > 1 && bp_hash.count(i->cname)) {
  //     i->somatic_score = -2;
  //   }
  
  // remove somatic SVs that overlap with germline svs
  if (sc.germline_svs.size()) {
    for (auto& i : bp_glob) {

      // don't filter normals or indels
      if (i->somatic_score <= 0 || i->isIndel())
	continue;

      GenomicRegion gr1 = i->BreakEndAsGenomicRegionLeft();
      GenomicRegion gr2 = i->BreakEndAsGenomicRegionRight();
      gr1.Pad(GERMLINE_CNV_PAD);
      gr2.Pad(GERMLINE_CNV_PAD);
      if (sc.germline_svs.OverlapSameInterval(gr1, gr2)) {
	i->somatic_score = -1;
      }
    }
    
  }

  // add the ref and alt tags
  for (auto& i : bp_glob)
    i->setRefAlt(unit.ref_genome.get(), sc.header); 
  
  // transfer local versions to thread store
  for (const auto& [_, a] : all_AlignedContigs_this_region) {
    if (a.hasVariant()) {
      unit.master_alc.push_back(a); //std::move(a));
    }
  }
  all_AlignedContigs_this_region.clear();
  //for (const auto& d : dmap)
  //  unit.m_disc_reads += d.second.reads.size();

  //for (const auto& a : alc)
  //  unit.m_bamreads_count += a.NumBamReads();
  for (auto& i : bp_glob) 
    if ( i->hasMinimal() )
      unit.m_bps.push_back(i);
  
  unit.st.bps_count = bp_glob.size();

  if (region.pos1 == 4287501) {
    for (const auto& bbb : unit.m_bps/*bp_glob*/) {
      std::cout << " FINAL BPGLOB " << bbb->printSimple(sc.header) << std::endl;
    }
  }
  
  // dump if getting to much memory
  unit.st.stop("pp");
  unit.ss << unit.st.logRuntime(sc.header) << "\n";
  unit.flush();
    
  // display the run time
  if (sc.total_regions_done % 1000 == 0) {
    sc.logger.log(sc.opts.verbose > 0, true,
		  svabaUtils::runTimeString(read_counts.first,
					    read_counts.second,
					    all_AlignedContigs_this_region.size(),
					    region,
					    sc.header,
					    unit.st, sc.start));
  }
  
  sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log, "...done with this region");
  return true;
}
