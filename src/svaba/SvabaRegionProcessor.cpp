#include "SvabaRegionProcessor.h"

extern "C" {
#include "fml.h"
}
#include <iomanip>  // For std::setw and std::right
#include <sstream>  // For std::ostringstream
#include <chrono>
#include <cstring>  // For memset

#include "SvabaDebug.h"
#include "SvabaUtils.h"
#include "SvabaOutputWriter.h"
#include "ContigAlignmentScore.h"

#include "SeqLib/BWAAligner.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/BFC.h"
#include "SeqLib/BamRecord.h"

#include "SvabaRead.h"
#include "SvabaAssemblerEngine.h"

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

#include "SvabaAssemblerConfig.h"

SvabaRegionProcessor::SvabaRegionProcessor(SvabaSharedConfig& sh_cf) : sc(sh_cf)
{ }

GRC SvabaRegionProcessor::runMateCollectionLoop(const GenomicRegion& region,
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
    return GRC();
  
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

  return somatic_mate_region_collection;
}


bool SvabaRegionProcessor::process(const SeqLib::GenomicRegion& region,
                                   svabaThreadUnit&             unit,
                                   size_t                       threadId)
{

  // SvABA2.0: are any dump-reads BAMs going to be written? Used to skip
  // bi:Z and bz:Z tag bookkeeping on the read hot path when no dump
  // output is requested — the tags feed only the weird/discordant/
  // corrected BAM writers in SvabaOutputWriter, so under the default
  // production command line (no --dump-reads) the tag work is pure
  // waste. Declared at function scope so both the r2c-alignment
  // (bz:Z) and BP-finalization (bi:Z) blocks see it.
  const bool need_read_tags =
      sc.opts.dump_weird_reads      ||   // compile-time constexpr
      sc.opts.dump_discordant_reads ||   // both flip together under
      sc.opts.dump_corrected_reads;      // --dump-reads

  // count for this unit
  unit.processed_count++;
  sc.total_regions_done++;
  if (sc.total_regions_done % 250 == 0) {
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
  
  // get the reference sequence of the local region
  string lregion;
  try {
    lregion = unit.ref_genome->QueryRegion(sc.header.IDtoName(region.chr), region.pos1, region.pos2);
  } catch (...) {
    std::cerr << " Caught exception for lregion with reg " << region.ToString(sc.header);
    lregion = "";
  }
  std::vector<std::pair<int,int>> repeats_to_avoid = svabaUtils::find_repeats(lregion, 20,10);
  SeqLib::GRC local_blacklist;
  for (auto& [key, walker] : unit.walkers) {  
    for (const auto& [st, sp] : repeats_to_avoid) {
      walker->local_blacklist.add(GenomicRegion(region.chr, region.pos1 + st, region.pos1 + sp));
    }
    walker->local_blacklist.MergeOverlappingIntervals();
    walker->local_blacklist.CreateTreeMap();
  }
  if (repeats_to_avoid.size())
    sc.logger.log(sc.opts.verbose > 1, sc.opts.verbose_log, "  found ref repeats at"); //debug
  for (const auto& [st,sp] : repeats_to_avoid) {
    sc.logger.log(sc.opts.verbose > 1, sc.opts.verbose_log, "  ",region.ChrName(sc.header),
		  ":",region.pos1 + st, "-",region.pos1+sp);    
  }
  
  
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
    
    sc.logger.log(sc.opts.verbose > 1, sc.opts.verbose_log, "---running svabaBamWalker", *walker);

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

    //debug
    // for (const auto& [walker_name,cm] : cigmap) {
    //   for (const auto& [key, ccount] : cm) {
    // 	std::cerr << "Walker: " << walker_name << " key: " <<
    // 	  key << " count " << ccount << std::endl;
    //   }
    // }
    
  } //end BAM loop

  sc.logger.log(sc.opts.verbose > 1, sc.opts.verbose_log,
		"...main region reads <case,control> ",
		"<", unit.st.weird_read_count.first, ",",
		unit.st.weird_read_count.second, ">");
  
  // adjust counts and timer
  unit.st.stop("r");

  // get the mate reads, if this is local assembly and has insert-size distro
  GRC mate_regions_for_native;
  if (!region.IsEmpty()) {
    mate_regions_for_native = runMateCollectionLoop(region, unit);
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
  
  // do kmer correction — pooled across all walkers (tumor + normal)
  // so BFC sees the combined k-mer spectrum for better error/signal
  // discrimination. Matches old svaba behavior. Single Train + single
  // ErrorCorrect instead of per-walker, halving the BFC cost.
  //
  // The BFC object lives on the thread unit (unit.pooled_bfc) so its
  // k-mer hash table (1M sub-tables at l_pre=20) persists across
  // regions. On the first Train(), fml_count allocates all 1M tables;
  // on subsequent calls, fml_count_into clears and reuses them via
  // kh_clear (memset of flags), eliminating ~2M malloc/free calls
  // per region — measured at 47% of worker time before this fix.
  if (sc.opts.ecCorrectType == "f") {

    // Lazy-init the per-thread BFC on first use
    if (!unit.pooled_bfc)
      unit.pooled_bfc = std::make_unique<BFC>();

    BFC& bfc = *unit.pooled_bfc;

    // Phase 1: pool all training reads + weird reads
    for (auto& [_, walker] : unit.walkers) {
      for (const auto& u : walker->train_reads)
        bfc.AddSequence(u.Seq, "", u.Name);
      for (const auto& r : walker->reads) {
        if (!r->to_assemble) continue;
        bfc.AddSequence(r->CorrectedSeq(), "", r->UniqueName());
      }
    }

    // Phase 2: train (reuses hash on 2nd+ call), then re-add
    // only correctable reads and run error correction
    bfc.Train();
    bfc.ClearReads();
    for (auto& [_, walker] : unit.walkers) {
      for (const auto& r : walker->reads) {
        if (!r->to_assemble) continue;
        bfc.AddSequence(r->CorrectedSeq(), "", r->UniqueName());
      }
    }
    bfc.ErrorCorrect();

    // Phase 3: drain corrected sequences and distribute back
    std::unordered_map<std::string, std::string> corrected_by_name;
    std::string s, nm;
    while (bfc.GetSequence(s, nm))
      corrected_by_name.emplace(std::move(nm), std::move(s));

    // Clear reads from the BFC (hash stays alive for next region)
    bfc.ClearReads();

    size_t num_reads_corrected = 0;
    for (auto& [_, walker] : unit.walkers) {
      for (auto& r : walker->reads) {
        if (!r->to_assemble) continue;
        auto it = corrected_by_name.find(r->UniqueName());
        if (it != corrected_by_name.end()) {
          r->SetCorrectedSeq(it->second);
          ++num_reads_corrected;
        }
      }
    }

    sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log,
		  "...BFC corrected ", SeqLib::AddCommas(num_reads_corrected));

    unit.st.stop("k");
  }
  
  // SvABA2.0: realign corrected read sequences to the reference so
  // BreakPoint::splitCoverage's "r2c better than native" gate has an
  // apples-to-apples comparison: same corrected sequence, same
  // svaba-internal BWA-MEM parameters, same scoring on both sides.
  //
  // Performance optimization: in production mode (no --dump-reads), we
  // DEFER this to after the r2c alignment loop and only realign reads
  // that actually hit a contig (HasR2C()). The full-reference BWA call
  // is the single most expensive operation (~39% of wall time) due to
  // cache misses on the ~3 GB FM-index, and most reads don't end up
  // supporting any variant. Reads without r2c hits never enter
  // splitCoverage, so their native alignment is never consulted.
  //
  // In --dump-reads mode we realign ALL reads here (before assembly)
  // so the corrected.bam output has full reference alignments for IGV.
  //
  // Reads that didn't go through correction (to_assemble == false)
  // keep their default sentinel (corrected_native_nm == -1) and the
  // gate falls back to the original BAM record for those.
  if (sc.opts.dump_corrected_reads) {
    for (const auto& [_, walker] : unit.walkers) {
      for (const auto& r : walker->reads) {
        if (!r->to_assemble) continue;

        BamRecordPtrVector aln;
        unit.bwa_aligner->alignSequence(r->CorrectedSeq(),
                                        r->UniqueName(),
                                        aln,
                                        false,
                                        0.6,
                                        0); // no secondary alignments

        // Pick the primary (non-supplementary, non-secondary, mapped)
        // record. If no usable record came back, leave the sentinel.
        for (const auto& a : aln) {
          if (!a->MappedFlag())       continue;
          if ( a->SecondaryFlag())    continue;
          if ( a->SupplementaryFlag()) continue;
          int32_t nm = -1;
          a->GetIntTag("NM", nm);
          r->corrected_native_cig = a->GetCigar();
          r->corrected_native_nm  = nm;
          break;
        }

        for (auto& a : aln) unit.all_corrected_reads.push_back(a);
      }
    }
    sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log,
                  "...realigned corrected reads to genome (dump mode, all reads)");
  }
  
  sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log,
		"...running assemblies for region ",
		region); 
  
  // set the contig prefix
  // Include the active assembler tag (fermi / sga) so that contig names
  // carry a record of which engine produced them, e.g.
  //   c_fermi_chr12_18000000_20000000_0
  //   c_sga_chr12_18000000_20000000_0
  string ctg_prefix = string("c_") + svaba::kAssemblerTag + "_" +
    region.ChrName(sc.header) + "_" +
    std::to_string(region.pos1) + "_" +
    std::to_string(region.pos2);
  
  // where to store the AlignedContigs for this region 
  std::unordered_map<std::string, AlignedContig> all_AlignedContigs_this_region;

  // setup the engine and peform assembly
  UnalignedSequenceVector all_unaligned_contigs_this_region;
#ifndef FERMI  
  {
    //    std::cerr << " SGA" << std::endl;
    svabaAssemblerEngine engine(ctg_prefix, sc.opts.sgaErrorRate,
				sc.opts.sgaMinOverlap,
				sc.readlen);

    size_t dbg = 0;
    for (const auto& [_, walker] : unit.walkers) {

      //debug
      // if (region.pos1 == 4287501) {
      // 	for (const auto& r : walker->reads)
      // 	  std::cout << " RSEQ " << r->CorrectedSeq() << std::endl;
      // }
      
      engine.fillReadTable(walker->reads);
      dbg += walker->reads.size();
    }
    
    sc.logger.log(sc.opts.verbose > 1, false, 
		  "...assembling reads with SGA");
    
    // do the actual assembly
    engine.performAssembly(1/*sc.opts.num_assembly_rounds*/);
    
    // retrieve contigs
    all_unaligned_contigs_this_region = engine.getContigs();

    // if (region.pos1 == 4287501) 
    //   std::cerr << " contigs " <<
    // 	all_unaligned_contigs_this_region.size() << std::endl;
  }
#else
    {
      //    std::cerr << "FERMI" << std::endl;      
      // build the reads structure
      size_t n_reads = 0;
      for (const auto& [_, walker] : unit.walkers) {
	for (const auto& r : walker->reads)
	  n_reads += static_cast<size_t>(r->to_assemble);
      }
      fseq1_t *fseq = (fseq1_t*)calloc(n_reads, sizeof(fseq1_t));

      size_t i = 0;
      for (const auto& [_, walker] : unit.walkers) {
	for (const auto& r : walker->reads) {

	  // skip reads that we want to track but not assemble
	  if (!r->to_assemble)
	    continue;
	  
	  const std::string& cseq = r->CorrectedSeq();
	  fseq[i].seq     = strdup(cseq.c_str());
	  fseq[i].l_seq   = cseq.length();

	  // fermi's fml_fltuniq reads qual even when ec_k=-1.
	  // Avoid constructing a std::string just to strdup it —
	  // malloc + memset directly.
	  fseq[i].qual    = (char*)malloc(fseq[i].l_seq + 1);
	  memset(fseq[i].qual, 'I', fseq[i].l_seq);
	  fseq[i].qual[fseq[i].l_seq] = '\0';
	  ++i;
	}
      }

      
      fml_opt_t opt;
      fml_opt_init(&opt);

      // SvABA2.0: loosen fermi-lite defaults for variant-supporting contigs.
      // Rationale (see fermi-lite/misc.c):
      //   - We already hand fermi CorrectedSeq(); a second EC pass with
      //     min_cnt=4 is prone to "correcting" low-VAF alt k-mers back to
      //     reference. Disable EC entirely (ec_k = -1 skips fml_correct).
      //   - fml_assemble clamps min_ensr/min_insr to [min_cnt, max_cnt].
      //     Lowering min_cnt lets tips/branches with 2 supporting reads
      //     survive graph cleaning.
      //   - MAG_F_POPOPEN enables aggressive trimming of "open" tips
      //     (dead-end branches). A heterozygous indel whose alt bubble
      //     is shorter than ~2.5x read length reads exactly like an
      //     open tip and gets trimmed. Clear the flag.
      //   - fml_opt_adjust sets min_elen = 2.5 * avg_read_len (~375bp for
      //     150bp reads). We can't override that through fml_assemble
      //     because it re-runs adjust internally, but disabling POPOPEN
      //     and lowering min_cnt already do most of the work.
      opt.ec_k            = -1;                  // skip fermi's EC (svaba pre-corrects)
      opt.min_cnt         = 2;                   // tolerate rare k-mers
      opt.max_cnt         = 8;                   // keep headroom for high-cov regions
      opt.mag_opt.flag   &= ~MAG_F_POPOPEN;      // don't aggressively trim variant tips
      // optional, very aggressive: opt.mag_opt.flag |= MAG_F_AGGRESSIVE is the
      // OPPOSITE of what we want (it POPS variant bubbles). Leave it cleared.
      
      sc.logger.log(sc.opts.verbose > 1, false, 
		  "...assembling reads with Fermi-lite");
    
      
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
#endif
    
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
    if ((int)i.Seq.length() < (sc.readlen * 1.2)) {
      SVABA_TRACE(i.Name, "TP1 SKIP contig too short: len=" << i.Seq.length()
                  << " threshold=" << (int)(sc.readlen * 1.2));
      continue;
    }
    
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
				    false, SECONDARY_FRAC,
				    0); // don' consider secondary alignments


    SVABA_TRACE(i.Name, "TP2 BWA aligned: " << human_alignments.size()
                << " alignments, contig_len=" << i.Seq.length());

    // provide a svaba-specific continuous "alignment score"
    for (auto& ba : human_alignments) {
      if (!ba) continue;
      if (!ba->MappedFlag()) continue;
      const auto s = svaba::scoreContigAlignment(*ba);
      svaba::tagContigAlignment(*ba, s);
    }
    
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
    
    // make the AlignedContig object for this contig
    AlignedContig ac(human_alignments, region, &sc);

    SVABA_TRACE(ac.getContigName(), "TP5 AlignedContig: frags=" << ac.getFragCount()
                << " global_bp=" << (ac.hasGlobalBP() ? "YES" : "NO")
                << " local_breaks=" << ac.getLocalBreakCount()
                << " hasVariant=" << (ac.hasVariant() ? "YES" : "NO"));

    // add this
    all_AlignedContigs_this_region.insert_or_assign(ac.getContigName(), std::move(ac));
  }
  unit.st.aligned_contig_count = all_aligned_contigs_this_region.size(); 

  sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log, "...contigs that meet size criteria: ",
		count_contigs_of_size);

  // if (region.pos1 == 4287501) {
  //   std::cerr <<" ALIGNED CONTIGS " << 
  //     all_aligned_contigs_this_region.size() << std::endl;
  //   for (const auto& r : all_unaligned_contigs_this_region)
  //     std::cout << "PASSED CONTIG " << r.Name << "\t" <<
  // 	r.Seq << std::endl;
  // }
	
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
  size_t n_ref_match_skipped = 0;

  for (auto& [_, walker] : unit.walkers) {
    for (auto& i : walker->reads) {
      ++seen;
      if (seen % 1000 == 0)
	sc.logger.log(sc.opts.verbose > 4, sc.opts.verbose_log, "- r2c ", SeqLib::AddCommas(seen));

      // get the sequence
      const std::string& seq = i->CorrectedSeq();

      // Fast pre-filter: if the BFC-corrected sequence is an exact
      // substring of the local reference, this read is a perfect
      // reference match and cannot support any variant. Skip the
      // expensive BWA-MEM r2c alignment entirely. This catches the
      // common case where high-NM reads (salvaged via the NM/len >
      // 0.02 path) had their sequencing errors corrected by BFC and
      // now match the reference perfectly. All reads align to the
      // forward strand so no RC check is needed.
      if (!seq.empty() &&
          lregion.find(seq) != std::string::npos) {
        ++n_ref_match_skipped;
        continue;
      }

      // do the read to contig alignment for a single read
      BamRecordPtrVector read2contigs;
      contig_bwa_aligner.alignSequence(seq,
				       i->Qname(),
				       read2contigs,
				       false,
				       0.60,
				       1000);

      // make sure we have only one alignment per contig.
      // Use a small vector + linear scan instead of unordered_set —
      // typical read hits 1-3 contigs, so hashing overhead dominates.
      // Store ChrID (int) instead of contig name strings.
      std::vector<int> seen_chrids;
      seen_chrids.reserve(4);

      // check which ones pass
      auto keep = [&](const BamRecordPtr& r) {
	int thisas = 0;
	r->GetIntTag("AS", thisas);

	// Inline match-base count from CIGAR to avoid GetCigar()
	// vector allocation inside NumMatchBases().
	uint32_t match_bases = 0;
	const uint32_t* cig_raw = bam_get_cigar(r->raw());
	for (uint32_t ci = 0; ci < r->raw()->core.n_cigar; ++ci) {
	  if (bam_cigar_op(cig_raw[ci]) == BAM_CMATCH)
	    match_bases += bam_cigar_oplen(cig_raw[ci]);
	}
	if ((double)match_bases * 0.5 > thisas)
	  return false;

	int chrid = r->ChrID();
	for (int id : seen_chrids)
	  if (id == chrid) return false;

	seen_chrids.push_back(chrid);
	return true;
      };
	
      read2contigs.erase(
			 std::remove_if(read2contigs.begin(), read2contigs.end(),
					[&](const BamRecordPtr& r) { return !keep(r); }),
			 read2contigs.end());

      // SvABA2.0: boundary-aware comma-append helper (matches the
      // dedupe semantics used for the bi:Z tag at BP finalization).
      // `need_read_tags` is declared at function scope above and
      // gates both this block (bz:Z) and the BP-finalization block
      // (bi:Z) further down.
      auto append_unique = [](std::string& cur, const std::string& id) {
	if (id.empty()) return false;
	if (cur.empty()) { cur = id; return true; }
	size_t pos = 0;
	while ((pos = cur.find(id, pos)) != std::string::npos) {
	  const bool left_ok  = (pos == 0) || cur[pos-1] == ',';
	  const bool right_ok = (pos + id.size() == cur.size()) ||
	                         cur[pos + id.size()] == ',';
	  if (left_ok && right_ok) return false;
	  pos += id.size();
	}
	cur += "," + id;
	return true;
      };

      // annotate the original read
      for (auto& r : read2contigs) {
	all_read2contigs.push_back(r); // add to total
	r2c this_r2c; // alignment of this read to this contig
	if (r->ReverseFlag())
	  this_r2c.rc = true;

	this_r2c.AddAlignment(r); // doesn't copy any memory
	const std::string& contig_name = all_unaligned_contigs_this_region[r->ChrID()].Name;
	i->AddR2C(contig_name, this_r2c);

	// SvABA2.0: stamp the bz:Z aux tag on the svabaRead pointer
	// (covers the weird/discordant BAMs, pointer-identity path) and
	// register for later stamping onto the corrected BAM via the
	// qname-keyed map. bz:Z is the "aligned to this contig"
	// superset; bi:Z (set later at BP finalization) is the
	// "supports a variant on this contig" subset.
	//
	// Gated on need_read_tags — the tags and the per-read map are
	// consumed only by the dump-reads BAM writers. Skipping this
	// under the default (--dump-reads off) is measurable on
	// high-coverage contigs; with --dump-reads on, nothing changes.
	if (need_read_tags) {
	  std::string cur;
	  i->GetZTag("bz", cur);
	  if (append_unique(cur, contig_name)) {
	    i->RemoveTag("bz");
	    i->AddZTag("bz", cur);
	  }
	  std::string& entry = unit.all_cnames_by_name[i->UniqueName()];
	  append_unique(entry, contig_name);
	}

	// if (region.pos1 == 4287501) {
	//   std::cout << "R2C " << r->Qname() << " - " <<
	//     r->Sequence() << " contig " << r->ChrName(sc.header) <<
	//     std::endl;
	// }

	// add the read to the right contig
	auto it = all_AlignedContigs_this_region.find(contig_name);
	if (it != all_AlignedContigs_this_region.end()) {
	  it->second.AddAlignedRead(i); // we are adding the svabaBamRead!! not the cread
	}
      } // read2contig alignment loop (per read)
    } // short read realignment loop (all reads)  
  } // short read realignment loop (all bam walkers)
  
  sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log,
		"...done aligning reads to contigs for ",
		SeqLib::AddCommas(all_read2contigs.size()),
		" r2c hits from ", SeqLib::AddCommas(seen), " reads (",
		SeqLib::AddCommas(n_ref_match_skipped),
		" skipped as perfect ref matches). Processing variants");


  // Deferred native realignment (production mode only): now that r2c is
  // done, realign only the reads that actually hit a contig. This avoids
  // the expensive full-reference BWA call for reads that will never enter
  // splitCoverage. In --dump-reads mode this was already done above for
  // ALL reads (before assembly).
  //
  // Additional optimization: when --always-realign-corrected is NOT set
  // (default), reads whose corrected sequence matches the original BAM
  // sequence skip BWA entirely — the input BAM's CIGAR/NM is directly
  // comparable because svaba uses the same BWA-MEM scoring internally.
  // Only reads that BFC actually modified (or that quality trimming
  // shortened) need re-alignment. The sentinel corrected_native_nm == -1
  // signals splitCoverage to fall back to the original BAM record, which
  // is exactly what we want for unchanged reads.
  if (!sc.opts.dump_corrected_reads) {
    // Build a local-reference BWA index for native realignment of
    // BFC-modified reads. The full genome FM-index is ~3 GB and every
    // bwt_sa/bwt_occ call hits main memory (~L3 miss). A local index
    // covering just this region's reference + all mate regions fits
    // in L1/L2 cache, making every BWA operation 10-50× cheaper per call.
    //
    // Pad each region by readlen on each side so reads near the
    // boundaries still find their full native alignment.
    // --always-realign-corrected forces the full genome index for
    // all reads (for non-BWA-aligned input BAMs).
    const int pad = std::max(sc.readlen, 200);

    // Collect reference sequences for the main region + all mate regions.
    // Each region becomes a separate "chromosome" in the local index so
    // reads from any source region can map back to their origin.
    UnalignedSequenceVector local_usv;
    bool local_build_ok = false;

    if (!sc.opts.alwaysRealignCorrected) {
      // Main region
      {
        const int32_t local_start = std::max(0, region.pos1 - pad);
        const int32_t local_end   = region.pos2 + pad;
        try {
          std::string seq = unit.ref_genome->QueryRegion(
              sc.header.IDtoName(region.chr), local_start, local_end);
          if (!seq.empty()) {
            // Name encodes chr:start-end so alignment coords can be
            // interpreted if needed, but for NM/CIGAR purposes we
            // only care about the edit distance, not the absolute
            // position on the genome.
            local_usv.push_back({
              sc.header.IDtoName(region.chr) + ":" +
                std::to_string(local_start) + "-" +
                std::to_string(local_end),
              std::move(seq), std::string()
            });
          }
        } catch (...) { /* skip this region */ }
      }

      // Mate regions — reads from here are in the walkers too and
      // need their own reference to align against.
      for (size_t i = 0; i < mate_regions_for_native.size(); ++i) {
        const auto& mr = mate_regions_for_native[i];
        const int32_t mr_start = std::max(0, mr.pos1 - pad);
        const int32_t mr_end   = mr.pos2 + pad;
        try {
          std::string seq = unit.ref_genome->QueryRegion(
              sc.header.IDtoName(mr.chr), mr_start, mr_end);
          if (!seq.empty()) {
            local_usv.push_back({
              sc.header.IDtoName(mr.chr) + ":" +
                std::to_string(mr_start) + "-" +
                std::to_string(mr_end),
              std::move(seq), std::string()
            });
          }
        } catch (...) { /* skip this mate region */ }
      }
    }

    // Build the local index. If it fails (e.g. empty ref), fall back
    // to the full genome aligner.
    BWAIndexPtr local_native_idx;
    std::unique_ptr<BWAAligner> local_native_aligner;
    if (!local_usv.empty()) {
      try {
        local_native_idx = std::make_shared<BWAIndex>();
        local_native_idx->ConstructIndex(local_usv);
        local_native_aligner = std::make_unique<BWAAligner>(local_native_idx);
        local_build_ok = true;
      } catch (...) {
        local_native_aligner.reset();
      }
    }

    // Choose which aligner to use: local (fast, cache-friendly) or
    // full genome (--always-realign-corrected, or local index failed).
    BWAAligner* native_aligner =
        local_native_aligner ? local_native_aligner.get()
                             : unit.bwa_aligner.get();

    size_t n_realigned = 0;
    size_t n_reused    = 0;
    size_t n_skipped   = 0;
    for (const auto& [_, walker] : unit.walkers) {
      for (const auto& r : walker->reads) {
        if (!r->to_assemble) continue;
        if (!r->HasR2C()) { ++n_skipped; continue; }

        // If --always-realign-corrected is NOT set and the corrected
        // sequence is identical to the original BAM sequence, the input
        // BAM's CIGAR/NM is already a valid native alignment (assuming
        // BWA was used for the input). Leave the sentinel
        // (corrected_native_nm == -1) so splitCoverage falls back to
        // GetCigar()/GetIntTag("NM",...) on the original BAM record.
        // CorrectedSeqChanged() compares against the BAM's 4-bit
        // encoding without allocating a Sequence() string.
        if (!sc.opts.alwaysRealignCorrected &&
            !r->CorrectedSeqChanged()) {
          ++n_reused;
          continue;
        }

        BamRecordPtrVector aln;
        native_aligner->alignSequence(r->CorrectedSeq(),
                                      r->UniqueName(),
                                      aln,
                                      false,
                                      0.6,
                                      0);

        for (const auto& a : aln) {
          if (!a->MappedFlag())       continue;
          if ( a->SecondaryFlag())    continue;
          if ( a->SupplementaryFlag()) continue;
          int32_t nm = -1;
          a->GetIntTag("NM", nm);
          r->corrected_native_cig = a->GetCigar();
          r->corrected_native_nm  = nm;
          break;
        }
        ++n_realigned;
      }
    }
    sc.logger.log(sc.opts.verbose > 0, sc.opts.verbose_log,
                  "...native realignment: ", SeqLib::AddCommas(n_realigned),
                  " re-aligned",
                  (local_build_ok
                    ? std::string(" (local ref, ") + std::to_string(local_usv.size()) + " seqs)"
                    : std::string(" (full genome)")),
                  ", ", SeqLib::AddCommas(n_reused),
                  " reused BAM, ", SeqLib::AddCommas(n_skipped),
                  " no r2c hit");
  }

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
    
    // add to the final structure
    //alc.push_back(std::move(a));
  }
  
  // if (region.pos1 == 4287501) {
  //   for (const auto& [_, bbb] : all_AlignedContigs_this_region) {
  //     std::cout << " ALC " << bbb.printToAlignmentsFile(sc.header) << std::endl;
  //   }
  // }
  
  unit.st.stop("as");

  // get the breakpoints
  BreakPointPtrVector bp_glob;
  for (auto& [_, i] : all_AlignedContigs_this_region) {
    BreakPointPtrVector tmp_allbreaks = i.getAllBreakPoints();
    bp_glob.insert(bp_glob.end(),
		   tmp_allbreaks.begin(),
		   tmp_allbreaks.end());
    //debig
    for (const auto& bp : bp_glob)
      assert(bp);
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
  for (auto& i : bp_glob) {
    i->checkBlacklist(sc.blacklist);

    // remove anything in a high-repeat blacklist
    for (const auto& [_, walker] : unit.walkers) {
      if (walker->local_blacklist.CountOverlaps(i->b1.gr) ||
	  walker->local_blacklist.CountOverlaps(i->b2.gr))
	i->confidence = "BLACKLIST";
    }
  }

  // add discordant reads
  for (auto& i : bp_glob) {
    i->CombineWithDiscordantClusterMap(dmap);
  }

  
  // if (region.pos1 == 4287501) {
  //   for (const auto& bbb : bp_glob) {
  //     std::cout << " BPGLOB " << bbb->printSimple(sc.header) << std::endl;
  //   }
  // }

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
    assert(tmpbp);
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
  bp_glob.erase(new_end, bp_glob.end()); // remove the duplicates after they are sorted to the back
  
  // if (region.pos1 == 4287501) {
  //   for (const auto& bbb : bp_glob) {
  //     std::cout << " AFTER BPGLOB " << bbb->printSimple(sc.header) << std::endl;
  //   }
  // }
  
  // add the coverage data to breaks for allelic fraction computation
  unordered_map<string, STCoverage*> covs;
  for (auto& [key, walker] : unit.walkers) {
    covs[key] = &walker->cov;
  }
  for (auto& i : bp_glob) {
    assert(i); 
    i->addCovs(covs);
  }

  // check the cigar matches
  for (auto& i : bp_glob) {
    i->indelCigarCheck(cigmap);
  }
  
  // score the breakpoints
  for (auto& i : bp_glob) {
    i->scoreBreakpoint(); 
  }
  
  // label somatic breakpoints that intersect directly with normal as NOT somatic
  unordered_set<string> norm_hash;
  for (auto& i : bp_glob) { // hash the normals
    if (i->somatic != SomaticState::SOMATIC_LOD &&
	i->confidence == "PASS") {
      for (const auto& h : i->getBreakEndHashes())
	norm_hash.insert(h);
    }
  }
  
  // find somatic that intersect with norm. Set somatic = 0;
  for (auto& i : bp_glob) {
    // don't need to check normal against normal
    if (i->somatic != SomaticState::SOMATIC_LOD)
      continue;
    
    // if overlaps then turn off somatic
    for (const auto& h : i->getBreakEndHashes())
      if (norm_hash.count(h))
	i->somatic = SomaticState::FAILED; // not somatic
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
      if (i->somatic != SomaticState::SOMATIC_LOD || i->isIndel())
	continue;

      GenomicRegion gr1 = i->BreakEndAsGenomicRegionLeft();
      GenomicRegion gr2 = i->BreakEndAsGenomicRegionRight();
      gr1.Pad(GERMLINE_CNV_PAD);
      gr2.Pad(GERMLINE_CNV_PAD);
      if (sc.germline_svs.OverlapSameInterval(gr1, gr2)) {
	i->somatic = SomaticState::FAILED;
      }
    }
    
  }

  // add the ref and alt tags
  for (auto& i : bp_glob)
    i->setRefAlt(unit.ref_genome.get(), sc.header); 
  
  // transfer local versions to thread store
  for (const auto& [_, a] : all_AlignedContigs_this_region) {
    if (a.hasVariant()) {
      SVABA_TRACE(a.getContigName(), "TP6 hasVariant=YES, pushing to master_alc");
      unit.master_alc.push_back(a); //std::move(a));
    } else {
      SVABA_TRACE(a.getContigName(), "TP6 hasVariant=NO — contig dropped from output"
                  << " frags=" << a.getFragCount()
                  << " indel_breaks=" << a.getIndelBreakCount());
    }
  }
  all_AlignedContigs_this_region.clear();
  //for (const auto& d : dmap)
  //  unit.m_disc_reads += d.second.reads.size();

  //for (const auto& a : alc)
  //  unit.m_bamreads_count += a.NumBamReads();
  for (auto& i : bp_glob) {
    SVABA_TRACE(i->cname, "TP23 hasMinimal check: dc.t=" << i->dc.tcount
                << " dc.n=" << i->dc.ncount << " t.split=" << i->t.split
                << " n.split=" << i->n.split << " local=" << (int)i->local
                << " confidence=" << i->confidence
                << " result=" << (i->hasMinimal() ? "PASS" : "FAIL"));
    if ( i->hasMinimal() ) {

      // SvABA2.0 (v3): assign this BP a unique, thread-stable ID
      // exactly once in its lifetime. Must happen BEFORE the bi:Z
      // stamping below, since the tag payload is now the bp_id (not
      // the cname). The `id.empty()` guard lets this loop be safely
      // reentered (e.g. across regions) without re-numbering BPs
      // that were finalized in a prior region. Because i is a
      // shared_ptr also referenced by the parent AlignedContig's
      // m_global_bp / m_local_breaks / fragment indel_breaks,
      // setting id here is immediately visible to printToR2CTsv
      // when it walks getAllBreakPoints() during writeUnit.
      if (i->id.empty())
        i->id = unit.next_bp_id();

      // SvABA2.0 (v3): tag every alt-supporting read with this
      // breakpoint's bp_id on the `bi:Z` aux tag. bp_id is a
      // thread-stable per-BP identifier (see SvabaThreadUnit::
      // next_bp_id) and is the same string emitted as col 52 of
      // bps.txt.gz and in the r2c.txt.gz `split_bps` / `disc_bps`
      // columns — so `samtools view | grep bi:Z:<bp_id>` returns
      // exactly the ALT-supporters for that specific variant row.
      // Pre-v3 the tag carried cnames (contig-level); a single
      // contig with multiple BPs (global + multi + indel) couldn't
      // be disambiguated. bz:Z still carries cnames (contig-level
      // "which contigs did this read r2c to") — see the r2c loop
      // above where that map is populated.
      //
      // Two parallel paths for bi:Z:
      //   (a) stamp the `bi:Z` aux tag on the original svabaRead
      //       shared_ptr. This covers weird/discordant BAM outputs,
      //       which share pointer identity with those records.
      //   (b) record UniqueName -> bp_ids in
      //       `unit.alt_bp_ids_by_name` so `writeUnit` can stamp
      //       the tag on all_corrected_reads (newly-aligned
      //       BamRecords from bwa that don't share pointer
      //       identity).
      // A read can support multiple BPs (split reads crossing
      // nearby calls; indel + flanking SV on one contig) so the tag
      // value is a comma-joined list, dedup'd on exact-match
      // boundaries.
      auto append_unique = [](std::string& cur, const std::string& id) {
        if (id.empty()) return false;
        if (cur.empty()) { cur = id; return true; }
        size_t pos = 0;
        while ((pos = cur.find(id, pos)) != std::string::npos) {
          const bool left_ok  = (pos == 0) || cur[pos-1] == ',';
          const bool right_ok = (pos + id.size() == cur.size()) ||
                                 cur[pos + id.size()] == ',';
          if (left_ok && right_ok) return false;
          pos += id.size();
        }
        cur += "," + id;
        return true;
      };

      const std::string& this_bp_id = i->id;

      auto tag_with_bp_id = [&](auto& r) {
        if (!r || this_bp_id.empty()) return;
        // (a) stamp the shared_ptr
        std::string cur;
        r->GetZTag("bi", cur);
        if (append_unique(cur, this_bp_id)) {
          r->RemoveTag("bi");
          r->AddZTag("bi", cur);
        }
        // (b) register in name->bp_ids map so the corrected BAM can
        // be tagged at write time
        std::string& entry = unit.alt_bp_ids_by_name[r->UniqueName()];
        append_unique(entry, this_bp_id);
      };

      // Same need_read_tags gate as the bz:Z stamping in the r2c loop
      // above — the bi:Z tag is consumed only by the dump-reads BAM
      // emitters in SvabaOutputWriter, so skipping the per-alt-read
      // GetZTag/RemoveTag/AddZTag triplets under the default (no
      // --dump-reads) path is a measurable CPU win on SVs with many
      // alt supporters. Under --dump-reads nothing changes.
      if (need_read_tags) {
        for (auto& r : i->reads)
          tag_with_bp_id(r);
        for (auto& [_, r] : i->dc.reads)
          tag_with_bp_id(r);
        for (auto& [_, r] : i->dc.mates)
          tag_with_bp_id(r);
      }

      unit.m_bps.push_back(i);
    }
  } // end for bp_glob

  unit.st.bps_count = bp_glob.size();

  // if (region.pos1 == 4287501) {
  //   for (const auto& bbb : unit.m_bps/*bp_glob*/) {
  //     std::cout << " FINAL BPGLOB " << bbb->printSimple(sc.header) << std::endl;
  //   }
  // }
  
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
