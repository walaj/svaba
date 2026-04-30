// LearnBamParams.cpp

#include "LearnBamParams.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/GenomicRegion.h"
#include "SvabaSharedConfig.h"
#include "SvabaOptions.h"
#include "SvabaLogger.h"
#include "gzstream.h"

#include <stdexcept>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <unordered_set>

namespace {

  // helper to extract read group ids from the header
  std::vector<std::string> extractReadGroups(const SeqLib::BamHeader& hdr) {
    std::vector<std::string> groups;
    std::istringstream lines(hdr.AsString());
    std::string line;
    while (std::getline(lines, line)) {
      if (line.rfind("@RG", 0) != 0) continue;
      std::istringstream fields(line);
      std::string field;
      while (std::getline(fields, field, '\t')) {
	if (field.rfind("ID:", 0) == 0) {
	  groups.push_back(field.substr(3));
	  break;
	}
      }
    }
    return groups;
  }

} // anonymous namespace

#include <ostream>

std::ostream& operator<<(std::ostream& os, const BamReadGroup& bg) {
  os << "r="  << bg.reads
     << " s="  << bg.supp
     << " u="  << bg.unmap
     << " qc=" << bg.qcfail
     << " d="  << bg.duplicate
     << " mu=" << bg.mate_unmap
     << " MQ=" << bg.mapq_max
     << " RL=" << bg.readlen_max
     << " mIS=" << bg.isize_median
     << " sdIS=" << bg.sd_isize
     << " n_isize=" << bg.isize_vec.size();
  return os;
}


LearnBamParams::LearnBamParams(SvabaSharedConfig& sc_,
			       const std::string& bamPath)
  : sc(sc_),
    bam_(bamPath)
{

  // make the new BamReader
  reader_ = std::make_shared<SeqLib::BamReader>();

  // open the BAM for learning
  if (!reader_->Open(bam_)) {
    throw std::runtime_error("Could not open bam file " + bam_ + " in LearnBamParams");
  }
}

size_t LearnBamParams::consumeReads(
    size_t max_reads,
    std::unordered_map<std::string, size_t>& rg_count,
    std::unordered_set<std::string>& satisfied_rgs,
    const std::vector<std::string>& groups)
{
  size_t consumed = 0;

  while (auto r = reader_->Next()) {
    if (consumed >= max_reads)
      break;

    // get read group tag
    std::string rg;
    if (!r->GetZTag("RG", rg))
      rg = "NA";

    ++consumed;

    // already saturated this RG — skip
    if (satisfied_rgs.count(rg))
      continue;

    // check if this read pushes the RG over the limit
    size_t& cnt = rg_count[rg];
    if (cnt >= static_cast<size_t>(sc.opts.perRgLearnLimit)) {
      satisfied_rgs.insert(rg);
      // check if all header-declared RGs are satisfied
      if (satisfied_rgs.size() >= groups.size())
	break;
      continue;
    }

    // refer to existing BamReadGroup or make new
    auto& bstats = bam_read_groups[rg];
    bstats.addRead(*r);
    ++cnt;
  }

  return consumed;
}

void LearnBamParams::learnParams() {

  // gather read group ids from the header
  auto hdr = reader_->Header();
  auto groups = extractReadGroups(hdr);
  int n_contigs = hdr.NumSequences();

  sc.logger.log(true, true, "......header declares ",
		groups.size(), " read groups across ",
		n_contigs, " contigs");

  // per-RG counts and satisfaction tracking
  std::unordered_map<std::string, size_t> rg_count;
  std::unordered_set<std::string> satisfied_rgs;
  size_t total_reads = 0;

  // --- Phase 1: multi-region sampling ---
  // Build sampling windows at the midpoint of each reference contig.
  // This ensures we see reads from every part of the genome, catching
  // read groups that only appear on later chromosomes.
  //
  // For standard contigs (chr1-chrY, ~24 contigs), the 1M-read-per-window
  // budget means we sample ~24M reads total in the worst case. For BAMs
  // with many small contigs (alt/decoy/random), we skip contigs shorter
  // than 1 Mb to avoid wasting time on regions with few reads.
  //
  // In human mode (default), we also skip non-standard chromosomes
  // (chrID > MAX_MATE_CHR_ID, i.e. chrM/alt/decoy/random) to focus
  // sampling on the primary assembly where isize/coverage stats are
  // representative. --non-human disables this gate.

  constexpr int32_t MIN_CONTIG_LEN = 1'000'000;   // skip tiny contigs
  constexpr size_t  READS_PER_WINDOW = 1'000'000;  // reads per sampling window
  constexpr int32_t WINDOW_HALF = 500'000;          // +/- 500kb around midpoint

  const int max_sample_chr = sc.opts.maxMateChrID; // 23 (through chrY) or -1 if --non-human

  for (int c = 0; c < n_contigs; ++c) {
    // early exit if all RGs are satisfied
    if (satisfied_rgs.size() >= groups.size())
      break;

    // skip non-standard chromosomes in human mode
    if (max_sample_chr >= 0 && c > max_sample_chr)
      continue;

    int32_t clen = hdr.GetSequenceLength(c);
    if (clen < MIN_CONTIG_LEN)
      continue;

    // sample from the midpoint of the contig
    int32_t mid = clen / 2;
    int32_t start = std::max(0, mid - WINDOW_HALF);
    int32_t end = std::min(clen, mid + WINDOW_HALF);

    SeqLib::GenomicRegion region(c, start, end);
    if (!reader_->SetRegion(region))
      continue; // skip if region query fails (no index?)

    size_t consumed = consumeReads(READS_PER_WINDOW, rg_count, satisfied_rgs, groups);
    total_reads += consumed;

    // log progress every 5 contigs
    if (c % 5 == 0 || c == n_contigs - 1) {
      sc.logger.log(sc.opts.verbose > 0, true,
		    "......sampling ", hdr.IDtoName(c),
		    " - ", SeqLib::AddCommas(total_reads), " reads, ",
		    satisfied_rgs.size(), "/", groups.size(), " RGs satisfied");
    }
  }

  // If some header-declared RGs weren't seen during multi-region
  // sampling, they're likely absent from the BAM entirely (phantom
  // header entries from upstream merges, etc.). Don't waste time on
  // a sequential scan — the multi-region pass already covered every
  // standard chromosome.
  if (satisfied_rgs.size() < groups.size()) {
    sc.logger.log(true, true,
		  "......", satisfied_rgs.size(), "/", groups.size(),
		  " RGs satisfied after multi-region sampling (",
		  SeqLib::AddCommas(total_reads), " reads); ",
		  "skipping sequential fallback — missing RGs likely absent from BAM");
  }

  sc.logger.log(true, true,
		"......learning complete: ", SeqLib::AddCommas(total_reads),
		" reads scanned, ",
		bam_read_groups.size(), " read groups found, ",
		satisfied_rgs.size(), "/", groups.size(), " header RGs satisfied");

  // log any RGs from the header that we never found reads for
  for (const auto& g : groups) {
    if (bam_read_groups.find(g) == bam_read_groups.end()) {
      sc.logger.log(true, true,
		    "......WARNING: read group '", g,
		    "' declared in header but no reads found");
    }
  }

  // dump raw isize data for R plotting BEFORE computeStats clears the vectors.
  // Use the BAM filename stem in the output name so multiple BAMs don't clobber.
  if (!sc.opts.analysisId.empty()) {
    // extract basename: /path/to/foo.bam -> foo
    std::string stem = bam_;
    auto slash = stem.rfind('/');
    if (slash != std::string::npos) stem = stem.substr(slash + 1);
    auto dot = stem.rfind('.');
    if (dot != std::string::npos) stem = stem.substr(0, dot);
    dumpLearnData(sc.opts.analysisId + "." + stem);
  }

  // compute isize stats
  // store the max readlen and mapq for this entire bam across RGs
  for (auto& br : bam_read_groups) {
    br.second.computeStats();
    readlen_max = std::max(readlen_max, br.second.readlen_max);
    mapq_max = std::max(mapq_max, br.second.mapq_max);
    isize_max = std::max(isize_max, br.second.isize_median);
  }

  // Log per-RG insert size stats and the derived discordant cutoff.
  //
  // Method:
  //   1. Sample up to perRgLearnLimit (default 1000) FR read pairs per
  //      read group from the midpoint of each standard chromosome.
  //   2. Discard the top 2% of |isize| values as outliers.
  //   3. Compute median and population SD on the remaining 98%.
  //   4. Discordant cutoff = median + SD * sdDiscCutoff
  //      (sdDiscCutoff defaults to 3.92 for tumor, 3.60 for normal).
  //   5. A read pair with |isize| > cutoff is classified as discordant.
  //
  // The cutoff is applied per-RG in svabaBamWalker::getIsizeCutoff(),
  // used by both readBam (read extraction) and TagDiscordant (read
  // tagging).  Normal uses a lower multiplier (sdDiscCutoffNormal)
  // so it's more sensitive to marginal discordant reads, preventing
  // germline events from being mis-called as somatic.
  sc.logger.log(true, true,
                "......insert size learning summary (per read group):");
  sc.logger.log(true, true,
                "......  method: median + SD * sdDiscCutoff on bottom 98% of FR pairs");
  sc.logger.log(true, true,
                "......  sdDiscCutoff=", sc.opts.sdDiscCutoff,
                " (tumor), sdDiscCutoffNormal=", sc.opts.sdDiscCutoffNormal,
                " (normal); up to ", sc.opts.perRgLearnLimit, " pairs/RG");
  for (const auto& [rg, bp] : bam_read_groups) {
    double cutoff = bp.isize_median + bp.sd_isize * sc.opts.sdDiscCutoff;
    sc.logger.log(true, true,
                  "......  RG='", rg,
                  "' n_pairs=", bp.n_isize_pairs,
                  " median=", static_cast<int>(bp.isize_median),
                  " sd=", static_cast<int>(bp.sd_isize),
                  " disc_cutoff(tumor)=", static_cast<int>(cutoff),
                  " readlen=", bp.readlen_max);
  }
}

void LearnBamParams::dumpLearnData(const std::string& prefix) const {

  std::string fn = prefix + ".learn.tsv.gz";

  ogzstream out(fn.c_str());
  if (!out.good()) {
    sc.logger.log(true, true, "WARNING: could not open ", fn, " for writing learn data");
    return;
  }

  // header
  out << "bam\trg\tisize\n";

  // one row per isize observation, per RG
  for (const auto& [rg, brg] : bam_read_groups) {
    for (uint32_t is : brg.isize_vec) {
      out << bam_ << '\t' << rg << '\t' << is << '\n';
    }
  }

  out.close();
  sc.logger.log(true, true, "......wrote learning data to ", fn);
}

void BamReadGroup::addRead(const SeqLib::BamRecord &r)
{

  // track the flags
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

  // track the mapq
  mapq_max = std::max(mapq_max, r.MapQuality());

  // track the insert size
  if (!r.PairMappedFlag() || r.Interchromosomal() || r.PairOrientation() != SeqLib::Orientation::FR)
    ;
  else if (!r.Interchromosomal())
    isize_vec.push_back(std::abs(r.InsertSize()));

  // track the read length
  readlen_max = std::max(readlen_max, r.Length());

}

 void BamReadGroup::computeStats() {

   if (isize_vec.empty()) {
     isize_median = 0.0;
     sd_isize     = 0.0;
     return;
   }

   // Remove the top 2% of values to filter out extreme outliers
   std::sort(isize_vec.begin(), isize_vec.end());
   size_t n = isize_vec.size();
   size_t keep = static_cast<size_t>(std::floor(n * 0.98));
   if (keep == 0) {
     // If 98% rounds down to 0, keep at least one element
     keep = 1;
   }
   isize_vec.resize(keep);

   n_isize_pairs = keep;

   // Calculate median (vector is already sorted)
   if (keep % 2 == 1) {
     isize_median = isize_vec[keep / 2];
   } else {
     isize_median = (isize_vec[keep / 2 - 1] + isize_vec[keep / 2]) / 2.0;
   }

    // Calculate population standard deviation around the median
    double sq_sum = 0.0;
    for (uint32_t val : isize_vec) {
      double diff = val - isize_median;
        sq_sum += diff * diff;
    }
    sd_isize = std::sqrt(sq_sum / keep);

    // Clean up memory
    isize_vec.clear();

 }
