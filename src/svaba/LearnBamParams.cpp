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
#include <fstream>
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

  // detect a recurrent 5' soft-clip "tag" per RG (untrimmed UMI/adapter/spacer).
  // Done before dumpLearnData so the learn CSV records it, and it's independent
  // of computeStats (which is isize-only and clears the isize vectors).
  for (auto& br : bam_read_groups)
    br.second.detectTag();

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
    tag_trim_5p = std::max(tag_trim_5p, br.second.tag_trim_5p);  // per-bam max tag
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
    // Visibility guard: a sane WGS discordant cutoff is well under ~10 kb. If
    // it is implausibly large the sample will miss discordant reads under the
    // cutoff (the false-somatic symptom this MAD bugfix addresses) — surface
    // it rather than failing silently.
    if (cutoff > 50000.0)
      sc.logger.log(true, true,
                    "......  WARNING: RG='", rg, "' discordant cutoff ~",
                    static_cast<long>(cutoff),
                    " bp is implausibly large — insert-size distribution may be"
                    " corrupted; discordant reads under this size will be MISSED");
  }

  // Flag any read group with a recurrent 5' soft-clip "tag" (untrimmed
  // UMI / inline barcode / adapter / library spacer). Such a tag mismatches
  // the reference at the read's 5' end on essentially every read, which sits
  // right at the leading edge of each read's overlap with its neighbor and so
  // breaks the overlap (string-graph) assembly genome-wide. svaba will trim it
  // off the read 5' ends before assembly (see sc.tag_trim_5p / TrimTag5p).
  for (const auto& [rg, bp] : bam_read_groups) {
    if (bp.tag_trim_5p <= 0) continue;
    sc.logger.log(true, true,
      "......  WARNING: RG='", rg, "' recurrent ", bp.tag_trim_5p,
      "-bp 5' soft-clip on ", static_cast<int>(bp.tag_frac * 100 + 0.5),
      "% of reads (", (bp.tag_seq.empty()
                         ? std::string("variable bases -> likely UMI")
                         : ("fixed '" + bp.tag_seq + "' -> adapter/primer/spacer")),
      ") -- looks like an untrimmed 5' tag. svaba will trim up to ",
      bp.tag_trim_5p, " bp from read 5' ends before assembly; for best results"
      " trim it upstream (fastp --trim_front1/2, cutadapt -u/-U, or UMI extraction).");
  }
}

void LearnBamParams::dumpLearnData(const std::string& prefix) const {

  std::string fn = prefix + ".learn.tsv.gz";

  ogzstream out(fn.c_str());
  if (!out.good()) {
    sc.logger.log(true, true, "WARNING: could not open ", fn, " for writing learn data");
    return;
  }

  // header (tag_trim_5p / tag_frac are per-RG constants repeated per row so the
  // recurrent-5'-tag finding is recorded in the dumped learning CSV)
  out << "bam\trg\tisize\ttag_trim_5p\ttag_frac\n";

  // one row per isize observation, per RG
  for (const auto& [rg, brg] : bam_read_groups) {
    for (uint32_t is : brg.isize_vec) {
      out << bam_ << '\t' << rg << '\t' << is << '\t'
          << brg.tag_trim_5p << '\t' << brg.tag_frac << '\n';
    }
  }

  out.close();
  sc.logger.log(true, true, "......wrote learning data to ", fn);
}

// ---- insert-size param cache: learn once, reuse across scatter shards ----------
void writeBamParams(const SvabaSharedConfig& sc, const std::string& file) {
  std::ofstream out(file);
  if (!out.good()) {
    sc.logger.log(true, true, "WARNING: cannot open bam-params file for write: ", file);
    return;
  }
  out << "bam\trg\tisize_median\tsd_isize\trg_readlen_max\trg_mapq_max\t"
         "bam_readlen_max\tbam_mapq_max\tbam_isize_max\treads\tn_isize_pairs\ttag_trim_5p\n";
  for (const auto& [prefix, lp] : sc.bamStats) {
    (void)prefix;
    for (const auto& [rg, brg] : lp.bam_read_groups) {
      out << lp.bamPath() << '\t' << rg << '\t'
          << brg.isize_median << '\t' << brg.sd_isize << '\t'
          << brg.readlen_max << '\t' << brg.mapq_max << '\t'
          << lp.readlen_max << '\t' << lp.mapq_max << '\t' << lp.isize_max << '\t'
          << brg.reads << '\t' << brg.n_isize_pairs << '\t' << brg.tag_trim_5p << '\n';
    }
  }
}

bool loadBamParams(SvabaSharedConfig& sc, const std::string& file) {
  std::ifstream in(file);
  if (!in.good()) return false;

  // bamStats is keyed by the run's sample prefix; map the file's PATH -> prefix
  std::unordered_map<std::string, std::string> path2prefix;
  for (const auto& b : sc.opts.bams) path2prefix[b.second] = b.first;

  std::string line;
  std::getline(in, line);  // header
  std::unordered_set<std::string> populated;
  size_t rows = 0;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    std::string bam, rg, med, sd, rrl, rmq, brl, bmq, bis, rds, np, tag;
    if (!std::getline(ss, bam, '\t') || !std::getline(ss, rg, '\t')) continue;
    std::getline(ss, med, '\t'); std::getline(ss, sd, '\t');
    std::getline(ss, rrl, '\t'); std::getline(ss, rmq, '\t');
    std::getline(ss, brl, '\t'); std::getline(ss, bmq, '\t'); std::getline(ss, bis, '\t');
    std::getline(ss, rds, '\t'); std::getline(ss, np, '\t');
    std::getline(ss, tag, '\t');  // tag_trim_5p (optional; absent in older caches)
    auto pit = path2prefix.find(bam);
    if (pit == path2prefix.end()) continue;   // a BAM not in this run
    const std::string& prefix = pit->second;
    auto it = sc.bamStats.find(prefix);
    if (it == sc.bamStats.end())
      it = sc.bamStats.emplace(prefix, LearnBamParams(sc, bam)).first;
    try {
      BamReadGroup& brg = it->second.bam_read_groups[rg];
      brg.isize_median = std::stod(med);
      brg.sd_isize     = std::stod(sd);
      brg.readlen_max  = std::stoi(rrl);
      brg.mapq_max     = std::stoi(rmq);
      if (!rds.empty()) brg.reads = std::stoull(rds);
      if (!np.empty())  brg.n_isize_pairs = std::stoull(np);
      if (!tag.empty()) brg.tag_trim_5p = std::stoi(tag);
      it->second.readlen_max = std::stoi(brl);
      it->second.mapq_max    = std::stoi(bmq);
      it->second.isize_max   = std::stod(bis);
      it->second.tag_trim_5p = std::max(it->second.tag_trim_5p, brg.tag_trim_5p);
    } catch (const std::exception&) { continue; }
    populated.insert(prefix);
    ++rows;
  }
  // require every BAM in the current run to be covered, else learn fresh
  for (const auto& b : sc.opts.bams) {
    if (!populated.count(b.first)) {
      sc.logger.log(true, true, "...bam-params file lacks an entry for ", b.second,
                    "; falling back to full insert-size learning");
      return false;
    }
  }
  return rows > 0;
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

  // accumulate the 5' soft-clip length distribution (for recurrent-tag
  // detection). The read's 5' end is the leading clip for forward reads and
  // the trailing clip for reverse reads (SEQ is stored reference-forward).
  if (r.MappedFlag()) {
    SeqLib::Cigar c = r.GetCigar();
    if (c.size() != 0) {
      ++clip5_total;
      const bool rev = r.ReverseFlag();
      const auto& edge = rev ? c.back() : c.front();
      if (edge.Type() == 'S') {
        int clip5 = static_cast<int>(edge.Length());
        if (clip5 >= 1 && clip5 < static_cast<int>(clip5_hist.size())) {
          ++clip5_hist[clip5];
          // sample the clipped bases from FORWARD reads only (leading bases,
          // no reverse-complement needed) to classify fixed-seq vs random/UMI
          if (!rev && clip5 <= 8) {
            std::string seq = r.Sequence();
            if (static_cast<int>(seq.size()) >= clip5)
              ++clip5_seq[seq.substr(0, clip5)];
          }
        }
      }
    }
  }

}

void BamReadGroup::detectTag() {
  tag_trim_5p = 0; tag_frac = 0.0; tag_seq.clear();
  if (clip5_total < 200) return;          // need a reasonable sample to be confident

  // modal 5' soft-clip length
  int best_len = 0; size_t best_cnt = 0;
  for (int L = 1; L < static_cast<int>(clip5_hist.size()); ++L)
    if (clip5_hist[L] > best_cnt) { best_cnt = clip5_hist[L]; best_len = L; }
  if (best_len == 0) return;

  // fraction of reads at the mode (±1 bp to absorb bwa boundary wobble)
  size_t cnt = clip5_hist[best_len];
  if (best_len - 1 >= 1) cnt += clip5_hist[best_len - 1];
  if (best_len + 1 < static_cast<int>(clip5_hist.size())) cnt += clip5_hist[best_len + 1];
  double frac = static_cast<double>(cnt) / static_cast<double>(clip5_total);

  constexpr double TAG_FRAC = 0.20;       // a clean BAM's 5' clips are sparse + length-spread
  if (frac < TAG_FRAC) return;            // no sharp short-clip spike -> no tag

  tag_trim_5p = best_len;
  tag_frac    = frac;

  // classify fixed-sequence (adapter/spacer) vs variable (UMI): is one clipped
  // sequence of the modal length dominant among the forward-read samples?
  size_t dom = 0, len_total = 0; std::string dom_seq;
  for (const auto& [s, n] : clip5_seq) {
    if (static_cast<int>(s.size()) != best_len) continue;
    len_total += n;
    if (n > dom) { dom = n; dom_seq = s; }
  }
  if (len_total > 0 && static_cast<double>(dom) / static_cast<double>(len_total) > 0.5)
    tag_seq = dom_seq;                    // fixed tag; leave empty => random/UMI
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

    // Robust scale estimate via MAD (median absolute deviation), scaled by
    // 1.4826 to be a Gaussian-equivalent SD.
    //
    // BUGFIX: the previous population SD over the 98%-trimmed values was NOT
    // robust. A normal BAM whose insert-size sample carries a heavy tail
    // (discordant / SV-spanning / mismapped pairs leaking into the midpoint
    // sampling windows) inflated sd_isize to ~150 kb, which pushed the
    // discordant cutoff (isize_median + sd_isize * sdDiscCutoff) to ~535 kb.
    // That silently blinded the normal to ALL discordant reads (even a 172 kb
    // deletion pair fell under the cutoff), so every discordant-supported SV
    // was mis-called somatic. MAD is unaffected by the tail: with the bulk of
    // pairs tight around the median, the median absolute deviation reflects the
    // concordant spread regardless of how large the outliers are.
    std::vector<uint32_t> dev;
    dev.reserve(keep);
    for (uint32_t val : isize_vec)
      dev.push_back(static_cast<uint32_t>(
          std::llround(std::abs(static_cast<double>(val) - isize_median))));
    std::sort(dev.begin(), dev.end());
    double mad = (keep % 2 == 1)
        ? static_cast<double>(dev[keep / 2])
        : (static_cast<double>(dev[keep / 2 - 1]) + dev[keep / 2]) / 2.0;
    sd_isize = 1.4826 * mad;

    // Clean up memory
    isize_vec.clear();

 }
