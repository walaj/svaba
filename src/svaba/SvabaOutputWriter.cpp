// svabaOutputWriter.cpp

#include "SvabaOutputWriter.h"
#include "SvabaThreadUnit.h"
#include "SvabaUtils.h"
#include "BreakPoint.h"
#include "DiscordantCluster.h"
#include "SvabaSharedConfig.h"

#include <htslib/sam.h>

#include <mutex>

using namespace std;

namespace {

  // one mutex to serialize ALL threads writes
  static mutex writeMutex_;

}

// Build a new BamHeader that prepends an @HD line (with the given sort order)
// to the @SQ block of `src`. Avoids SeqLib::BamHeader::AsString(), which
// blindly dereferences h->text; in modern htslib that field is lazy and is
// often NULL until sam_hdr_str() materializes it, causing a segfault.
static SeqLib::BamHeader make_header_with_hd(const SeqLib::BamHeader& src,
                                             const std::string& so) {
  std::string txt = "@HD\tVN:1.6\tSO:" + so + "\n";

  const bam_hdr_t* h = src.get();
  if (h) {
    // sam_hdr_str forces htslib to regenerate the textual form from its
    // internal representation and caches it in h->text.
    const char* cur = sam_hdr_str(const_cast<bam_hdr_t*>(h));
    if (cur) txt += cur;
  }

  return SeqLib::BamHeader(txt);
}

SvabaOutputWriter::SvabaOutputWriter(SvabaLogger& logger_, SvabaOptions& opts_)
  : logger(logger_), opts(opts_)
{ }

// init(): call once, from run_svaba.cpp immediately after you've parsed
// your opts and before you shoot off any threads.
void SvabaOutputWriter::init(const string& analysis_id,
			     const SeqLib::BamHeader& b_header) {

  // keep a local reference to the header
  bam_header_ = b_header;

  // open our gzipped text outputs
  //
  // Per-thread r2c.txt.gz streams are opened by svabaThreadUnit itself
  // (see its ctor) — nothing to open here for that output. What's left
  // on this shared-writer path is only the small coordination outputs
  // where mutex contention is a non-issue.
  svabaUtils::fopen(analysis_id + ".bps.txt.gz",        os_allbps_);
  svabaUtils::fopen(analysis_id + ".discordant.txt.gz",  os_discordant_);
  svabaUtils::fopen(analysis_id + ".runtime.txt",  os_runtime_);

  // write the header line for breakpoints file:
  os_allbps_ << BreakPoint::header();
  for (auto &p : opts.bams) 
    os_allbps_ << "\t" << p.first << "_" << p.second;
  os_allbps_ << "\n";

  // write the header line for runtime file
  os_runtime_ <<
    svabaUtils::svabaTimer::header << "\n";
  
  // write the header line for discordant clusters:
  os_discordant_ << DiscordantCluster::header() << "\n";

  // BAM output for aligned contigs
  std::string aligned_contigs_bam_path = analysis_id;
  aligned_contigs_bam_path.append(".contigs.bam");
  b_contig_writer_ = SeqLib::BamWriter(SeqLib::BAM);
  auto hdr_unsorted = make_header_with_hd(b_header, "unsorted");
  b_contig_writer_.SetHeader(hdr_unsorted);  
  if (!b_contig_writer_.Open(aligned_contigs_bam_path)) {    
    std::cerr << "ERROR: could not open aligned contig writer " <<
      aligned_contigs_bam_path << std::endl;
    exit(EXIT_FAILURE);
  }
  b_contig_writer_.WriteHeader();
  
  
}

void SvabaOutputWriter::writeUnit(svabaThreadUnit& unit,
				  SvabaSharedConfig& sc) {

  // SvABA2.0: stamp a comma-joined aux tag on any record whose
  // UniqueName (svabaRead) or Qname (bwa-aligned corrected record)
  // appears in one of the thread-unit maps. The svabaRead pointer
  // path already tagged weird/discordant records during the r2c and
  // BP-finalization passes, but the corrected records are newly-
  // aligned by bwa and need tagging here. We apply uniformly so the
  // tag is present in every tagged BAM.
  //
  //   bi:Z — bp_ids of variants this read supports as ALT (subset,
  //          v3 — pre-v3 this was cnames; switched to bp_ids so the
  //          tag matches the per-BP resolution of r2c.txt.gz's
  //          split_bps / disc_bps columns and bps.txt.gz's col 52
  //          bp_id. `samtools view | grep bi:Z:<bp_id>` pulls the
  //          ALT-supporters for exactly that variant row.)
  //   bz:Z — cnames of contigs this read r2c'd to (superset; the
  //          contig-level grouping key. A read may r2c to a contig
  //          yet not support any of its BPs.)
  //
  // A read can appear in both maps at the same time (every
  // ALT-supporter has an r2c alignment to its contig). The two
  // tags use different identifier namespaces now: bp_id for bi:Z,
  // cname for bz:Z. Join back to bps.txt.gz via col 30 (cname) or
  // col 52 (bp_id) as appropriate.
  auto stamp_tag = [&](auto& rec, const char* tag,
                       const std::unordered_map<std::string,std::string>& m,
                       const std::string& key) {
    auto it = m.find(key);
    if (it == m.end() || it->second.empty()) return;
    std::string cur;
    rec.GetZTag(tag, cur);
    const std::string& to_add = it->second;
    if (cur.empty()) {
      rec.AddZTag(tag, to_add);
      return;
    }
    // boundary-aware merge/dedupe
    std::string merged = cur;
    size_t s = 0;
    while (s <= to_add.size()) {
      size_t e = to_add.find(',', s);
      if (e == std::string::npos) e = to_add.size();
      std::string tok = to_add.substr(s, e - s);
      if (!tok.empty()) {
        bool present = false;
        size_t pos = 0;
        while ((pos = merged.find(tok, pos)) != std::string::npos) {
          const bool left_ok  = (pos == 0) || merged[pos-1] == ',';
          const bool right_ok = (pos + tok.size() == merged.size()) ||
                                 merged[pos + tok.size()] == ',';
          if (left_ok && right_ok) { present = true; break; }
          pos += tok.size();
        }
        if (!present) merged += "," + tok;
      }
      if (e == to_add.size()) break;
      s = e + 1;
    }
    if (merged != cur) {
      rec.RemoveTag(tag);
      rec.AddZTag(tag, merged);
    }
  };
  auto stamp_bi = [&](auto& rec, const std::string& key) {
    stamp_tag(rec, "bi", unit.alt_bp_ids_by_name, key);
  };
  auto stamp_bz = [&](auto& rec, const std::string& key) {
    stamp_tag(rec, "bz", unit.all_cnames_by_name, key);
  };

  // write the weird reads
  if (sc.opts.dump_weird_reads) {

    auto it = unit.writers.find("w");
    if (it == unit.writers.end()) {
      std::cerr << " BAM writer for weird reads not found" << std::endl;
      exit(EXIT_FAILURE);
    }

    for (const auto& r : unit.all_weird_reads) {

      // stamp bi:Z / bz:Z BEFORE SetQname, while UniqueName is still
      // the map key we indexed under
      stamp_bi(*r, r->UniqueName());
      stamp_bz(*r, r->UniqueName());

      // switch out the qname so it indicates also which BAM it was from
      r->SetQname(r->UniqueName());

      r->RemoveTag("RG");
      r->RemoveTag("PG");
      bool ok = it->second->WriteRecord(*r);
      if (!ok)
	std::cerr << "...unable to write weird read record" << std::endl;
    }
  }

  // write the discordant reads
  if (sc.opts.dump_discordant_reads) {
    auto it = unit.writers.find("d");
    if (it == unit.writers.end()) {
      std::cerr << " BAM writer for discordant reads not found" << std::endl;
    }

    for (const auto& r : unit.all_weird_reads) {
      if (r->dd > 0) {

	// stamp bi:Z / bz:Z before qname mangling
	stamp_bi(*r, r->UniqueName());
	stamp_bz(*r, r->UniqueName());

	// switch out the qname so it indicates also which BAM it was from
	r->SetQname(r->UniqueName());

	// write it
	r->RemoveTag("RG");
	r->RemoveTag("PG");
	bool ok = it->second->WriteRecord(*r);
	if (!ok)
	  std::cerr << "...unable to write discordant read record" << std::endl;
      }
    }

  }

  // write the corrected reads
  if (sc.opts.dump_corrected_reads) {
    auto it = unit.writers.find("c");
    if (it == unit.writers.end()) {
      std::cerr << " read found with prefix " << " c " <<
	" that is not in the weird read writers\n";
      exit(EXIT_FAILURE);
    }

    for (const auto& r : unit.all_corrected_reads) {
      // corrected records are bwa-aligned BamRecords whose Qname was
      // set to the original svabaRead UniqueName at alignSequence()
      // time. That's the key both maps are indexed by.
      stamp_bi(*r, r->Qname());
      stamp_bz(*r, r->Qname());

      r->RemoveTag("RG");
      r->RemoveTag("PG");
      bool ok = it->second->WriteRecord(*r);
      if (!ok)
	std::cerr << "...unable to write corrected read record" << std::endl;
    }
  }

  // SvABA2.0: per-thread r2c.txt.gz write path.
  //
  // Happens BEFORE writeMutex_ is acquired because unit.r2c_out_ is a
  // thread-local gzip stream (not shared across threads) — no coordination
  // needed, and keeping the gzip deflate out from under the shared mutex
  // lets all 16 threads compress in parallel, which is the whole point of
  // the per-thread split (see CLAUDE.md for the full rationale). Mirrors
  // the per-thread BAM writes above, which also run outside the lock.
  //
  // Gated on opts.dump_alignments (--dump-reads) — without it no stream
  // was opened and we skip both the serialization and the write.
  if (opts.dump_alignments) {
    for (const auto& alc : unit.master_alc) {
      if (alc.hasVariant()) {
        // printToR2CTsv emits a contig row followed by one row per
        // r2c-aligned read (already newline-terminated internally),
        // so no trailing newline here.
        // r2c_out_ is unique_ptr<ogzstream>, so dereference. It's non-null
        // here because opts.dump_alignments implies ctor already allocated.
        *unit.r2c_out_ << alc.printToR2CTsv(bam_header_);
      }
    }
  }

  lock_guard<mutex> guard(writeMutex_); // lock the shared writers

  sc.total_regions_done += unit.processed_since_memory_dump;
  unit.processed_since_memory_dump = 0;

  // discordant clusters
  // hardcoding "false" for readtracking for simplicity
  for (const auto& kv : unit.m_disc) {
    const auto& dc = kv.second;
    if (dc.valid())
      os_discordant_ << dc.toFileString(bam_header_, false) << "\n";
  }

  // write contig alignments to BAM
  for (auto& i : unit.master_contigs) {
    bool ok = b_contig_writer_.WriteRecord(*i);
    assert(ok);
  }

  // breakpoints
  for (auto& bp : unit.m_bps) {
    
    if ( bp->hasMinimal() )
      //	 && (bp.confidence != "NOLOCAL" || bp.complex_local) )
      {
	os_allbps_ << bp->toFileString(sc.header) << "\n";
      }
  }
  
  // discordant reads
  //std::cerr << " DC " << unit.all_discordant_reads.size() << std::endl;    
  // if (opts.dump_discordant_reads) {
  //   for (const auto& r : unit.all_discordant_reads) {
  //     bool ok = b_discordant_read_writer_.WriteRecord(*r);
  //     assert(ok);
  //   }
  // }

  // weird reads
  // if (opts.dump_weird_reads) {
  //   for (const auto& r : unit.all_weird_reads) {
  //     bool ok = b_weird_read_writer_.WriteRecord(*r);
  //     assert(ok);
  //   }
  // }

  // // corrected reads
  // if (opts.dump_corrected_reads) {
  //   for (const auto& r : unit.all_corrected_reads) {
  //     bool ok = b_corrected_read_writer_.WriteRecord(*r);
  //     assert(ok);
  //   }
  // }

  // runtime
  os_runtime_ << unit.ss.str();
}

void SvabaOutputWriter::close() {

  // SvABA2.0: alignments.txt.gz is gone entirely; r2c.txt.gz is per-thread
  // (closed in svabaThreadUnit's destructor). Nothing to close here for
  // those outputs.
  os_allbps_.close();
  os_discordant_.close();
  
  // if (opts.dump_discordant_reads) {
  //   if (!b_discordant_read_writer_.Close()) {
  //     std::cerr << "Unable to close discordant read writer" << std::endl;
  //   }
  // }
  
  // if (opts.dump_weird_reads) {
  //   if (!b_weird_read_writer_.Close()) {
  //     std::cerr << "Unable to close weird read writer" << std::endl;
  //   }
  // }
  
  // if (opts.dump_corrected_reads) {
  //   if (!b_corrected_read_writer_.Close()) {
  //     std::cerr << "Unable to close corrected read writer" << std::endl;
  //   }
  // }
  
  if (!b_contig_writer_.Close()) {
    std::cerr << "Unable to close contigs bam writer" << std::endl;
  }
  
}
