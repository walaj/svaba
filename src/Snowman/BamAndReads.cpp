#include "BamAndReads.h"
#include "SnowUtils.h"

using namespace std;

#define BUFF 500

// values above this are discarded
#define MATE_LOOKUP_LIM 10
#define MATE_LOOKUP_WID 100

// divide the larger BAM chunk into smaller reigons for smaller assemblies
void BamAndReads::divideIntoRegions() {

  GenomicRegionVector grv = interval.divideWithOverlaps(littlechunk, window_pad);
  
  tree = GenomicRegion::createTreeMap(grv);
  
  for (auto i : grv) {
    AssemblyRegionUP ar(new AssemblyRegion(i));
    arvec.push_back(move(ar));
  }
  
  return;
}

// read in the BAM file and assign reads to AssemblyRegions
void BamAndReads::readBam() {

  // start the timer
  clock_t startr;
  startr = clock();

  ReadVec this_reads;
  GenomicRegionVector black = _read_bam(this_reads, -1);

  blacklist.insert(blacklist.end(), black.begin(), black.end());

  // transfer to assembly regions, place on heap
  for (auto& i : this_reads)
    addRead(i);

#ifdef HAVE_HTSLIB
  if (hts_itr) {
    hts_itr_destroy(hts_itr);
    hts_itr = NULL;
  }
#endif

  // end the timer
  read_time = clock() - startr;
  
}

// use the interval tree to see if a read should be added to different regions
void BamAndReads::addRead(Read &r) {

  // add to all reads
  string sr;
  r_get_SR(r, sr);
  m_allreads[sr] = r;
  m_all_non_mate_reads[sr] = r;

  unique_reads++;

  GenomicIntervalVector giv; 

  // check read itself
  tree[r_id(r)].findOverlapping(r_pos(r), r_pos(r), giv);
  size_t num_read_intervals = giv.size();

  // check mate to see if we should add this read to mates region too
  // if mate is in an interval, don't add to "partner windows" because we will read this already
  // "partner" windows is for pairmate regions to retreive separately.
  if (r_is_mapped(r) && r_mpos(r) >= interval.pos1 && r_mpos(r) <= interval.pos2 && r_mid(r) == interval.chr)
    tree[r_mid(r)].findOverlapping(r_mpos(r), r_mpos(r), giv);  
  bool mate_in_interval = giv.size() > num_read_intervals;

  // loop through assembly regions, and then place the read there if it or its mate hits
  // TODO break if already found all
  for (auto& i : arvec) 
    for (auto& j : giv)
      if (i->region.pos1 == j.start) { // read or mate hits this assembly region
        reads++;
	i->reads.push_back(r);
	
	// add the partner window to this assembly region
	if (!mate_in_interval && r_is_mmapped(r) && (r_mapq(r) > 0)) // the mate is NOT in an assembly window. It's a partner to look up
	  i->partner_windows.push_back(GenomicRegion(r_mid(r), r_mpos(r)-BUFF, r_mpos(r)+BUFF));
	break;
      }
}

/*
// use the interval tree to see if a read should be added to different regions
void BamAndReads::addRead(BamAlignment &a) {

  BamAlignmentUP ba(new BamAlignment(a));

  unique_reads++;

  GenomicIntervalVector giv; 

  // check read itself
  tree[ba->RefID].findOverlapping(ba->Position, ba->Position, giv);
  size_t num_read_intervals = giv.size();

  // check mate to see if we should add this read to mates region too
  if (ba->IsMateMapped() && ba->MatePosition >= interval.pos1 && ba->MatePosition <= interval.pos2 && ba->MateRefID == interval.chr)
    tree[ba->MateRefID].findOverlapping(ba->MatePosition, ba->MatePosition, giv);  
  bool mate_in_interval = giv.size() > num_read_intervals;

  // loop through assembly regions, and then hits
  // TODO break if already found all
  for (auto& i : arvec) 
    for (auto& j : giv)
      if (i->region.pos1 == j.start) {
        reads++;
	i->reads.push_back(ba);
	
	// add the partner window to this assembly region
	if (!mate_in_interval && ba->IsMateMapped() && ba->MapQuality > 0) {
	  i->partner_windows.push_back(GenomicRegion(ba->MateRefID, ba->MatePosition-BUFF, ba->MatePosition+BUFF));
	}
	break;
      }
}
*/

/*
// use the interval tree to see if a read should be added to different regions
void BamAndReads::addMateRead(BamAlignment &a) {

  mate_unique_reads++;

  BamAlignmentUP ba(new BamAlignment(a));
  // try the min read structure
  //MinReadSP ms(new MinRead(ba->QueryBases, ba->Position, ba->RefID, ba->Position, ba->Position + ba->Length));

  // loop through assembly regions
  // if a read overlaps a partner window for this AssemblyRegion, add it
  for (auto& i : arvec) {
    GenomicIntervalVector giv;
    i->tree_pw[ba->RefID].findOverlapping(ba->Position, ba->Position, giv);

    if (giv.size() > 0) {
      i->reads.push_back(ba);
      mate_reads++;
    }
  }

}
*/

// use the interval tree to see if a read should be added to different regions
void BamAndReads::addMateRead(Read &r) {

  string sr;
  r_get_SR(r, sr);
  m_allreads[sr] = r;
  m_all_mate_reads[sr] = r;

  mate_unique_reads++;

  // loop through assembly regions
  // if a read overlaps a partner window for this AssemblyRegion, add it
  for (auto& i : arvec) {
    GenomicIntervalVector giv;
    i->tree_pw[r_id(r)].findOverlapping(r_pos(r), r_pos(r), giv);

    if (giv.size() > 0) {
      i->reads.push_back(r);
      mate_reads++;
    }
  }

}


// start a new BamReader, initial a new set of assembly intervals
BamAndReads::BamAndReads(GenomicRegion gr, MiniRulesCollection *tmr, int tverb, string tbam, string tprefix, int tlittle, int tpad) : 
  interval(gr), mr(tmr), verbose(tverb), bam(tbam), prefix(tprefix), littlechunk(tlittle), window_pad(tpad) {

  // open the HTS reader
#ifdef HAVE_HTSLIB
  fp = bgzf_open(bam.c_str(), "r"); 

  if (!fp) {
    cerr << "Error using HTS reader on opening " << bam << endl;
    exit(EXIT_FAILURE);
  }
  br = bam_hdr_read(fp);
  // open the header with HTS
  if (!br) {
    cerr << "Error using HTS reader on opening " << bam << endl;
    exit(EXIT_FAILURE);
  }

  //HTS set region
  if (!idx)
    idx = hts_idx_load(bam.c_str(), HTS_FMT_BAI);
  hts_itr = sam_itr_queryi(idx, interval.chr, interval.pos1, interval.pos2);
  if (!hts_itr) {
    std::cerr << "Error: Failed to set region: " << interval << endl; 
    exit(EXIT_FAILURE);
  }
#endif

#ifdef HAVE_BAMTOOLS
  // put the reader on the heap
  reader = new BamReader();

  // open the reader
  if (!reader->Open(bam)) {
    cerr << "Error: Cannot open " << bam << " for reading" << endl;
    exit(EXIT_FAILURE);
  }

  // get the index
  if (!reader->LocateIndex()) {
    
    // try finding it manually
    string bai = bam;
    if (!reader->OpenIndex(bai + ".bai")) {
      bai = SnowUtils::scrubString(bai, ".bam");
      bai += ".bai";
      if (!reader->OpenIndex(bai)) {
	cerr << "Error: Cannot locate index file for " << bam << endl;
	exit(EXIT_FAILURE);
      }
    }

  }

  // set the region
  if (!reader->SetRegion(interval.chr, interval.pos1, interval.chr, interval.pos2)) {
    cerr << "Failed to set the region " << interval << endl;
    exit(EXIT_FAILURE);
  }
#endif
  
  // set set of assembly intervals
  divideIntoRegions();

}

// define how to print this
ostream& operator<<(ostream &out, BamAndReads &bar) {

  out << "Unique reads: " << bar.unique_reads << endl;
  out << "Reads: "        << bar.reads << endl;
  out << "Read time "     << bar.read_time << endl;
  out << "Mate unique reads: " << bar.mate_unique_reads << endl;
  out << "Mate Reads: "        << bar.mate_reads << endl;
  out << "Mate read time "     << bar.mate_read_time << endl;

  return out; 
}

void BamAndReads::calculateMateRegions() {

  GenomicRegionVector grv;

  for (auto& i : arvec) {
    
    // merge the partner windows to find what regions to read
    i->partner_windows = GenomicRegion::mergeOverlappingIntervals(i->partner_windows);
    // make the tree
    i->tree_pw = GenomicRegion::createTreeMap(i->partner_windows);
    
    // make a growing list all partner regions
    for (auto& q : i->partner_windows) {
      if (q.width() > (BUFF * 2+1)) { // BUFF * 2 + 1 is width if only one read. 
	q.pad(500);
	grv.push_back(q);
      }
    }
  }

  // merge them down
  grv = GenomicRegion::mergeOverlappingIntervals(grv);
  
  if (verbose > 3) {
    cout << "Calculated mate regions" << endl;
    for (auto& i : grv)
      cout << i << endl;
  }

  // find out how many of each type go to each
  /*  GenomicIntervalTreeMap this_tree = GenomicRegion::createTreeMap(grv);
  for (auto& i : arvec) {
    for (auto& r : i->reads) {
      
    }
    }*/

  mate_regions.insert(mate_regions.begin(), grv.begin(), grv.end());

}

void BamAndReads::readMateBam() {

  // start the timer
  clock_t startr;
  startr = clock();

  // loop through the mated regions and grab the reads
  for (auto& i : mate_regions) {

#ifdef HAVE_BAMTOOLS
    if (!reader->SetRegion(i.chr, i.pos1, i.chr, i.pos2)) {
      cerr << "Failed to set MATE region at " << i << " on bam: " << bam << endl;
      exit(EXIT_FAILURE);
    }
#endif

#ifdef HAVE_HTSLIB
    hts_itr = sam_itr_queryi(idx, i.chr, i.pos1, i.pos2);
  if (!hts_itr) {
    std::cerr << "Error: Failed to set region: " << i << endl; 
    exit(EXIT_FAILURE);
  }
#endif

    ReadVec this_reads;
    _read_bam(this_reads, 3000); // dont read more than 3000 reads

#ifdef HAVE_HTSLIB    
    if (hts_itr) {
      hts_itr_destroy(hts_itr);
      hts_itr = NULL;
    }
#endif

    if (verbose > 1)
      cout << "READ BAM: mate reads size " << this_reads.size() << " on region " << i << endl;

    for (auto& r : this_reads) 
      addMateRead(r);

  }

  // end the timer
  mate_read_time = clock() - startr;

}

// 
GenomicRegionVector BamAndReads::_read_bam(ReadVec &reads, int limit) {

  GenomicRegionVector black;
  
  ReadVec bam_buffer;

  // prepare to deduplicate reads, based on name and seq
  unordered_map<string, bool> name_map;
  unordered_map<string, bool> seq_map;

  int pileup = 0;
#ifdef HAVE_HTSLIB
  void *dum; // needed by hts_itr_next, for some reason...
#endif 

  for (;;) {

    Read r;
#ifdef HAVE_HTSLIB
    bam1_t* b = bam_init1(); 
    if (hts_itr_next(fp, hts_itr, b, dum) <= 0) {
      bam_destroy1(b);
      break; 
    }
    r = std::shared_ptr<bam1_t> (b, free_delete());
#else
    GET_READ(r);
#endif

    // immediately add cigar
    if (r_cig_size(r) > 1) {

      stringstream ss; 
      int pos = r_pos(r); 
      
      for (int i = 0; i < r_cig_size(r); i++) {
	if (r_cig_type(r,i) == 'D' || r_cig_type(r,i) == 'I') {	
	  //ss << r_id(r) << "_" << pos << "_" << /*r_cig_len(r,i) <<*/ r_cig_type(r, i);
	  ss << r_id(r) << "_" << pos << "_" << r_cig_len(r,i) << r_cig_type(r, i);
	  cigmap[ss.str()]++;
	  ss.str("");

	}
	if (!(r_cig_type(r, i) == 'I') && !(r_cig_type(r, i) == 'S') && !(r_cig_type(r,i) == 'H'))
	  pos += r_cig_len(r, i);
      }

    }
    
    // check if read passes rules. 
    string rule_pass = mr->isValid(r);

    // if it makes it, check that it's not a duplicate
    bool pass_rule_and_dup = false;
    if (rule_pass != "") {

      // deduplicate by query-bases / position
      string sname = to_string(r_id(r)) + "_" + to_string(r_pos(r)) + "_" + to_string(r_mid(r)) + "_" + to_string(r_mpos(r));    
      // deduplicate by Name
      string uname = r_qname(r) + "_" + to_string(r_is_first(r));
      
      // its not already add, insert
      bool uname_pass = false, sname_pass = false;
      if (name_map.count(uname) == 0) { // && seq_map.count(sname) == 0) {  
	uname_pass = true;
	name_map.insert(pair<string, int>(uname, true));
      } 
      if (seq_map.count(sname) == 0) {
	sname_pass = true;
	seq_map.insert(pair<string, int>(sname, true));
      }

      pass_rule_and_dup = uname_pass && sname_pass;

    }

    // read passes rule AND is not a duplicate
    if ( pass_rule_and_dup ) {

      // keep track of pile
      int32_t nm = 0; 
      r_get_int32_tag(r, "NM", nm);
      if (r_mapq(r) == 0 || nm >= 4) 
	pileup++;

      // clear it out
      r_remove_tag(r, "R2");
      r_remove_tag(r, "Q2");
      r_remove_tag(r, "OQ");
      
      r_add_Z_tag(r, "RL", rule_pass);
      string srn =  prefix+to_string(r_flag(r)) + "_" + r_qname(r);
      r_add_Z_tag(r, "SR", srn);

      bam_buffer.push_back(r);

      size_t buffer_lim = 100;
      // deal with bam buff
      if (bam_buffer.size() >= buffer_lim) {

        int pos_width = r_pos(bam_buffer.back()) - r_pos(bam_buffer[0]); //bam_buffer.back().Position - bam_buffer[0].Position;

	// check if it has too many discordant reads in different directions
	GenomicRegionVector grv_tmp;
	for (auto& i : bam_buffer) {
	  if (pos_width <= MATE_LOOKUP_WID) { // if the buffer have ~10x weird read coverage or higher, check for weird discordance
	    if (r_is_pmapped(i))
	      grv_tmp.push_back(GenomicRegion(r_mid(i), r_mpos(i) - 5000, r_mpos(i) + 5000));
	    grv_tmp = GenomicRegion::mergeOverlappingIntervals(grv_tmp);
	  }
	}
	  
	// check if bad region
	if (pileup >= buffer_lim * 0.8 && pos_width <= 20) {
	  for (auto& i : bam_buffer)
	    if (r_mapq(i) > 0) 
	      reads.push_back(i);
	  if (verbose > 4)
	    cout << "Detected mapq 0 pileup of " << pileup << " at " << (r_id(r)+1) << ":" << r_pos(bam_buffer[0]) << "-" << r_pos(bam_buffer.back()) << endl;
	} else if (grv_tmp.size() >= MATE_LOOKUP_LIM) {
	  if (verbose > 4) {
	    cout << "Detected bad discordant pileup with " << grv_tmp.size() << " regions at " << (r_id(r)+1) << ":" << r_pos(bam_buffer[0]) << "-" << r_pos(bam_buffer.back()) << endl;
	    for (auto &y : grv_tmp)
	      cout << "      " << y << endl;
	  }
	  GenomicRegion bl(interval.chr, r_pos(bam_buffer[0]), r_pos(bam_buffer.back()));
	  bl.pad(10);
	  black.push_back(bl);
	  for (auto& it : bam_buffer) {
	    reads.push_back(it);
	  }
	}
	// it's OK 
	else {
	  for (auto& it : bam_buffer) {
	    reads.push_back(it);
	  }
	  //	    addRead(it);
	}

	bam_buffer.clear();
	pileup = 0;

	// check the size. If it's too high, kill the reader. -1 means never kill
	if (reads.size() > limit && limit >= 0)
	  return GenomicRegionVector();

      } // end buffer check
      
    } // end save read checking

  } // end read while loop

  // write the final buffer
  for (auto& it : bam_buffer) { 
    reads.push_back(it);
    //addRead(it);
  }

  return black;
}

void BamAndReads::removeBlacklist(GenomicIntervalTreeMap &bt) {

  size_t before = 0, after = 0;
  for (auto& v : arvec) {
    before += v->reads.size();
    v->removeBlacklist(bt);
    after += v->reads.size();
  }

  if (verbose > 1)
    cout << "before reads: " << before << " after " << after << endl;

}

void AssemblyRegion::removeBlacklist(GenomicIntervalTreeMap &bt) {

  ReadVec new_reads;
  GenomicRegionVector new_partner_windows;

  for (auto& r : reads) {
    GenomicIntervalVector giv;
    bt[r_id(r)].findOverlapping(r_pos(r), r_pos(r), giv);
    bt[r_mid(r)].findOverlapping(r_mpos(r), r_mpos(r), giv);
    if (giv.size() == 0) {
      new_reads.push_back(r);
    }
  }

  reads = new_reads;
  new_reads.clear();

  for (auto& p : partner_windows) {
    GenomicIntervalVector giv;
    bt[p.chr].findOverlapping(p.pos1, p.pos2, giv);
    if (giv.size() == 0)
      new_partner_windows.push_back(p);
  }

  partner_windows = new_partner_windows;
  
}
