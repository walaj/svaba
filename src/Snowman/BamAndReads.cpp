#include "BamAndReads.h"
#include "SnowUtils.h"

using namespace std;

#define BUFF 500
// values above this are discarded
#define MATE_LOOKUP_LIM 20


// divide the larger BAM chunk into smaller reigons for smaller assemblies
void BamAndReads::divideIntoRegions() {

  GenomicRegionVector grv = interval.divideWithOverlaps(5000, 500);
  
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

  BamAlignmentVector reads;
  _read_bam(reads, -1);

  // transfer to assembly regions, place on heap
  for (auto& i : reads)
    addRead(i);

  // end the timer
  read_time = clock() - startr;

}

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

  // if mate is outside, add to discordant region pile
  //if (!mate_in_interval && ba->IsMateMapped()) {
  //  disc.push_back(a);

  //}

  // try the min read structure
  //MinReadSP ms(new MinRead(ba->QueryBases, a->Position, a->RefID, a->Position, a->Position + a->Length));

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
	
	//i->mrv.push_back(ms);
	
	break;
      }
    
}

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


// start a new BamReader, initial a new set of assembly intervals
BamAndReads::BamAndReads(GenomicRegion gr, MiniRulesCollection *tmr, int tverb, string tbam, string tprefix) : 
  interval(gr), mr(tmr), verbose(tverb), bam(tbam), prefix(tprefix) {

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
      if (q.width() > (BUFF * 2+1)) // BUFF * 2 + 1 is width if only one read. 
      grv.push_back(q);
    }
  }

  // merge them down
  grv = GenomicRegion::mergeOverlappingIntervals(grv);
  
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
    if (!reader->SetRegion(i.chr, i.pos1, i.chr, i.pos2)) {
      cerr << "Failed to set MATE region at " << i << " on bam: " << bam << endl;
      exit(EXIT_FAILURE);
    }
    BamAlignmentVector reads;
    _read_bam(reads, 3000); // dont read more than 3000 reads
    
    //debug
    if (verbose > 1)
      cout << "mate reads size " << reads.size() << " on region " << i << endl;

    for (auto& r : reads) 
      addMateRead(r);
  }

  // end the timer
  mate_read_time = clock() - startr;

}

// 
void BamAndReads::_read_bam(BamAlignmentVector &reads, int limit) {

  BamAlignmentVector bam_buffer;
  vector<int> mapq_buffer;

  // prepare to deduplicate reads, based on name and seq
  unordered_map<string, bool> name_map;
  unordered_map<string, bool> seq_map;

  int pileup = 0;

  BamTools::BamAlignment a;
  while (reader->GetNextAlignmentCore(a)) {

    // clear the name because it's just a relic?
    a.Name = ""; 

    // check if read passes rules. 
    string rule_pass = mr->isValid(a);

    // if it makes it, check that it's not a duplicate
    bool pass_rule_and_dup = false;
    if (rule_pass != "") {

      // build it if we haven't
      if (a.Name == "")
	a.BuildCharData();

      // deduplicate by query-bases / position
      string sname = to_string(a.RefID) + "_" + to_string(a.Position) + "_" + to_string(a.MateRefID) + "_" + to_string(a.MatePosition) + a.QueryBases;    
      // deduplicate by Name
      string uname = a.Name + "_" + to_string(a.IsFirstMate());
      
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
      int32_t nm; 
      if (!a.GetTag("NM", nm)) nm = 0;
      if (a.MapQuality == 0 || nm >= 4) 
	pileup++;

      // clear it out
      a.RemoveTag("R2");
      a.RemoveTag("Q2");
      a.RemoveTag("OQ");

      // add a tag to say which region/rule it passes
      a.AddTag("RL","Z",rule_pass);

      // add a tag to give it a unique read name
      a.AddTag("SR", "Z", prefix + to_string(a.AlignmentFlag) + "_" + a.Name);
      
      bam_buffer.push_back(move(a));

      size_t buffer_lim = 100;
      // deal with bam buff
      if (bam_buffer.size() >= buffer_lim) {

	// check if it has too many discordant reads in different directions
	GenomicRegionVector grv_tmp;
	for (auto& it : bam_buffer) {
	  if (it.IsMateMapped() && it.IsMapped())
	    grv_tmp.push_back(GenomicRegion(it.MateRefID, it.MatePosition - 5000, it.MatePosition + 5000));
	  grv_tmp = GenomicRegion::mergeOverlappingIntervals(grv_tmp);
	}
	  
	// check if bad region
	int buf_width = bam_buffer.back().Position - bam_buffer[0].Position;
	if (pileup >= buffer_lim * 0.8 && buf_width <= 40) {
	  /*	  for (auto& it : bam_buffer)
	    if (it.MapQuality > 0) 
	      reads.push_back(move(it));
	  */
		//      addRead(it);
	  if (verbose > 4)
	    cout << "Detected mapq 0 pileup of " << pileup << " at " << a.RefID+1 << ":" << bam_buffer[0].Position << "-" << bam_buffer.back().Position << endl;
	} else if (grv_tmp.size() >= MATE_LOOKUP_LIM) {
	  if (verbose > 4) {
	    cout << "Detected bad discordant pileup with " << grv_tmp.size() << " regions at " << a.RefID+1 << ":" << bam_buffer[0].Position << "-" << bam_buffer.back().Position << endl;
	    for (auto &y : grv_tmp)
	      cout << "      " << y << endl;
	  }
	  // its bad, skip the whole thing
	}
	// it's OK 
	else {
	  for (auto& it : bam_buffer) 
	    reads.push_back(move(it));
	  //	    addRead(it);
	}

	bam_buffer.clear();
	pileup = 0;

	// check the size. If it's too high, kill the reader. -1 means never kill
	if (reads.size() > limit && limit >= 0)
	  return;

      } // end buffer check
      
    } // end save read checking

  } // end read while loop

  // write the final buffer
  for (auto& it : bam_buffer) { 
    reads.push_back(move(it));
    //addRead(it);
  }

}
