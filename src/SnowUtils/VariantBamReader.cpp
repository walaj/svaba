#include "VariantBamReader.h"

using namespace std;
using namespace BamTools;

// Trim the sequence by removing low quality bases from either end
int32_t VariantBamReader::qualityTrimRead(int qualTrim, int32_t &startpoint, shared_ptr<bam1_t> &b) {

    int endpoint = -1; //seq.length();
    startpoint = 0;
    int i = 0; 
    
    uint8_t * qual = bam_get_qual(b);
    
    // get the start point (loop forward)
    while(i < b->core.l_qseq) {
      if (qual[i] >= qualTrim) {
          startpoint = i;
          break;
	}
	i++;
    }

    // get the end point (loop backwards)
    i = b->core.l_qseq - 1; //seq.length() - 1;
    while(i >= 0) {
        if (qual[i] >= qualTrim) { //ps >= qualTrim) {
	  endpoint = i + 1; // endpoint is one past edge
          break;
	}
	i--;
    }

    // check that they aren't all bad
    if (startpoint == 0 && endpoint == -1) {
      //trimmed_seq = 0; //trimmed_seq = "";
      //seq = "";
      //qual = "";
      return 0;
    }


    // Clip the read
    /*    uint8_t * query_seq = bam_get_seq(b);
    for (int i = startpoint; i < (endpoint - startpoint); i++) {
      cout << bam_seqi(query_seq, i) << endl;
      trimmed_seq[i] = BASES[bam_seqi(query_seq, i)];
      }*/
    //    seq =   seq.substr(startpoint, endpoint);
    //qual = qual.substr(startpoint, endpoint)

    return (endpoint - startpoint);

}

// Trim the sequence by removing low quality bases from either end
int32_t VariantBamReader::qualityTrimRead(int qualTrim, int32_t &startpoint, shared_ptr<BamTools::BamAlignment> &b) {

    int endpoint = -1; //seq.length();
    startpoint = 0;
    int i = 0; 
    
    // get the start point (loop forward)
    while(i < b->Qualities.length()) {
        int ps = char2phred(b->Qualities[i]);
      if (ps >= qualTrim) {
          startpoint = i;
          break;
	}
	i++;
    }

    // get the end point (loop backwards)
    i = b->Qualities.length() - 1; //core.l_qseq - 1; //seq.length() - 1;
    while(i >= 0) {
        int ps = char2phred(b->Qualities[i]);
        if (ps >= qualTrim) { //ps >= qualTrim) {
	  endpoint = i + 1; // endpoint is one past edge
          break;
	}
	i--;
    }

    // check that they aren't all bad
    if (startpoint == 0 && endpoint == -1) {
      //trimmed_seq = 0; //trimmed_seq = "";
      //seq = "";
      //qual = "";
      return 0;
    }


    // Clip the read
    /*    uint8_t * query_seq = bam_get_seq(b);
    for (int i = startpoint; i < (endpoint - startpoint); i++) {
      cout << bam_seqi(query_seq, i) << endl;
      trimmed_seq[i] = BASES[bam_seqi(query_seq, i)];
      }*/
    //    seq =   seq.substr(startpoint, endpoint);
    //qual = qual.substr(startpoint, endpoint)

    return (endpoint - startpoint);

}


// Trim the sequence by removing low quality bases from either end
void VariantBamReader::qualityTrimRead(int qualTrim, std::string &seq, std::string &qual) {

    assert(seq.size() == qual.size());

    int endpoint = -1; //seq.length();
    int startpoint = 0;
    int i = 0; 
 
    // get the start point (loop forward)
    while(i < (int)seq.length()) {
        int ps = char2phred(qual[i]);
        if (ps >= qualTrim) {
          startpoint = i;
          break;
	}
	i++;
    }

    // get the end point (loop backwards)
    i = seq.length() - 1;
    while(i >= 0) {
        int ps = char2phred(qual[i]);
        if (ps >= qualTrim) {
          endpoint = i + 1; // endpoint is one past edge
          break;
	}
	i--;
    }
    // check that they aren't all bad
    if (startpoint == 0 && endpoint == -1) {
      seq = "";
      qual = "";
      return;
    }

    // Clip the read
    seq =   seq.substr(startpoint, endpoint);
    qual = qual.substr(startpoint, endpoint);

    return;

}

// obtain the clipped count
unsigned VariantBamReader::getClipCount(BamAlignment &a) {
  
  std::vector<int> clipSize;
  std::vector<int> readPos;
  std::vector<int> genPos;
  a.GetSoftClips(clipSize, readPos, genPos, false);
  
  // get the clip number
  unsigned clipnum = 0;
  for(std::vector<int>::iterator j=clipSize.begin(); j != clipSize.end(); ++j)
    clipnum += *j;
  return clipnum;
}

// read in reads from a BAM in order. If m_writer != NULL write the bam, otherwise
// store in &bav.
/*bool VariantBamReader::writeVariantBam(BamQC &qc, BamAlignmentVectorP bav) {

  ReadCount rc_main;
  ReadCount rc_this;

  int pileup = 0;
  
  BamTools::BamAlignment a * = new BamTools::BamAlignment;

  BamAlignmentVectorP bam_buffer;
  vector<int> mapq_buffer;

  while (m_reader->GetNextAlignmentCore(*a)) {

    // clear the name because it's just a relic?
    a->Name = "";

    rc_this.total++;
    rc_main.total++;
        
    // print the message every 500k reads and check 
    // if we should continue or quit
    if (m_verbose > 0 && rc_main.total % 500000 == 0) {
      printMessage(rc_main, a);

      // zero the counters
      rc_this = ReadCount();

      // kill if seen 25m reads, and it's looking bad
      int perclimit = 50;
      if (rc_main.percent() >= perclimit && rc_main.total > 25000000) { 
	cerr << "This is a a really bad BAM after checking out 25m+ reads. Killing job. Percent weird reads: " << rc_main.percent() << " is above limit of " << perclimit << endl;
	cerr << "Reading in region" << m_region << endl;
	exit(EXIT_FAILURE);
      }
    }

    // check if read passes rules. 
    string rule_pass = m_mr->isValid(a);

    // build the qc
    if (qc.use)
      qc.addRead(a);

    if ( rule_pass != "" && !qc_only ) {

      // build it if we haven't
      if (a->Name == "")
	a->BuildCharData();

      // keep track of pile
      if (a->MapQuality == 0) 
	pileup++;

      // add a tag to say which region/rule it passes
      a->AddTag("RL","Z",rule_pass);

      bam_buffer.push_back(a);

      rc_this.keep++;
      rc_main.keep++;
      
      size_t buffer_lim = 100;
      // deal with bam buff
      if (bam_buffer.size() >= buffer_lim) {
	// check if bad region
	int buf_width = bam_buffer.back().Position - bam_buffer[0].Position;
	if (pileup >= buffer_lim * 0.8 && buf_width <= 40) {
	  for (auto it = bam_buffer.begin(); it != bam_buffer.end(); it++) 
	    if (it->MapQuality > 0) {
	      if (m_writer)
		m_writer->SaveAlignment(*it);
	      else 
		bav.push_back(*it);
	    }
	  if (m_verbose > 2)
	    cout << "Detected mapq 0 pileup of " << pileup << " at " << a->RefID+1 << ":" << bam_buffer[0].Position << "-" << bam_buffer.back().Position << endl;
	} 
	// it's OK or its in full region
	else if (bam_buffer.size() >= buffer_lim) {
	  for (auto it = bam_buffer.begin(); it != bam_buffer.end(); it++) {
	    if (m_writer)
	      m_writer->SaveAlignment(*it);
	    else 
	      bav.push_back(*it);
		  
	  }
	}

	bam_buffer.clear();
	pileup = 0;

      } // end buffer check
      
    } // end save read checking

  } // end read while loop

  // write the final buffer
  for (auto it = bam_buffer.begin(); it != bam_buffer.end(); it++)
    if (m_writer)
      m_writer->SaveAlignment(*it);
    else 
      bav.push_back(*it);
  
  
  // print the final message
  if (m_verbose > 0)
    printMessage(rc_main, a);
  
  return true;

}
*/

bool VariantBamReader::setBamRegion(GenomicRegion gp) {
  m_region = gp;

  // set the region
  if (!m_reader->SetRegion(m_region.chr, m_region.pos1, m_region.chr, m_region.pos2)) {
    std::cerr << "Error: Failed to set region: " << gp << endl; 
    exit(EXIT_FAILURE);
  }

  return true;
}

// closes the BamWriter and makes an index file
void VariantBamReader::MakeIndex() {

  m_writer->Close();
  
  // open the file 
  BamReader reader;
  if (!reader.Open(m_out)) {
    cerr << "Error: Could not open the output BAM to create index " << m_out << endl;
    exit(EXIT_FAILURE);
  }

  // create the index
  if (!reader.CreateIndex()) {
    cerr << "Error: Could not create the output BAM index for " << m_out << endl;
    exit(EXIT_FAILURE);
  }

  reader.Close();

}

// make a new object and put the reader and writer on the heap.
// this is also opens the files for reading/writing and checks
// that they are readable/writable
VariantBamReader::VariantBamReader(string in, string out, MiniRulesCollection *mr, int verbose) {

  m_mr = mr;
  m_bam = in;
  m_out = out;
  m_verbose = verbose;

  // open the reader
  m_reader = new BamReader();
  if (!m_reader->Open(in)) {
    cerr << "Error: Cannot open " << in << " for reading" << endl;
    exit(EXIT_FAILURE);
  }

  // get the index
  if (!m_reader->LocateIndex()) {
    
    // try finding it manually
    string bai = in;
    if (!m_reader->OpenIndex(bai + ".bai")) {
      bai = SnowUtils::scrubString(bai, ".bam");
      bai += ".bai";
      if (!m_reader->OpenIndex(bai)) {
	cerr << "Error: Cannot locate index file for " << in << endl;
	exit(EXIT_FAILURE);
      }
    }

  }

  // open the writer
  m_writer = new BamWriter();
  if (!m_writer->Open(out, m_reader->GetHeaderText(), m_reader->GetReferenceData())) {
    cerr << "Error: Cannot open BAM for writing " << out << endl;
    exit(EXIT_FAILURE);
  }


}

// print position and number of reads and percentage kept directly to stdout
void VariantBamReader::printMessage(const ReadCount &rc_main, const BamAlignment * a) const {

  char buffer[100];
  string posstring = SnowUtils::AddCommas<int>(a->Position);
  sprintf (buffer, "Reading read %11s at position %2s:%-11s. Kept %11s (%2d%%) [running count across whole BAM]",  
	   rc_main.totalString().c_str(), GenomicRegion::chrToString(a->RefID).c_str(), posstring.c_str(),  
	   rc_main.keepString().c_str(), rc_main.percent());
  printf ("%s\n",buffer);

}


// remove duplicate reads by name and alignment flag
void VariantBamReader::deduplicateReads(const BamAlignmentVector &inbav, BamAlignmentVector &outbav) {
  
  unordered_map<string, bool> name_map;
  unordered_map<string, bool> seq_map;

  for (BamAlignmentVector::const_iterator it = inbav.begin(); it != inbav.end(); it++) {

    // deduplicate by query-bases / position
    string sname = to_string(it->RefID) + "_" + to_string(it->Position) + "_" + to_string(it->MateRefID) + "_" + to_string(it->MatePosition) + it->QueryBases;    
    // deduplicate by Name
    string uname = it->Name + "_" + to_string(it->IsFirstMate());

    // its not already add, insert
    if (name_map.count(uname) == 1 && seq_map.count(sname) == 1) {  
      name_map.insert(pair<string, int>(uname, true));
      seq_map.insert(pair<string, int>(sname, true));
      outbav.push_back(*it);
    } 

  }
}
