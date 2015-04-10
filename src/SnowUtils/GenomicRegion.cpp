#include "GenomicRegion.h"
#include "gzstream.h"

static const GenomicRegionVector centromeres =
{
  {0 ,121350000,138000000},
  {1 ,90500000, 96800000},
  {2 ,87900000, 93900000},
  {3 ,48200000, 52700000},
  {4 ,46100000, 50700000},
  {5 ,58700000, 63300000},
  {6 ,57880000, 61980000},
  {7 ,43100000, 48100000},
  {8 ,47300000, 50700000},
  {9 ,38000000, 42700000},
  {10,51600000, 55700000},
  {11,33300000, 38200000},
  {12,1        ,19500000},
  {13,1        ,19100000},
  {14,1        ,20700000},
  {15,34600000 ,46436000},
  {16,22200000 ,25800000},
  {17,15400000, 19000000},
  {18,24400000, 28600000},
  {19,25600000, 29400000},
  {20,10900000, 14300000},
  {21,1        ,17900000},
  {22,58100000 ,63000000},
  {23,11600000 ,13400000},
};

GenomicRegionVector GenomicRegion::non_centromeres = 
{
{0,1,121500000},
{0,138000000,CHR_LEN[0]},
{1,1,90500000},
{1,96800000,CHR_LEN[1]},
{2,1,87900000},
{2,93900000,CHR_LEN[2]},
{3,1,48200000},
{3,52700000,CHR_LEN[3]},
{4,1,46100000},
{4,50700000,CHR_LEN[4]},
{5,1,58700000},
{5,63300000,CHR_LEN[5]},
{6,1,57880000},
{6,61700000,CHR_LEN[6]},
{7,1,43760000},
{7,48100000,CHR_LEN[7]},
{8,1,47300000},
{8,50700000,CHR_LEN[8]},
{9,1,38000000},
{9,42700000,CHR_LEN[9]},
{10,1,51600000},
{10,55700000,CHR_LEN[10]},
{11,1,33300000},
{11,38200000,CHR_LEN[11]},
{12,19500000,CHR_LEN[12]},
{13,19100000,CHR_LEN[13]},
{14,20700000,CHR_LEN[14]},
{15,1,34600000},
{15,46436000,CHR_LEN[15]},
{16,1,22200000},
{16,25800000,CHR_LEN[16]},
{17,1,15400000},
{17,19000000,CHR_LEN[17]},
{18,1,24600000},
{18,28600000,CHR_LEN[18]},
{19,1,25600000},
{19,29400000,CHR_LEN[19]},
{20,1,10900000},
{20,14300000,CHR_LEN[20]},
{21,17900000,CHR_LEN[21]},
{22,1,58100000},
{22,63000000,CHR_LEN[22]}}; // exclude Y below
//{23,1,11600000},
//{23,13400000,250000000}};
  
// return the width of the genomic region
int GenomicRegion::width() const {
  return pos2 - pos1 + 1;
}

// returns 0 for no overlaps, 1 for partial and 2 for complete
int GenomicRegion::getOverlap(const GenomicRegion gr) const {

  if (gr.chr != chr)
    return 0;
  
  bool gr1_in = gr.pos1 >= pos1 && gr.pos1 <= pos2;
  bool gr2_in = gr.pos2 >= pos1 && gr.pos2 <= pos2;
  bool pos1_in = pos1 >= gr.pos1 && pos1 <= gr.pos2;
  bool pos2_in = pos2 >= gr.pos1 && pos2 <= gr.pos2;

  if ( (gr1_in && gr2_in) || (pos1_in && pos2_in) )
    return 2;

  if (gr1_in || gr2_in || pos1_in || pos2_in)
    return 1;

  return 0;

}

// write genomic region to a string
string GenomicRegion::toString() const {
  stringstream out;
  out << chrToString(chr)  << ":" << SnowUtils::AddCommas<int>(pos1) << "-" << SnowUtils::AddCommas<int>(pos2) << "(" << strand << ")"; 
  return out.str();
}

void GenomicRegion::pad(int pad) {
  pos1 = max(1, pos1-pad);
  pos2 = min(pos2+pad, 250000000); // 2500000000 is dummy for now. should be chr end
}


// determine if something overlaps with centromere 
int GenomicRegion::centromereOverlap() const {
  for (GenomicRegionVector::const_iterator it = centromeres.begin(); it != centromeres.end(); it++) 
    if (this->getOverlap(*it) > 0)
      return this->getOverlap(*it);
  return 0;
}

bool GenomicRegion::operator<(const GenomicRegion& b) const {
  return (chr < b.chr) || (chr == b.chr && pos1 < b.pos1) || (chr==b.chr && pos1 == b.pos1 && pos2 < b.pos2);
}

bool GenomicRegion::operator==(const GenomicRegion &b) const {
  return (chr == b.chr && pos1 == b.pos1 && b.pos2 == pos2);
}

bool GenomicRegion::operator<=(const GenomicRegion &b) const {
  return (*this < b || *this == b);
}


//bool GenomicIntervalLessThan(const GenomicInterval& a, const GenomicInterval &b) {
//  return (abs(a.start) < abs(b.start)) || (abs(a.start) == abs(b.start) && abs(a.stop) < abs(b.stop));
//}

std::ostream& operator<<(std::ostream& out, const GenomicInterval& gi) {
  string strand;
  if (gi.start < 0)
    strand = "(+)";
  else
    strand = "(-)";
  out << abs(gi.start) << "-" << abs(gi.stop) << strand;
  return out;
}

std::ostream& operator<<(std::ostream& out, const GenomicRegion& gr) {
  out << gr.toString();
  return out;
}

// constructor for GenomicRegion that takes strings. Assumes chr string is in 
// natural (1, ..., X) or (chr1, ..., chrX) format. That is, it converts to
// BamTools format with a -1 operation.
GenomicRegion::GenomicRegion(string t_chr, string t_pos1, string t_pos2) {

  chr = GenomicRegion::chrToNumber(t_chr);
  try {
    t_pos1 = SnowUtils::scrubString(t_pos1, ",");
    t_pos2 = SnowUtils::scrubString(t_pos2, ",");
    pos1 = stoi(t_pos1);
    pos2 = stoi(t_pos2);
  } catch (...) { 
    cerr << "stoi failed in GenomicRegion constructor. Tried: " << t_pos1 << " " << t_pos2 << endl;
  }
}

// constructor to take a pair of coordinates to define the genomic interval
GenomicRegion::GenomicRegion(int t_chr, int t_pos1, int t_pos2) {
  chr = t_chr;
  pos1 = t_pos1;
  pos2 = t_pos2;
}

// convert a chromosome string into a number
int GenomicRegion::chrToNumber(string ref) {

  // remove the chr identifier if it is there
  if (ref.find("chr") != string::npos)
    ref = ref.substr(3, ref.size() - 3);

  string ref_id = ref;;
  if (ref_id == "X")
    ref_id = "23";
  else if (ref_id == "Y")
    ref_id = "24";
  else if (ref_id == "M" || ref_id == "MT")
    ref_id = "25";
  
  int out = -1;
  try {
    out = stoi(ref_id);
  } catch (...) {
    //cerr << "Caught error trying to convert " << ref << " to number" << endl;
  }

  //assert(out > 0);
  return (out-1); // offset by one becuase chr1 = 0 in BamAlignment coords
}

// convert a chromosome number to a string. Assumes 
// a natural ordering (1, ...), not BamTools ordering (0, ...)
string GenomicRegion::chrToString(int ref) {
  string ref_id;
  if (ref == 22)
    ref_id = "X";
  else if (ref == 23)
    ref_id = "Y";
  else if (ref == 24)
    ref_id = "M";
  else
    ref_id = to_string(ref+1);
  assert(ref_id != "23");
  return ref_id;
}

// return regions containing the whole genome
GenomicRegionVector GenomicRegion::getWholeGenome() {
  GenomicRegionVector gv;
  for (int i = 0; i < 24; i++)
    gv.push_back(GenomicRegion(i, 1, CHR_LEN[i]));
  return gv;
}

// checks whether a GenomicRegion is empty
bool GenomicRegion::isEmpty() const {
  return chr == 0 && pos1 == 0 && pos2 == 0;
}

// parse a region file
void GenomicRegion::regionFileToGRV(string file, int pad, GenomicRegionVector * grv) {

  igzstream iss(file.c_str());
  if (!iss || file.length() == 0) { 
    cerr << "Region file does not exist: " << file << endl;
    exit(EXIT_FAILURE);
  }

  string line;

  // get the header line to check format
  string header, header2;;
  igzstream iss_h(file.c_str());
  getline(iss_h, header, '\n');
  getline(iss_h, header2, '\n');

  ////////////////////////////////////
  // MUTECT CALL STATS
  ////////////////////////////////////
  if (header.find("MuTect") != string::npos) { 

    cout << "Reading MuTect CallStats"  << endl;
    string curr_chr = "dum";
    
    while (getline(iss, line, '\n')) {
      size_t counter = 0;
      string chr, pos, judge;
      istringstream iss_line(line);
      string val;
      if (line.find("KEEP") != string::npos) {
	while(getline(iss_line, val, '\t')) {
	  switch (counter) { 
	  case 0 : chr = val; break; 
	  case 1 : pos = val; break;
	  }
	  if (counter >= 1)
	    break;
	  counter++;
	  
	  if (curr_chr != chr) {
	    cout << "...reading chr" << chr << endl;
	    curr_chr = chr;
	  }

	}
	if (chrToNumber(chr) >= 0) {
	  GenomicRegion gr(chr, pos, pos);
	  gr.pad(pad);
	  grv->push_back(gr);
	}
      } // end "keep" conditional
    } // end main while

  }
  ////////////////////////////////////
  // MuTect2 Indel BED
  ////////////////////////////////////
  else if ( (header.find("track") != string::npos) || (header2.find("ActiveRegions") != string::npos) ) {
    cout << "Reading Mutect2 BED" << endl;
      string curr_chr = "dum";
    while (getline(iss, line, '\n')) {
      size_t counter = 0;
      string chr, pos1, pos2;
      istringstream iss_line(line);
      string val;

      if (line.find("\t1.0") != string::npos) {
	while(getline(iss_line, val, '\t')) {
	  switch (counter) { 
	  case 0 : chr = val; break; 
	  case 1 : pos1 = val; break;
	  case 2 : pos2 = val; break;
	  }
	  if (counter >= 2)
	    break;
	  counter++;
	  
	  if (chr != curr_chr) {
	    cout << "...reading chr" << chr << endl;
	    curr_chr = chr;
	  }

	}
	if (chrToNumber(chr) >= 0) {
	  GenomicRegion gr(chr, pos1, pos2);
	  gr.pad(pad);
	  grv->push_back(gr);
	}
      } // end "keep" conditional
    } // end main while
    
  }
  ////////////////////////////////////
  // VCF file
  ////////////////////////////////////
  else if (header.find("vcf") != string::npos || header.find("VCF") != string::npos) {
    cout << "Parsing VCF file "  << endl;
    while (getline(iss, line, '\n')) {
      if (line.length() > 0) {
	if (line.at(0) != '#') { // its a valid line
	  istringstream iss_this(line);
	  int count = 0;
	  string val, chr, pos;
	  
	  while (getline(iss_this, val, '\t')) {
	    switch (count) {
	    case 0 : chr = val;
	    case 1 : pos = val;
	    }
	    count++;
	  }
	  if (count < 3) {
	    cerr << "Didn't parse VCF line properly: " << line << endl;
	  } else {
	    GenomicRegion gr(chr, pos, pos);
	    gr.pad(pad);
	    grv->push_back(gr);
	  }
	    
	}
      }
    }
  }
  ///////////////////////////////////
  // regular BED file 
  ///////////////////////////////////
  else if (file.find(".bed")) {
    cout << "Reading normal BED" << endl;
    string curr_chr = "dum";
    while (getline(iss, line, '\n')) {
      size_t counter = 0;
      string chr, pos1, pos2;
      istringstream iss_line(line);
      string val;

      if (line.find("#") == string::npos) {
	while(getline(iss_line, val, '\t')) {
	  switch (counter) { 
	  case 0 : chr = val; break; 
	  case 1 : pos1 = val; break;
	  case 2 : pos2 = val; break;
	  }
	  if (counter >= 2)
	    break;
	  counter++;
	  
	  if (chr != curr_chr) {
	    cout << "...reading chr" << chr << endl;
	    curr_chr = chr;
	  }

	}
	if (chrToNumber(chr) >= 0) {
	  GenomicRegion gr(chr, pos1, pos2);
	  gr.pad(pad);
	  grv->push_back(gr);
	}
      } // end "keep" conditional
    } // end main while
    
 
  }
  ////////////////////////////////////
  // csv region file
  ////////////////////////////////////
  else {
    while (getline(iss, line, '\n')) {
      
      size_t counter = 0;
      string chr, pos1, pos2, val;
      
      istringstream iss_line(line);
      while(getline(iss_line, val, ',')) {
	switch (counter) {
	case 0 : chr  = val;  break;
	case 1 : pos1 = val; break;
	case 2 : pos2 = val; break;
	}
	counter++;
      }

      if (chrToNumber(chr) >= 0) {
	GenomicRegion gr(chr, pos1, pos2);
	gr.pad(pad);
	grv->push_back(gr);
      }
    }
  }
    ////////////////////////////////////
    
  }

// reduce a set of GenomicRegions into the minium overlapping set (same as GenomicRanges "reduce")
GenomicRegionVector GenomicRegion::mergeOverlappingIntervals(const GenomicRegionVector &grv) {

  // make the list
  GenomicRegionList intervals(grv.begin(), grv.end());

  intervals.sort();
  GenomicRegionList::iterator inext(intervals.begin());
  ++inext;
  for (GenomicRegionList::iterator i(intervals.begin()), iend(intervals.end()); inext != iend;) {
    if((i->pos2 > inext->pos1) && (i->chr == inext->chr))
      {
	if(i->pos2 >= inext->pos2) intervals.erase(inext++);
	else if(i->pos2 < inext->pos2)
	  { i->pos2 = inext->pos2; intervals.erase(inext++); }
      }
    else { ++i; ++inext; }
  }

  // move it over to a grv
  GenomicRegionVector v{ std::make_move_iterator(std::begin(intervals)), 
      std::make_move_iterator(std::end(intervals)) };
  return v;

}

// convert a GRV into an interval tree (map), where each chrom is its own tree
GenomicIntervalTreeMap GenomicRegion::createTreeMap(const GenomicRegionVector &grv) {

  GenomicRegionVector gg = grv;
  sort(gg.begin(), gg.end());

  GenomicIntervalMap map;
  for (auto it : gg) {
    map[it.chr].push_back(GenomicInterval(it.pos1, it.pos2, it));
  }

  GenomicIntervalTreeMap tree;
  for (auto it : map) 
    tree[it.first] = GenomicIntervalTree(it.second);
  return tree;
}

// send a grv to BED format
void GenomicRegion::sendToBed(const GenomicRegionVector &grv, const string file) {
  
  if (grv.size() ==  0)
    return; 
  ofstream ofile(file.c_str(), ios::out);
  for (auto it : grv)
    ofile << chrToString(it.chr) << "\t" << it.pos1 << "\t" << it.pos2 << endl;
  ofile.close();

}


// divide a region into pieces of width and overlaps
GenomicRegionVector GenomicRegion::divideWithOverlaps(int width, int ovlp) {

  // undefined otherwise
  assert(width > ovlp);

  int start = pos1;
  int end = pos1 + width;

  GenomicRegionVector grv;
  
  // region is smaller than width
  if ( end >= pos2 ) {
    grv.push_back(*this);
    return grv; // return a GenomicRegionVector of itself
  }

  // loop through the sizes until done
  while (end <= pos2) {
    grv.push_back(GenomicRegion(chr, start, end));
    end += width - ovlp; // make the new one
    start += width - ovlp;
  }
  assert(grv.size() > 0);
  
  // finish the last one
  start = grv.back().pos2 - width;
  end = pos2;
  grv.push_back(GenomicRegion(chr, start, end));

  return grv;

}

//
/*uslong GenomicRegion::posToBigPos(int refid, int pos) {
  
  if (refid < 25)
    return 0;
  
  return CHR_CLEN[refid] + pos;
  
  }*/

/*Genomic2DInterval GenomicRegion::makeGenomic2DInterval(GenomicRegion &r1, GenomicRegion &r2) {

  Genomic2DInterval g2i;

  uslong bp1 = GenomicRegion::posToBigPos(r1.refID, r1.pos1);
  uslong bp2 = GenomicRegion::posToBigPos(r1.refID, r1.pos1);

  g2i.start = bp1 * GENOME_LENGTH + bp1;
  g2i.stop  = g2GenomicRegion::posToBigPos(r1.refID, r1.pos1);

  }*/
