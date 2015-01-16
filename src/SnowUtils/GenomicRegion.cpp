#include "GenomicRegion.h"

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
  
static const GenomicRegionVector blacklist = 
{{0,564449,    570371}, 
{ 0,724136,    727043},
{ 0,825006,    825115},
{ 0,2583334,   2634374},
{ 0,4363064,   4363242},
{ 0,5725866,   5736651},
{ 0,16839923,  16841396},
{ 0,38077347,  38077423},
{ 0,91852785,  91853147},
{ 0,104163724, 104163860},
{ 0,108112972, 108113707},
{ 0,121351474, 121487059},
{ 0,142535434, 142543081},
{ 0,142723256, 142723968},
{ 0,142792613, 142793303},
{ 0,142835822, 142837333},
{ 0,143274490, 143284340},
{ 0,145277108, 145277572},
{ 0,149033183, 149035829},
{ 0,156186169, 156186712},
{ 0,224199390, 224204260},
{ 0,233318475, 233318498},
{ 0,236260366, 236260821},
{ 0,237766308, 237766764},
{ 0,238105345, 238105511},
{ 0,238108025, 238108378},
{ 0,238108645, 238109697},
{ 9,18841533, 18862467},
{ 9,20035661, 20037171},
{ 9,36722282, 36723650},
{ 9,38772277, 38819357},
{ 9,38868892, 38889025},
{ 9,39076515, 39155771},
{ 9,42354835, 42548642},
{ 9,42596676, 42602082},
{ 9,42596700, 42602110},
{ 9,42661264, 42667623},
{ 9,42790522, 42818398},
{ 9,135498649,135502716},
{ 10,6831669,  6831838},
{ 10,10529403, 10531969},
{ 10,48671444, 48902406},
{ 10,48931242, 48964015},
{ 10,50318471, 50784078},
{ 10,51090700, 51374066},
{ 10,51567242, 51594226},
{ 10,54694046, 55027975},
{ 10,73221660, 73221946},
{ 10,85194913, 85195322},
{ 10,87524468, 87525005},
{ 10,103275584,103281729},
{ 10,122874443,122874443},
{ 11,20704285, 20704583},
{ 11,34372315, 34372825},
{ 11,34432130, 34857010},
{ 11,37989447, 38441828},
{ 11,38531376, 38531930},
{ 11,41757383, 41757545},
{ 11,127650407,127651075},
{ 11,132061320,132062046},
{ 12,56545728, 56545925},
{ 12,110076444,110076782},
{ 13,18999935, 19056900},
{ 13,32953263, 32954381},
{ 13,84637832, 84639038},
{ 13,90341302, 90341516},
{ 14,19999941, 20044132},
{ 15,32493036, 32570826},
{ 15,32590063, 32598801},
{ 15,33237130, 33241330},
{ 15,33864355, 34023306},
{ 15,34180542, 34197081},
{ 15,34530115, 34542632},
{ 15,35193580, 35285885},
{ 15,46385718, 46456668},
{ 15,46497639, 46500515},
{ 15,47538629, 47539297},
{ 16,19355538, 19356096},
{ 16,19502495, 19506773},
{ 16,21905167, 21906712},
{ 16,22018524, 22032049},
{ 16,22221073, 22263006},
{ 16,25263010, 25268059},
{ 16,25415551, 25417559},
{ 16,31149365, 31149981},
{ 16,33478114, 33478372},
{ 16,41381502, 41382591},
{ 16,41463538, 41464075},
{ 16,41464478, 41465015},
{ 16,41465562, 41467288},
{ 16,51183038, 51183763},
{ 16,55868618, 55868752},
{ 16,75158031, 75158430},
{ 17,96416,    97552},
{ 17,105658,   112233},
{ 17,2842252,  2842356},
{ 17,15393801, 15393992},
{ 17,18510894, 18520356},
{ 17,44126235, 44126593},
{ 17,45379603, 45379864},
{ 17,50319086, 50319301},
{ 17,77772846, 77773065},
{ 18,246006,   247844},
{ 18,22877614, 22877696},
{ 18,23235030, 23235504},
{ 18,24182398, 24186210},
{ 18,24385474, 24633168},
{ 18,27730611, 28262682},
{ 18,36066445, 36066810},
{ 18,36756398, 36800948},
{ 18,37759473, 37797722},
{ 18,44914313, 44916340},
{ 18,44960681, 44962681},
{ 1,739925,    740994},
{ 1,49456729,  49457067},
{ 1,88124390,  88124903},
{ 1,89830421,  89880514},
{ 1,90371401,  90394776},
{ 1,90443001,  90545431},
{ 1,91595080,  91616015},
{ 1,92267428,  92326280},
{ 1,115695017, 115695281},
{ 1,117781085, 117781300},
{ 1,132966248, 132989300},
{ 1,132994855, 133007983},
{ 1,133011824, 133013298},
{ 1,133036250, 133040042},
{ 1,133044095,133045945},
{ 1,143848503,143848792},
{ 1,148022736,148022878},
{ 1,149639207,149639515},
{ 1,156120500,156120610},
{ 1,162135000,162139241},
{ 1,230045426,230045796},
{ 19,26257032,26320267},
{ 19,29517710,29521147},
{ 19,29803876,29833334},
{ 19,55932703,55936114},
{ 19,62916702,62918053},
{ 20,9647205, 9648529},
{ 20,9694896, 9704962},
{ 20,9825451, 9827612},
{ 20,9827612, 9845233},
{ 20,9881895, 9882569},
{ 20,10084922,10088004},
{ 20,10492876,10493049},
{ 20,10599428,10599915},
{ 20,10697886,10860890},
{ 20,11186054,11188131},
{ 20,14338127,14369791},
{ 20,18800575,18800997},
{ 20,27228003,27228242},
{ 20,46796081,46796336},
{ 21,16847814,16862659},
{ 21,18876789,18884510},
{  2,25508897, 25509131},
{  2,73159606, 73161131},
{  2,75696297, 75699304},
{  2,75717841, 75720426},
{  2,80995858, 81014459},
{  2,90311686, 90507410},
{  2,93504815, 93519133},
{  2,96335934, 96337436},
{  2,160665423,160665642},
{  2,196625514,196625860},
{  2,197825427,197834080},
{  3,9987,     12694},
{  3,12276463, 12292424},
{  3,12641862, 12642305},
{  3,21583630, 21583719},
{  3,27732004, 27732240},
{  3,47774268, 47774416},
{  3,49085372, 49342114},
{  3,49488472, 49662085},
{  3,52659961, 52688986},
{  3,56194229, 56194584},
{  3,65473858, 65473941},
{  3,68264186, 68266830},
{  3,70296565, 70296841},
{  3,76807083, 76807320},
{  3,78929660, 78929920},
{  3,156374749,156377226},
{  3,156384860,156387314},
{  3,163342479,163342744},
{  3,190190746,190203442},
{  3,190801869,190802909},
{  3,190943802,190943962},
{  3,190987268,190990949},
{  3,191026302,191044344},
{  4,17517177, 17600940},
{  4,21477364, 21497415},
{  4,34177882, 34197574},
{  4,45908253, 46411114},
{  4,49405493, 49554574},
{  4,71146650, 71146996},
{  4,79945807, 79948223},
{  4,93903068, 93906726},
{  4,97746524, 97746679},
{  4,99381556, 99390873},
{  4,105889063,105889263},
{  4,123095972,123097432},
{  4,134258949,134264271},
{  4,174541634,174542177},
{  5,58735349, 58739031},
{  5,58745955, 58780547},
{  5,61880095, 61944008},
{  5,62189892, 62206612},
{  5,62207809, 62230644},
{  5,62283965, 62284581},
{  5,133593944,133594201},
{  5,137059142,137059326},
{  5,150665074,150665281},
{  5,157731310,157735525},
{  6,43878355, 43878530},
{  6,45291516, 45291740},
{  6,56437808, 56442977},
{  6,57253980, 57254183},
{  6,57255310, 57255444},
{  6,57261829, 57261998},
{  6,57544726, 57556913},
{  6,57811488, 57836990},
{  6,57939184, 58055539},
{  6,61054285, 62454680},
{  6,64059157, 64066183},
{  6,64951348, 64956223},
{  8,44070617, 44070871},
{  8,44873123, 44902307},
{  8,45355954, 45357644},
{  8,45435109, 45443517},
{  8,66494170, 66494805},
{  8,66767710, 66864329},
{  8,66970914, 67005594},
{  8,67315122, 67321036},
{  8,67789868, 67792893},
{  8,68410775, 68435115},
{  8,69677073, 69687998},
{8,69689770,   69711497},
{8,69947961,   70011196},
{ 8,70076144,  70076855},
{ 8,70318723,  70327683},
{8,72653073,   72653572},
{8,78790077,   78790255},
{8,79186574,   79187026},
{8,141019938,  141021783},
{22,55206111,  55206740},
{ 22,55207753, 55208152},
{22,55208300,  55208643},
{ 22,55208980, 55209208},
{22,55209655,  55210006},
{ 22,58330488, 58330843},
{22,58373806,  58373962},
{22,58377680,  58377864},
{22,58415350,  58416387},
{22,58432411,  58432680},
{22,58485887,  58486241},
{22,58488898,  58494528},
{22,58499466,  58504235},
{22,58506076,  58528214},
{22,58528184,  58536883},
{22,58544061,  58582415},
{22,61681834,  61919683},
{22,62003205,  62041580},
{22,83658929,  83659019},
{22,108297348, 108297886},
{22,114959057, 115006437},
{22,125605623, 125607351},
{22,125714985, 125715338},
{22,125864844, 125864980},
{22,125865719, 125865874},
{11,66451267,  66451604}, // JEREMIAH START
{2,196625362,  196626363},
{1,109815579,  109816597},
{1,33141000,   33142700}};

GenomicRegion::GenomicRegion(const string s) {

  // if it's a 1:1234567-2342342342(-) format, parse it this way. Otherwise,
  // parse as if it is a c_1:1234567_asaf string
  if ( (s.find('(') != string::npos) && (s.find(')') != string::npos) ) { // it s a (-) string

    string schr, spos1, spos2, sstrand;
    stringstream iss(s);
    getline(iss, schr, ':');
    getline(iss, spos1, '-');
    getline(iss, spos2, '(');
    getline(iss, sstrand, ')');
    
    chr = stoi(schr)-1;
    pos1 = stoi(spos1);
    pos2 = stoi(spos2);
    strand = sstrand.at(0);
    return;
  }
  

  // set the chromosome
  chr = -1;
  try {
    std::regex reg("^c_([0-9XYM]+):*");
    std::smatch match;
    if (std::regex_search(s.begin(), s.end(), match, reg))
      chr = std::stoi(match[1].str()) - 1;
  } catch (...) {
    std::cerr << "Caught error with contig: " << s << std::endl;
    chr = 0;
  }
  
  // set pos1
  pos1 = -1; 
  try { 
    std::regex creg("^c_[0-9XYM]+:([0-9]+)-*");
    std::smatch cmatch;
    if (std::regex_search(s.begin(), s.end(), cmatch, creg))
      pos1 = std::stoi(cmatch[1].str());
  } catch (...) {
    std::cerr << "Caught error with contig: " << s << std::endl;
    pos1 = 1;
  }
  
  // set pos2
  pos2 = -1;
  try {
    std::regex creg2("^c_[0-9XYM]+:[0-9]+-([0-9]+).*");
    std::smatch cmatch2;
    if (std::regex_search(s.begin(), s.end(), cmatch2, creg2))
      pos2 = std::stoi(cmatch2[1].str());
  } catch (...) {
    std::cerr << "Caught error with contig: " << s << std::endl;
    pos2 = 1;
  }
}

// determine if something overlaps with blacklist regions
int GenomicRegion::blacklistOverlap() const {
  for (GenomicRegionVector::const_iterator it = blacklist.begin(); it != blacklist.end(); it++) 
    if (this->getOverlap(*it) > 0)
      return this->getOverlap(*it);
  return 0;
}

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

//
string GenomicRegion::toString() const {
  stringstream out;
  out << chr  << ":" << pos1 << "-" << pos2;
  return out.str();
}

string GenomicRegion::toStringOffset() const {
  stringstream out;
  out << chrToString(chr)  << ":" << pos1 << "-" << pos2 << "(" << strand << ")"; 
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

  /*
  if (chr < b.chr ||)
    return true;
  if (chr > b.chr)
    return false;

  // ref ids are same
  if (pos1 < b.pos1)
    return true;
  if (pos1 > b.pos1)
    return false;

  //pos1 are the same
  if (pos2 < b.pos2)
    return true;

  // they are the exact same
  return false;
  */
}

std::ostream& operator<<(std::ostream& out, const GenomicRegion& gr) {
  out << gr.toStringOffset();
  return out;
}

// constructor for GenomicRegion that takes strings. Assumes chr string is in 
// natural (1, ..., X) or (chr1, ..., chrX) format. That is, it converts to
// BamTools format with a -1 operation.
GenomicRegion::GenomicRegion(string t_chr, string t_pos1, string t_pos2) {

  chr = GenomicRegion::chrToNumber(t_chr);
  pos1 = stoi(t_pos1);
  pos2 = stoi(t_pos2);
}

// constructor to take a pair of coordinates to define the genomic interval
GenomicRegion::GenomicRegion(int t_chr, int t_pos1, int t_pos2) {
  chr = t_chr;
  pos1 = t_pos1;
  pos2 = t_pos2;
  //abspos1 = GenomicRegion::convertPos(t_chr, t_pos1);
  //abspos2 = GenomicRegion::convertPos(t_chr, t_pos2);
}

// constructor to store just one genomic point
GenomicRegion::GenomicRegion(int t_chr, int t_pos1) {
  chr = t_chr;
  pos1 = t_pos1;
  pos2 = t_pos1;
  //abspos1 = GenomicRegion::convertPos(t_chr, t_pos1);
  //abspos2 = abspos1;
}


// read a row from a BED file
void BEDRow::readNextRow(std::istream & stream) {

  string chr, pos1, pos2;
  stream >> chr >> pos1 >> pos2;

  // initialize the GRanges object
  gr = GenomicRegion(chr, pos1, pos2);

  // Skip until the end of the line.
  stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

// convert a chromosome string into a number
int GenomicRegion::chrToNumber(string ref) {

  // remove the chr identifier if it is there
  if (ref.find("chr") != string::npos)
    ref = ref.substr(3, ref.size() - 3);

  string ref_id = ref;;
  if (ref_id == "X")
    ref_id = "24";
  else if (ref_id == "Y")
    ref_id = "25";
  else if (ref_id == "M")
    ref_id = "26";
  
  int out = -1;
  try {
    out = stoi(ref_id);
  } catch (...) {
    cerr << "Caught error trying to convert " << ref << " to number" << endl;
  }

  assert(out > 0);
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

// convert a chr and pos into an absolute genomic position
slong GenomicRegion::convertPos(unsigned refid, unsigned pos, bool revstrand /* false */) {
  if (refid > 24)
    refid = 24;
  cout << "CHR_CLEN: "<< CHR_CLEN[refid] << " refid " << refid << " pos " << pos << endl;
  slong out = CHR_CLEN[refid] + pos;
  cout << "VAL: " << out << endl;
  return (revstrand ? -out : out);
}

// return regions containing the whole genome
GenomicRegionVector GenomicRegion::getWholeGenome() {
  GenomicRegionVector gv;
  for (int i = 0; i < 24; i++)
    gv.push_back(GenomicRegion(i, 1, CHR_LEN[i]));
  return gv;
}

// returns the intersection between genomic region and vector
/*
GenomicRegionVector GenomicRegion::intersection(GenomicRegionVector &grv) {

  // make the interval tree from the regions
  GenomicIntervalVector vec;
  for (GenomicRegionVector::const_iterator it = grv.begin(); it != grv.end(); it++)
    vec.push_back(GenomicInterval(it->abspos1, it->abspos2, *it));
  GenomicTree tree = GenomicTree(vec);

  // find the overlaps
  GenomicIntervalVector result;
  tree.findOverlapping(abspos1, abspos2, result);

  // fill out the resulting genomic region vector
  GenomicRegionVector out;
  for (GenomicIntervalVector::const_iterator it = result.begin(); it != result.end(); it++)
    out.push_back(it->value);

  return out;
}
*/

/*slong GenomicRegion::RPtoNum(string rp) { 
  
  slong out = 0;

  return out;
  }*/

void DiscordantCluster::addRead(string name) {
  unordered_map<string, bool>::iterator ff = qnames.find(name);
  if (ff == qnames.end())
    qnames.insert(pair<string, bool>(name, true));
  return;
}

double DiscordantCluster::getMeanMapq() const {
  double mean = 0;
  if (mapq.size() > 0)
    mean = accumulate(mapq.begin(), mapq.end(), 0.0) / mapq.size();
  return mean;
}

DiscordantCluster::DiscordantCluster(string tcluster) {
  cluster = tcluster;
  
  // parse out the GenomicRegions
  stringstream iss(tcluster);
  string region_string1, region_string2;
  getline(iss, region_string1, '_');
  getline(iss, region_string2, '_');
  reg1 = GenomicRegion(region_string1);
  reg2 = GenomicRegion(region_string2);
}


// define how to print this to stdout
std::ostream& operator<<(std::ostream& out, const DiscordantCluster& dc) {
  out << dc.reg1.chr+1 << ":" << dc.reg1.pos1 << "(" << dc.reg1.strand << ")" << "--" << 
    dc.reg2.chr+1 << ":" << dc.reg2.pos1 << "(" << dc.reg2.strand << ")" << " Tcount: " << dc.tcount << 
    " Ncount: "  << dc.ncount;
  return out;
}


// define how to print to file
string DiscordantCluster::toFileString() const { 
  string sep = ",";
  stringstream ss;

  // take the position closest to the break
  int pos1 = (reg1.strand=='+') ? reg1.pos2 : reg1.pos1;
  int pos2 = (reg2.strand=='+') ? reg2.pos2 : reg2.pos1;
  double meanmapq = getMeanMapq();

  //chr1, pos1, str1, chr2, pos2, str2, tum, norm, contig
  ss << reg1.chr+1 << sep << pos1 << sep << reg1.strand << sep <<
        reg2.chr+1 << sep << pos2 << sep << reg2.strand << sep <<
    tcount << sep << ncount << sep << contig << sep << meanmapq;
  return ss.str(); 
}

// define how to sort theses
bool DiscordantCluster::operator<(const DiscordantCluster &b) const {

  if (reg1.chr < b.reg1.chr)
    return true;
  if (reg1.pos1 < b.reg1.pos1)
    return true;
  return false;
  
}

// checks whether a GenomicRegion is empty
bool GenomicRegion::isEmpty() const {
  return chr == 0 && pos1 == 0 && pos2 == 0;
}
