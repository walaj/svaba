#include "KmerFilter.h"

#include "ReadTable.h"

int KmerFilter::correctReads(SeqLib::BamRecordVector& vec, SeqLib::BamRecordVector& ref_reads) {

  // first you have to make the index if not there
  if (!pBWT)
    __makeIndex(ref_reads);
  if (!pBWT)
    return 0;

  if (!vec.size())
    return 0;

  int corrected_reads = 0;

  //int intervalCacheLength = 10; // SGA defaul
  //int intervalCacheLength = 1; //vec[0].Length();
  BWTIntervalCache* pIntervalCache = nullptr; //new BWTIntervalCache(intervalCacheLength, pBWT);
  BWTIndexSet indices; //(pBWT, pRBWT, nullptr, pIntervalCache);
  indices.pBWT = pBWT;
  indices.pCache = pIntervalCache;

  KmerCountMap kmerCache;

  for (auto& r : vec) {

    // only correct valid reads
    //if (!r.GetIntTag("VR"))
    //  continue;

    // non-clipped mapped reads with no mismatches are OK (nothing to correct)
    if (r.GetIntTag("NM") == 0 && r.NumClip() == 0 && r.MappedFlag()) 
      continue;

    std::string readSequence = r.QualitySequence(); //QualityTrimmedSequence(4, dum);

    std::string origSequence = readSequence;
    int n = readSequence.length();
    if (n < m_kmer_len) // can't correct, too short
      continue;
    int nk = n - m_kmer_len + 1;
    std::vector<int> minPhredVector(nk, 25); // 25 is a dummy value

    // Are all kmers in the read well-represented?
    bool allSolid = false;
    bool done = false;
    int rounds = 0;
    int maxAttempts = 3; //m_params.numKmerRounds;

    while (!done && nk > 0) {
      // Compute the kmer counts across the read
      // and determine the positions in the read that are not covered by any solid kmers
      // These are the candidate incorrect bases
      std::vector<int> countVector(nk, 0);
      std::vector<int> solidVector(n, 0);

      for(int i = 0; i < nk; ++i)
        {
	  std::string kmer;
	  try { 
	    kmer = readSequence.substr(i, 31);
	  } catch (...) {
	    std::cerr << "KmerFilter substr out of bounds. seqlen " << readSequence.length() << 
	      " stat " << i << " length " << 31 << std::endl;
	  }
	  
	  int count = 0;
	  KmerCountMap::iterator iter = kmerCache.find(kmer);
	  if (iter != kmerCache.end()) {
	    count = iter->second; 
	  } else {
	    count = BWTAlgorithms::countSequenceOccurrences(kmer, indices);
	    kmerCache.insert(std::make_pair(kmer, count));
	  }
	  
	  // Get the phred score for the last base of the kmer
	  int phred = minPhredVector[i];
	  countVector[i] = count;
	  
	  // Determine whether the base is solid or not based on phred scores
	  int threshold = CorrectionThresholds::Instance().getRequiredSupport(phred);
	  if(count >= threshold)
            {
	      for(int j = i; j < i + 31; ++j)
		solidVector[j] = 1;
            }
	}
  

      allSolid = true;
      for(int i = 0; i < n; ++i) {
	if(solidVector[i] != 1)
	  allSolid = false;
      }
      
      // Stop if all kmers are well represented or we have exceeded the number of correction rounds
      if(allSolid || rounds++ > maxAttempts)
	break; 
      
      // Attempt to correct the leftmost potentially incorrect base
      bool corrected = false;
      for(int i = 0; i < n; ++i)
        {
	  if(solidVector[i] != 1)
            {
	      // Attempt to correct the base using the leftmost covering kmer
	      int phred = 25; //workItem.read.getPhredScore(i);
	      int threshold = CorrectionThresholds::Instance().getRequiredSupport(phred);
	      
	      int left_k_idx = (i + 1 >= 31 ? i + 1 - 31 : 0);
	      corrected = attemptKmerCorrection(i, left_k_idx, std::max(countVector[left_k_idx], threshold), readSequence, indices);
	      if(corrected)
		break;
	      
	      // base was not corrected, try using the rightmost covering kmer
	      size_t right_k_idx = std::min(i, n - 31);
	      corrected = attemptKmerCorrection(i, right_k_idx, std::max(countVector[right_k_idx], threshold), readSequence, indices);
	      if(corrected)
		break;
            }
        }

	
      // If no base in the read was corrected, stop the correction process
      if(!corrected)
        {
	  assert(!allSolid);
	  done = true;
        }

    } // end while    
    
    // if allsolid
    if( readSequence != origSequence)
      {
	++corrected_reads;
	assert(readSequence.length());
	r.AddZTag("KC", readSequence);
	//std::cerr << ssi.id << std::endl;
	//std::cerr << "**************Read corrected from\to " << std::endl << "\t" << ssi.seq.toString() << std::endl << "\t" << readSequence << std::endl;
        //result.correctSequence = readSequence;
        //result.kmerQC = true;
      }
    else
      {
	//std::cerr << "Read NOT corrected from\to " << std::endl << "\t" << ssi.seq.toString() << std::endl << "\t" << readSequence << std::endl;
        //result.correctSequence = workItem.read.seq.toString();
        //result.kmerQC = false;
      }
  }

  delete pIntervalCache;

  return corrected_reads;

  return 0;
}


// directly from SGA, Jared Simpson
bool KmerFilter::attemptKmerCorrection(size_t i, size_t k_idx, size_t minCount, std::string& readSequence, BWTIndexSet& inds)
{
  assert(i >= k_idx && i < k_idx + m_kmer_len);
  size_t base_idx = i - k_idx;
  char originalBase = readSequence[i];
  std::string kmer;
  try {
    kmer = readSequence.substr(k_idx, m_kmer_len);
  } catch(...) {
    std::cerr << "KmerFilter::attemptKmerCorrection substr out of bounds. seqlen " << readSequence.length() << 
      " start " << k_idx << " length " << m_kmer_len << std::endl;

  }
  size_t bestCount = 0;
  char bestBase = '$';

#if KMER_TESTING
  std::cout << "i: " << i << " k-idx: " << k_idx << " " << kmer << " " << reverseComplement(kmer) << "\n";
#endif

  for(int j = 0; j < DNA_ALPHABET::size; ++j)
    {
      char currBase = ALPHABET[j];
      if(currBase == originalBase)
	continue;
      kmer[base_idx] = currBase;
      size_t count = BWTAlgorithms::countSequenceOccurrences(kmer, inds);

#if KMER_TESTING
      printf("%c %zu\n", currBase, count);
#endif
      if(count >= minCount)
        {
	  // Multiple corrections exist, do not correct
	  if(bestBase != '$')
	    return false;

	  bestCount = count;
	  bestBase = currBase;
        }
    }

  if(bestCount >= minCount)
    {
      assert(bestBase != '$');
      readSequence[i] = bestBase;
      return true;
    }
  return false;
}

void KmerFilter::makeIndex(const std::vector<char*>& v) {

  ReadTable pRT;
  pRT.setZero();
  
  int dd = 0;
  // make the reads tables
  for (auto& i : v) {
    
    SeqItem si;
    std::string seq(i);

    // if the read is good, add it to the table so we can use for kmer index
    if (seq.length() >= 40 && seq.find("N") == std::string::npos) {
      si.id = std::to_string(dd);
      si.seq = seq;
      pRT.addRead(si);
      ++dd;
    }

  }

  if (pRT.getCount() == 0)
    return;
  
  // make suffix array
  pSAf = new SuffixArray(&pRT, 1, false);
  // make BWT
  pBWT= new RLBWT(pSAf, &pRT);


}

//
void KmerFilter::__makeIndex(SeqLib::BamRecordVector& vec) {

  ReadTable pRT;
  pRT.setZero();

  int dd = 0;
  // make the reads tables
  for (auto& i : vec) {

    SeqItem si;
    std::string sr, seq = "";

    //sr = i.GetZTag("SR");
    seq = i.QualitySequence(); //i.QualityTrimmedSequence(4, dum);
    //seq = i.Sequence();

    // if the read is good, add it to the table so we can use for kmer index
    if (seq.length() >= 40 && seq.find("N") == std::string::npos) {
      si.id = std::to_string(dd);
      si.seq = seq;
      pRT.addRead(si);
      ++dd;
    }

  }

  if (pRT.getCount() == 0)
    return;
  
  // make suffix array
  pSAf = new SuffixArray(&pRT, 1, false);
  // make BWT
  pBWT= new RLBWT(pSAf, &pRT);


}

