//-----------------------------------------------
// Copyright 201- Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
//        and Kees Albers (caa@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// DindelHMM Realigns read against haplotype using 
// the homopolymer aware HMM described in 
// Albers et al. (2010) Genome Research
//
#include <iostream>
#include <string>
#include <cmath>
#include <cassert>

#include "DindelRealignWindow.h"
#include "DindelHMM.h"
const int DINDEL_DEBUG=0;

#ifdef DINDELHMMSTANDALONE
class DindelRead
{
	public:
		DindelRead(const std::string & seq) { m_seq = seq; };

		int length() const { return (int) m_seq.length(); };
		double getQual(int b) const { return 20.0; };
		char getBase(int b) const { return m_seq[b]; };
		std::string m_seq;
};

class DindelHaplotype
{
	public:
		DindelHaplotype(const std::string & seq) { m_seq = seq; };
		int length() const { return (int) m_seq.length(); };
		const std::string & getSequence() const { return m_seq; }
		int getHomopolymerLength(int b) const { return 1; }
 		std::string m_seq;

};
char complement(char b) { return b;}
#endif


DindelHMM::DindelHMM(DindelRead & read, const DindelMultiHaplotype & haplotype) : m_pRead(&read), m_pHaplotype(& haplotype)
{
        
	if (DINDEL_DEBUG) std::cerr << "CALLED DindelHMM::DindelHMM " << std::endl;
}

void DindelHMM::getSeedPositions(std::set<int> & positions)
{
    positions.clear();

    const std::vector<unsigned int> & readKeys = m_pRead->getHashKeys();

    const DindelSequenceHash & hapHash = m_pHaplotype->getHash();

    HashMap<int, int> pos_to_freq;
    std::map<int, std::list<int> > freq_to_pos;
    for (int x=0;x<int(readKeys.size());x++)
    {
        std::list<int> hpos_set = hapHash.lookupList(readKeys[x]);
        if (!hpos_set.empty())
        {
            int numHits=0;
            for (std::list<int>::const_iterator hpos = hpos_set.begin();hpos!=hpos_set.end() && numHits<4;++hpos,++numHits) {
                int rpos = *hpos-x;
                HashMap<int, int>::iterator it = pos_to_freq.find(rpos);
                if (it == pos_to_freq.end()) pos_to_freq[rpos]=1; else it->second++;
            }
        }
    }

    // convert

    for (HashMap<int, int>::const_iterator it = pos_to_freq.begin();it!=pos_to_freq.end();it++)
	    freq_to_pos[it->second].push_back(it->first);

    // output two best hits
    if (DINDEL_DEBUG) std::cerr << " SeedPositionsREAD: " << m_pRead->getID();
    if (DINDEL_DEBUG) std::cerr <<  " seedpositions " << m_pRead->getID();
    int c=0;
    for (std::map<int, std::list<int> >::const_reverse_iterator it = freq_to_pos.rbegin();it != freq_to_pos.rend(); it++)
    {
	    if (DINDEL_DEBUG) std::cerr << "freq: " << it->first << " pos: ";
	    for (std::list<int>::const_iterator it2 = it->second.begin();it2!=it->second.end();it2++)
	    {
		if (DINDEL_DEBUG) 
		{
			std::cerr << " (" << *it2 << " " << *it2+m_pHaplotype->getSingleMappingHaplotype(0).getRefStart() << ")\n";
			std::cerr << m_pHaplotype->getSequence() << " isREF " << m_pHaplotype->isReference() << std::endl;
			int start = *it2;
			  std::string pad;
			  std::string rseq = m_pRead->getSequence();
			  std::cerr << "rseq.size(): " << rseq.size() << std::endl;
			  if (start>0)
			  {
			     pad = std::string(start-1, ' ');
			  }
			  else 
			  {
			      rseq = m_pRead->getSequence().substr(-start, m_pRead->length());
			  }

          		std::cerr << pad << rseq << std::endl;
		}

	    	positions.insert(*it2);
		c++;
		if (c>=3) break;
	    }
    	    if (DINDEL_DEBUG) std::cerr << " done\n";

	    if (c>=3) break;
    }

    /*
     * doesn't work any more if reads don't come from a BAM file
    if (positions.empty())
    {
	// try the mapper coordinates
	if (!m_pRead->isUnmapped())
	{
    		int start = m_pRead->getBamStartPosAdjusted();
        	int end = m_pRead->getBAMEndAdjusted();
		
		positions.insert(start-m_pHaplotype->getRefStart());
		positions.insert(end-m_pRead->length()-m_pHaplotype->getRefStart());
	} 
    }
    */
}

ReadHaplotypeAlignment DindelHMM::getAlignment()
{
    bool rcRead = m_pRead->getRCRead();

    std::set<int> positions;

    getSeedPositions(positions);

    if (DINDEL_DEBUG && positions.empty()) std::cerr << "WARNING: HMM seed positions empty" << std::endl;
	
    // std::cerr << "seed positions: " << positions.size() << std::endl;

    ReadHaplotypeAlignment best(-1000.0, -1);
    double max_ll = -1000.0;

    for (std::set<int>::const_iterator iter = positions.begin(); iter != positions.end(); iter++)
    {
        int hFirstBase = *iter-DINDEL_HMM_BANDWIDTH/2;
        ReadHaplotypeAlignment rha = DindelHMMForward<DINDEL_HMM_BANDWIDTH>(m_pRead, m_pHaplotype, hFirstBase, rcRead);
        if (rha.logLik>max_ll)
        {
            best = rha;
            max_ll = best.logLik;
        }
    }

    return best;
}

/*
ReadHaplotypeAlignment DindelHMMForward(const DindelRead & read, const DindelHaplotype & haplotype, int hFirstBase, bool rcRead, int BandWidth)
{
}
*/


#ifdef DINDELHMMSTANDALONE
int main()
{
/*
	//                                     0         0         0         0         0
	DindelHaplotype hap = DindelHaplotype("ATTCCATCTACTTGGATATACGGATCAGCATGACTCCCTATCTTTATTC");
//	DindelRead read =          DindelRead("ATTCCATCTACTTGGATATACGGATCAGCATG");
	DindelRead read =          DindelRead("ATTCCATCTACTTGGATATCGGATCAGCATG");
*/


	



	/*
	DindelHaplotype hap = DindelHaplotype("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789");
	DindelRead read =          DindelRead("ABCDEFGHIJKLMMNOPQRST");

	for (int x=0;x<100000;x++)
		HMMForward(read, hap, -5, false); 
	*/

	DindelHaplotype hap = DindelHaplotype("GTCAGAGGCCAGACGCCACCGAGACTCCACCTCCCAGCATGCTAGCTCCAATTCCACCTCTACAGCAGCCTAGTCCTGAATCCACACCACAGCAGCCTAGCCCTGAATCCACACCACAGCA");
	DindelRead read = DindelRead("CGAGACTCCACCTCCCAGCATGCTAGCTCCAATTCCACCTCTCAGCAGCCTAGTCCTGAATCCACACCACAGCAGC");
	HMMForward(read, hap, -5, false);





}
#endif
/*
# l           2    3    4    5    6    7    8    9    10     NO_INS P( NO_INS(l)| NO_INS(l+1) )
#             ^          \
#             |           \
#             | normal     \ del 
# l+1         3    4    5    6    7    8    9   10    11     NO_INS
#
#
# l+2         4    5    6    7    8    9   10   11    12


# you can go from 

#             


# l           2    3    4    5    6    7    8    9    10        NO_INS
#                       |             
#                       |             
#                       |  
# l+1         3    4    5    6    7    8    9   10    11    INS
#                     /
#                   /                       X    X     X      P(NO_INS|NO_INS)=1 for x >= BandWidth-3
# l+2         4    5    6    7    8    9   10   11    12    NO_INS

# l           2    3    4    5    6    7    8    9    10        INS
#                                /               
#                               /              
#                             /                                                                                            
# l+1         3    4    5    6    7    8    9   10    11    INS                  and states as long as P(INS|NO_INS) is correct
#                     /
#                   /                       X    X     X      P(NO_INS|NO_INS)=1 for x >= BandWidth-3
# l+2         4    5    6    7    8    9   10   11    12    NO_INS

*/
