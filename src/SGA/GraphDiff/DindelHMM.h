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
#ifndef DINDELHMM_H
#define	DINDELHMM_H
#include "DindelRealignWindow.h"

const int DEBUGDINDELHMM=0;

typedef double Float;

class DindelHMM
{
    public:

        // Constructor
        DindelHMM(DindelRead & read, const DindelMultiHaplotype & haplotype);

        // Functions
        ReadHaplotypeAlignment getAlignment();
        
    private:

        //
        void getSeedPositions(std::set<int> & positions);

        // Data
        DindelRead* m_pRead;
        const DindelMultiHaplotype * m_pHaplotype;
};

template<int BandWidth> ReadHaplotypeAlignment DindelHMMForward(const DindelRead * pRead, 
                                                                const DindelMultiHaplotype * pHaplotype, 
                                                                int hFirstBase, 
                                                                bool rcRead)
{
	// code for aligning assuming last base of read is aligned to the haplotype
	// hFirstBase gives relative haplotype base of first base in read
	
	if (DEBUGDINDELHMM) 
        std::cerr << " _hmm bandwidth: " << BandWidth << " hFirstBase: " << hFirstBase << std::endl;
	
	double hp[] = { 2.9e-5, 2.9e-5,2.9e-5, 2.9e-5, 4.3e-5, 1.1e-4, 2.4e-4, 5.7e-4, 1.0e-3, 1.4e-3 };
	const int MAXHP = 10;

	//const int BandWidth=32;
	const Float GAP_EXT = 0.5;
	const Float GAP_STOP = 1.0 - GAP_EXT;

	int rlen = pRead->length();
	int hlen = pHaplotype->length();

	const std::string & hapSeq = pHaplotype->getSequence();

	Float *curr = new Float[BandWidth*2];
	Float *next = new Float[BandWidth*2];
	Float *obs = new Float[BandWidth];
	Float *gap_prob = new Float[BandWidth];

	// initialize curr to ones
	for (int x=0;x<BandWidth*2;x++) 
	{
		curr[x]=1.0; // NOTE NOT in log domain
		next[x]=0.0;
	}
	if (DEBUGDINDELHMM) std::cerr << "RLEN: " << rlen << std::endl;

	double norm=0.0;
	for (int l=0;l<rlen;l++)
	{
		//if (DEBUGDINDELHMM) std::cerr << "\n\n**************************************************************" << std::endl;
		int readBaseIndex = (!rcRead)?l:(rlen-1-l);
		char rb = (!rcRead)?pRead->getBase(readBaseIndex):complement(pRead->getBase(readBaseIndex));

		// probabilities of correctly and incorrectly observing the read base

		Float pr=1.0 - exp ( (-2.3026/10.0)*double(pRead->getQual(readBaseIndex) ));
        	Float p_base_correct = (.25+.75*pr);
        	Float p_base_incorrect =(.75+1e-10-.75*pr);

		// lower and upper haplotype position for this read base

		int sb = 0;		// start in band
		int eb = BandWidth-1;   // end in band (inclusive)
		int lh = hFirstBase+l; // start base on haplotype
		if (lh<0)
		{
			sb=-lh;
			lh=0;
		}

		if (sb>=BandWidth && l<rlen-1) 
		{
			// this is left of the pHaplotype->
			// We don't need to explicitly do a forward pass, as the base is assumed to match by default.
			norm += log(p_base_correct);
			continue;
			//sb = BandWidth;
		}

		int uh = lh+BandWidth;
		if (uh>=hlen)
		{
			eb=BandWidth-(1+uh-hlen);
			if (eb<-1) eb=-1;
			uh=hlen-1;              // both are inclusive
		}
		if (DEBUGDINDELHMM && 1)
		{
			std::cerr << "l: " << l << " lh: " << lh << " uh: " << uh << " sb: " << sb << " eb: " << eb << " "; // << std::endl;
		}

		if (sb>=BandWidth) sb=BandWidth;
		// set observations and gap_prob
		for (int x=0;x<sb;++x)
		{
			obs[x] = p_base_correct;
			gap_prob[x] = hp[0]/2.0; // divide by two to get the insertion/deletion prob

		}
		for (int x=eb+1;x<BandWidth;++x)
		{
			obs[x] = p_base_correct;
			gap_prob[x] = hp[0]/2.0;
		}

		int h=lh;
		for (int x=sb;x<=eb;++x,++h)
		{
			obs[x] = (hapSeq[h]==rb)?p_base_correct:p_base_incorrect;
			int hplen = pHaplotype->getHomopolymerLength(h);
			Float prob;
			if (hplen<MAXHP) prob=hp[hplen]; else
			{
				prob=hp[9]+4.3e-4*double(hplen-10);
				if (prob>0.95) prob=0.95;
			}
			gap_prob[x]=prob/2.0;
		}

		if (DEBUGDINDELHMM && 1)
		{
			//std::cerr << "l: " << l << std::endl;
			for (int x=0;x<BandWidth;x++) std::cerr << "(" << obs[x] << " gp " << log(gap_prob[x]) << ")" << std::endl;
			//std::cerr << "norm: " << norm << std::endl;
		}

		// UPDATES
		if (l<rlen-1)
		{
			// do forward pass

		    	// P( INS(l) | INS(l+1) ) P (l+1)

		    	// INSERTION <= INSERTION
			for (int x=1;x<BandWidth;++x) next[BandWidth+x-1] += GAP_EXT*curr[BandWidth+x]*p_base_correct;

			// INSERTION <= NO_INSERTION
			for (int x=1;x<BandWidth;++x) next[x-1] += gap_prob[x]*curr[x+BandWidth]*p_base_correct;

			// PREMULTIPLY curr and obs
			for (int x=0;x<BandWidth;++x) curr[x]*=obs[x];

			// NO_INSERTION <= INSERTION
			for (int x=0;x<BandWidth;++x) next[BandWidth+x] += GAP_STOP*curr[x];

			// NO_INSERTION <= NO_INSERTION
			for (int x=0;x<BandWidth;++x) next[x] += (1.0-2.0*gap_prob[x])*curr[x]; // this is the no-indel transition
			for (int x=0;x<BandWidth-1;++x) next[x+1] += 0.632333*gap_prob[x]*curr[x]; // deletion of 1
    			for (int x=0;x<BandWidth-2;++x) next[x+2] += 0.232622*gap_prob[x]*curr[x];
    			for (int x=0;x<BandWidth-3;++x) next[x+3] += 0.085577*gap_prob[x]*curr[x];
    			for (int x=0;x<BandWidth-4;++x) next[x+4] += 0.031482*gap_prob[x]*curr[x];
			
    			for (int x=0;x<BandWidth-5;++x) next[x+5] += 0.011582*gap_prob[x]*curr[x];
    			for (int x=0;x<BandWidth-6;++x) next[x+6] += 0.004261*gap_prob[x]*curr[x];
    			for (int x=0;x<BandWidth-7;++x) next[x+7] += 0.001567*gap_prob[x]*curr[x];
    			for (int x=0;x<BandWidth-8;++x) next[x+8] += 0.000577*gap_prob[x]*curr[x];
			
		}
		else
		{
			if (DEBUGDINDELHMM) std::cerr << "normalizing" << std::endl;
			// add prior and observations for last locus
			for (int x=0;x<BandWidth;x++) next[x] = curr[x]*obs[x]*(1.0-2.0*gap_prob[x])/Float(BandWidth);
			for (int x=0;x<BandWidth;x++) next[x+BandWidth] = curr[x+BandWidth]*p_base_correct*gap_prob[x]/Float(BandWidth);
		}

		// normalize
		Float sum=0.0;
		for (int x=0;x<BandWidth*2;x++) sum += next[x];
		norm += log(sum);


		for (int x=0;x<BandWidth*2;x++)
		{
			next[x]/=sum;
			if (next[x]<1e-10) next[x]=1e-10; // introduces small rounding error....
			curr[x]=0.0; // note curr will be swapped with next below
		}

		if (DEBUGDINDELHMM)
		{
			//std::cerr << "xyz l: " << l << std::endl;
			if (1) for (int x=0;x<BandWidth;x++) std::cerr << "x: " << x << "\t" << next[x] << " " << next[x+BandWidth] << std::endl;
			std::cerr << "[ " << p_base_correct << " l: " << l << " " << " norm: " << norm << " sum: " << sum << "]" << std::endl;
		}


		// END UPDATES

		// flip next and curr
		Float *tmp = next;
		next = curr;
		curr = tmp;
	} // end forward passes

	// norm should give the log-likelihood


	// check if there is a state that has posterior >0.95
	Float postProb = -1.0;
	int state = -1;
	for (int x=0;x<BandWidth;x++) 
	{
		if (curr[x]>postProb)
		{
			postProb=curr[x];
			state=x;
		}
	}
	int hapPosLastReadBase = (state!=-1)?(hFirstBase+rlen-1+state):-1;

	// cleanup
	delete[] curr;
	delete[] next;
	delete[] obs;
	delete[] gap_prob;

	//std::cerr << "NORM POS " << norm << " hlen: " << hlen << " check: " << (hapPosLastReadBase-rlen+1>=hlen) << " check2: " << (hFirstBase+BandWidth>=hlen) << std::endl;
	assert(norm<0.0);

	ReadHaplotypeAlignment rha;
	//rha.logLik = addLogs(norm, -pRead->getMappingQual()*.23026-log(double(BandWidth))); // cap based on mapping quality
	rha.logLik = norm;
	rha.postProbLastReadBase = postProb;
	rha.hapPosLastReadBase=hapPosLastReadBase;
	return rha;
}

ReadHaplotypeAlignment DindelHMMForward(const DindelRead & read, 
                                        const DindelHaplotype & haplotype, 
                                        int hFirstBase, 
                                        bool rcRead, 
                                        int BandWidth);

#endif	/* DINDELHMM_H */

