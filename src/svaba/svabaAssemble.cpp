//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// assemble - convert read overlaps into contigs
//

#include <iostream>
#include <fstream>
#include "Util.h"
#include "svabaAssemble.h"
#include "SGUtil.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"
#include "EncodedString.h"
#include <unistd.h>
#include <string>
#include "SGSearch.h"

//#define DEBUG_ASSEMBLY 1

void walkExtra(StringGraph * pGraph, SGWalkVector& outWalks);

StringGraph* assemble(std::stringstream& asqg_stream, int minOverlap, int maxEdges, bool bExact, 
	      int trimLengthThreshold, bool bPerformTR, bool bValidate, int numTrimRounds, 
              int resolveSmallRepeatLen, int numBubbleRounds, double maxBubbleGapDivergence, 
              double maxBubbleDivergence, int maxIndelLength, int cutoff, std::string prefix, 
		      SeqLib::UnalignedSequenceVector &contigs, bool get_components)
{

  AssemblyOptions ao;

  StringGraph * pGraph = SGUtil::loadASQG(asqg_stream, minOverlap, true, maxEdges);
  pGraph->m_get_components = get_components;

  if(bExact)
    pGraph->setExactMode(true);
  
  // Visitor functors
  SGTransitiveReductionVisitor trVisit;
  SGGraphStatsVisitor statsVisit;
  SGTrimVisitor trimVisit(trimLengthThreshold);
  SGContainRemoveVisitor containVisit;
  SGValidateStructureVisitor validationVisit;
  
  while(pGraph->hasContainment())
    pGraph->visit(containVisit);

  // Remove any extraneous transitive edges that may remain in the graph
  if(bPerformTR)
    {
      std::cout << "Removing transitive edges\n";
      pGraph->visit(trVisit);
    }
  
  // Compact together unbranched chains of vertices
  pGraph->simplify();
  
  if(bValidate)
    {
      pGraph->visit(validationVisit);
    }

  // Remove dead-end branches from the graph
  if(numTrimRounds > 0) {
      int numTrims = numTrimRounds;
      while(numTrims-- > 0)
	pGraph->visit(trimVisit);
    }
  
  // Resolve small repeats
  if(resolveSmallRepeatLen > 0 && false) {
      SGSmallRepeatResolveVisitor smallRepeatVisit(resolveSmallRepeatLen);
    }
  
  if(numBubbleRounds > 0)
    {
      SGSmoothingVisitor smoothingVisit(ao.outVariantsFile, maxBubbleGapDivergence, maxBubbleDivergence, maxIndelLength);
      int numSmooth = numBubbleRounds;
      while(numSmooth-- > 0)
	pGraph->visit(smoothingVisit);
      pGraph->simplify(); 
    }

  pGraph->renameVertices(prefix);

  SGVisitorContig av;
  pGraph->visit(av);
  
  SeqLib::UnalignedSequenceVector tmp = av.m_ct;
  for (SeqLib::UnalignedSequenceVector::const_iterator it = tmp.begin(); it != tmp.end(); it++) {
    if ((int)(it->Seq.length()) >= cutoff) {
      contigs.push_back({it->Name + "C", it->Seq, std::string()}); //postpend with character to distribugish _2 from _22
    }
  }
  
  /*  if (walk_all) {
      SGWalkVector outWalks;
      walkExtra(pGraph, outWalks);
      for (auto& i : outWalks) {
      std::string seqr = i.getString(SGWT_START_TO_END);
      if ((int)seqr.length() >= cutoff) 
      contigs.push_back(Contig(i.pathSignature(), seqr));
      }
      }*/
  

  return pGraph;
  //    delete pGraph;
   
}

void walkExtra(StringGraph * pGraph, SGWalkVector& outWalks) {
  
  typedef std::vector<VertexPtrVec> ComponentVector;
    VertexPtrVec allVertices = pGraph->getAllVertices();
    ComponentVector components;
#ifdef DEBUG_ASSEMBLY
    std::cerr << "selecting connected componsnet" << std::endl;
#endif
    
    SGSearchTree::connectedComponents(allVertices, components);
    
    // Select the largest component
    int selectedIdx = -1;
    size_t largestSize = 0;
    
    for(size_t i = 0; i < components.size(); ++i)
      {
	if(components[i].size() > largestSize)
	  {
	    selectedIdx = i;
	    largestSize = components[i].size();
	  }
      }
    
#ifdef DEBUG_ASSEMBLY
    std::cerr << "looping vertices" << std::endl;
#endif
    
    assert(selectedIdx != -1);    
    VertexPtrVec selectedComponent = components[selectedIdx];
    // Build a vector of the terminal vertices
    VertexPtrVec terminals;
    for(size_t i = 0; i < selectedComponent.size(); ++i)
      {
	Vertex* pVertex = selectedComponent[i];
	size_t asCount = pVertex->getEdges(ED_ANTISENSE).size();
	size_t sCount = pVertex->getEdges(ED_SENSE).size();
	
	if(asCount == 0 || sCount == 0)
	  terminals.push_back(pVertex);
      }
    
#ifdef DEBUG_ASSEMBLY
    std::cerr << "getting walks on " << terminals.size() << " vertices "  << std::endl;
#endif
    
    // Find walks between all-pairs of terminal vertices
    SGWalkVector tempWalks;
    if (terminals.size() > 2 && terminals.size() < 10) {
      for(size_t i = 0; i < terminals.size(); ++i)
	{
	  for(size_t j = i + 1; j < terminals.size(); j++)
	    {
	      Vertex* pX = terminals[i];
	      Vertex* pY = terminals[j];
	      int maxDistance = 2;
	      SGSearch::findWalks(pX, pY, ED_SENSE, maxDistance, 1000000, false, tempWalks);
	      SGSearch::findWalks(pX, pY, ED_ANTISENSE, maxDistance, 1000000, false, tempWalks);   
	    }
	}
    }
    
#ifdef DEBUG_ASSEMBLY
    std::cerr << "deduping walks on " << terminals.size() << " vertices "  << std::endl;
#endif
    
    // Remove duplicate walks
    std::map<std::string, SGWalk> walkMap;
    for(size_t i = 0; i < tempWalks.size(); ++i)
      {
	std::string walkString = tempWalks[i].getString(SGWT_START_TO_END);
	walkMap.insert(std::make_pair(walkString, tempWalks[i]));
      }
    // Copy unique walks to the output
    for(std::map<std::string, SGWalk>::iterator mapIter = walkMap.begin(); mapIter != walkMap.end(); ++mapIter)
      outWalks.push_back(mapIter->second);
    
    // Sort the walks by string length
    std::sort(outWalks.begin(), outWalks.end(), SGWalk::compareByTotalLength);
    
    int numOutputWalks = 10;
    if(numOutputWalks > 0 && numOutputWalks < (int)outWalks.size())
      {
	assert(outWalks.begin() + numOutputWalks < outWalks.end());
	outWalks.erase(outWalks.begin() + numOutputWalks, outWalks.end());
      }
  
}
