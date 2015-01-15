///----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// PairedDeBruijnHaplotypeBuilder - Build haplotypes
// using a de Bruijn graph and read pair constraints
//
#include "PairedDeBruijnHaplotypeBuilder.h"
#include "BWTAlgorithms.h"
#include "SGSearch.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"
#include "Profiler.h"
#include "HapgenUtil.h"
#include "multiple_alignment.h"
#include "KmerOverlaps.h"
#include "Verbosity.h"

//
//
//
PairedDeBruijnHaplotypeBuilder::PairedDeBruijnHaplotypeBuilder(const GraphCompareParameters& params) : m_parameters(params)
{

}

//
PairedDeBruijnHaplotypeBuilder::~PairedDeBruijnHaplotypeBuilder()
{

}

//
void PairedDeBruijnHaplotypeBuilder::setInitialHaplotype(const std::string& str)
{
    m_startingKmer = str;
}

// Run the bubble construction process
HaplotypeBuilderReturnCode PairedDeBruijnHaplotypeBuilder::run(StringVector& out_haplotypes)
{
    PROFILE_FUNC("GraphCompare::buildVariantStringGraph")
    assert(!m_startingKmer.empty());
    
    // We search until we find the first common vertex in each direction
    //size_t MIN_TARGET_COUNT = m_parameters.bReferenceMode ? 1 : 2;
    size_t MAX_ITERATIONS = 5000;
    size_t MAX_SIMULTANEOUS_BRANCHES = 100;
    size_t MAX_TOTAL_BRANCHES = 500;
    int MIN_TARGET_DISTANCE = 50;

    // Tracking stats
    size_t max_simul_branches_used = 0;
    size_t total_branches = 0;
    size_t iterations = 0;

    // To simplify the search we use the full read sequences
    // to guide the de Bruijn graph assembly. In addition, we
    // require the search to end on a kmer of a read pair of
    // one of the reads containing the initial kmer. 
    // Set up these data structures now
    DBGPathGuide guide(m_startingKmer.size());
    std::set<std::string> target_set;
    selectGuideAndTargetKmers(m_startingKmer, false, guide, target_set);
    selectGuideAndTargetKmers(reverseComplement(m_startingKmer), true, guide, target_set);

    if(Verbosity::Instance().getPrintLevel() > 3)
        printf("PairedDBGHaplotype: found %zu targets\n", target_set.size());

    if(target_set.empty())
        return HBRC_OK;

    // Initialize the graph
    StringGraph* pGraph = new StringGraph;
    BuilderExtensionQueue queue;

    Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(m_startingKmer, m_startingKmer);
    pVertex->setColor(GC_BLACK);
    pGraph->addVertex(pVertex);

    // Add the vertex to the extension queue
    queue.push(BuilderExtensionNode(pVertex, ED_SENSE, 0));
    queue.push(BuilderExtensionNode(pVertex, ED_ANTISENSE, 0));

    std::vector<Vertex*> sense_join_vector;
    std::vector<Vertex*> antisense_join_vector;

    bool extension_allowed[ED_COUNT] = { true, true };

    // Perform the extension. The while conditions are heuristics to avoid searching
    // the graph too much 
    while(!queue.empty() && iterations++ < MAX_ITERATIONS && queue.size() < MAX_SIMULTANEOUS_BRANCHES && total_branches < MAX_TOTAL_BRANCHES)
    {
        if(queue.size() > max_simul_branches_used)
            max_simul_branches_used = queue.size();

        BuilderExtensionNode curr = queue.front();
        queue.pop();

        if(!extension_allowed[curr.direction])
            continue;

        // Calculate de Bruijn extensions for this node
        std::string vertStr = curr.pVertex->getSeq().toString();
        AlphaCount64 extensionCounts = BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(vertStr, m_parameters.variantIndex.pBWT, curr.direction);

        // Check whether to accept this edge into the graph
        // We currently only use the counts and not the guide kmers
        bool count_passed[4] = { false, false, false, false };

        for(size_t i = 0; i < DNA_ALPHABET::size; ++i)
        {
            char b = DNA_ALPHABET::getBase(i);
            size_t count = extensionCounts.get(b);
            count_passed[i] = count >= m_parameters.minDBGCount;
        }

        // If there is any guide p-mer found, only use the guide pmers to extend the graph. Otherwise
        // use all the pmers that passed the count threshold.
        std::string extensions;
        for(size_t i = 0; i < DNA_ALPHABET::size; ++i)
        {
            char b = DNA_ALPHABET::getBase(i);
            if( count_passed[i] ) 
                extensions.push_back(b);
        }

        std::string extensionsUsed;
        for(size_t i = 0; i < extensions.size(); ++i)
        {
            char b = extensions[i];
            extensionsUsed.push_back(b);
            std::string newStr = VariationBuilderCommon::makeDeBruijnVertex(vertStr, b, curr.direction);

            // Create the new vertex and edge in the graph
            // Skip if the vertex already exists
            if(pGraph->getVertex(newStr) != NULL)
                continue;
            
            // Allocate the new vertex and add it to the graph
            Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(newStr, newStr);
            pVertex->setColor(GC_BLACK);
            pGraph->addVertex(pVertex);

            // Add edges
            VariationBuilderCommon::addSameStrandDeBruijnEdges(pGraph, curr.pVertex, pVertex, curr.direction);
            size_t ref_count = BWTAlgorithms::countSequenceOccurrences(newStr, m_parameters.referenceIndex);
            
            //if(target_set.find(newStr) != target_set.end() && curr.distance >= MIN_TARGET_DISTANCE)
            if((ref_count == 1 || target_set.find(newStr) != target_set.end()) && curr.distance >= MIN_TARGET_DISTANCE)
            {
                if(curr.direction == ED_SENSE)
                    sense_join_vector.push_back(pVertex);
                else
                    antisense_join_vector.push_back(pVertex);
                extension_allowed[curr.direction] = false;
            }
            else
            {
                // Add the vertex to the extension queue
                queue.push(BuilderExtensionNode(pVertex, curr.direction, curr.distance + 1));
            }
        }
        
        // Update the total number of times we branches the search
        if(!extensionsUsed.empty())
            total_branches += extensionsUsed.size() - 1;
    }
    
    if(Verbosity::Instance().getPrintLevel() > 2)
    {
        printf("[PairedDBG] iterations: %zu total_branches: %zu left found? %d right found? %d\n", 
                 iterations, total_branches, !extension_allowed[ED_ANTISENSE], !extension_allowed[ED_SENSE]);
    }

    // If the graph construction was successful, walk the graph
    // between the endpoints to make a string
    // Generate haplotypes between every pair of antisense/sense join vertices
    for(size_t i = 0; i < antisense_join_vector.size(); ++i) {
        for(size_t j = 0; j < sense_join_vector.size(); ++j) {
            SGWalkVector outWalks;
            SGSearch::findWalks(antisense_join_vector[i],
                                sense_join_vector[j],
                                ED_SENSE,
                                100000, // max distance to search
                                10000, // max nodes to search
                                true, // exhaustive search
                                outWalks);

            for(size_t k = 0; k < outWalks.size(); ++k)
            {

                out_haplotypes.push_back(outWalks[k].getString(SGWT_START_TO_END));
            }
        }
    }

    delete pGraph;
    return HBRC_OK;
}
 
//
void PairedDeBruijnHaplotypeBuilder::selectGuideAndTargetKmers(const std::string& kmer, 
                                                               bool is_kmer_reversed, 
                                                               DBGPathGuide& guide,
                                                               std::set<std::string>& target_set) const
{
    size_t k = kmer.size();
    std::vector<size_t> ids = getReadIDs(kmer);

    // Set up guide
    for(size_t i = 0; i < ids.size(); ++i)
    {
        std::string read = BWTAlgorithms::extractString(m_parameters.variantIndex.pBWT, ids[i]);
        if(is_kmer_reversed)
            read = reverseComplement(read);
        guide.addSequence(read);
    }

    // Set up targets
    for(size_t i = 0; i < ids.size(); i++)
    {
        size_t pair_id = ids[i] % 2 == 0 ? ids[i] + 1 : ids[i] - 1;
        std::string pair = BWTAlgorithms::extractString(m_parameters.variantIndex.pBWT, pair_id);

        for(size_t j = 0; j < pair.size() - k + 1; ++j)
        {
            std::string test_kmer = pair.substr(j, k);
            size_t ref_count = BWTAlgorithms::countSequenceOccurrences(test_kmer, m_parameters.referenceIndex);
            size_t base_count = BWTAlgorithms::countSequenceOccurrences(test_kmer, m_parameters.baseIndex);
            if(ref_count == 1 && base_count > 0)
                target_set.insert(is_kmer_reversed ? test_kmer : reverseComplement(test_kmer));
        }
    }
}

//       
std::vector<size_t> PairedDeBruijnHaplotypeBuilder::getReadIDs(const std::string& kmer) const
{
    BWTInterval interval = BWTAlgorithms::findInterval(m_parameters.variantIndex, kmer);
    std::vector<size_t> out;
    for(int64_t i = interval.lower; i <= interval.upper; ++i)
        out.push_back(m_parameters.variantIndex.pSSA->calcSA(i, m_parameters.variantIndex.pBWT).getID());
    return out;
}
