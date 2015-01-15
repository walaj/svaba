///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// StringHaplotypeBuilder - Build read coherent haplotypes
// using read overlaps.
//
#include "StringHaplotypeBuilder.h"
#include "BWTAlgorithms.h"
#include "SGSearch.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"
#include "Profiler.h"
#include "HapgenUtil.h"
#include "multiple_alignment.h"
#include "KmerOverlaps.h"
#include "CorrectionThresholds.h"

//#define SHOW_MULTIPLE_ALIGNMENT 1
//#define SHOW_GRAPH 1
#define OVERLAP_HAP_DEBUG 1

//
//
//
StringHaplotypeBuilder::StringHaplotypeBuilder(const GraphCompareParameters& params) : m_parameters(params), 
                                                                                       m_numReads(0), 
                                                                                       m_overlapLengthFrac(0.75)
{
    assert(m_parameters.minOverlap > 25);

    ErrorCorrectParameters correction_params;
    correction_params.pOverlapper = NULL;
    correction_params.indices = m_parameters.variantIndex;
    correction_params.algorithm = ECA_KMER;

    // Overlap-based corrector params
    correction_params.minOverlap = 51;
    correction_params.numOverlapRounds = 1;
    correction_params.minIdentity = 0.95f;
    correction_params.conflictCutoff = 5;
    correction_params.depthFilter = 100;
    correction_params.printOverlaps = false;

    // k-mer based corrector params
    correction_params.numKmerRounds = 10;
    correction_params.kmerLength = 31;
    CorrectionThresholds::Instance().setBaseMinSupport(3);

    m_graph = new StringGraph;
    m_corrector = new ErrorCorrectProcess(correction_params);
    m_extractor = new OverlapExtractorWithCorrection(31, m_parameters.variantIndex, m_corrector, m_parameters.minOverlap, 0.95);
}

//
StringHaplotypeBuilder::~StringHaplotypeBuilder()
{
    delete m_corrector;
    delete m_graph;
    delete m_extractor;
}

// The source string is the string the bubble starts from
void StringHaplotypeBuilder::setInitialHaplotype(const std::string& sequence)
{
#ifdef OVERLAP_HAP_DEBUG
    printf("\n\n***** Starting new haplotype %s\n", sequence.c_str());
#endif
    m_initial_kmer_string = sequence;
}

// Run the bubble construction process
HaplotypeBuilderReturnCode StringHaplotypeBuilder::run(StringVector& out_haplotypes)
{
    PROFILE_FUNC("StringHaplotypeBuilder::run")
    assert(m_parameters.variantIndex.pBWT != NULL);
    assert(!m_initial_kmer_string.empty());

    // Start the haplotype generation process by finding reads containing the initial kmer
    StringVector query_kmers(1, m_initial_kmer_string);
    StringVector reads;
    getReadsForKmers(query_kmers, m_parameters.kmer, &reads);

#ifdef OVERLAP_HAP_DEBUG
    printf("Found %zu initial reads\n", reads.size());
#endif

    // Correct sequencing errors in the reads
    correctReads(&reads);
    if(reads.empty())
        return HBRC_OK;

#ifdef OVERLAP_HAP_DEBUG
    printf("Corrected %zu initial reads\n", reads.size());
#endif

    // Select the seed read for extension
    size_t seed_idx = 0;
    while(seed_idx < reads.size() && reads[seed_idx].empty())
        seed_idx++;

    if(seed_idx >= reads.size())
        return HBRC_OK;

    std::string seed = reads[seed_idx];
    Vertex* seed_vertex = new(m_graph->getVertexAllocator()) Vertex(seed, seed);
    m_graph->addVertex(seed_vertex);

    Vertex* right_endpoint = NULL;
    Vertex* left_endpoint = NULL;

    right_endpoint = extendSeed(seed_vertex, ED_SENSE);
    if(right_endpoint != NULL)
        left_endpoint = extendSeed(seed_vertex, ED_ANTISENSE);

    if(left_endpoint != NULL && right_endpoint != NULL)
    {
        SGWalkVector walks;
        SGSearch::findWalks(left_endpoint, right_endpoint, ED_SENSE, 2000, 10000, true, walks);
        for(size_t i = 0; i < walks.size(); ++i)
        {
            out_haplotypes.push_back(walks[i].getString(SGWT_START_TO_END));
            std::cout << ">h" << i << "\n" << out_haplotypes.back() << "\n";
        }
    }

    /*
    std::stringstream graph_name;
    graph_name << "graph.r";
    m_graph->writeDot(graph_name.str() + ".dot");
    m_graph->writeASQG(graph_name.str() + ".asqg");
    exit(1);
    */

    return HBRC_OK;
}

Vertex* StringHaplotypeBuilder::extendSeed(Vertex* seed_vertex, EdgeDir direction)
{
    BuilderExtensionNode node(seed_vertex, direction, 0);
    BuilderExtensionQueue queue;
    queue.push(node);

    size_t limit = 100;
    while(!queue.empty() && limit-- > 0)
    {
        printf("CoreExtension. Queue size: %zu\n", queue.size());

        // Check if we are extending a unipath
        if(queue.size() == 1)
        {
            if(queue.front().distance > 10 || (queue.front().distance > 5 && isJoinSequence(queue.front().pVertex->getStr(), queue.front().direction)))
            {
                std::cout << "Single path, ending\n";
                return queue.front().pVertex;
            }
        }

        BuilderExtensionNode current_node = queue.front();
        queue.pop();

        std::cout << "\textending d: " << current_node.distance << " " << current_node.pVertex->getStr() << "\n";
        SequenceOverlapPairVector extensions = computeExtensions(current_node.pVertex->getStr(), current_node.direction);
        if(extensions.size() == 0)
        {
            std::cout << "  No extensions\n";
            trimTip(current_node.pVertex, current_node.direction);
            continue;
        }
        
        for(size_t i = 0; i < extensions.size(); ++i)
        {
            std::string incoming_sequence = extensions[i].sequence[1];
            Vertex* incoming_vertex = m_graph->getVertex(incoming_sequence);
            if(incoming_vertex == NULL)
            {
                // we need to add a new vertex in the graph for this sequence
                incoming_vertex = new(m_graph->getVertexAllocator()) Vertex(incoming_sequence, incoming_sequence);
                m_graph->addVertex(incoming_vertex);
                BuilderExtensionNode incoming_node(incoming_vertex, current_node.direction, current_node.distance + 1);
                queue.push(incoming_node);
            }

            // add edges to the graph
            addEdge(current_node.pVertex, incoming_vertex, extensions[i].overlap);
        }
    }

    return NULL;
}

//
SequenceOverlapPairVector StringHaplotypeBuilder::computeExtensions(std::string sequence, EdgeDir direction)
{
    SequenceOverlapPairVector overlaps = getCorrectedOverlaps(sequence, direction);
    std::sort(overlaps.begin(), overlaps.end(), SequenceOverlapPair::sortByOverlapLengthDesc);
    if(overlaps.empty())
        return overlaps;

    // Remove transitive overlaps
    size_t longest_idx = 0;
    double longest_length = overlaps[longest_idx].overlap.getOverlapLength();
    while(1) 
    {
        // Get the label of this irreducible edge
        std::string irr_label = getEdgeLabel(overlaps[longest_idx].sequence[1], overlaps[longest_idx].overlap);
        std::cout << "Irreducible: " << overlaps[longest_idx].sequence[1] << "\n";
        std::cout << "Extension: " << irr_label << "\n";

        // Check whether the extension of the shorter overlaps have this sequence as an extension prefix
        for(size_t i = longest_idx + 1; i < overlaps.size(); ++i)
        {
            // Skip eliminated sequences
            if(overlaps[i].sequence[1].empty())
                continue;
            
            std::string label = getEdgeLabel(overlaps[i].sequence[1], overlaps[i].overlap);
            std::cout << "TransExtension?: " << label << "\n";
            
            bool is_transitive = false;
            if(direction == ED_SENSE)
                is_transitive = label.substr(0, irr_label.size()) == irr_label;
            else
                is_transitive = label.substr(label.size() - irr_label.size()) == irr_label;

            std::cout << "  isTrans?: " << is_transitive << "\n";

            // Eliminate this sequence
            if(is_transitive)
                overlaps[i].sequence[1] = "";
        }

        // Done elimination of transitive edges through this string, move to the next longest non-eliminated sequence
        while(++longest_idx < overlaps.size())
        {
            if(!overlaps[longest_idx].sequence[1].empty())
                break;
        }

        // Check if all the reads have been classified 
        if(longest_idx == overlaps.size())
            break;
    }

    SequenceOverlapPairVector out_pairs;
    std::cout << "Irreducible overlaps for read: " << sequence << "\n";
    for(size_t i = 0; i < overlaps.size(); ++i)
    {
        double ol = overlaps[i].overlap.getOverlapLength();
        if(!overlaps[i].sequence[1].empty() && ol / longest_length >= m_overlapLengthFrac)
        {
            overlaps[i].overlap.printAlignment(sequence, overlaps[i].sequence[1]);
            out_pairs.push_back(overlaps[i]);
        }
    }

    return out_pairs;
}
        
//
std::string StringHaplotypeBuilder::getEdgeLabel(const std::string& s2, const SequenceOverlap& overlap)
{
    // Get the unmatched part of s2
    int ms = overlap.match[1].start;
    int me = overlap.match[1].end;

    if(ms == 0)
        return s2.substr(me + 1);
    else
        return s2.substr(0, ms);
}

//
void StringHaplotypeBuilder::getReadsForKmers(const StringVector& kmer_vector, size_t k, StringVector* reads)
{
    PROFILE_FUNC("StringHaplotypeBuilder::getReadsForKmers")
    SeqRecordVector fwd_si;
    SeqRecordVector rev_si;

    // Forward reads
    HapgenUtil::extractHaplotypeReads(kmer_vector, m_parameters.variantIndex, k, false, m_parameters.maxReads, m_parameters.maxExtractionIntervalSize, &fwd_si, NULL);

    // Reverse reads
    HapgenUtil::extractHaplotypeReads(kmer_vector, m_parameters.variantIndex, k, true, m_parameters.maxReads, m_parameters.maxExtractionIntervalSize, &rev_si, NULL);

    // Copy reads into the positioned read vector, initially with unset positions
    for(size_t i = 0; i < fwd_si.size(); ++i)
    {
        std::string read_sequence = fwd_si[i].seq.toString();
        reads->push_back(read_sequence);
    }

    // reversed reads
    for(size_t i = 0; i < rev_si.size(); ++i) 
    {
        std::string read_sequence = reverseComplement(rev_si[i].seq.toString());
        reads->push_back(read_sequence);
    }
}
SequenceOverlapPairVector StringHaplotypeBuilder::getCorrectedOverlaps(const std::string& sequence, EdgeDir direction)
{
    SequenceOverlapPairVector initial_overlaps = m_extractor->getExactOverlaps(sequence);

    // Filter out the overlaps that aren't in the requested direction
    SequenceOverlapPairVector directional_overlaps;
    for(size_t i = 0; i < initial_overlaps.size(); ++i)
    {
        // Recompute the overlap, using the previous overlap as a guide
        const SequenceOverlap& overlap = initial_overlaps[i].overlap;
        bool is_prefix_overlap = overlap.match[0].start == 0;
        bool is_suffix_overlap = overlap.match[0].end == ((int)sequence.size() - 1);

        bool valid_direction = (is_prefix_overlap && direction == ED_ANTISENSE) ||
                               (is_suffix_overlap && direction == ED_SENSE);

        if(valid_direction)
            directional_overlaps.push_back(initial_overlaps[i]);
    }
    return directional_overlaps;
}

//
SequenceOverlapPairVector StringHaplotypeBuilder::getCorrectedOverlaps2(const std::string& sequence, EdgeDir direction)
{
    PROFILE_FUNC("StringHaplotypeBuilder::getCorrectedOverlaps")
    
    // Extract reads that share a short kmer with the input sequence
    
    size_t k = 31;
    SequenceOverlapPairVector overlap_vector = KmerOverlaps::retrieveMatches(sequence, k, m_parameters.minOverlap,
                                                                             0.95, 2, m_parameters.variantIndex);
    
    //SequenceOverlapPairVector overlap_vector = m_extractor->queryOverlaps(sequence);
    
    // Copy out the read sequences so they can be corrected
    // We use the original sequencing strand of the read here - if
    // it was reversed when we extracted it, we flip it back
    StringVector reads;
    for(size_t i = 0; i < overlap_vector.size(); ++i)
    {
        if(!overlap_vector[i].is_reversed)
            reads.push_back(overlap_vector[i].sequence[1]);
        else
            reads.push_back(reverseComplement(overlap_vector[i].sequence[1]));
    }

#ifdef OVERLAP_HAP_DEBUG
    printf("Extracted %zu reads\n", overlap_vector.size());
#endif
    
    //m_extractor->queryOverlaps(sequence);

    // Correct the reads
    correctReads(&reads);
    assert(reads.size() == overlap_vector.size());

#ifdef SHOW_MULTIPLE_ALIGNMENT
    MultipleAlignment corrected_ma;
    corrected_ma.addBaseSequence("query", sequence, "");
#endif

    // Overlap each corrected read with the input sequence
    // Discard any read that does perfectly match with a long enough overlap
    SequenceOverlapPairVector out;
    
    for(size_t i = 0; i < reads.size(); ++i) 
    {
        // Skip uncorrected reads and those identical to the query
        if(reads[i].empty())
            continue;

        // Change the strand of the corrected read back to the same as the query sequence
        std::string incoming = reads[i];
        if(overlap_vector[i].is_reversed)
            incoming = reverseComplement(incoming);
        
        // Skip identical matches
        if(incoming == sequence)
            continue;

        // Recompute the overlap, using the previous overlap as a guide
        const SequenceOverlap& initial_overlap = overlap_vector[i].overlap;
        SequenceOverlap overlap = Overlapper::extendMatch(sequence, incoming, initial_overlap.match[0].start, initial_overlap.match[1].start, 2);

        bool is_prefix_overlap = overlap.match[0].start == 0;
        bool is_suffix_overlap = overlap.match[0].end == ((int)sequence.size() - 1);

        bool valid_direction = (is_prefix_overlap && direction == ED_ANTISENSE) ||
                               (is_suffix_overlap && direction == ED_SENSE);

#ifdef SHOW_MULTIPLE_ALIGNMENT
        corrected_ma.addOverlap("crct_in", incoming, "", overlap);
#endif

        if(overlap.edit_distance == 0 && overlap.getOverlapLength() >= m_parameters.minOverlap && valid_direction)
        {
            SequenceOverlapPair sop;
            sop.sequence[1] = incoming;
            sop.overlap = overlap;
            out.push_back(sop);
        }
    }

#ifdef SHOW_MULTIPLE_ALIGNMENT
        corrected_ma.print(200);
#endif
    return out;
}

//
void StringHaplotypeBuilder::correctReads(StringVector* reads)
{
    PROFILE_FUNC("StringHaplotypeBuilder::correctReads")
    StringVector out_reads;
    for(size_t i = 0; i < reads->size(); ++i)
    {
        // Perform correction
        std::string corrected_sequence;
        DNAString dna(reads->at(i));
        SeqRecord record = { "null", dna, "" };
        //SeqRecord record("null", dna.toString(), ""); //JEREMIAH
        SequenceWorkItem wi(0, record);
        ErrorCorrectResult r = m_corrector->correct(wi);
        if(r.kmerQC || r.overlapQC)
            corrected_sequence = r.correctSequence.toString();
        // We allow the sequence to be null if the read could not be corrected
        out_reads.push_back(corrected_sequence);
    }
    reads->swap(out_reads);
}

//
void StringHaplotypeBuilder::trimTip(Vertex* x, EdgeDir direction)
{
    EdgePtrVec x_edges = x->getEdges(direction);
    if(x_edges.size() > 0)
        return; // not a tip

    // Check if we should recurse to the neighbors of x
    EdgePtrVec x_opp_edges = x->getEdges(!direction);
    
    Vertex* x_neighbor = NULL;
    if(x_opp_edges.size() == 1)
        x_neighbor = x_opp_edges.front()->getEnd();

    // Remove x from the graph
    m_graph->removeConnectedVertex(x);

    // Recurse to x's neighbors if it is now a tip
    if(x_neighbor != NULL && x_neighbor != x)
        trimTip(x_neighbor, direction);
}

//
bool StringHaplotypeBuilder::isJoinSequence(const std::string& sequence, EdgeDir dir)
{
    // Check if this sequence forms a complete k-mer path in the base/reference sequence
    size_t k = m_parameters.kmer;
    if(sequence.size() < k)
        return false;

    // Check if the start/end kmer of the sequence is present in the reference
    // depending on the extension direction
    std::string kmer_seq;
    if(dir == ED_ANTISENSE)
        kmer_seq = sequence.substr(0, k);
    else
        kmer_seq = sequence.substr(sequence.size() - k);

    // Count occurrences of the terminal kmer in the reference
    size_t base_count = BWTAlgorithms::countSequenceOccurrences(kmer_seq, m_parameters.baseIndex); 
    return base_count > 0;
}

// 
void StringHaplotypeBuilder::addEdge(Vertex* vertex1, Vertex* vertex2, SequenceOverlap overlap)
{
    Overlap sga_overlap(vertex1->getID(), overlap.match[0].start, overlap.match[0].end, vertex1->getSeqLen(),
                        vertex2->getID(), overlap.match[1].start, overlap.match[1].end, vertex2->getSeqLen(), false, 0);
    SGAlgorithms::createEdgesFromOverlap(m_graph, sga_overlap, true);
}
