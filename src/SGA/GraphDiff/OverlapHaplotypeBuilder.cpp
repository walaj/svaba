///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapHaplotypeBuilder - Build read coherent haplotypes
// using read overlaps.
//
#include "OverlapHaplotypeBuilder.h"
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
//#define OVERLAP_HAP_DEBUG 1

//
//
//
OverlapHaplotypeBuilder::OverlapHaplotypeBuilder(const GraphCompareParameters& params) : m_parameters(params), m_numReads(0)
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
}

//
OverlapHaplotypeBuilder::~OverlapHaplotypeBuilder()
{
    delete m_corrector;
    delete m_graph;
}

// The source string is the string the bubble starts from
void OverlapHaplotypeBuilder::setInitialHaplotype(const std::string& sequence)
{
#ifdef OVERLAP_HAP_DEBUG
    printf("\n\n***** Starting new haplotype %s\n", sequence.c_str());
#endif
    m_initial_kmer_string = sequence;
}

// Run the bubble construction process
HaplotypeBuilderReturnCode OverlapHaplotypeBuilder::run(StringVector& out_haplotypes)
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::run")
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

#ifdef OVERLAP_HAP_DEBUG
    printf("Corrected %zu initial reads\n", reads.size());
#endif

    // Start the graph using the corrected reads that contain the initial kmer
    if(!buildInitialGraph(reads))
        return HBRC_OK;
    
    bool done = false;

    int MAX_ROUNDS = 10;
    size_t MAX_TIPS = 20;
    size_t MAX_GRAPH_SIZE = 500;
    int round = 0;

    while(!done) 
    {
        // Extend the graph by finding new overlaps for the reads at the tips
        extendGraph();

        // Clean up the graph
        SGIdenticalRemoveVisitor dupVisit;
        m_graph->visit(dupVisit);
        m_graph->setContainmentFlag(false);

        // Remove transitive edges
        SGTransitiveReductionVisitor trVisit;
        m_graph->visit(trVisit);

        // Check for walks that cover seed vertices
        checkWalks(&out_haplotypes);
        
        size_t num_vertices = m_graph->getNumVertices();
        size_t num_tips = findTips().size();

        done = round++ >= MAX_ROUNDS || !out_haplotypes.empty() || num_tips > MAX_TIPS || num_vertices >= MAX_GRAPH_SIZE;
        
#ifdef SHOW_GRAPH
        std::stringstream graph_name;
        graph_name << "graph.r" << round;
        m_graph->writeDot(graph_name.str() + ".dot");
        m_graph->writeASQG(graph_name.str() + ".asqg");
        printf("graph size: %zu tips: %zu round: %d done: %d\n", num_vertices, num_tips, round, done);
#endif
    }

    return HBRC_OK;
}

//
bool OverlapHaplotypeBuilder::buildInitialGraph(const StringVector& reads)
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::buildInitialGraph")
    // Compute initial ordering of reads based on the position of the
    // starting kmer sequence. If the starting kmer was corrected out
    // of a read, it is discarded.
    StringVector ordered_reads;
    orderReadsInitial(m_initial_kmer_string, reads, &ordered_reads);

    if(ordered_reads.size() < m_parameters.minDiscoveryCount)
        return false;

#ifdef SHOW_MULTIPLE_ALIGNMENT
    //DEBUG print MA
    MultipleAlignment ma = buildMultipleAlignment(ordered_reads);
    ma.print(200);
#endif
    // Insert initial reads into graph
    for(size_t i = 0; i < ordered_reads.size(); ++i)
        insertVertexIntoGraph("seed-", ordered_reads[i]);
    return true;
}

//
void OverlapHaplotypeBuilder::extendGraph()
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::extendGraph")
    
    //
    ExtendableTipVector vertices_to_trim;

    // Find tip vertices
    ExtendableTipVector tips = findTips();

    for(size_t i = 0; i < tips.size(); ++i)
    {
        Vertex* x = m_graph->getVertex(tips[i].id);
        assert(x != NULL);
        EdgeDir dir = tips[i].direction;

#ifdef OVERLAP_HAP_DEBUG
        printf("Vertex %s can be extended\n", x->getID().c_str());
#endif

        // Find corrected reads that perfectly overlap this vertex
        StringVector overlapping_reads = getCorrectedOverlaps(x->getSeq().toString(), dir);

#ifdef OVERLAP_HAP_DEBUG
        printf("Found %zu overlaps\n", overlapping_reads.size());
#endif

        // Insert the new reads into the graph
        for(size_t j = 0; j < overlapping_reads.size(); ++j)
        {
            // Determine whether the incoming read is possible join point
            bool is_join = isJoinSequence(overlapping_reads[j], dir);

            std::stringstream label;
            label << (is_join ? "join-" : "extend-");
            label << (dir == ED_ANTISENSE ? "left-" : "right-");
            insertVertexIntoGraph(label.str(), overlapping_reads[j]);
        }

        // Check if x is still a tip in this direction. If so, we trim it and its branch from the graph
        EdgePtrVec x_edges = x->getEdges(dir);
        if(x_edges.size() == 0)
            vertices_to_trim.push_back(tips[i]);
    }
    
    // Trim non-extended vertices
    for(size_t i = 0; i < vertices_to_trim.size(); ++i)
    {
        Vertex* x = m_graph->getVertex(vertices_to_trim[i].id);

        // This vertex may have been removed in a previous iteration
        if(x != NULL)
            trimTip(x, vertices_to_trim[i].direction);
    }
}

// 
void OverlapHaplotypeBuilder::insertVertexIntoGraph(const std::string& prefix, const std::string& sequence)
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::insertVertexIntoGraph")
    // Check if a vertex with this sequence already exists, if so
    // we do not create a new vertex
    if(m_used_reads.find(sequence) != m_used_reads.end())
        return;

    // Create the vertex
    std::stringstream id_ss;
    id_ss << prefix << m_numReads++;
    Vertex* pVertex = new(m_graph->getVertexAllocator()) Vertex(id_ss.str(), sequence);
    m_graph->addVertex(pVertex);

#ifdef OVERLAP_HAP_DEBUG
    std::cout << "Inserting vertex " << id_ss.str() << "\n";
#endif

    // Get a list of vertices that have a k-mer match to this sequence.
    SharedVertexKmerVector candidate_vertices = getCandidateOverlaps(pVertex);

    // Align the sequence against candidate reads
    for(size_t i = 0; i < candidate_vertices.size(); ++i) 
    {
        Vertex* existing_vertex = candidate_vertices[i].vertex;

        // Skip self-edge
        if(existing_vertex == pVertex)
            continue;

        std::string existing_sequence = existing_vertex->getSeq().toString();

        // Try to compute the overlap using the shared kmer as a seed match
        SequenceOverlap overlap;
        std::string shared_kmer = sequence.substr(candidate_vertices[i].kmer_index, m_vertex_map_kmer);
        size_t pos_0 = existing_sequence.find(shared_kmer);
        size_t pos_1 = sequence.find(shared_kmer);
        assert(pos_0 != std::string::npos && pos_1 != std::string::npos);

        // Check for secondary occurrences
        if(existing_sequence.find(shared_kmer, pos_0 + 1) != std::string::npos || 
                    sequence.find(shared_kmer, pos_1 + 1) != std::string::npos) 
        {
            // One of the reads has a second occurrence of the kmer. Use
            // the slow overlapper.
            overlap = Overlapper::computeOverlap(existing_sequence, sequence);
        } 
        else 
        {
            // Seed the match using kmer position
            overlap = Overlapper::extendMatch(existing_sequence, sequence, pos_0, pos_1, 1);
        }

        /*
        printf("Overlap: %s - %s\n", existing_vertex->getID().c_str(), pVertex->getID().c_str());
        overlap.printAlignment(existing_sequence, sequence);
        */

        if(overlap.edit_distance == 0 && overlap.getOverlapLength() >= m_parameters.minOverlap)
        {
            // Add an overlap to the graph
            // Translate the sequence overlap struture into an SGA overlap
            Overlap sga_overlap(existing_vertex->getID(), overlap.match[0].start, overlap.match[0].end, existing_sequence.size(),
                                             id_ss.str(), overlap.match[1].start, overlap.match[1].end, sequence.size(), false, 0);
            SGAlgorithms::createEdgesFromOverlap(m_graph, sga_overlap, true);
        }
    }

    // Insert the sequence into the used reads set
    m_used_reads.insert(sequence);

    // Update kmer map too
    updateKmerVertexMap(pVertex);
}

//
SharedVertexKmerVector OverlapHaplotypeBuilder::getCandidateOverlaps(const Vertex* incoming_vertex)
{
    // Get a list of vertices that have a k-mer match to this sequence.
    SharedVertexKmerVector candidate_vertices;
    const std::string sequence = incoming_vertex->getSeq().toString();
    size_t nk = sequence.size() - m_vertex_map_kmer + 1;
    for(size_t i = 0; i < nk; ++i) 
    {
        std::string query_kmer = sequence.substr(i, m_vertex_map_kmer);
        HashMap<std::string, StringVector>::iterator find_iter = m_kmer_vertex_id_cache.find(query_kmer);
        if(find_iter != m_kmer_vertex_id_cache.end())
        {
            // Get the vertices for these ids
            for(StringVector::iterator iter = find_iter->second.begin(); iter != find_iter->second.end(); ++iter)
            {
                Vertex* candidate = m_graph->getVertex(*iter);
                if(candidate != NULL)
                {
                    SharedVertexKmer svk = { i, candidate };
                    candidate_vertices.push_back(svk);
                }
            }
        }
    }

    // Sort the vector by the vertex pointer and remove duplicates
    // Remove duplicates from the candidate list
    std::sort(candidate_vertices.begin(), candidate_vertices.end(), SharedVertexKmer::sortByVertex);
    SharedVertexKmerVector::iterator new_end = std::unique(candidate_vertices.begin(), candidate_vertices.end(), SharedVertexKmer::equalByVertex);
    candidate_vertices.resize(new_end - candidate_vertices.begin());
    return candidate_vertices;
}

// 
void OverlapHaplotypeBuilder::updateKmerVertexMap(const Vertex* incoming_vertex)
{
    const std::string sequence = incoming_vertex->getSeq().toString();
    size_t nk = sequence.size() - m_vertex_map_kmer + 1;

    // Insert kmers into the kmer pointer cache
    for(size_t i = 0; i < nk; ++i) 
    {
        std::string query_kmer = sequence.substr(i, m_vertex_map_kmer);
        m_kmer_vertex_id_cache[query_kmer].push_back(incoming_vertex->getID());
    }
}

// 
void OverlapHaplotypeBuilder::checkWalks(StringVector* walk_strings)
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::checkWalks")
    // Make a vector of the join sequences
    VertexPtrVec seed_vertices;
    
    VertexPtrVec left_join_vertices;
    VertexPtrVec right_join_vertices;

    VertexPtrVec vertices = m_graph->getAllVertices();
    for(size_t i = 0; i < vertices.size(); ++i)
    {
        // Classify vertex
        Vertex* x = vertices[i];

        if(x->getID().find("join") != std::string::npos)
        {
            if(x->getID().find("left") != std::string::npos && isUniqueJoin(x, ED_ANTISENSE))
                left_join_vertices.push_back(x);

            if(x->getID().find("right") != std::string::npos && isUniqueJoin(x, ED_SENSE))
                right_join_vertices.push_back(x);
        }

        if(x->getID().find("seed") != std::string::npos)
            seed_vertices.push_back(vertices[i]);
    }

#ifdef OVERLAP_HAP_DEBUG
    printf("[%zu %zu join candidates]\n", left_join_vertices.size(), right_join_vertices.size());
#endif

    //
    for(size_t i = 0; i < left_join_vertices.size(); ++i)
    {
        for(size_t j = 0; j < right_join_vertices.size(); ++j)
        {
            /*
            printf("Trying %s -> %s\n", left_join_vertices[i]->getID().c_str(), 
                                        right_join_vertices[j]->getID().c_str());
            */

            // Try to find a walk between this pair of join vertices
            SGWalkVector walks;
            SGSearch::findWalks(left_join_vertices[i], right_join_vertices[j], ED_SENSE, 2000, 10000, true, walks);

            if(!walks.empty())
            {
                // Check if the walk set is self-contained
                bool self_contained = areWalksValid(walks);
                if(!self_contained)
                    continue;

                // Check how many seed vertices this walk covers
                size_t seeds_covered = countCoveredVertices(seed_vertices, walks);
                //printf("    walk covers %zu of %zu seeds\n", seeds_covered, seed_vertices.size());
                double fraction_covered = (double)seeds_covered / seed_vertices.size();
                if(fraction_covered > 0.5)
                {
                    // Add the sequence of the completed walks to the output vector
                    for(size_t i = 0; i < walks.size(); ++i)
                        walk_strings->push_back(walks[i].getString(SGWT_START_TO_END));
                    return;
                }
            }
        }
    }
}

//
bool OverlapHaplotypeBuilder::areWalksValid(const SGWalkVector& walks)
{
    // Build ID set
    std::set<Vertex*> vertex_set;
    for(size_t i = 0; i < walks.size(); ++i)
    {
        VertexPtrVec verts = walks[i].getVertices();
        for(size_t j = 0; j < verts.size(); ++j)
            vertex_set.insert(verts[j]);
    }

    // For every internal vertex of each walk, check that it only has edges to vertices in the set
    for(size_t i = 0; i < walks.size(); ++i)
    {
        VertexPtrVec verts = walks[i].getVertices();
        for(size_t j = 1; j < verts.size() - 1; ++j)
        {
            EdgePtrVec edges = verts[j]->getEdges();
            for(size_t k = 0; k < edges.size(); ++k)
            {
                if(vertex_set.find(edges[k]->getEnd()) == vertex_set.end())
                    return false; // edge to a non-contained vertex
            }
        }
    }
    return true;
}

//
bool OverlapHaplotypeBuilder::isUniqueJoin(Vertex* x, EdgeDir direction)
{
    // To avoid generating multiple paths over the same part of the graph
    // we only use vertices as join points where the graph splits, ends
    // or goes to non-join sequences. This avoids the case
    // where we have chains of unambiguous join vertices
    // and we try to build paths using all pairs of them
    EdgePtrVec dir_edges = x->getEdges(direction);
    EdgePtrVec opp_edges = x->getEdges(!direction);
    if(dir_edges.size() != 1 || opp_edges.size() > 1)
    {
        return true;
    }
    else
    {
        // Single edge, only add if the edge
        // is not to a join vertex
        assert(dir_edges.size() == 1);
        Vertex* y = dir_edges.front()->getEnd();
        return y->getID().find("join") == std::string::npos;
    }
}

//
ExtendableTipVector OverlapHaplotypeBuilder::findTips() const
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::findTips")
    ExtendableTipVector tips;
    VertexPtrVec vertices = m_graph->getAllVertices();
    for(size_t i = 0; i < vertices.size(); ++i)
    {
        size_t edge_count[ED_COUNT];
        edge_count[ED_SENSE] = 0;
        edge_count[ED_ANTISENSE] = 0;

        Vertex* x = vertices[i];
        EdgePtrVec x_edges = x->getEdges();
        for(size_t j = 0; j < x_edges.size(); ++j)
        {
            Edge* xy = x_edges[j];

            // Do not count edges to duplicate vertices
            Overlap o = xy->getOverlap();
            if(!o.isContainment())
                edge_count[xy->getDir()]++;
        }

        // Skip completely disconnected vertices
        if(edge_count[ED_SENSE] == 0 && edge_count[ED_ANTISENSE] == 0)
            continue;

        if(edge_count[ED_SENSE] == 0)
        {
            ExtendableTip tip = { x->getID(), ED_SENSE };
            tips.push_back(tip);
        }

        if(edge_count[ED_ANTISENSE] == 0)
        {
            ExtendableTip tip = { x->getID(), ED_ANTISENSE };
            tips.push_back(tip);
        }
    }
    return tips;
}

//
void OverlapHaplotypeBuilder::trimTip(Vertex* x, EdgeDir direction)
{
    // Check if we should recurse to the neighbors of x
    EdgePtrVec x_opp_edges = x->getEdges(!direction);
    
    /*
    Vertex* x_neighbor = NULL;
    if(x_opp_edges.size() == 1)
        x_neighbor = x_opp_edges.front()->getEnd();
    */

    // Remove x from the graph
    m_graph->removeConnectedVertex(x);

    /*
    // Recurse to x's neighbors if it is now a tip
    if(x_neighbor != NULL)
        trimTip(x_neighbor, direction);
        */
}       

//
bool OverlapHaplotypeBuilder::isJoinSequence(const std::string& sequence, EdgeDir dir)
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
StringVector OverlapHaplotypeBuilder::getCorrectedOverlaps(const std::string& sequence, EdgeDir direction)
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::getCorrectedOverlaps")

    // Extract reads that share a short kmer with the input sequence
    size_t k = 31;
    SequenceOverlapPairVector overlap_vector = KmerOverlaps::retrieveMatches(sequence, k, m_parameters.minOverlap,
                                                                             0.95, 2, m_parameters.variantIndex);

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

    // Correct the reads
    correctReads(&reads);
    assert(reads.size() == overlap_vector.size());

#ifdef SHOW_MULTIPLE_ALIGNMENT
    MultipleAlignment corrected_ma;
    corrected_ma.addBaseSequence("query", sequence, "");
#endif

    // Overlap each corrected read with the input sequence
    // Discard any read that does perfectly match with a long enough overlap
    StringVector out_reads;
    for(size_t i = 0; i < reads.size(); ++i) 
    {
        // Skip uncorrected reads
        if(reads[i].empty())
            continue;

        // Change the strand of the corrected read back to the same as the query sequence
        std::string incoming = reads[i];
        if(overlap_vector[i].is_reversed)
            incoming = reverseComplement(incoming);

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
            out_reads.push_back(incoming);
    }

#ifdef SHOW_MULTIPLE_ALIGNMENT
        corrected_ma.print(200);
#endif
    return out_reads;
}

//
void OverlapHaplotypeBuilder::getReadsForKmers(const StringVector& kmer_vector, size_t k, StringVector* reads)
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::getReadsForKmers")
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

//
size_t OverlapHaplotypeBuilder::countCoveredVertices(const VertexPtrVec& vertices, const SGWalkVector& walks)
{
    size_t count = 0;
    std::set<Vertex*> walk_vertices;
    for(size_t i = 0; i < walks.size(); ++i)
    {
        for(size_t j = 0; j < walks[i].getNumVertices(); ++j)
            walk_vertices.insert(walks[i].getVertex(j));
    }

    for(size_t i = 0; i < vertices.size(); ++i)
    {
        if(walk_vertices.find(vertices[i]) != walk_vertices.end())
            count += 1;
    }

    return count;
}

//
void OverlapHaplotypeBuilder::correctReads(StringVector* reads)
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::correctReads")
    StringVector out_reads;
    for(size_t i = 0; i < reads->size(); ++i)
    {
        // Check if this read has been corrected before
        std::string corrected_sequence;
        HashMap<std::string, std::string>::iterator find_iter = m_correction_cache.find(reads->at(i));
        if(find_iter != m_correction_cache.end())
        {
            corrected_sequence = find_iter->second;
        }
        else
        {
            // Perform correction
            DNAString dna(reads->at(i));
            SeqRecord record = { "null", dna, "" };
            //SeqRecord record("null", dna.toString(), ""); //JEREMIAH
            SequenceWorkItem wi(0, record);
            ErrorCorrectResult r = m_corrector->correct(wi);
            if(r.kmerQC || r.overlapQC)
                corrected_sequence = r.correctSequence.toString();
            
            // Insert the sequence into the cache
            // We allow it to be empty 
            m_correction_cache.insert(std::make_pair(reads->at(i), corrected_sequence));
        }

        // We allow the sequence to be null if the read could not be corrected
        out_reads.push_back(corrected_sequence);
    }
    reads->swap(out_reads);
}

static bool sortPairSecondAscending(const std::pair<size_t, size_t>& a, const std::pair<size_t, size_t>& b)
{
    return a.second > b.second;
}

// Order the list of reads based on the fact that they all share the same kmer
void OverlapHaplotypeBuilder::orderReadsInitial(const std::string& initial_kmer, const StringVector& reads, StringVector* ordered_vector)
{
    std::vector<std::pair<size_t, size_t> > read_kmer_vector;;
    for(size_t i = 0; i < reads.size(); ++i)
    {
        size_t pos = reads[i].find(initial_kmer);
        if(pos != std::string::npos)
            read_kmer_vector.push_back(std::make_pair(i, pos));
    }

    // 
    if(read_kmer_vector.empty())
        return;

    // Sort the vector by the second element
    std::sort(read_kmer_vector.begin(), read_kmer_vector.end(), sortPairSecondAscending);

    // Insert the reads into the ordered list
    for(size_t i = 0; i < read_kmer_vector.size(); ++i)
        ordered_vector->push_back(reads[read_kmer_vector[i].first]);
}

//
MultipleAlignment OverlapHaplotypeBuilder::buildMultipleAlignment(const StringVector& ordered_vector)
{
    MultipleAlignment multiple_alignment;

    // Insert first sequence
    StringVector::const_iterator prev_iterator = ordered_vector.begin();
    multiple_alignment.addBaseSequence("base", *prev_iterator, "");

    // Insert remaining sequences
    StringVector::const_iterator read_iterator = prev_iterator;
    read_iterator++;
    while(read_iterator != ordered_vector.end())
    {
        SequenceOverlap overlap = Overlapper::computeOverlap(*prev_iterator, *read_iterator, ungapped_params);
        multiple_alignment.addExtension("seq", *read_iterator, "", overlap);
        prev_iterator = read_iterator;
        read_iterator++;
    }

    return multiple_alignment;
}

