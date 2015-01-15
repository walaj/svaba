//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// LRAlignment - Collection of algorithms for performing
// long read matches against an FM-index
//
#include "LRAlignment.h"
#include "QuickBWT.h"
#include "MultiAlignment.h"
#include "StdAlnTools.h"

//#define BWA_COMPAT_DEBUG 1
//#define BWA_COMPAT_DEBUG_RESOLVE 1

namespace LRAlignment
{
static const int MINUS_INF = -0x3fffffff;
#define MASK_LEVEL 0.90f

// Initialize the LRCell with default values
void LRCell::initializeDefault()
{
    interval.lower = 0;
    interval.upper = 0;
    I = MINUS_INF;
    D = MINUS_INF;
    G = MINUS_INF;

    parent_cidx = 0;
    q_len = 0;
    t_len = 0;
    parent_idx = -1;
    u_idx = -1;

    clearChildren();
}

// Reset the child indices for the cell
void LRCell::clearChildren()
{
    children_idx[0] = -1;
    children_idx[1] = -1;
    children_idx[2] = -1;
    children_idx[3] = -1;
}

// Returns true if any child of the cell is uninitialized
bool LRCell::hasUninitializedChild() const
{
    return children_idx[0] == -1 || children_idx[1] == -1 || children_idx[2] == -1 ||  children_idx[3] == -1;
}

// Implementation of bwa-sw algorithm.
// Roughly follows Heng Li's implementation 
void bwaswAlignment(const std::string& query, const BWT* pTargetBWT, const SampledSuffixArray* pTargetSSA, const LRParams& params, LRHitVector& outHits)
{
    // Construct an FM-index of the query sequence
    BWT* pQueryBWT = NULL;
    SuffixArray* pQuerySA = NULL;
    createQuickBWT(query, pQueryBWT, pQuerySA);
    
    // Interval hash table used to mark cells as being duplicate
    LRHash dupHash;

    // Initialize the hash table of DAWG nodes
    LRHash dawgHash;
    initializeDAWGHash(pQueryBWT, dawgHash);
    
    // Initialize a stack of elements with a single entry for the root node
    // of the query DAWG
    LRStack stack;
    
    // High scoring alignments are stored as LRHits in these vectors
    // positionHitsVector stores up to 2 hits starting at every base of the query sequence
    // terminalHitsVector stores hits to sequence prefixes 
    LRHitVector positionHitsVector(2*query.size());
    LRHitVector terminalHitsVector;

    // Each dawg node is added to the pendingVector initially
    // Once all predecessors of the node have been visited,
    // the node is moved to the stack
    LRPendingVector pendingVector;
    size_t num_pending = 0;

    // Initialize a stack entry for the root node with an empty scoring cell
    LRStackEntry* pInitial = new LRStackEntry;
    pInitial->interval.lower = 0;
    pInitial->interval.upper = pQueryBWT->getBWLen() - 1;

    //
    LRCell x;
    x.initializeDefault();
    x.G = 0;
    x.interval.lower = 0;
    x.interval.upper = pTargetBWT->getBWLen() - 1;

    pInitial->cells.push_back(x);
    stack.push(pInitial);

    // The main loop of the algorithm - traverse the DAWG
    while(!stack.empty())
    {
        LRStackEntry* v = stack.top();
        stack.pop();
        
        size_t old_n = v->cells.size();

        // TODO: bandwidth test and max depth ?
        
        // Descend into the children of the current dawg node
        // If the interval update succeeds, calculate scores between
        // the child node and all the LRCells of the current node
        for(int qci = 0; qci < DNA_ALPHABET::size; ++qci)
        {
            char query_child_base = DNA_ALPHABET::getBase(qci);
            BWTInterval child_interval = v->interval;
            BWTAlgorithms::updateInterval(child_interval, query_child_base, pQueryBWT);

            if(!child_interval.isValid())
                continue;
    
            // Create an new array of cells for the scores between the child node
            // and all the nodes in the prefix tree
            LRStackEntry* u = new LRStackEntry;
            u->interval = child_interval;

            // Loop over the nodes in v
            for(size_t i = 0; i < v->cells.size(); ++i)
            {
                LRCell* p = &v->cells[i];

                if(p->interval.upper == -1)
                    continue; // duplicate that has been deleted
                LRCell x; // the cell being calculated
                LRCell* c[4]; // pointers to cells required to calculate x

                c[0] = &x; 
                x.G = MINUS_INF;
                x.u_idx = p->u_idx = -1;

                bool add_x_to_u = false;
                if(p->parent_idx >= 0) 
                {
                    int upos = v->cells[p->parent_idx].u_idx;

                    // Set pointers to the cell
                    c[1] = upos >= 0 ? &u->cells[upos] : NULL;
                    c[2] = p;
                    c[3] = &v->cells[p->parent_idx];

                    int match_score = qci == p->parent_cidx ? params.alnParams.match : -params.alnParams.mismatch;
                    int score = fillCells(params, match_score, c);

                    if(score > 0)
                    {
                        // this cell has a positive score
                        // set the x's parent position in u and set the
                        // parent's child position to x
                        x.parent_idx = upos;
                        p->u_idx = u->cells.size(); // x will be added to u in this position
                        
                        if(x.parent_idx >= 0)
                            u->cells[upos].children_idx[p->parent_cidx] = p->u_idx;
                        add_x_to_u = true;
                    }
                }
                else
                {
                    if(p->D > p->G - params.alnParams.gap_open)
                        x.D = p->D - params.alnParams.gap_ext; // extend gap
                    else
                        x.D = p->G - params.alnParams.gap_open_extend; // open a new gap

                    if(x.D > 0)
                    {
                        x.G = x.D;
                        x.I = MINUS_INF;
                        x.parent_idx = -1;
                        p->u_idx = u->cells.size();
                        add_x_to_u = true;
                    }
                }

                if(add_x_to_u)
                {
                    // set the remaining fields in the current cell
                    x.clearChildren();
                    x.parent_cidx = p->parent_cidx;
                    x.interval = p->interval;
                    x.q_len = p->q_len + 1;
                    x.t_len = p->t_len;
                    u->cells.push_back(x);

                    // TODO: heap heuristics from bwa-sw?
                }

                // Check if we should descend into another node of the prefix trie of the target
                if( (x.G > params.alnParams.gap_open_extend /* && heap test */) || i < old_n)
                {
                    if(p->hasUninitializedChild())
                    {
                        for(int tci = 0; tci < DNA_ALPHABET::size; ++tci)
                        {
                            if(p->children_idx[tci] != -1)
                                continue; // already added
                            char target_child_base = DNA_ALPHABET::getBase(tci);
                            BWTInterval target_child_interval = p->interval;
                            BWTAlgorithms::updateInterval(target_child_interval, target_child_base, pTargetBWT);
                            if(!target_child_interval.isValid()) // child with this extension base does not exist
                            {
                                p->children_idx[tci] = -2;
                                continue;
                            }

                            // Create new entry on array v
                            LRCell y;
                            y.G = y.I = y.D = MINUS_INF;
                            y.interval = target_child_interval;
                            y.parent_cidx = tci;
                            y.q_len = p->q_len;
                            y.t_len = p->t_len + 1;


                            y.parent_idx = i;
                            y.clearChildren();

                            p->children_idx[tci] = v->cells.size();
                            v->cells.push_back(y);

                            p = &v->cells[i]; // p may have been invalidated in the push, update
                        } // for tcl
                    } // if has uninitialized
                } // if X.h
            } // for all v
		
            // Save high-scoring cells as hits
            if(!u->cells.empty())
            {
                saveBestPositionHits(pQuerySA, u, params.alnParams.threshold, positionHitsVector);
                //saveAllHits(pQuerySA, pTargetSSA, pTargetBWT, u, params.threshold, terminalHitsVector);
                //saveTerminalHits(pQuerySA, pTargetSSA, pTargetBWT, u, params.threshold, terminalHitsVector);
            }

            // Update the stack by adding u or pushing it to the pending vector
            num_pending += updateStack(&stack, u, &pendingVector, &dawgHash, params);

        } // for qci
        
        // done with v
        delete v;
    } // for all stack

    assert(num_pending == 0);

    // Append the two hits vectors together
    positionHitsVector.insert(positionHitsVector.end(), terminalHitsVector.begin(), terminalHitsVector.end());
    terminalHitsVector.clear();
    resolveDuplicateHits(pTargetBWT, pTargetSSA, positionHitsVector, 2);

    // set the output hits
    outHits.insert(outHits.end(), positionHitsVector.begin(), positionHitsVector.end());

    delete pQueryBWT;
    delete pQuerySA;
}

// Update the stack to contain entry u, after considering any possible merges
// with StackEntries from the pending vector
int updateStack(LRStack* pStack, 
                LRStackEntry* u, 
                LRPendingVector* pPendingVector, 
                LRHash* pDawgHash, 
                const LRParams& params)
{
    // Find the iterator for u in the dawgHash
    LRHash dupHash;

    uint64_t key = u->interval.lower << 32 | u->interval.upper;
    LRHash::iterator hashIter = pDawgHash->find(key);
    assert(hashIter != pDawgHash->end() && (uint32_t)hashIter->second > 0);
    --hashIter->second;
    uint32_t count = (uint32_t)hashIter->second;
    uint32_t position = hashIter->second >> 32;

    // Check if an entry in the pending array exists for this DAWG node
    int change = 0;
    if(position > 0)
    {
        LRStackEntry* w = (*pPendingVector)[position - 1];

        // An entry in the pending array has been created for this query substing
        // Merge u into the interval
        if(!u->cells.empty())
        {
            // Swap so w is the StackEntry wth more cells
            if(w->cells.size() < u->cells.size())
            {
                w = u;
                u = (*pPendingVector)[position - 1];
                (*pPendingVector)[position - 1] = w;
            }
            mergeStackEntries(w, u);
        }

        if(count == 0)
        {
            // this node in the dawg will not be visited again
            // move the stack entry from the pending list to the stack
            removeDuplicateCells(w, dupHash);
            cutTail(w, params);
            pStack->push(w);
            (*pPendingVector)[position - 1] = 0;
            change -= 1;
        }

        // u is empty or merged, it is no longer needed
        delete u;

    }
    else if(count > 0)
    {
        // Create an entry in the pending queue for the current node of the DAWG
        if(!u->cells.empty())
        {
            pPendingVector->push_back(u);
            change += 1;

            // Save the position of u in the pending vector into the hash
            // all subsequent traversals that visit this node of the dawg
            // will get merged into this position. index + 1 is stored
            // so that position == 0 indicates the empty case
            hashIter->second = (uint64_t)pPendingVector->size() << 32 | count;
        }
        else
        {
            // u has no cells to calculate, discard it
            delete u;
        }
    }
    else // count == 0, pos == 0
    {
        // This substring is unique, push u straight onto the stack
        cutTail(u, params);
        pStack->push(u);
    }

    return change;
}

// Process a list of cells and save a best hit for each position
void saveBestPositionHits(const SuffixArray* pQuerySA, LRStackEntry* u, int threshold, LRHitVector& hits)
{
    for(size_t i = 0; i < u->cells.size(); ++i)
    {
        LRCell* p = &u->cells[i];
        if(p->G < threshold)
            continue;
        
        for(int64_t k = u->interval.lower; k <= u->interval.upper; ++k)
        {
            // Calculate the beginning of the alignment using the suffix array
            // of the query
            SAElem e = pQuerySA->get(k);
            int beg = e.getPos();
            int end = beg + p->q_len;

            // Save the best hit for alignments starting at beg in positions hits[2*beg]
            // and the second best hit in hits[2*beg+1]
            LRHit* q = NULL;
            if(p->G > hits[beg*2].G)
            {
                // move the previous best to the second best slot
                hits[beg*2+1] = hits[beg*2];
                q = &hits[2*beg];
            }
            else if(p->G > hits[beg*2+1].G)
            {
                // the current hit is the second best at this position
                q = &hits[2*beg+1];
            }

            if(q)
            {
                q->interval = p->interval;
                q->length = p->t_len;
                q->G = p->G;
                q->q_start = beg;
                q->q_end = end;
                q->flag = 0;
                q->num_seeds = 0;
            }
        }
    }

}


// Process a list of cells and save alignment hits meeting the threshold
void saveAllHits(const SuffixArray* pQuerySA, const SampledSuffixArray* /*pTargetSSA*/,
                 const BWT* /*pTargetBWT*/, LRStackEntry* u, 
                 int threshold, LRHitVector& hits)
{
    for(size_t i = 0; i < u->cells.size(); ++i)
    {
        LRCell* p = &u->cells[i];
        if(p->G < threshold)
            continue;
        
        for(int64_t k = u->interval.lower; k <= u->interval.upper; ++k)
        {
            // Calculate the beginning of the alignment using the suffix array
            // of the query
            SAElem e = pQuerySA->get(k);
            int beg = e.getPos();
            int end = beg + p->q_len;

            // Save the best hit for alignments starting at beg in positions hits[2*beg]
            // and the second best hit in hits[2*beg+1]
            LRHit q;
            q.interval = p->interval;
            q.length = p->t_len;
            q.G = p->G;
            q.q_start = beg;
            q.q_end = end;
            q.flag = 0;
            q.num_seeds = 0;
            hits.push_back(q);
        }
    }
}

// Process all the cells in the stack entry and save hits to cells that represent
// terminal intervals (those intervals containing the start of some reads)
void saveTerminalHits(const SuffixArray* pQuerySA, const SampledSuffixArray* /*pTargetSSA*/, 
                      const BWT* pTargetBWT, LRStackEntry* u, 
                      int threshold, LRHitVector& hits)
{
    for(size_t i = 0; i < u->cells.size(); ++i)
    {
        LRCell* p = &u->cells[i];
        if(p->G < threshold)
            continue;
        
        // If the target bwt interval contains a '$' symbol, push the hit
        AlphaCount64 extCount = BWTAlgorithms::getExtCount(p->interval, pTargetBWT);
        if(extCount.get('$') > 0)
        {
            // Create a hit for every beginning position in the query
            for(int64_t k = u->interval.lower; k <= u->interval.upper; ++k)
            {
                // Calculate the beginning of the alignment using the suffix array
                // of the query
                SAElem e = pQuerySA->get(k);
                int beg = e.getPos();
                int end = beg + p->q_len;

                LRHit q;
                q.interval = p->interval;
                q.length = p->t_len;
                q.G = p->G;
                q.q_start = beg;
                q.q_end = end;
                q.flag = 0;
                q.num_seeds = 0;
                hits.push_back(q);
            }
        }
    }
}

// Initialize the hash table of DAWG nodes by inserting an element
// into the hash table for each distinct substring of the query
void initializeDAWGHash(BWT* pQueryBWT, LRHash& hashTable)
{
    std::stack<uint64_t> stack;
    stack.push(pQueryBWT->getBWLen() - 1);
    while(!stack.empty())
    {
        uint64_t x = stack.top();
        stack.pop();

        BWTInterval interval(x >> 32, (uint32_t)x);
        AlphaCount64 lower = pQueryBWT->getFullOcc(interval.lower - 1);
        AlphaCount64 upper = pQueryBWT->getFullOcc(interval.upper);
        for(int ci = 0; ci < DNA_ALPHABET::size; ++ci)
        {
            char b = DNA_ALPHABET::getBase(ci);    
            size_t pb = pQueryBWT->getPC(b);
            interval.lower = pb + lower.get(b);
            interval.upper = pb + upper.get(b) - 1;

            if(!interval.isValid())
                continue;
            x = (uint64_t)interval.lower << 32 | interval.upper;

            // insert (or update) this key in the hash
            LRHash::iterator iter = hashTable.find(x);
            if(iter == hashTable.end())
            {
                hashTable.insert(std::make_pair(x,1));
                stack.push(x);
            }
            else
            {
                iter->second++;
            }
        }
    }
}

// Merge two stack entries which represent the same node in the DAWG
void mergeStackEntries(LRStackEntry* u, LRStackEntry* v)
{
    // Update parent index and child entries for v cell's to reflect
    // their new positions at the end of array of u
    size_t num_u = u->cells.size();
	for (size_t i = 0; i != v->cells.size(); ++i) 
    {
        LRCell* p = &v->cells[i];
		if(p->parent_idx >= 0)
            p->parent_idx += num_u;
        for(size_t j = 0; j < 4; ++j)
    		if (p->children_idx[j] >= 0) 
                p->children_idx[j] += num_u;
	}
    u->cells.insert(u->cells.end(), v->cells.begin(), v->cells.end());
}

// Remove duplicated cells from the stack entry
void removeDuplicateCells(LRStackEntry* u, LRHash& hash)
{
    hash.clear();
    int n = 0;
    for(size_t i = 0; i < u->cells.size(); ++i)
    {
        int j = -1; // element to be deleted
        LRCell* p = &u->cells[i];

        if(p->interval.upper == -1)
            continue; // already deleted

        uint64_t key = p->interval.lower << 32 | p->interval.upper;
        LRHash::iterator iter = hash.find(key);
        if(iter != hash.end())
        {
            // Check if the score of p is greater than the stored cell
            // with the same key. Set j to be the index of u with the lower
            // score.
            uint32_t currScore = iter->second;
            if((int)currScore >= p->G)
            {
                j = i;
            }
            else
            {
                j = iter->second >> 32;
                iter->second = i << 32 | p->G;
            }
        }
        else
        {
            // initialize hash for cell i
            hash[key] = i << 32 | p->G;
        }

        if(j >= 0)
        {
#ifdef BWA_COMPAT_DEBUG
            printf("marking p [%zu %zu] as deleted cidx: %d [DUP]\n", p->interval.lower, p->interval.upper, p->parent_cidx);
#endif
            p = &u->cells[j];
            p->interval.lower = -1;
            p->interval.upper = -1;
            p->G = 0;
            if(p->parent_idx >= 0)
                u->cells[p->parent_idx].children_idx[p->parent_cidx] = -3;
            n += 1;
        }
    }
#ifdef BWA_COMPAT_DEBUG
    printf("removed %d duplicate entries\n", n);
#endif
}

// Fill the values of C[0] depending on the values in the other 3 cells
int fillCells(const LRParams& params, int match_score, LRCell* c[4])
{
	int G = c[3] ? c[3]->G + match_score : MINUS_INF;
	if(c[1]) 
    {
		if(c[1]->I > c[1]->G - params.alnParams.gap_open)
            c[0]->I = c[1]->I - params.alnParams.gap_ext; // extend gap
        else
            c[0]->I = c[1]->G - params.alnParams.gap_open_extend; // open new gap
		if (c[0]->I > G)
            G = c[0]->I; // new best score
	} 
    else
    {
        c[0]->I = MINUS_INF;
    }

	if(c[2])
    {
		if(c[2]->D > c[2]->G - params.alnParams.gap_open)
            c[0]->D = c[2]->D - params.alnParams.gap_ext; // extend gap
        else
            c[0]->D = c[2]->G - params.alnParams.gap_open_extend; // open new gap

		if (c[0]->D > G) 
            G = c[0]->D; // new best score
	} 
    else
    {
        c[0]->D = MINUS_INF;
    }
    
	return(c[0]->G = G);
}

// Remove duplicate hits to the same target sequence from the hits vector
int resolveDuplicateHitsByID(const BWT* pTargetBWT, const SampledSuffixArray* pTargetSSA, LRHitVector& hits, int /*IS*/)
{
    if(hits.empty())
        return 0;
    assert(pTargetBWT != NULL);
    assert(pTargetSSA != NULL);

    // Convert each hit interval to target coordinates
    LRHitVector newHits;
    for(size_t i = 0; i < hits.size(); ++i)
    {
        LRHit* p = &hits[i];
#ifdef BWA_COMPAT_DEBUG_RESOLVE
        printf("Total hits: %zu\n", hits.size());
        printf("Reallocating hits array to size: %d\n", n);
#endif
    
        if(!p->interval.isValid() || p->G <= 0)
            continue;

        for(int64_t k = p->interval.lower; k <= p->interval.upper; ++k)
        {
            LRHit tmp = *p;
            SAElem elem = pTargetSSA->calcSA(k, pTargetBWT);
            tmp.targetID = elem.getID();
            tmp.t_start = elem.getPos();
            tmp.interval.lower = 0;
            tmp.interval.upper = -1;
#ifdef BWA_COMPAT_DEBUG_RESOLVE
            printf("Created new hit at position (%d, %d)\n", tmp.targetID, tmp.position);
#endif
            newHits.push_back(tmp);
        }
    }
    
    // Swap the new hits structure with the old hits
    hits.swap(newHits);

    // sort hits by targetID, then score
    std::sort(hits.begin(), hits.end(), LRHit::compareIDandG);
    
    // Loop over the hits, zeroing all but the highest scoring hit for each target
    uint64_t prevID = std::numeric_limits<uint64_t>::max();
    size_t new_n = 0;
    for(size_t i = 0; i < hits.size(); ++i)
    {
//        std::cout << "ID: " << hits[i].targetID << " score: " << hits[i].G << "\n";
        if(prevID == hits[i].targetID)
            hits[i].G = 0;
        else
            new_n += 1;
        prevID = hits[i].targetID;
    }

    std::sort(hits.begin(), hits.end(), LRHit::compareG);
    hits.resize(new_n);
    std::cout << "FINAL HITS: " << new_n << "\n";
#ifdef BWA_COMPAT_DEBUG_RESOLVE
    printf("Final number of hits after duplicate removal: %zu\n", hits.size());
#endif
    return hits.size();
}
// Remove duplicated hits from the hits vector
int resolveDuplicateHits(const BWT* pTargetBWT, const SampledSuffixArray* pTargetSSA, LRHitVector& hits, int IS)
{
    if(hits.empty())
        return 0;
    assert(pTargetBWT != NULL);

    if(pTargetBWT != NULL && pTargetSSA != NULL)
    {
        // Convert each hit to target coordinates
        int old_n = hits.size();
        int n = 0;
        
        // Count the number of valid hits 
        LRHitVector newHits;
        for(size_t i = 0; i < hits.size(); ++i)
        {
            LRHit* p = &hits[i];
            if(p->interval.isValid() && p->interval.size() <= IS)
                n += p->interval.size();
            else if(!p->interval.isValid()) // bwa compatibility hack
                n += 1; 
            else if(p->G > 0)
                n += 1;
        }
#ifdef BWA_COMPAT_DEBUG_RESOLVE
        printf("Total hits: %zu\n", hits.size());
        printf("Reallocating hits array to size: %d\n", n);
#endif
        newHits.resize(n);
    
        // Populate the newHits vector with the hits to keep
        int j = 0;
        for(int i = 0; i < old_n; ++i)
        {
            LRHit* p = &hits[i];
            if(p->interval.isValid() && p->interval.size() <= IS)
            {
                for(int64_t k = p->interval.lower; k <= p->interval.upper; ++k)
                {
                    newHits[j] = *p;
                    SAElem elem = pTargetSSA->calcSA(k, pTargetBWT);
                    newHits[j].targetID = elem.getID();
                    newHits[j].t_start = elem.getPos();
                    newHits[j].interval.lower = 0;
                    newHits[j].interval.upper = -1;
#ifdef BWA_COMPAT_DEBUG_RESOLVE
                    printf("Created new hit at position (%d, %d)\n", newHits[j].targetID, newHits[j].position);
#endif
                    ++j;
                }
            }
            else if(p->G > 0)
            {
                newHits[j] = *p;
                SAElem elem = pTargetSSA->calcSA(p->interval.lower, pTargetBWT);
                newHits[j].targetID = elem.getID();
                newHits[j].t_start = elem.getPos();
                newHits[j].interval.lower = 0;
                newHits[j].interval.upper = -1;
                newHits[j].flag |= 1;
#ifdef BWA_COMPAT_DEBUG_RESOLVE
                printf("Created new hit at position (%d, %d)\n", newHits[j].targetID, newHits[j].position);
#endif
                ++j;
            }
        }

        // Swap the new hits structure with the old hits
        hits.swap(newHits);
    }

    // sort hits by score
    std::sort(hits.begin(), hits.end(), LRHit::compareG);
    
    int i,j;
    // Remove hits with significant overlaps on the query or target sequences
    for(i = 1; i < (int)hits.size(); ++i)
    {
        LRHit* p = &hits[i];
        if(p->G == 0)
            break;
        for(j = 0; j < i; ++j)
        {
            LRHit* q = &hits[j];
            bool compatible = true;

            // hit q has already been removed
            if(q->G == 0)
                continue;

            // if the hits are to different targets they can't overlap
            if(p->targetID != q->targetID)
                continue;

            if(p->interval.upper == -1 && q->interval.upper == -1)
            {
                // Calculate the overlap length on the query
                int qol = (p->q_end < q->q_end? p->q_end : q->q_end) - (p->q_start > q->q_start? p->q_start : q->q_start);
                if(qol < 0)
                    qol = 0;

                if((float)qol / (p->q_end - p->q_start) > MASK_LEVEL || (float)qol / (q->q_end - q->q_start) > MASK_LEVEL) 
                {
                    int64_t tol = 0;
                    if(p->t_start + p->length < q->t_start + q->length)
                        tol = p->t_start + p->length;
                    else
                        tol = q->t_start + q->length;

                    if(p->t_start > q->t_start)
                        tol -= p->t_start;
                    else
                        tol -= q->t_start;
                    
					if((double)tol / p->length > MASK_LEVEL || (double)tol / q->length > MASK_LEVEL)
						compatible = false;

#ifdef BWA_COMPAT_DEBUG_RESOLVE
                    printf("    idx = (%d,%d) id=(%d,%d) G=(%d,%d) qol: %d tol: %zu compat: %d\n", 
                            i,j, (int)p->targetID, (int)q->targetID, p->G, q->G, qol, tol, compatible);
#endif
                    if(!compatible)
                    {
                        p->G = 0;
                        break;
                    }
                }
            }
        }
    }

    int new_n = i;

    for(i = j = 0; i < new_n; ++i)
    {
        if(hits[i].G == 0)
            continue;
        if(i != j)
            hits[j++] = hits[i];
        else
            ++j;
    }
    hits.resize(j);
#ifdef BWA_COMPAT_DEBUG_RESOLVE
    printf("Final number of hits after duplicate removal: %zu\n", hits.size());
#endif
    return hits.size();
}

// Remove cells from a stack entry based on some criteria
void cutTail(LRStackEntry* u, const LRParams& params)
{
    switch(params.cutTailAlgorithm)
    {
        case LRCA_Z_BEST:
            cutTailByZBest(u, params);
            break;
        case LRCA_Z_BEST_STRATA:
            cutTailByStratifiedZBest(u, params);
            break;
        case LRCA_SCORE_FRAC:
            cutTailByScorePercent(u, params);
            break;
        case LRCA_NONE:
            return;
        default:
            assert(false);
    }
}

// Remove cells from a stack entry (query dawg node) by their score
// as expressed as a fraction of the maximum possible score for the
// query
void cutTailByScorePercent(LRStackEntry* u, const LRParams& params)
{
    //printf("cutTail starting split score: %2.lf\n", split);
    for(size_t i = 0; i < u->cells.size(); ++i)
    {
        LRCell* p = &u->cells[i];
        double cellScore = (double)p->G * 100.0f / (p->q_len * 1);

        if(cellScore < params.percentCutoff)
        {
#ifdef BWA_COMPAT_DEBUG
            if(p->interval.upper > -1)
                printf("marking p [%d %d] as deleted [CT]\n", (int)p->interval.lower, (int)p->interval.upper);
#endif
            //printf("cutTail cutting cell with qlen: %d %d %2.1lf\n", p->q_len, p->G, cellScore);
            p->interval.lower = 0;
            p->interval.upper = -1;
            p->G = 0;
            if(p->parent_idx >= 0)
                u->cells[p->parent_idx].children_idx[p->parent_cidx] = -1;
        }
//        else
//            printf("cutTail keeping cell with qlen: %d %d %2.1lf\n", p->q_len, p->G, cellScore);

    }
}

// Remove cells from a stack entry, only keeping the top zBest scoring nodes
void cutTailByZBest(LRStackEntry* u, const LRParams& params)
{
    assert(params.zBest > 0);
    // Save an int vector of scores
    IntVector scores;
    for(size_t i = 0; i < u->cells.size(); ++i)
    {
        LRCell* p = &u->cells[i];
        if(p->interval.upper != -1 && p->G > 0)
            scores.push_back(-p->G);
    }

    if((int)scores.size() <= params.zBest)
        return;

    // Partially sort the scores to select the T-th best score
    std::nth_element(scores.begin(), scores.begin() + params.zBest, scores.end());
    double split = -scores[params.zBest];

#ifdef BWA_COMPAT_DEBUG
    printf("[CT] split score: %2.lf\n", split);
#endif

    int n = 0;
    //printf("cutTail starting split score: %2.lf\n", split);
    for(size_t i = 0; i < u->cells.size(); ++i)
    {
        LRCell* p = &u->cells[i];
        if(p->G == split)
            ++n;
        if(p->G < split || (p->G == split && n >= params.zBest))
        {
#ifdef BWA_COMPAT_DEBUG
            if(p->interval.upper > -1)
                printf("marking p [%d %d] as deleted [CT]\n", (int)p->interval.lower, (int)p->interval.upper);
#endif
            //printf("cutTail cutting cell with qlen: %d %d %2.1lf\n", p->q_len, p->G, cellScore);
            p->interval.lower = 0;
            p->interval.upper = -1;
            p->G = 0;
            if(p->parent_idx >= 0)
                u->cells[p->parent_idx].children_idx[p->parent_cidx] = -1;
        }
        /*else
            printf("cutTail keeping cell with qlen: %d %d %2.1lf\n", p->q_len, p->G, cellScore);
        */

    }
}

// Remove cells from a stack entry, only keeping the top zBest scoring nodes for each query length
void cutTailByStratifiedZBest(LRStackEntry* u, const LRParams& params)
{
    assert(params.zBest > 0);
    if(u->cells.empty())
        return;

    // Calculate the range of query lengths
    int minQ = std::numeric_limits<int>::max();
    int maxQ = 0;
    for(size_t i = 0; i < u->cells.size(); ++i)
    {
        LRCell* p = &u->cells[i];
        if(p->q_len < minQ)
            minQ = p->q_len;
        if(p->q_len > maxQ)
            maxQ = p->q_len;
    }

    // Calculate cut points for each strata
    int span = maxQ - minQ + 1;
    assert(span > 0);

    //
    std::vector<IntVector> scores(span);
    IntVector cutoffs(span);
    IntVector counts(span);

    for(size_t i = 0; i < u->cells.size(); ++i)
    {
        LRCell* p = &u->cells[i];
        int strata = p->q_len - minQ;
        assert(strata >= 0 && strata < (int)scores.size());

        if(p->interval.upper != -1 && p->G > 0)
            scores[strata].push_back(-p->G);
    }

    // Perform a partial sort of each strata and determine cutoff points
    for(int i = 0; i < span; ++i)
    {
        if((int)scores[i].size() <= params.zBest)
        {
            cutoffs[i] = 0; // keep all positively-scored nodes
            continue;
        }

        // Partially sort the scores to select the T-th best score
        std::sort(scores[i].begin(), scores[i].end());
        cutoffs[i] = -scores[i][params.zBest];
    }

#ifdef BWA_COMPAT_DEBUG
    printf("[CT] split score: %2.lf\n", split);
#endif

    //printf("cutTail starting split score: %2.lf\n", split);
    for(size_t i = 0; i < u->cells.size(); ++i)
    {
        LRCell* p = &u->cells[i];
        int strata = p->q_len - minQ;
        int c = cutoffs[strata];

        if(p->G == c)
            ++counts[strata];
        if(p->G < c || (p->G == c && counts[strata] >= params.zBest))
        {
#ifdef BWA_COMPAT_DEBUG
            if(p->interval.upper > -1)
                printf("marking p [%d %d] as deleted [CT]\n", (int)p->interval.lower, (int)p->interval.upper);
#endif
            //printf("cutTail cutting cell with qlen: %d %d %2.1lf\n", p->q_len, p->G, cellScore);
            p->interval.lower = 0;
            p->interval.upper = -1;
            p->G = 0;
            if(p->parent_idx >= 0)
                u->cells[p->parent_idx].children_idx[p->parent_cidx] = -1;
        }
        /*else
            printf("cutTail keeping cell with qlen: %d %d %2.1lf\n", p->q_len, p->G, cellScore);
        */

    }
}


// Construct a multiple alignment from a vector of hits
MultiAlignment convertHitsToMultiAlignment(const std::string& query, 
                                           const BWT* pTargetBWT, 
                                           const SampledSuffixArray* /*pTargetSSA*/,
                                           const LRParams& params,
                                           const LRHitVector& hits)
{
    // Set up stdaln data structures
	size_t max_target = StdAlnTools::calculateMaxTargetLength(query.size(), params.alnParams);

    path_t* path = (path_t*)calloc(max_target + query.size(), sizeof(path_t));
    
    AlnParam par;
    int matrix[25];
    StdAlnTools::setAlnParam(par, matrix, params.alnParams);

    MAlignDataVector mAlignVec;
    uint8_t* pQueryT = StdAlnTools::createPacked(query);

    for(size_t i = 0; i < hits.size(); ++i)
    {
        LRHit tempHit = hits[i];
        assert(tempHit.targetID != (uint64_t)-1);
        std::string target = BWTAlgorithms::extractString(pTargetBWT, tempHit.targetID);
        uint8_t* pTargetT = StdAlnTools::createPacked(target);

        //extendHitFullLength(tempHit, pQueryT, pTargetT, query.size(), target.size(), &par);

        int q_start_pos, q_end_pos, t_start_pos, t_end_pos;
        q_start_pos = tempHit.q_start;
        q_end_pos = tempHit.q_end; // exclusive
        t_start_pos = tempHit.t_start;
        t_end_pos = t_start_pos + tempHit.length; // exclusive
        
        // Get pointers to the substrings of interest
        int ql = q_end_pos - q_start_pos;
        uint8_t* pQuerySub = pQueryT + q_start_pos; 
        
        int tl = t_end_pos - t_start_pos;
        uint8_t* pTargetSub = pTargetT + t_start_pos;
        
        int path_len = 0;
		/*int score =*/ aln_global_core(pTargetSub, tl, pQuerySub, ql, &par, path, &path_len);
        int cigarLen = 0;
		uint32_t* pCigar = aln_path2cigar32(path, path_len, &cigarLen);
        
        //
        MAlignData maData;
        maData.str = target.substr(t_start_pos, t_end_pos - t_start_pos);
        maData.position = q_start_pos;

        // Add initial padding to cigar
        maData.expandedCigar.append(maData.position, 'S');
        
        // Add alignment symbols
        for (int j = 0; j != cigarLen; ++j)
            maData.expandedCigar.append(pCigar[j]>>4, "MID"[pCigar[j]&0xf]);

        /*        
        printf("AlignmentDP (%zu)\n", hits[i].targetID);
        printf("Coords: q[%d, %d] t[%d, %d]\n", q_start_pos, q_end_pos, t_start_pos, t_end_pos);
        printf("C: %s\n", maData.expandedCigar.c_str());
        */
        mAlignVec.push_back(maData);

        // Cleanup
        delete [] pTargetT;
        free(pCigar);
    }

    // Cleanup
    free(path);
    delete [] pQueryT;

    return MultiAlignment(query, mAlignVec);
}

// Attempt to extend a hit to the left and right using aln_extend_core
// from stdaln
void extendHitFullLength(LRHit& hit, uint8_t* pQueryPacked, uint8_t* pTargetPacked, 
                         int q_length, int t_length, AlnParam* pStdAlnPar)
{
    int q_start = hit.q_start;
    int q_end = hit.q_end; // exclusive
    int t_start = hit.t_start;
    int t_end = t_start + hit.length; // exclusive

    // Attempt left-extension
    if(t_start != 0 && q_start != 0)
    {
        int q_diff = q_start;
        int t_diff = t_start;
        int d = std::min(q_diff, t_diff);
        q_start -= d;
        t_start -= d;
    }

    // Attempt right-extension
    if(q_end < q_length && t_end < t_length)
    {
        // extension match lengths
        int tel = t_length - t_start;
        int qel = q_length - q_start;
        uint8_t* pTargetPackedSub = pTargetPacked + t_start;
        uint8_t* pQueryPackedSub = pQueryPacked + q_start;
        path_t path;
        int score = aln_extend_core(pTargetPackedSub, tel, pQueryPackedSub, qel, pStdAlnPar, &path, 0, hit.G, NULL);
        if(score > hit.G)
        {
            hit.G = score;
            hit.length = path.i;
            hit.q_end = path.j + hit.q_start;
        }
    }
}



};
