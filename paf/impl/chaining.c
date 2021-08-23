#include "paf.h"
#include "../inc/paf.h"

/*
 * Functions for chaining together pafs
 */

static int intcmp(int64_t i, int64_t j) {
    return i > j ? 1 : (i < j ? -1 : 0);
}

/*
 * Compare two pafs by query start coordinates - Note: doesn't care about query sequence name
 */
static int paf_cmp_by_query_location(const void *a, const void *b) {
    Paf *p1 = (Paf *)a, *p2 = (Paf *)b;
    return intcmp(p1->query_start, p2->query_start);
}

/*
 * Represents a chain of paf alignments
 */
typedef struct _chain Chain;
struct _chain {
    Paf *paf;
    int64_t score;
    Chain *pChain;
};

/*
 * Compare chains first by query sequence name, then target sequence name,
 * then target end coordinate,
 * and finally by query end coordinate.
 */
static int chain_cmp_by_location(const void *a, const void *b) {
    Paf *p1 = ((Chain *)a)->paf, *p2 = ((Chain *)b)->paf;
    int i = strcmp(p1->query_name, p2->query_name);
    if(i == 0) {
        i = strcmp(p1->target_name, p2->target_name);
        if (i == 0) {
            i = intcmp(p1->target_end, p2->target_end);
            if (i == 0) {
                i = intcmp(p1->query_end, p2->query_end);
            }
        }
    }
    return i;
}

/*
 * Compare two chains by score.
 */
static int chain_cmp_by_score(const void *a, const void *b) {
    Chain *c1 = (Chain *)a, *c2 = (Chain *)b;
    return intcmp(c1->score, c2->score);
}

/*
 * Get chains that the paf contained in "Chain" could be chained with.
 */
static stSortedSetIterator *get_predecessor_chains(stSortedSet *active_chained_alignments, Chain *chain) {
    int64_t i=chain->paf->query_end, j=chain->paf->target_end;
    chain->paf->query_end = chain->paf->query_start;
    chain->paf->target_end = chain->paf->target_start;
    Chain *chain2 = stSortedSet_searchLessThanOrEqual(active_chained_alignments, chain);
    chain->paf->query_end = i;
    chain->paf->target_end = j;
    if(chain2 == NULL) {
        return stSortedSet_getIterator(active_chained_alignments);
    }
    stSortedSetIterator *it = stSortedSet_getIteratorFrom(active_chained_alignments, chain2);
    stSortedSet_getNext(it);
    stSortedSet_getNext(it);
    //assert(stSortedSet_getPrevious(it) == chain2);
    return it;
}

/*
 * Joins together the pafs in a chain into a single paf
 */
static Paf *chain_to_paf(Chain *chain) {
    Paf *p = chain->paf;

    while(chain->pChain != NULL) { // Join together links in the chain
        Paf *q = p;
        p = chain->pChain->paf; // p is now the previous paf in the chain to q

        // Checks that we can chain these together
        assert(strcmp(p->target_name, q->target_name) == 0);
        assert(strcmp(p->query_name, q->query_name) == 0);
        assert(p->query_end <= q->query_start);
        assert(p->target_end <= q->target_start);
        assert(p->same_strand == q->same_strand);

        // If there is a cigar then join them
        if(p->cigar != NULL) {
            // Find the last op in p
            Cigar *cigar = p->cigar;
            while(cigar->next != NULL) {
                cigar = cigar->next;
            }
            // Add in any needed indels to account for gaps between p and q
            if(p->query_end < q->query_start) {
                cigar->next = st_calloc(1, sizeof(Cigar)); // Allocate a cigar record
                cigar->next->length = q->query_start - p->query_end;
                cigar->next->op = query_insert;
                cigar = cigar->next;
            }
            if(p->target_end < q->target_start) {
                cigar->next = st_calloc(1, sizeof(Cigar)); // Allocate a cigar record
                cigar->next->length = q->target_start - p->target_end;
                cigar->next->op = query_delete;
                cigar = cigar->next;
            }
            // Link them together
            cigar->next = q->cigar;
        }

        // Update the coordinates
        p->query_end = q->query_end;
        p->target_end = q->target_end;
        p->score = q->score;

        // Shift back to the prior link in the chain and cleanup
        Chain *c = chain;
        chain = chain->pChain;
        free(q);
        free(c);
    }

    free(chain); // Now free the final wrapper in the chain

    return p;
}

/*
 * Chains together the input pafs. Ignores strand.
 */
stList *paf_chain_ignore_strand(stList *pafs, int64_t (*gap_cost)(int64_t, void *),
                                void *gap_cost_params, int64_t max_gap_length) {
    stList_sort(pafs, paf_cmp_by_query_location); // Sort alignments by query start coordinate

    stSortedSet *active_chained_alignments = stSortedSet_construct3(chain_cmp_by_location, NULL); // The set of
    // alignments being chained,
    // sorted by chromosome, target end coordinate, then query end coordinate. Each record is of type Chain

    stSortedSet *chains = stSortedSet_construct3(chain_cmp_by_score, free); // The set of all chains, sorted by score

    stList *to_remove = stList_construct(); // List of chained alignments to remove from the active_chained_alignments
    // array at the end of each loop

    // For each alignment
    for(int64_t i=0; i<stList_length(pafs); i++) {
        Paf *paf = stList_get(pafs, i);
        Chain *chain = st_calloc(1, sizeof(Chain));
        chain->paf = paf;
        chain->score = paf->score;

        // Find highest scoring chains that alignment could be chained with:
        stSortedSetIterator *it = get_predecessor_chains(active_chained_alignments, chain); // This is an iterator over chains that the alignment could be joined to

        //stSortedSetIterator *it = stSortedSet_size(active_chained_alignments) > 0 ?
        //        stSortedSet_getIteratorFrom(active_chained_alignments, stSortedSet_getLast(active_chained_alignments)) :
        //                          stSortedSet_getIterator(active_chained_alignments);
        //stSortedSet_getNext(it);
        //stSortedSet_getNext(it);

        Chain *pChain;
        while((pChain = stSortedSet_getPrevious(it)) != NULL) {

            if(strcmp(paf->query_name, pChain->paf->query_name) != 0 ||
               strcmp(paf->target_name, pChain->paf->target_name) != 0 ||
               paf->same_strand != pChain->paf->same_strand) {
                break; // Can not chain, and no further predecessors can exist
            }

            if(paf->query_start < pChain->paf->query_end) { // Chain ends of the query after the current paf starts,
                // so can not chain, but further predecessors may exist
                continue;
            }

            // If the query gap is larger than max_gap_length then we can remove it from the set that can be chained
            if(paf->query_start - pChain->paf->query_end > max_gap_length) {
                stList_append(to_remove, pChain);
                continue; // further predecessors may exist
            }

            if (paf->target_start < pChain->paf->target_end) {
                continue; // Can not chain, but further predecessors may exist
            }
            if (paf->target_start - pChain->paf->target_end > max_gap_length) { // If the target gap is longer
                // than we can chain to then there are no more alignments we can chain to
                break;
            } else { // We can chain to this alignment
                int64_t chain_score = paf->score + pChain->score -
                        gap_cost(paf->query_start - pChain->paf->query_end, gap_cost_params) -
                        gap_cost(paf->target_start - pChain->paf->target_end, gap_cost_params);
                if (chain_score > chain->score) {
                    chain->score = chain_score;
                    chain->pChain = pChain;
                }
            }
        }
        stSortedSet_destructIterator(it);

        // Add the paf to the chained alignments
        stSortedSet_insert(active_chained_alignments, chain);

        // Add the chain to the final set of chains
        stSortedSet_insert(chains, chain);

        // Remove any chainable alignments
        while(stList_length(to_remove) > 0) {
            stSortedSet_remove(active_chained_alignments, stList_pop(to_remove));
        }
    }

    // Get chains, from highest scoring to lowest
    stList *outputChains = stList_construct();
    while(stSortedSet_size(chains) > 0) { // For each alignment, in order of score (high to low)
        Chain *chain = stSortedSet_remove(chains, stSortedSet_getLast(chains));
        stList_append(outputChains, chain);
        while(1) {
            if(chain->pChain == NULL) {
                break; // End of the chain, break
            }
            if (stSortedSet_search(chains, chain->pChain) == NULL) { // If the previous link in the chain is already in a chain
                chain->pChain = NULL; // Break the current chain
                break; // As we're now at the end of the chain, break
            }
            chain = chain->pChain;
            assert(stSortedSet_search(chains, chain) == chain);
            stSortedSet_remove(chains, chain);
        }
    }

    // Now finally convert chains back to single pafs
    for(int64_t i=0; i<stList_length(outputChains); i++) {
        stList_set(outputChains, i, chain_to_paf(stList_get(outputChains, i)));
    }

    // Cleanup
    assert(stList_length(to_remove) == 0);
    stList_destruct(to_remove);
    stSortedSet_destruct(active_chained_alignments);
    assert(stSortedSet_size(chains) == 0);
    stSortedSet_destruct(chains);

    return outputChains;
}

/*
 * Makes it so that a reverse strand alignment can be chained by "mirroring" the target sequence coordinates
 */
void invert_target_strand(Paf *p) {
    int64_t i = p->target_start;
    p->target_start = -p->target_end;
    p->target_end = -i;
}

stList *paf_chain(stList *pafs, int64_t (*gap_cost)(int64_t, void *), void *gap_cost_params, int64_t max_gap_length) {
    // Split into forward and reverse strand alignments
    stList *positive_strand_pafs = stList_construct();
    stList *negative_strand_pafs = stList_construct();
    for(int64_t i=0; i<stList_length(pafs); i++) {
        Paf *p = stList_get(pafs, i);
        if(p->same_strand) {
            stList_append(positive_strand_pafs, p);
        }
        else {
            invert_target_strand(p);
            stList_append(negative_strand_pafs, p);
        }
    }

    stList *positive_chained_pafs = paf_chain_ignore_strand(positive_strand_pafs, gap_cost, gap_cost_params, max_gap_length);
    stList *negative_chained_pafs = paf_chain_ignore_strand(negative_strand_pafs, gap_cost, gap_cost_params, max_gap_length);

    // Correct negative strand coordinates
    for(int64_t i=0; i<stList_length(negative_chained_pafs); i++) {
        invert_target_strand(stList_get(negative_chained_pafs, i));
    }

    stList_appendAll(positive_chained_pafs, negative_chained_pafs);

    // Cleanup
    stList_setDestructor(negative_chained_pafs, NULL);
    stList_destruct(negative_chained_pafs);

    // Check
    for(int64_t i=0; i<stList_length(positive_chained_pafs); i++) {
        paf_check(stList_get(positive_chained_pafs, i));
    }

    return positive_chained_pafs;
}
