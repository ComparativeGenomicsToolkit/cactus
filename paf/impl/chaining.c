#include "paf.h"
#include "../inc/paf.h"

/*
 * Functions for chaining together pafs
 */

static int intcmp(int64_t i, int64_t j) {
    return i > j ? 1 : (i < j ? -1 : 0);
}

static int paf_cmp_by_location(const void *a, const void *b) {
    Paf *p1 = (Paf *)a, *p2 = (Paf *)b;
    int i = strcmp(p1->query_name, p2->query_name);
    if(i == 0) {
        i = strcmp(p1->target_name, p2->target_name);
        if(i == 0) {
            i = intcmp(p1->query_start, p2->query_start);
            if(i == 0) {
                i = intcmp(p1->target_start, p2->target_start);
            }
        }
    }
    return i;
}

typedef struct _chain Chain;
struct _chain {
    Paf *paf;
    int64_t score;
    Chain *pChain;
};

static int chain_cmp_by_location(const void *a, const void *b) {
    Paf *p1 = ((Chain *)a)->paf, *p2 = ((Chain *)b)->paf;
    int i = strcmp(p1->query_name, p2->query_name);
    if(i == 0) {
        i = strcmp(p1->target_name, p2->target_name);
        if(i == 0) {
            i = intcmp(p1->target_end, p2->target_end);
            if(i == 0) {
                i = intcmp(p1->query_end, p2->query_end);
            }
        }
    }
    return i;
}

static int chain_cmp_by_score(const void *a, const void *b) {
    Chain *c1 = (Chain *)a, *c2 = (Chain *)b;
    return intcmp(c1->score, c2->score);
}

static stSortedSetIterator *get_chains(stSortedSet *chained_alignments, Chain *chain) {
    int64_t i=chain->paf->query_end, j=chain->paf->target_end;
    chain->paf->query_end = chain->paf->query_start;
    chain->paf->target_end = chain->paf->target_start;
    Chain *chain2 = stSortedSet_searchLessThan(chained_alignments, chain);
    chain->paf->query_end = i;
    chain->paf->target_end = j;
    if(chain2 == NULL) {
        return stSortedSet_getIterator(chained_alignments);
    }
    return stSortedSet_getIteratorFrom(chained_alignments, chain2);
}

static Paf *chain_to_paf(Chain *chain) {
    Paf *p = chain->paf;

    while(chain->pChain != NULL) { // Join together links in the chain
        Paf *q = p;
        p = chain->pChain->paf;
        // p is not the previous paf in the chain and q the subsequent one
        p->query_end = q->query_end;
        p->target_end = q->target_end;
        p->score = q->score;

        // If there is a cigar then join them
        if(p->cigar != NULL) {
            Cigar *cigar = p->cigar;
            while(cigar->next != NULL) {
                cigar = cigar->next;
            }
            cigar->next = p->cigar;
        }
        Chain *c = chain;
        chain = chain->pChain;
        free(q);
        free(c);
    }

    free(chain); // Now free the final wrapper in the chain

    paf_check(p);

    return p;
}

stList *paf_chain(stList *pafs, int64_t (*gap_cost)(int64_t, void *), void *gap_cost_params, int64_t max_gap_length) {
    stList_sort(pafs, paf_cmp_by_location); // Sort alignments by chromosome, forward-or-reverse strand alignment,
    // query start coordinate, then target query coordinate

    stSortedSet *chained_alignments = stSortedSet_construct3(chain_cmp_by_location, free); // The set of chained alignments,
    // sorted by chromosome, forward-or-reverse strand alignment,
    // target end coordinate, then query end coordinate. Each record is of type Chain

    stList *to_remove = stList_construct(); // List of chained alignments to remove from the chain_alignments array at the end of each loop

    // For each alignment
    for(int64_t i=0; i<stList_length(pafs); i++) {
        Paf *paf = stList_get(pafs, i);
        Chain *chain = st_calloc(1, sizeof(Chain));
        chain->paf = paf;
        chain->score = paf->score;

        // Find highest scoring chains that alignment could be chained with:
        stSortedSetIterator *it = get_chains(chained_alignments, chain); // This is an iterator over chains that the alignment could be joined to
        Chain *pChain;
        while((pChain = stSortedSet_getPrevious(it)) != NULL) {
            // If the query gap is larger than max_gap_length then we can remove it from the set that can be chained
            if(paf->query_start - pChain->paf->query_end > max_gap_length) {
                stList_append(to_remove, pChain);
            }
            else if(paf->target_start - pChain->paf->target_end > max_gap_length) { // If the target gap is longer than we can chain to
                // then there are no more alignments we can chain to
                break;
            }
            else { // We can chain to this alignment
                int64_t chain_score = pChain->score + gap_cost(paf->query_start - pChain->paf->query_end, gap_cost_params) +
                                      gap_cost(paf->target_start - pChain->paf->target_end, gap_cost_params);
                if(chain_score > chain->score) {
                    chain->score = chain_score;
                    chain->pChain = pChain;
                }
            }
        }

        // Add the paf to the chained alignments
        stSortedSet_insert(chained_alignments, chain);

        // Remove any chainable alignments
        while(stList_length(to_remove) > 0) {
            stSortedSet_remove(chained_alignments, stList_pop(to_remove));
        }
    }

    // Get chains, from highest scoring to lowest
    stList *chains = stSortedSet_getList(chained_alignments);
    stList_sort(chains, chain_cmp_by_score);
    stList *outputChains = stList_construct();
    while(stList_length(chains) > 0) { // For each alignment, in order of score (high to low)
        Chain *chain = stList_pop(chains);
        if(stSortedSet_search(chained_alignments, chain) != NULL) { // If not already in a chain
            stList_append(outputChains, chain);
            while(1) {
                stSortedSet_remove(chained_alignments, chain); // Indicate the alignment is in a chain
                if(chain->pChain == NULL) {
                    break; // End of the chain, break
                }
                if (stSortedSet_search(chained_alignments, chain->pChain) == NULL) { // If the previous link in the chain is already in a chain
                    chain->pChain = NULL; // Break the current chain
                    break; // As we're now at the end of the chain, break
                }
                chain = chain->pChain;
            }
        }
    }

    // Now finally convert chains back to single pafs
    for(int64_t i=0; i<stList_length(outputChains); i++) {
        stList_set(outputChains, i, chain_to_paf(stList_get(outputChains, i)));
    }

    // Cleanup
    assert(stList_length(to_remove) == 0);
    stList_destruct(to_remove);
    assert(stSortedSet_size(chained_alignments) == 0);
    stSortedSet_destruct(chained_alignments);
    assert(stList_length(chains) == 0);
    stList_destruct(chains);

    return outputChains;
}

