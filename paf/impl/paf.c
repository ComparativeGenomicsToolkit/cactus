#include "paf.h"
#include <ctype.h>
#include "../inc/paf.h"

/*
 * Library functions for manipulating paf files.
 */

void paf_destruct(Paf *paf) {
    Cigar *c = paf->cigar;
    while(c != NULL) { // Cleanup the individual cigar records
        Cigar *c2 = c->next;
        free(c);
        c = c2;
    }
    free(paf);
}

static Cigar *parse_cigar_record(char **c) {
    Cigar *cigar = st_calloc(1, sizeof(Cigar)); // Allocate a cigar record
    // Calculate the number of characters representing the length of the record
    int64_t i=0;
    while(isdigit((*c)[i])) {
        i++;
    }
    char t = (*c)[i]; // The type of the cigar operation
    if(t != 'M' && t != 'I' && t != 'D') {
        st_errAbort("Got an unexpected character paf cigar string: %c\n", t);
    }
    cigar->op = t == 'M' ? match : (t == 'I' ? query_insert : query_delete);
    (*c)[i] = ' ';
    cigar->length = atoi(*c);
    *c = &((*c)[i+1]);
    return cigar;
}

static Cigar *parse_cigar(char *cigar_string) {
    if(cigar_string[0] == '\0') { // If is the empty string
        return NULL;
    }
    Cigar *cigar = parse_cigar_record(&cigar_string);
    Cigar *pCigar = cigar;
    while(cigar_string[0] != '\0') {
        Cigar *nCigar = parse_cigar_record(&cigar_string);
        pCigar->next = nCigar;
        pCigar = nCigar;
    }
    return cigar;
}

Paf *paf_parse(char *paf_string) {
    Paf *paf = st_calloc(1, sizeof(Paf));

    stList *tokens = stString_split(paf_string); // Tokenize the record

    // Get query coordinates
    paf->query_name = stString_copy(stList_get(tokens, 0));
    paf->query_length = atoi(stList_get(tokens, 1));
    paf->query_start = atoi(stList_get(tokens, 2));
    paf->query_end = atoi(stList_get(tokens, 3));

    // Is the alignment forward or reverse
    char strand = ((char *)stList_get(tokens, 4))[0];
    if(strand != '+' && strand != '-') {
        st_errAbort("Got an unexpected strand character (%c) in a paf string: %s\n", strand, paf_string);
    }
    paf->same_strand = strand == '+';

    // Get target coordinates
    paf->target_name = stString_copy(stList_get(tokens, 5));
    paf->target_length = atoi(stList_get(tokens, 6));
    paf->target_start = atoi(stList_get(tokens, 7));
    paf->target_end = atoi(stList_get(tokens, 8));

    // Get the core alignment metric attributes of the record
    paf->num_matches = atoi(stList_get(tokens, 9));
    paf->num_bases = atoi(stList_get(tokens, 10));
    paf->mapping_quality = atoi(stList_get(tokens, 11));

    // Parse the remaining optional tags
    for(int64_t i=12; i<stList_length(tokens); i++) {
        stList *tag = stString_splitByString(stList_get(tokens, i), ":");
        char *type = stList_get(tag, 0);
        if(strcmp(type, "tp") == 0) {
            paf->type = ((char *)stList_get(tag, 2))[0];
            assert(paf->type == 'P' || paf->type == 'S' || paf->type == 'I');
        } else if (strcmp(type, "AS") == 0) {
            paf->score = atoi(stList_get(tag, 2));
        } else if(strcmp(type, "cg") == 0) {
            paf->cigar = parse_cigar(stList_get(tag, 2));
        }
        stList_destruct(tag);
    }

    // Cleanup
    stList_destruct(tokens);

    return paf;
}

Paf *paf_read(FILE *fh) {
    char *c = stFile_getLineFromFile(fh);
    //st_uglyf("Got the following line: %s\n", c);
    if(c == NULL) {
        return NULL;
    }
    Paf *paf = paf_parse(c);
    //st_uglyf("Got the following record: %s\n", paf_print(paf));
    free(c);
    return paf;
}

int64_t cigar_number_of_records(Paf *paf) {
    int64_t i=0;
    Cigar *c = paf->cigar;
    while(c != NULL) {
        i++;
        c = c->next;
    }
    return i;
}

char *paf_print(Paf *paf) {
    // Generous estimate of size needed for each paf record.
    int64_t buf_size = 12 * cigar_number_of_records(paf) + 130 + strlen(paf->query_name) + strlen(paf->target_name);
    char *buffer = st_malloc(sizeof(char) * buf_size); // Giving a generous
    int64_t i = sprintf(buffer, "%s %" PRIi64 " %" PRIi64" %" PRIi64" %c %s %" PRIi64" %" PRIi64" %" PRIi64
                                " %" PRIi64 " %" PRIi64 " %" PRIi64,
                        paf->query_name, paf->query_length, paf->query_start, paf->query_end,
                        paf->same_strand ? '+' : '-',
                        paf->target_name, paf->target_length, paf->target_start, paf->target_end,
                        paf->num_matches, paf->num_bases, paf->mapping_quality);
    if(paf->type != '\0') {
        i += sprintf(buffer+i, " tp:A:%c", paf->type);
    }
    if(paf->score != INT_MAX) {
        i += sprintf(buffer+i, " AS:i:%" PRIi64, paf->score);
    }
    if(i > buf_size) {
        st_errAbort("Size of paf record exceeded buffer size\n");
    }
    if(paf->cigar != NULL) {
        i += sprintf(buffer+i, " cg:Z:");
        Cigar *c = paf->cigar;
        while(c != NULL) {
            i += sprintf(buffer+i, "%" PRIi64 "%c", c->length, c->op == match ? 'M' : (c->op == query_insert ? 'I' : 'D'));
            c = c->next;
            if(i > buf_size) {
                st_errAbort("Size of paf record exceeded buffer size\n");
            }
        }
    }
    if(i > buf_size) {
        st_errAbort("Size of paf record exceeded buffer size\n");
    }
    return buffer;
}

void paf_write(Paf *paf, FILE *fh) {
    char *c = paf_print(paf);
    fprintf(fh, "%s\n", c);
    free(c);
}

void paf_check(Paf *paf) {
    if(paf->query_start < 0 || paf->query_start >= paf->query_length) {
        st_errAbort("Paf query start coordinates are invalid, %s", paf_print(paf));
    }
    if(paf->query_start > paf->query_end || paf->query_end >= paf->query_length) {
        st_errAbort("Paf query end coordinates are invalid, %s", paf_print(paf));
    }
    if(paf->target_start < 0 || paf->target_start >= paf->target_length) {
        st_errAbort("Paf target start coordinates are invalid, %s", paf_print(paf));
    }
    if(paf->target_start > paf->target_end || paf->target_end >= paf->target_length) {
        st_errAbort("Paf target end coordinates are invalid, %s", paf_print(paf));
    }
    if(paf->cigar != NULL) {
        // Check that cigar alignment, if present, matches the alignment:
        int64_t i=0, j=0;
        Cigar *cigar = paf->cigar;
        while(cigar != NULL) {
            if(cigar->op == match || cigar->op == query_insert) {
                i += cigar->length;
            }
            if(cigar->op == match || cigar->op == query_delete) {
                j += cigar->length;
            }
            cigar = cigar->next;
        }
        if(i != paf->query_end - paf->query_start) {
            st_errAbort("Paf cigar alignment does not match query length: %" PRIi64 " vs. %" PRIi64 " %s", i,
                        paf->query_end - paf->query_start, paf_print(paf));
        }
        if(j != paf->target_end - paf->target_start) {
            st_errAbort("Paf cigar alignment does not match target length: %" PRIi64 " vs. %" PRIi64 " %s", j,
                        paf->target_end - paf->target_start, paf_print(paf));
        }
    }
}

static void swap(void **a, void **b) {
    void *c = *a;
    *a = *b;
    *b = c;
}

void paf_invert(Paf *paf) {
    swap((void **)&paf->query_start, (void **)&paf->target_start);
    swap((void **)&paf->query_end, (void **)&paf->target_end);
    swap((void **)&paf->query_length, (void **)&paf->target_length);
    swap((void **)&paf->query_name, (void **)&paf->target_name);

    Cigar *c = paf->cigar;
    while(c != NULL) {
        if(c->op == query_insert) {
            c->op = query_delete;
        }
        else if(c->op == query_delete) {
            c->op = query_insert;
        }
        c = c->next;
    }
}

stList *read_pafs(FILE *fh) {
    stList *pafs = stList_construct3(0, (void (*)(void *))paf_destruct);
    Paf *paf;
    while((paf = paf_read(fh)) != NULL) {
        paf_check(paf);
        stList_append(pafs, paf);
    }
    return pafs;
}

void write_pafs(FILE *fh, stList *pafs) {
    for(int64_t i=0; i<stList_length(pafs); i++) {
        paf_write(stList_get(pafs, i), fh);
    }
}

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

stSortedSetIterator *get_chains(stSortedSet *chained_alignments, Chain *chain) {
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

    // Cleanup
    assert(stList_length(to_remove) == 0);
    stList_destruct(to_remove);
    assert(stSortedSet_size(chained_alignments) == 0);
    stSortedSet_destruct(chained_alignments);
    assert(stList_length(chains) == 0);
    stList_destruct(chains);

    return outputChains;
}

