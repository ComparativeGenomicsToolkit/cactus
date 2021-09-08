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
    cigar->length = atoll(*c);
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
    paf->query_length = atoll(stList_get(tokens, 1));
    paf->query_start = atoll(stList_get(tokens, 2));
    paf->query_end = atoll(stList_get(tokens, 3));

    // Is the alignment forward or reverse
    char strand = ((char *)stList_get(tokens, 4))[0];
    if(strand != '+' && strand != '-') {
        st_errAbort("Got an unexpected strand character (%c) in a paf string: %s\n", strand, paf_string);
    }
    paf->same_strand = strand == '+';

    // Get target coordinates
    paf->target_name = stString_copy(stList_get(tokens, 5));
    paf->target_length = atoll(stList_get(tokens, 6));
    paf->target_start = atoll(stList_get(tokens, 7));
    paf->target_end = atoll(stList_get(tokens, 8));

    // Get the core alignment metric attributes of the record
    paf->num_matches = atoll(stList_get(tokens, 9));
    paf->num_bases = atoll(stList_get(tokens, 10));
    paf->mapping_quality = atoll(stList_get(tokens, 11));

    // Parse the remaining optional tags
    for(int64_t i=12; i<stList_length(tokens); i++) {
        stList *tag = stString_splitByString(stList_get(tokens, i), ":");
        char *type = stList_get(tag, 0);
        if(strcmp(type, "tp") == 0) {
            paf->type = ((char *)stList_get(tag, 2))[0];
            assert(paf->type == 'P' || paf->type == 'S' || paf->type == 'I');
        } else if (strcmp(type, "AS") == 0) {
            paf->score = atoll(stList_get(tag, 2));
        } else if(strcmp(type, "cg") == 0) {
            paf->cigar = parse_cigar(stList_get(tag, 2));
        }
        else if(strcmp(type, "tl") == 0) {
            paf->tile_level = atoll(stList_get(tag, 2));
        }
        stList_destruct(tag);
    }

    // Cleanup
    stList_destruct(tokens);

    return paf;
}

Paf *paf_read(FILE *fh) {
    char *c = stFile_getLineFromFile(fh);
    if(c == NULL) {
        return NULL;
    }
    Paf *paf = paf_parse(c);
    free(c);

    paf_check(paf);

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
    int64_t i = sprintf(buffer, "%s\t%" PRIi64 "\t%" PRIi64"\t%" PRIi64"\t%c\t%s\t%" PRIi64"\t%" PRIi64"\t%" PRIi64
                                "\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64,
                        paf->query_name, paf->query_length, paf->query_start, paf->query_end,
                        paf->same_strand ? '+' : '-',
                        paf->target_name, paf->target_length, paf->target_start, paf->target_end,
                        paf->num_matches, paf->num_bases, paf->mapping_quality);
    if(paf->type != '\0') {
        i += sprintf(buffer+i, "\ttp:A:%c", paf->type);
    }
    if(paf->score != INT_MAX) {
        i += sprintf(buffer+i, "\tAS:i:%" PRIi64, paf->score);
    }
    if(paf->tile_level != -1) {
        i += sprintf(buffer+i, "\ttl:i:%" PRIi64, paf->tile_level);
    }
    if(i > buf_size) {
        st_errAbort("Size of paf record exceeded buffer size\n");
    }
    if(paf->cigar != NULL) {
        i += sprintf(buffer+i, "\tcg:Z:");
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
    if(paf->query_start > paf->query_end || paf->query_end > paf->query_length) {
        st_errAbort("Paf query end coordinates are invalid, %s", paf_print(paf));
    }
    if(paf->target_start < 0 || paf->target_start >= paf->target_length) {
        st_errAbort("Paf target start coordinates are invalid, %s", paf_print(paf));
    }
    if(paf->target_start > paf->target_end || paf->target_end > paf->target_length) {
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

int64_t paf_get_number_of_aligned_bases(Paf *paf) {
    int64_t aligned_bases = 0;
    Cigar *c = paf->cigar;
    while(c != NULL) {
        if(c->op == match) {
            aligned_bases += c->length;
        }
        c = c->next;
    }
    return aligned_bases;
}

static Cigar *cigar_reverse(Cigar *c) {
    if(c == NULL) {
        return NULL;
    }
    Cigar *p = NULL;
    // p -> c -> n -> q
    while(c->next != NULL) {
        Cigar *n = c->next;
        c->next = p; // p <- c , n -> q
        p = c;
        c = n;
    }
    c->next = p; // p <- c <- n
    return c;
}

static Cigar *cigar_trim(int64_t *query_c, int64_t *target_c, Cigar *c, int64_t end_bases_to_trim, int sign) {
    int64_t bases_trimmed = 0;
    while(c != NULL && (c->op != match || bases_trimmed < end_bases_to_trim)) {
        if(c->op == match) { // can trim this alignment
            if(bases_trimmed + c->length > end_bases_to_trim) {
                int64_t i = end_bases_to_trim - bases_trimmed;
                c->length -= i;
                *query_c += sign*i;
                *target_c += sign*i;
                assert(c->length > 0);
                break;
            }
            bases_trimmed += c->length;
            *query_c += sign*c->length;
            *target_c += sign*c->length;
        }
        else if(c->op == query_insert) {
            *query_c += sign*c->length;
        }
        else {
            assert(c->op == query_delete);
            *target_c += sign*c->length;
        }
        Cigar *c2 = c;
        c = c->next;
        free(c2);
    }
    return c;
}

void paf_trim_ends(Paf *paf, int64_t end_bases_to_trim) {
    // Trim front end
    paf->cigar = cigar_trim(&(paf->query_start), &(paf->target_start), paf->cigar, end_bases_to_trim, 1);
    // Trim back end
    paf->cigar = cigar_reverse(cigar_trim(&(paf->query_end), &(paf->target_end), cigar_reverse(paf->cigar), end_bases_to_trim, -1));
}

void paf_trim_end_fraction(Paf *p, float percentage) {
    // trim
    assert(percentage >= 0 && percentage <= 1.0);
    int64_t aligned_bases = paf_get_number_of_aligned_bases(p);
    int64_t end_bases_to_trim = (aligned_bases * percentage)/2.0;
    st_logDebug("For alignment of %" PRIi64 " query bases, %"
    PRIi64 " target bases and %" PRIi64 " aligned bases trimming %" PRIi64 " bases from each paf end\n",
            p->query_end - p->query_start, p->target_end - p->target_start, aligned_bases, end_bases_to_trim);
    paf_trim_ends(p, end_bases_to_trim);
}
