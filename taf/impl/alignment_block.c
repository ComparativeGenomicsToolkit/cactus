#include "taf.h"
#include "ond.h"

void alignment_destruct(Alignment *alignment) {
    Alignment_Row *row = alignment->row;
    while(row != NULL) {
        Alignment_Row *r = row;
        row = row->n_row;
        free(r);
    }
    free(alignment);
}

static stList *get_rows_in_a_list(Alignment_Row *row) {
    stList *l = stList_construct();
    while(row != NULL) {
        stList_append(l, row);
        row = row->n_row;
    }
    return l;
}

static bool rows_equal(Alignment_Row *left_row, Alignment_Row *right_row) {
    // Do the rows match
    return strcmp(left_row->sequence_name, right_row->sequence_name) == 0 && left_row->strand == right_row->strand &&
            left_row->start + left_row->length <= right_row->start;
}

void alignment_link_adjacent(Alignment *left_alignment, Alignment *right_alignment) {
    stList *left_rows = get_rows_in_a_list(left_alignment->row);
    stList *right_rows = get_rows_in_a_list(right_alignment->row);
    // get the alignment of the rows
    WFA *wfa = WFA_construct(left_rows, right_rows, (bool (*)(void *, void *))rows_equal, 1, 1); // Use unit gap and mismatch costs for the diff
    stList *aligned_rows = WFA_get_alignment(wfa);
    // connect up the rows
    assert(stList_length(aligned_rows) % 2 == 0); // must be even length
    for(int64_t i=0; i<stList_length(aligned_rows); i+=2) {
        Alignment_Row *left_row = stList_get(aligned_rows, i);
        Alignment_Row *right_row = stList_get(aligned_rows, i+1);
        left_row->r_row = right_row;
        right_row->l_row = left_row;
    }
    // clean up
    stList_destruct(left_rows);
    stList_destruct(right_rows);
    stList_destruct(aligned_rows);
    WFA_destruct(wfa);
}

