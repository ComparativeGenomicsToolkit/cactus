#include "taf.h"

void alignment_destruct(Alignment *alignment) {
    Alignment_Row *row = alignment->row;
    while(row != NULL) {
        Alignment_Row *r = row;
        row = row->n_row;
        free(r);
    }
    free(alignment);
}

void alignment_link_adjacent(Alignment *left_alignment, Alignment *right_alignment) {
    Alignment_Row *left_row = left_alignment->row;
    Alignment_Row *right_row = right_alignment->row;

    while(left_row != NULL && right_row != NULL) {
        // Do the rows match
        if(strcmp(left_row->sequence_name, right_row->sequence_name) == 0) {
            if(left_row->strand == right_row->strand) {
                if(left_row->start + left_row->length <= right_row->start) {
                    right_row->p_row = left_row;
                    left_row = left_row->n_row;
                    right_row = right_row->n_row;
                    continue;
                }
            }
        }
        left_row = left_row->n_row;
        right_row = right_row->n_row;
    }
}

