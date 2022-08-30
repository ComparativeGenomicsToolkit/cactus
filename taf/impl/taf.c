#include "taf.h"

Alignment *taf_read_column(FILE *fh, Alignment *p_column) {
    return NULL;
}

void write_column(Alignment_Row *row, int64_t column, FILE *fh, bool run_length_encode_bases) {
    char base = '\0';
    int64_t base_count = 0;
    while(row != NULL) {
        if(row->bases[column] == base) {
            base_count++;
        }
        else {
            if(base != '\0') {
                if(run_length_encode_bases) {
                    fprintf(fh, "%c %" PRIi64 " ", base, base_count);
                }
                else {
                    for (int64_t i = 0; i < base_count; i++) {
                        fprintf(fh, "%c", base);
                    }
                }
            }
            base = row->bases[column];
            base_count = 1;
        }
        row = row->n_row;
    }
    if(base != '\0') {
        if(run_length_encode_bases) {
            fprintf(fh, "%c %" PRIi64 " ", base, base_count);
        }
        else {
            for (int64_t i = 0; i < base_count; i++) {
                fprintf(fh, "%c", base);
            }
        }
    }
}

void write_coordinates(Alignment_Row *p_row, Alignment_Row *row, FILE *fh) {
    int64_t i = 0;
    fprintf(fh, " ;");
    while(p_row != NULL) { // Write any row deletions
        if(p_row->r_row == NULL) { // if the row is deleted
            fprintf(fh, " d %" PRIi64 "", i);
        }
        p_row = p_row->n_row; i++;
    }
    i = 0;
    while(row != NULL) { // Now write the new rows
        if(row->l_row == NULL) { // if the row is inserted
            fprintf(fh, " i %" PRIi64 " %s %" PRIi64 " %c %" PRIi64 "",
                    i, row->sequence_name, row->start, row->strand ? '+' : '-', row->sequence_length);
        }
        else {
            if(strcmp(row->l_row->sequence_name, row->sequence_name) == 0 &&
                row->l_row->strand == row->strand &&
                row->l_row->start + row->l_row->length <= row->start) {
                if(row->l_row->start + row->l_row->length < row->start) { // if there is an indel
                    fprintf(fh, " g %" PRIi64 " %" PRIi64 "",
                            i, row->start - (row->l_row->start + row->l_row->length));
                }
            }
            else { // Substitute one row for another
                fprintf(fh, " s %" PRIi64 " %s %" PRIi64 " %c %" PRIi64 "",
                        i, row->sequence_name, row->start, row->strand ? '+' : '-', row->sequence_length);
            }
        }
        row = row->n_row; i++;
    }
}

void taf_write_block(Alignment *p_alignment, Alignment *alignment, FILE *fh, bool run_length_encode_bases) {
    Alignment_Row *row = alignment->row;
    if(row != NULL) {
        int64_t column_no = strlen(row->bases);
        assert(column_no > 0);
        write_column(row, 0, fh, run_length_encode_bases);
        write_coordinates(p_alignment != NULL ? p_alignment->row : NULL, row, fh);
        fprintf(fh, "\n");
        for(int64_t i=1; i<column_no; i++) {
            write_column(row, i, fh, run_length_encode_bases);
            fprintf(fh, "\n");
        }
    }
}
