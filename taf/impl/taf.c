#include "taf.h"

Alignment *taf_read_column(FILE *fh, Alignment *p_column) {

}

void write_column(Alignment_Row *row, int64_t column, FILE *fh) {
    char base = '\0';
    int64_t base_count = 0;
    while(row != NULL) {
        if(row->bases[column] == base) {
            base_count++;
        }
        else {
            if(base != '\0') {
                fprintf(fh, "%c %" PRIi64 " ", base, base_count);
            }
            base = row->bases[column];
            base_count = 1;
        }
        row = row->n_row;
    }
    if(base != '\0') {
        fprintf(fh, "%c %" PRIi64 " ", base, base_count);
    }
}

void write_coordinates(Alignment_Row *row, FILE *fh) {
    while(row != NULL) {
        if(row->p_row == NULL) {
            fprintf(fh, "a %s %" PRIi64 " %c %" PRIi64 " ",
                    row->sequence_name, row->start, row->strand ? '+' : '-', row->sequence_length);
        }
        row = row->n_row;
    }
}

void taf_write_block(Alignment *alignment, FILE *fh) {
    Alignment_Row *row = alignment->row;
    if(row != NULL) {
        int64_t column_no = strlen(row->bases);
        assert(column_no > 0);
        write_column(row, 0, fh);
        write_coordinates(row, fh);
        fprintf(fh, "\n");
        for(int64_t i=1; i<column_no; i++) {
            write_column(row, i, fh);
            fprintf(fh, "\n");
        }
    }
}
