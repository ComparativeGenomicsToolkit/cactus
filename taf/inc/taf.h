#ifndef STTAF_H_
#define STTAF_H_

#include "sonLib.h"

/*
 * Structures to represent blocks of an alignment
 */
typedef struct _row Alignment_Row;

typedef struct _alignment {
    Alignment_Row *row; // An alignment is just a sequence of rows
} Alignment;

struct _row { // Each row encodes the information about an aligned sequence
    char *sequence_name; // name of sequence
    int64_t start, length, sequence_length; // zero based, half open coordinates
    bool strand; // nonzero is "+" else "-"
    char *bases; // [A-Za-z*+]* string
    Alignment_Row *l_row;  // connection to a left row (may be NULL)
    Alignment_Row *r_row;  // connection to a right row (may be NULL)
    Alignment_Row *n_row;  // the next row in the alignment
};

/*
 * Clean up the memory for an alignment
 */
void alignment_destruct(Alignment *alignment);

/*
 * Use the O(ND) alignment to diff the rows between two alignments and connect together their rows
 * so that we can determine which rows in the right_alignment are a continuation of rows in the
 * left_alignment. We use this for efficiently outputting TAF.
 */
void alignment_link_adjacent(Alignment *left_alignment, Alignment *right_alignment);


/*
 * Read a maf alignment block from the file stream. Return NULL if none available
 */
Alignment *maf_read_block(FILE *fh);

/*
 * Write a maf block
 */
void maf_write_block(Alignment *alignment, FILE *fh);

/*
 * Read a taf column.
 */
Alignment *taf_read_column(FILE *fh, Alignment *p_column);

/*
 * Write a taf block.
 */
void taf_write_block(Alignment *p_alignment, Alignment *alignment, FILE *fh, bool run_length_encode_bases);

#endif /* STTAF_H_ */

