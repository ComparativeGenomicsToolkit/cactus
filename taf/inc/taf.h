#ifndef STTAF_H_
#define STTAF_H_

#include "sonLib.h"

/*
 * a
 * s       simDog.chr6     437451  11      +       593897  CCCGTCAGTGT
 * s       simHuman.chr6   446327  11      +       601863  TCCGCCAAGGT
 * s       simMouse.chr6   460751  11      +       636262  TTCATCAGAGT
 * s       simRat.chr6     470339  11      +       647215  TTCATTAGGGT
 * a
 * s       simCow.chr6     445326  8       +       602619  TTTTCCCA
 * s       simDog.chr6     437462  8       +       593897  TT-TTCCG
 * s       simHuman.chr6   446338  8       +       601863  TTCTTCCG
 * s       simMouse.chr6   460762  8       +       636262  TTTTACCG
 * s       simRat.chr6     470350  8       +       647215  TTTTACCG
 */

typedef struct _row Alignment_Row;

typedef struct _alignment {
    Alignment_Row *row;
} Alignment;

void alignment_destruct(Alignment *alignment);

void alignment_link_adjacent(Alignment *left_alignment, Alignment *right_alignment);

struct _row {
    char *sequence_name; // name of sequence
    int64_t start, length, sequence_length; // zero based, half open coordinates
    bool strand; // nonzero is "+" else "-"
    char *bases; // [A-Za-z*+]* string
    Alignment_Row *p_row;  // any connection to a previous row
    Alignment_Row *n_row;  // the next row in the alignment
};

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
void taf_write_block(Alignment *alignment, FILE *fh);

#endif /* STTAF_H_ */

