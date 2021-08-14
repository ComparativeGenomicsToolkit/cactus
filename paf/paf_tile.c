#include "paf.h"

/*
 * paf_tile: Greedily assign a "level" to each alignment, so that alignments with the lowest level are tesslating, high-scoring alignments.
 * (1) Load local alignments files (PAF)
 * (2) Sort alignments by score, from best-to-worst
 * (3) Create integer array representing counts of alignments to bases in the genome, setting values initially to 0.
 * (4) Iterate through alignments, from best-to-worst, find maximum?? count, q, of aligned bases to a base covered by the alignment,
 * set the "level" of the alignment to q+1, increase by one the aligned bases count of each base covered by the alignment.
 * (5) Output local alignments file, adding the alignment level tag to each alignment, sorted by score from best-to-worst
 */

int main(int argc, char *argv[]) {

}

