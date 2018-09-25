/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLib.h"
#include "pairwiseAlignment.h"

/*
 * Script takes a set of pairwise alignments using the lastz cigar format and returns a modified
 * set such that all alignments are reported with respect to the positive strand of the first sequence
 * and such that all alignments are mirrored, so that they are additionally reporting after flipping the first
 * sequence for the second sequence.
 */

void invertStrands(struct PairwiseAlignment *pairwiseAlignment) {
	/*
	 * Inverts the strands of the alignment.
	 */
	// Flips the strands of first sequence

	if(pairwiseAlignment->start1 != pairwiseAlignment->end1) { // If alignment has non zero length on the first sequence
		int64_t start = pairwiseAlignment->start1;
		pairwiseAlignment->start1 = pairwiseAlignment->end1;
		pairwiseAlignment->end1 = start;
	}
	pairwiseAlignment->strand1 = pairwiseAlignment->strand1 ? 0 : 1;

	if(pairwiseAlignment->start1 != pairwiseAlignment->end1) { // If alignment has non zero length on the second sequence
		int64_t start = pairwiseAlignment->start2;
		pairwiseAlignment->start2 = pairwiseAlignment->end2;
		pairwiseAlignment->end2 = start;
	}
	pairwiseAlignment->strand2 = pairwiseAlignment->strand2 ? 0 : 1;

	// Invert the order of the operations
	listReverse(pairwiseAlignment->operationList);
}

void cigarReverse(struct PairwiseAlignment *pairwiseAlignment) {
	/*
	 * Flips the query and target sequences
	 */

	// Swap the 1s and 2s
	char *contig1 = pairwiseAlignment->contig1;
	int64_t start1 = pairwiseAlignment->start1;
	int64_t end1 = pairwiseAlignment->end1;
	int64_t strand1 = pairwiseAlignment->strand1;

	pairwiseAlignment->contig1 = pairwiseAlignment->contig2;
	pairwiseAlignment->start1 = pairwiseAlignment->start2;
	pairwiseAlignment->end1 = pairwiseAlignment->end2;
	pairwiseAlignment->strand1 = pairwiseAlignment->strand2;

	pairwiseAlignment->contig2 = contig1;
	pairwiseAlignment->start2 = start1;
	pairwiseAlignment->end2 = end1;
	pairwiseAlignment->strand2 = strand1;

	// Invert the operations
	struct AlignmentOperation *op;
	for(int64_t i=0; i<pairwiseAlignment->operationList->length; i++) {
		op = pairwiseAlignment->operationList->list[i];
		assert(op->length >= 0);
		if(op->opType == PAIRWISE_INDEL_Y) {
			op->opType = PAIRWISE_INDEL_X;
		}
		else if(op->opType == PAIRWISE_INDEL_X) {
			op->opType = PAIRWISE_INDEL_Y;
		}
	}
}

int main(int argc, char *argv[]) {
	/*
	 * For each alignment in the input file copy the alignment to the output file and additionally
	 * write out the alignment with the first and second sequences reversed. For each alignment written out
	 * we ensure the alignment is reported with respect to the positive strand of the first reported sequence.
	 */
	st_setLogLevelFromString(argv[1]);

    FILE *fileHandleIn = stdin;
	FILE *fileHandleOut = stdout;

	if(argc == 4) {
		fileHandleIn = fopen(argv[2], "r");
		fileHandleOut = fopen(argv[3], "w");
	}
	else {
		assert(argc == 2);
	}

    struct PairwiseAlignment *pairwiseAlignment;

    while ((pairwiseAlignment = cigarRead(fileHandleIn)) != NULL) {

        // Write out original cigar
    	if(!pairwiseAlignment->strand1) {
    		invertStrands(pairwiseAlignment);
    	}
    	checkPairwiseAlignment(pairwiseAlignment);
        cigarWrite(fileHandleOut, pairwiseAlignment, 0);

        // Write out mirror cigar (with query and target reversed)
        cigarReverse(pairwiseAlignment);
        if(!pairwiseAlignment->strand1) {
        	invertStrands(pairwiseAlignment);
        }
        checkPairwiseAlignment(pairwiseAlignment);
        cigarWrite(fileHandleOut, pairwiseAlignment, 0);

        // Cleanup
        destructPairwiseAlignment(pairwiseAlignment);
    }
    fclose(fileHandleIn);
    fclose(fileHandleOut);

    //while(1);

    return 0;
}
